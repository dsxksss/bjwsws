"""
SFE 流程 (4) 的准备：为每个单点突变生成 forward/reverse 的 Flex DDG mut 文件 +
tasks_sfe.csv（共 2N 个 Flex DDG 作业），submit_batch.py 可直接吃。

按文档协议：
  - Forward：取 WT 复合物的 N 个构象（帧），按单点突变算 ddG。
      输入 = WT 帧 PDB；mut 文件内容 = 突变名本身（如 DM31A，D→A）
  - Reverse：取突变体复合物的 N 个构象，虚拟恢复为 WT 算 ddG。
      输入 = 突变体帧 PDB；mut 文件内容 = 反向（如 AM31D，A→D，首尾字母对调）

帧文件命名约定（来自 sfe_split_frames.py）：
  - WT 帧：       <wt-frames-dir>/<wt-tag>_NN.pdb          （所有突变共享同一套）
  - 突变体帧：    <mut-frames-root>/<mut>/<mut>_NN.pdb      （每个突变一套）

task_name 编码方向与帧号（sfe_<mut>_fwd_NN / sfe_<mut>_rev_NN），供 sfe_fetch.py 归位。

用法：
    python sfe_prepare.py \\
        --mut-list prod/single_mut_list.txt \\
        --wt-frames-dir prod/frames_WT --wt-tag WT \\
        --mut-frames-root prod/frames_mut \\
        --n-frames 21 \\
        --out-dir prod
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

SINGLE_MUT_RE = re.compile(r"^([A-Z])([A-Za-z0-9])(\d+)([A-Z])$")


def parse_single_mut(name: str) -> tuple[str, str, int, str] | None:
    m = SINGLE_MUT_RE.match(name.strip())
    if not m:
        return None
    return m.group(1), m.group(2), int(m.group(3)), m.group(4)


def read_mut_list(path: Path) -> list[str]:
    out: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        out.append(line)
    return out


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="生成 SFE 的 forward/reverse Flex DDG mut 文件 + tasks_sfe.csv",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--mut-list", required=True, help="单点突变列表（每行 DM31A）")
    p.add_argument("--wt-frames-dir", required=True,
                   help="WT 帧目录（<wt-tag>_NN.pdb，所有突变共享）")
    p.add_argument("--wt-tag", default="WT", help="WT 帧文件名前缀（默认 WT）")
    p.add_argument("--mut-frames-root", required=True,
                   help="突变体帧根目录（每突变子目录 <root>/<mut>/<mut>_NN.pdb）")
    p.add_argument("--n-frames", type=int, default=21,
                   help="每组构象帧数（默认 21；测试可调小如 3）")
    p.add_argument("--out-dir", default="prod", help="输出目录（默认 ./prod）")
    p.add_argument("--limit", type=int, default=0,
                   help="只取前 N 个突变（0=全部）")
    p.add_argument("--strict", action="store_true",
                   help="缺帧文件时报错退出（默认仅告警，仍写入任务）")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    mut_list_path = Path(args.mut_list)
    if not mut_list_path.exists():
        print(f"❌ 找不到突变列表: {mut_list_path}", file=sys.stderr)
        return 2

    wt_dir = Path(args.wt_frames_dir)
    mut_root = Path(args.mut_frames_root)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    mutfile_dir = out_dir / "sfe_mutfiles"
    mutfile_dir.mkdir(exist_ok=True)

    muts = read_mut_list(mut_list_path)
    if args.limit > 0:
        muts = muts[: args.limit]
    n = args.n_frames

    rows: list[dict[str, str]] = []
    n_missing = 0
    skipped_multi = 0

    for mut in muts:
        parsed = parse_single_mut(mut)
        if not parsed:
            skipped_multi += 1
            print(f"⚠ 跳过非单点突变: {mut}")
            continue
        orig, chain, pos, new = parsed

        # forward / reverse mut 文件（强制 LF，避免 Rosetta 解析 CRLF 出错）
        fwd_content = f"{orig}{chain}{pos}{new}"   # 突变名本身
        rev_content = f"{new}{chain}{pos}{orig}"   # 首尾对调，恢复 WT
        fwd_file = mutfile_dir / f"fwd_{mut}.txt"
        rev_file = mutfile_dir / f"rev_{mut}.txt"
        fwd_file.write_text(fwd_content + "\n", encoding="utf-8", newline="\n")
        rev_file.write_text(rev_content + "\n", encoding="utf-8", newline="\n")

        mut_frames_dir = mut_root / mut
        for i in range(1, n + 1):
            wt_frame = wt_dir / f"{args.wt_tag}_{i:02d}.pdb"
            mut_frame = mut_frames_dir / f"{mut}_{i:02d}.pdb"
            if not wt_frame.exists():
                n_missing += 1
                msg = f"缺 WT 帧: {wt_frame}"
                if args.strict:
                    print(f"❌ {msg}", file=sys.stderr); return 2
                print(f"⚠ {msg}")
            if not mut_frame.exists():
                n_missing += 1
                msg = f"缺突变体帧: {mut_frame}"
                if args.strict:
                    print(f"❌ {msg}", file=sys.stderr); return 2
                print(f"⚠ {msg}")

            rows.append({
                "task_name":   f"sfe_{mut}_fwd_{i:02d}",
                "mutation":    mut,
                "module_name": "Flex DDG",
                "module_id":   "",
                "method_name": "",
                "flow_name":   "",
                "flow_id":     "",
                "params":      f"PDB File={wt_frame.as_posix()};Mut File={fwd_file.as_posix()}",
            })
            rows.append({
                "task_name":   f"sfe_{mut}_rev_{i:02d}",
                "mutation":    mut,
                "module_name": "Flex DDG",
                "module_id":   "",
                "method_name": "",
                "flow_name":   "",
                "flow_id":     "",
                "params":      f"PDB File={mut_frame.as_posix()};Mut File={rev_file.as_posix()}",
            })

    tasks_csv = out_dir / "tasks_sfe.csv"
    fieldnames = ["task_name", "mutation", "module_name", "module_id",
                  "method_name", "flow_name", "flow_id", "params"]
    with tasks_csv.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)

    n_muts = len({r["mutation"] for r in rows})
    print()
    print("=" * 50)
    print(f"SFE 任务清单已生成: {tasks_csv}")
    print("=" * 50)
    print(f"  突变数         : {n_muts}（每个 {n} 帧 × 正反 = {2*n} 个 Flex DDG 作业）")
    print(f"  作业总数       : {len(rows)}")
    if skipped_multi:
        print(f"  跳过多点突变   : {skipped_multi}")
    if n_missing:
        print(f"  ⚠ 缺失帧文件   : {n_missing}（MD/拆帧未完成？确认后再提交）")
    print()
    print("下一步：")
    print(f"  python submit_batch.py --config config.json --tasks {tasks_csv} -o sfe_jobs.jsonl")
    return 0


if __name__ == "__main__":
    sys.exit(main())
