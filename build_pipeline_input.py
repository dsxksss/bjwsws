"""
桥接 stage2(打分汇总) → stage3(CRPCA run_pipeline)。

把 fetch_batch.py 产出的 scores.csv 转成 run_pipeline.py 能直接吃的
single_point_scores.csv：
  1. 修正 location —— fetch_batch 里 location 是 PDB resSeq；run_pipeline 要的是
     master 序列（抗体链按 --antibody-chains 顺序拼接）的 1-based 线性位置。
     用 CRPCA 自己的 extract_master_from_pdb 计算，零复制风险。
  2. 校验 original_aa 与 master 一致（提前抓错，免得 run_pipeline 才报）。
  3. （可选）合并 SFE —— 把 step02c 产出的 sfe_scores.csv 的分值并进 sfe 列。

依赖：biopython（CRPCA 也依赖它）。不需要 torch/gpytorch。

示例：
    python build_pipeline_input.py \\
        --scores prod/scores.csv \\
        --pdb test_files/7l7e-ab-ag-complex.pdb \\
        --antibody-chains M,N \\
        --sfe prod/sfe_scores.csv \\
        --output prod/single_point_scores.csv
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

# 五工具列（与 CRPCA scoring.SCORE_COLUMNS 一致）
SCORE_COLUMNS = ["sfe", "fep", "rosetta_flex", "foldx", "abnativ"]
OUT_FIELDS = ["mutation", "chain", "location", "original_aa", "mutant_aa"] + SCORE_COLUMNS


def _crpca_default_root() -> Path:
    return Path(__file__).resolve().parent / "_crpca_download" / "CRPCA"


def build_residue_map(pdb: str, antibody_chains: list[str], crpca_root: Path):
    """
    用 CRPCA 的 extract_master_from_pdb 建 (chain, resSeq) → (linear_index, master_aa) 映射。
    返回 (mapping_dict, master_len)。
    """
    if str(crpca_root) not in sys.path:
        sys.path.insert(0, str(crpca_root))
    try:
        from local_pipeline.common.pdb_structure import extract_master_from_pdb
    except ImportError as e:
        print(f"❌ 无法导入 CRPCA pdb_structure（检查 --crpca-root 与 biopython）: {e}",
              file=sys.stderr)
        raise SystemExit(2)

    master, mapping = extract_master_from_pdb(pdb, antibody_chains)
    by_key: dict[tuple[str, int], tuple[int, str]] = {}
    for r in mapping:
        by_key[(r["chain"], int(r["resseq"]))] = (int(r["linear_index"]), r["aa"])
    return by_key, len(master)


def load_sfe(sfe_path: Path) -> dict[str, float]:
    """读 step02c 产出的 sfe_scores.csv，返回 {mutation: sfe}。优先 sfe 列，退回 sfe_ddg。"""
    out: dict[str, float] = {}
    with sfe_path.open(encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        col = "sfe" if "sfe" in (reader.fieldnames or []) else "sfe_ddg"
        if col not in (reader.fieldnames or []):
            raise ValueError(f"{sfe_path} 既无 sfe 也无 sfe_ddg 列")
        for row in reader:
            mut = (row.get("mutation") or "").strip()
            val = (row.get(col) or "").strip()
            if not mut or not val:
                continue
            try:
                out[mut] = float(val)
            except ValueError:
                pass
    return out


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="把 scores.csv 转成 run_pipeline 就绪的 single_point_scores.csv")
    p.add_argument("--scores", required=True, help="fetch_batch.py 产出的 scores.csv")
    p.add_argument("--pdb", required=True, help="PDB 结构（与 run_pipeline 用同一个）")
    p.add_argument("--antibody-chains", required=True,
                   help="抗体链，按 master 拼接顺序，如 M,N（重链在前）")
    p.add_argument("--sfe", default=None, help="（可选）step02c 产出的 sfe_scores.csv，合并进 sfe 列")
    p.add_argument("--output", "-o", default="single_point_scores.csv", help="输出 CSV")
    p.add_argument("--crpca-root", default=None,
                   help="CRPCA 仓库根目录（默认 ./_crpca_download/CRPCA）")
    p.add_argument("--strict", action="store_true",
                   help="original_aa 与 master 不一致时直接报错（默认仅警告并跳过该行）")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    crpca_root = Path(args.crpca_root) if args.crpca_root else _crpca_default_root()
    ab_chains = [c.strip() for c in args.antibody_chains.split(",") if c.strip()]

    res_map, master_len = build_residue_map(args.pdb, ab_chains, crpca_root)
    print(f"master 长度: {master_len}；抗体链顺序: {ab_chains}")

    sfe_map = load_sfe(Path(args.sfe)) if args.sfe else {}
    if args.sfe:
        print(f"读入 SFE: {len(sfe_map)} 条 ← {args.sfe}")

    in_rows = list(csv.DictReader(Path(args.scores).open(encoding="utf-8-sig", newline="")))
    out_rows: list[dict] = []
    n_fixed = n_skip = n_sfe_filled = 0

    for row in in_rows:
        mut = (row.get("mutation") or "").strip()
        chain = (row.get("chain") or "").strip()
        orig = (row.get("original_aa") or "").strip()
        try:
            resseq = int(str(row.get("location")).strip())
        except (ValueError, TypeError):
            print(f"⚠ {mut}: location 非整数，跳过", file=sys.stderr)
            n_skip += 1
            continue

        key = (chain, resseq)
        if key not in res_map:
            msg = f"{mut}: (chain={chain}, resSeq={resseq}) 不在 master 映射中"
            if args.strict:
                print(f"❌ {msg}", file=sys.stderr); return 2
            print(f"⚠ {msg}，跳过", file=sys.stderr)
            n_skip += 1
            continue

        linear_index, master_aa = res_map[key]
        if orig and orig != master_aa:
            msg = f"{mut}: original_aa={orig} 与 master[{linear_index}]={master_aa} 不一致"
            if args.strict:
                print(f"❌ {msg}", file=sys.stderr); return 2
            print(f"⚠ {msg}，跳过", file=sys.stderr)
            n_skip += 1
            continue

        out = {k: (row.get(k, "") or "") for k in OUT_FIELDS}
        out["location"] = linear_index
        out["original_aa"] = master_aa  # 以 master 为准
        # 合并 SFE（仅当原 sfe 为空且 sfe_map 有值）
        if (not out.get("sfe")) and mut in sfe_map:
            out["sfe"] = sfe_map[mut]
            n_sfe_filled += 1
        out_rows.append(out)
        n_fixed += 1

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=OUT_FIELDS)
        w.writeheader()
        w.writerows(out_rows)

    print()
    print("=" * 50)
    print(f"已生成 run_pipeline 输入: {out_path}")
    print("=" * 50)
    print(f"  写出行数        : {n_fixed}")
    print(f"  跳过(映射/校验) : {n_skip}")
    print(f"  SFE 填入        : {n_sfe_filled}")
    print()
    print("下一步（在 CRPCA conda 环境 / 容器内）：")
    print(f"  python run_pipeline.py --single-point-scores {out_path} \\")
    print(f"      --pdb {args.pdb} --antibody-chains {args.antibody_chains} --output-dir work/run1")
    return 0


if __name__ == "__main__":
    sys.exit(main())
