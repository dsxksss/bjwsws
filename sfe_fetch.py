"""
SFE 流程 (4)→(5)：回收 forward/reverse 的 Flex DDG 作业结果，汇总成 42 列宽表
（mutation,chain,forward_01..N,reverse_01..N），可选直接调 step02c 出 sfe_scores.csv。

逻辑与 fetch_batch.py 一致（查状态 → Done 则下载 → 就地更新 jsonl），但聚合方式不同：
按 task_name 里编码的 方向(fwd/rev) 和 帧号(NN)，把每个作业的单个 ddG 归位到宽表单元格。

只有「正反 2N 个 ddG 全齐」的突变才写进宽表（step02c 不接受缺值）；不齐的会被列出。

用法：
    python sfe_fetch.py --config config.json --input sfe_jobs.jsonl \\
        --n-frames 21 --wide-output prod/sfe_flex_ddg_wide.csv \\
        --run-step02c --sfe-output prod/sfe_scores.csv
"""
from __future__ import annotations

import argparse
import csv
import json
import logging
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

from _common import (
    build_session,
    find_result_dir,
    log,
    now_iso,
    read_jsonl,
    rewrite_jsonl,
    setup_console_logging,
)

TERMINAL_STATUS = {"Done", "Abort", "Cancel", "Failed"}
TASK_RE = re.compile(r"^sfe_(?P<mut>[A-Za-z0-9]+)_(?P<dir>fwd|rev)_(?P<nn>\d+)$")
SINGLE_MUT_RE = re.compile(r"^([A-Z])([A-Za-z0-9])(\d+)([A-Z])$")

# step02c 在 vendored CRPCA 里的相对位置
_HERE = Path(__file__).resolve().parent
_CRPCA_ROOT = _HERE / "_crpca_download" / "CRPCA"
_STEP02C = _CRPCA_ROOT / "local_pipeline" / "step02c_compute_sfe_from_flex_ddg.py"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="回收 SFE 的 Flex DDG 作业并汇总成宽表")
    p.add_argument("--config", required=True, help="登录配置 JSON 路径")
    p.add_argument("--input", "-i", required=True, help="sfe submit 生成的 jsonl")
    p.add_argument("--n-frames", type=int, default=21, help="每组帧数（默认 21）")
    p.add_argument("--wide-output", default="sfe_flex_ddg_wide.csv",
                   help="42 列宽表输出路径")
    p.add_argument("--no-fetch", action="store_true",
                   help="跳过查询/下载，只用已下载结果汇总宽表")
    p.add_argument("--run-step02c", action="store_true",
                   help="宽表写好后直接调 step02c 出 sfe_scores.csv")
    p.add_argument("--sfe-output", default="sfe_scores.csv",
                   help="--run-step02c 时的 sfe_scores.csv 输出路径")
    p.add_argument("--step02c-python", default=None,
                   help="跑 step02c 用的 python（需有 numpy/pandas）。"
                        "默认用当前解释器；若 SDK 环境没装 numpy/pandas，"
                        "可指向另一个装了的 python")
    p.add_argument("--invert-for-sampling", action="store_true",
                   help="传给 step02c：输出 sfe=-sfe_ddg（越大越好）")
    return p.parse_args()


def _fetch_pending(session, records: list[dict], in_path: Path) -> Counter:
    """查状态 + 下载 Done 结果 + 就地更新 jsonl（与 fetch_batch 同逻辑）。"""
    stats: Counter = Counter()
    updated = False
    for rec in records:
        name = rec.get("task_name", "?")
        prev = rec.get("status")
        job_id = rec.get("job_id", -1)
        if prev in TERMINAL_STATUS:
            stats[f"Skip_{prev}"] += 1
            continue
        if job_id == -1:
            stats["Skip_Failed"] += 1
            continue
        try:
            job = session.query_job_by_id(job_id)
        except Exception as e:
            log.warning(f"? {name} (job {job_id}) 查询失败: {e}")
            stats["QueryError"] += 1
            continue
        if job is None:
            stats["NotFound"] += 1
            continue
        job_status = getattr(job, "status", None)
        if job_status in ("Done", "Abort", "Cancel"):
            try:
                ok = session.fetch_and_write_job_result(job_id, use_task_name=True)
            except Exception as e:
                log.error(f"✗ {name} (job {job_id}) 下载异常: {e}")
                stats["FetchError"] += 1
                continue
            rec["status"] = job_status
            rec["fetch_time"] = now_iso()
            rec["result_dir"] = f"{job_id}-{rec.get('module_name') or 'job'}"
            updated = True
            stats[job_status if ok else f"FetchFailed_{job_status}"] += 1
            log.info(f"{'✓' if job_status=='Done' else '⚠'} {name} (job {job_id}) -> {job_status}")
        else:
            stats["Doing"] += 1
            log.info(f"… {name} (job {job_id}) 仍在运行 (status={job_status})")
    if updated:
        rewrite_jsonl(in_path, records)
        log.info(f"已更新 {in_path}")
    return stats


def _read_single_ddg(result_dir: Path) -> float | None:
    """读 Flex DDG 结果里的 ddG.csv，取其单个 ddG 值。"""
    files = list(result_dir.rglob("ddG.csv"))
    if not files:
        return None
    with files[0].open(encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f):
            v = (row.get("ddG") or "").strip()
            if v:
                try:
                    return float(v)
                except ValueError:
                    log.warning(f"[SFE] {files[0]} ddG 非数值: {v!r}")
    return None


def _chain_of(mut: str) -> str:
    m = SINGLE_MUT_RE.match(mut)
    return m.group(2) if m else ""


def aggregate_wide(
    records: list[dict],
    result_base: Path,
    n_frames: int,
    wide_output: Path,
) -> dict[str, int]:
    """把每个 fwd/rev 帧作业的 ddG 归位到宽表，仅写正反全齐的突变。"""
    # cells[mut]["forward"|"reverse"][frame_idx] = ddg
    cells: dict[str, dict[str, dict[int, float]]] = defaultdict(
        lambda: {"forward": {}, "reverse": {}})

    n_parsed = n_no_dir = n_no_ddg = 0
    for rec in records:
        if rec.get("status") != "Done":
            continue
        m = TASK_RE.match(rec.get("task_name", ""))
        if not m:
            continue
        mut = m.group("mut")
        direction = "forward" if m.group("dir") == "fwd" else "reverse"
        idx = int(m.group("nn"))

        rdir = find_result_dir(result_base, rec.get("job_id"))
        if not rdir:
            n_no_dir += 1
            continue
        ddg = _read_single_ddg(rdir)
        if ddg is None:
            n_no_ddg += 1
            continue
        cells[mut][direction][idx] = ddg
        n_parsed += 1

    fwd_cols = [f"forward_{i:02d}" for i in range(1, n_frames + 1)]
    rev_cols = [f"reverse_{i:02d}" for i in range(1, n_frames + 1)]
    header = ["mutation", "chain"] + fwd_cols + rev_cols

    complete: list[str] = []
    incomplete: list[tuple[str, int, int]] = []
    rows: list[list] = []
    for mut in sorted(cells.keys()):
        fwd = cells[mut]["forward"]
        rev = cells[mut]["reverse"]
        nf = sum(1 for i in range(1, n_frames + 1) if i in fwd)
        nr = sum(1 for i in range(1, n_frames + 1) if i in rev)
        if nf == n_frames and nr == n_frames:
            row = [mut, _chain_of(mut)]
            row += [fwd[i] for i in range(1, n_frames + 1)]
            row += [rev[i] for i in range(1, n_frames + 1)]
            rows.append(row)
            complete.append(mut)
        else:
            incomplete.append((mut, nf, nr))

    wide_output.parent.mkdir(parents=True, exist_ok=True)
    with wide_output.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)

    if incomplete:
        print(f"  ⚠ {len(incomplete)} 个突变正反未齐（需各 {n_frames} 个），未写入宽表：")
        for mut, nf, nr in incomplete[:10]:
            print(f"      {mut}: forward {nf}/{n_frames}, reverse {nr}/{n_frames}")
    return {
        "ddg_parsed": n_parsed, "complete_muts": len(complete),
        "incomplete_muts": len(incomplete), "no_result_dir": n_no_dir, "no_ddg": n_no_ddg,
    }


def _load_result_dir(config_path: str) -> str:
    cfg = json.loads(Path(config_path).read_text(encoding="utf-8"))
    return cfg.get("result_dir", "wemol_results")


def main() -> int:
    setup_console_logging(level=logging.INFO)
    args = parse_args()

    in_path = Path(args.input)
    if not in_path.exists():
        log.error(f"找不到输入文件: {in_path}")
        return 2

    if not args.no_fetch:
        session = build_session(args.config)
        records = read_jsonl(in_path)
        log.info(f"读到 {len(records)} 条记录，开始查状态/下载")
        st = _fetch_pending(session, records, in_path)
        print()
        print("=" * 50)
        print("回收汇总")
        print("=" * 50)
        print(f"  ✓ Done : {st['Done']}   ⚠ Abort/Cancel : {st['Abort']+st['Cancel']}")
        print(f"  … Doing(下次再扫) : {st['Doing']}   ? 异常 : {st['QueryError']+st['NotFound']+st['FetchError']}")
        if st["Doing"]:
            print(f"\n还有 {st['Doing']} 个未完成；等一会儿重跑本脚本。宽表只汇总已 Done 的。")

    latest = read_jsonl(in_path)
    result_base = Path(_load_result_dir(args.config))
    wide_out = Path(args.wide_output)

    print()
    print("=" * 50)
    print(f"汇总宽表 → {wide_out}")
    print("=" * 50)
    agg = aggregate_wide(latest, result_base, args.n_frames, wide_out)
    print(f"  归位 ddG 数      : {agg['ddg_parsed']}")
    print(f"  正反齐全的突变   : {agg['complete_muts']}")
    print(f"  未齐全的突变     : {agg['incomplete_muts']}")
    print(f"  找不到结果目录   : {agg['no_result_dir']}   缺 ddG.csv : {agg['no_ddg']}")
    print(f"  宽表             : {wide_out}")

    if args.run_step02c:
        if agg["complete_muts"] == 0:
            print("\n⚠ 没有正反齐全的突变，跳过 step02c。")
            return 0
        if not _STEP02C.exists():
            print(f"\n❌ 找不到 step02c: {_STEP02C}", file=sys.stderr)
            return 2
        sfe_out = Path(args.sfe_output)
        step02c_py = args.step02c_python or sys.executable
        cmd = [step02c_py, str(_STEP02C),
               "--input", str(wide_out), "--output", str(sfe_out),
               "--n-conformations", str(args.n_frames)]
        if args.invert_for_sampling:
            cmd.append("--invert-for-sampling")
        env_note = f"PYTHONPATH={_CRPCA_ROOT}"
        print(f"\n运行 step02c（{env_note}）：\n  {' '.join(cmd)}")
        import os
        env = dict(os.environ)
        env["PYTHONPATH"] = str(_CRPCA_ROOT) + os.pathsep + env.get("PYTHONPATH", "")
        r = subprocess.run(cmd, env=env)
        if r.returncode != 0:
            print("⚠ step02c 返回非 0，请检查上面的报错。", file=sys.stderr)
            return r.returncode
        print(f"✓ sfe_scores.csv → {sfe_out}")
    else:
        print(f"\n下一步（出 sfe 分数）：")
        print(f"  PYTHONPATH={_CRPCA_ROOT} python {_STEP02C} \\")
        print(f"      --input {wide_out} --output {args.sfe_output} --n-conformations {args.n_frames}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
