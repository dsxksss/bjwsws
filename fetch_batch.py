"""
批量回收 Wemol 作业结果 + 汇总成总表。

逻辑：
- 读 submit_batch.py 生成的 jsonl
- 对每条非终态记录，查询作业状态
- 状态为 Done/Abort/Cancel 时下载结果，并就地更新 jsonl
- 仍在 Doing 的保持不动，下次再扫
- 最后扫所有已完成的结果目录，按模块解析打分，输出总表 scores.csv

可重复运行：每次都会用当前可用的 Done 记录覆盖 scores.csv。
跑完一次还有 Doing 的话，等一段时间再跑本脚本就行。

示例：
    python fetch_batch.py --config config.json \\
                          --input submitted_jobs.jsonl \\
                          --output scores.csv
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from collections import Counter
from pathlib import Path

from _common import (
    aggregate_to_csv,
    build_session,
    log,
    now_iso,
    read_jsonl,
    rewrite_jsonl,
    setup_console_logging,
)


# 终态：这些状态不再处理
TERMINAL_STATUS = {"Done", "Abort", "Cancel", "Failed"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="批量回收 Wemol 作业结果")
    p.add_argument("--config", required=True, help="登录配置 JSON 路径")
    p.add_argument(
        "--input", "-i",
        required=True,
        help="submit_batch.py 生成的 jsonl 文件路径",
    )
    p.add_argument(
        "--output", "-o",
        default="scores.csv",
        help="汇总结果 CSV 路径（默认 scores.csv）",
    )
    p.add_argument(
        "--include-multipoint",
        action="store_true",
        help="保留多点突变（默认仅保留单点）",
    )
    p.add_argument(
        "--no-aggregate",
        action="store_true",
        help="跳过最终汇总步骤，只做 fetch",
    )
    p.add_argument(
        "--foldx-interface",
        default=None,
        help="FoldX interface to use for the run_pipeline foldx column, e.g. M_N.",
    )
    p.add_argument(
        "--legacy-output",
        action="store_true",
        help="Write legacy raw module columns instead of the run_pipeline input schema.",
    )
    return p.parse_args()


def main() -> int:
    setup_console_logging(level=logging.INFO)
    args = parse_args()

    in_path = Path(args.input)
    if not in_path.exists():
        log.error(f"找不到输入文件: {in_path}")
        return 2

    session = build_session(args.config)
    records = read_jsonl(in_path)
    log.info(f"读到 {len(records)} 条记录，开始扫描")

    stats: Counter = Counter()
    updated = False

    for rec in records:
        name = rec.get("task_name", "?")
        prev_status = rec.get("status")
        job_id = rec.get("job_id", -1)

        # 已是终态 / 提交时就失败的，跳过
        if prev_status in TERMINAL_STATUS:
            stats[f"Skip_{prev_status}"] += 1
            continue
        if job_id == -1:
            stats["Skip_Failed"] += 1
            continue

        # 查作业状态
        try:
            job = session.query_job_by_id(job_id)
        except Exception as e:
            log.warning(f"? {name} (job {job_id}) 查询失败: {e}")
            stats["QueryError"] += 1
            continue

        if job is None:
            log.warning(f"? {name} (job {job_id}) 查询不到，跳过")
            stats["NotFound"] += 1
            continue

        job_status = getattr(job, "status", None)

        if job_status in ("Done", "Abort", "Cancel"):
            ok = False
            try:
                ok = session.fetch_and_write_job_result(job_id, use_task_name=True)
            except Exception as e:
                log.error(f"✗ {name} (job {job_id}) 下载结果异常: {e}")
                stats["FetchError"] += 1
                continue

            rec["status"]     = job_status
            rec["fetch_time"] = now_iso()
            # SDK 把结果落到 RESULT_FILE_SAVE_DIR/<job_id>-<module_name>/...
            # 这里只记录顶层目录提示用户去哪儿找；具体子目录由 SDK 决定
            rec["result_dir"] = f"{job_id}-{rec.get('module_name') or rec.get('flow_name') or 'job'}"
            updated = True

            if ok:
                stats[job_status] += 1
                symbol = "✓" if job_status == "Done" else "⚠"
                log.info(f"{symbol} {name} (job {job_id}) -> {job_status}")
            else:
                stats[f"FetchFailed_{job_status}"] += 1
                log.warning(f"⚠ {name} (job {job_id}) 状态={job_status} 但下载失败")
        else:
            # Doing / 其他中间态：保持原状，下次再扫
            stats["Doing"] += 1
            log.info(f"… {name} (job {job_id}) 仍在运行 (status={job_status})")

    if updated:
        rewrite_jsonl(in_path, records)
        log.info(f"已更新 {in_path}")

    # ─── 汇总 ───────────────────────────────────────────────
    print()
    print("=" * 50)
    print("回收汇总")
    print("=" * 50)
    print(f"  ✓ Done                 : {stats['Done']}")
    print(f"  ⚠ Abort                : {stats['Abort']}")
    print(f"  ⚠ Cancel               : {stats['Cancel']}")
    print(f"  … Doing (下次再扫)     : {stats['Doing']}")
    print(f"  - 已 Done 跳过         : {stats['Skip_Done']}")
    print(f"  - 已 Abort/Cancel 跳过 : {stats['Skip_Abort'] + stats['Skip_Cancel']}")
    print(f"  - 提交即失败           : {stats['Skip_Failed']}")
    print(f"  ? 查询/下载异常        : {stats['QueryError'] + stats['NotFound'] + stats['FetchError']}")
    print()

    if stats["Doing"] > 0:
        print(f"还有 {stats['Doing']} 条作业未完成；等一段时间后重新运行本脚本即可。")
    elif (stats["Done"] + stats["Skip_Done"]) == len(records):
        print("所有作业已完成并下载。")

    # ─── 汇总打分到总表 ──────────────────────────────────────
    if args.no_aggregate:
        return 0

    # 读回最新的 records（fetch 可能更新过 status）
    latest = read_jsonl(in_path)
    result_base = Path(_load_result_dir_from_config(args.config))
    out_csv = Path(args.output)

    log.info(f"开始汇总打分到 {out_csv}（结果目录: {result_base}）")
    agg_stats = aggregate_to_csv(
        records=latest,
        result_dir_base=result_base,
        output_csv=out_csv,
        include_multipoint=args.include_multipoint,
        pipeline_format=not args.legacy_output,
        foldx_interface=args.foldx_interface,
    )

    print()
    print("=" * 50)
    print("打分汇总")
    print("=" * 50)
    print(f"  解析的模块作业数 : {agg_stats['parsers_called']}")
    print(f"  缺少解析器       : {agg_stats['parser_missing']}")
    print(f"  找不到结果目录   : {agg_stats['no_result_dir']}")
    print(f"  突变行数         : {agg_stats['muts_parsed']}")
    print(f"  打分列数         : {agg_stats['columns']}")
    print(f"  输出文件         : {out_csv}")
    print(f"  FoldX ambiguous blanks : {agg_stats.get('foldx_ambiguous', 0)}")
    return 0


def _load_result_dir_from_config(config_path: str) -> str:
    """从 config.json 读 result_dir 字段（与 build_session 内部默认保持一致）。"""
    cfg = json.loads(Path(config_path).read_text(encoding="utf-8"))
    return cfg.get("result_dir", "wemol_results")


if __name__ == "__main__":
    sys.exit(main())
