"""
批量提交 Wemol 作业脚本。

特性：
- 读取 CSV / JSON 任务清单
- 串行提交（use_delay=False，立即拿到 job_id）
- job_id 实时落盘到 jsonl —— 进程中断可断点续跑
- 自动捕获 JobRunMaxNumLimit 并 sleep 重试
- 已成功提交过的 task_name 自动跳过（幂等）

示例：
    python submit_batch.py \\
        --config config.json \\
        --tasks tasks.csv \\
        --output submitted_jobs.jsonl \\
        --delay 10 \\
        --retry-delay 120
"""
from __future__ import annotations

import argparse
import logging
import sys
from collections import Counter
from pathlib import Path
from time import sleep

from _common import (
    already_submitted,
    append_jsonl,
    build_session,
    load_tasks,
    log,
    now_iso,
    read_jsonl,
    setup_console_logging,
    submit_one,
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="批量提交 Wemol 作业")
    p.add_argument("--config", required=True, help="登录配置 JSON 路径")
    p.add_argument("--tasks", required=True, help="任务清单文件（.csv 或 .json）")
    p.add_argument(
        "-o", "--output",
        default="submitted_jobs.jsonl",
        help="提交记录 jsonl 输出路径（默认: submitted_jobs.jsonl，已存在则追加并断点续跑）",
    )
    p.add_argument(
        "--delay",
        type=float, default=10.0,
        help="每提交一次作业后的固定等待秒数（默认 10）",
    )
    p.add_argument(
        "--retry-delay",
        type=float, default=120.0,
        help="遇到 JobRunMaxNumLimit 后等待多少秒重试（默认 120）",
    )
    return p.parse_args()


def main() -> int:
    setup_console_logging(level=logging.INFO)
    args = parse_args()

    session = build_session(args.config)
    tasks = load_tasks(args.tasks)
    output = Path(args.output)

    log.info(f"任务清单共 {len(tasks)} 条；记录文件 {output}")

    stats = Counter()
    for task in tasks:
        name = task["task_name"]

        # 幂等：已成功提交过的跳过
        if already_submitted(output, name):
            log.info(f"⤳ 跳过已提交: {name}")
            stats["Skip"] += 1
            continue

        # 提交循环：JobRunMaxNumLimit 重试，其他异常即时记录失败
        while True:
            try:
                job_id = submit_one(session, task)

                record = {
                    "task_name":   name,
                    "mutation":    task.get("mutation"),
                    "module_name": task.get("module_name"),
                    "module_id":   task.get("module_id"),
                    "flow_name":   task.get("flow_name"),
                    "flow_id":     task.get("flow_id"),
                    "method_name": task.get("method_name"),
                    "params":      task.get("params", {}),
                    "job_id":      job_id,
                    "submit_time": now_iso(),
                    "status":      "Submitted" if job_id != -1 else "Failed",
                    "fetch_time":  None,
                    "result_dir":  None,
                    "error":       None if job_id != -1 else "run_job 返回 -1",
                }
                append_jsonl(output, record)

                if job_id != -1:
                    log.info(f"✓ {name} -> job_id={job_id}")
                    stats["Submitted"] += 1
                else:
                    log.error(f"✗ {name}: SDK 返回 job_id=-1")
                    stats["Failed"] += 1

                sleep(args.delay)
                break

            except Exception as e:
                err_msg = str(e.args[0]) if e.args else str(e)
                err_meta = e.args[1] if len(e.args) > 1 else {}

                if "JobRunMaxNumLimit" in err_msg:
                    running = err_meta.get("User.JobRunNum", "?")
                    limit = err_meta.get("Limit.JobRunMaxNum", "?")
                    log.warning(
                        f"⏸ 平台达到并发上限 ({running}/{limit})，"
                        f"sleep {args.retry_delay}s 后重试 {name}"
                    )
                    sleep(args.retry_delay)
                    continue
                else:
                    record = {
                        "task_name":   name,
                        "module_name": task.get("module_name"),
                        "module_id":   task.get("module_id"),
                        "flow_name":   task.get("flow_name"),
                        "flow_id":     task.get("flow_id"),
                        "method_name": task.get("method_name"),
                        "params":      task.get("params", {}),
                        "job_id":      -1,
                        "submit_time": now_iso(),
                        "status":      "Failed",
                        "fetch_time":  None,
                        "result_dir":  None,
                        "error":       f"{type(e).__name__}: {err_msg}",
                    }
                    append_jsonl(output, record)
                    log.error(f"✗ {name}: {err_msg}")
                    stats["Failed"] += 1
                    break

    # ─── 汇总 ───────────────────────────────────────────────
    total_in_file = len(read_jsonl(output))
    print()
    print("=" * 50)
    print("提交汇总")
    print("=" * 50)
    print(f"  本次提交成功 (Submitted) : {stats['Submitted']}")
    print(f"  本次提交失败 (Failed)    : {stats['Failed']}")
    print(f"  本次跳过 (Skip 已提交)   : {stats['Skip']}")
    print(f"  记录文件中累计条数       : {total_in_file}  ({output})")
    print()
    print(f"下一步：等作业跑完后运行")
    print(f"    python fetch_batch.py --config {args.config} --input {output}")
    return 0 if stats["Failed"] == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
