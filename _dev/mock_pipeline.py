"""
本地 mock 工具：不连 Wemol，给 tasks.csv 生成假的 jsonl + 假的结果目录文件。
用于验证 submit / fetch / aggregate 整条 pipeline 的连通性。

工作流：
    python prepare_inputs.py --mut-list single_mut_list.txt ... -o prep/
    python mock_pipeline.py  --tasks prep/tasks.csv -o submitted_jobs.jsonl \\
                              --result-dir wemol_results
    python aggregate_results.py -i submitted_jobs.jsonl --result-dir wemol_results

打分值是随机数（FlexDDG / FoldX 用 -5~+5，FEP 用 -10~+10），仅用于打通流程。
"""
from __future__ import annotations

import argparse
import csv
import json
import logging
import random
import re
import sys
from datetime import datetime
from pathlib import Path


log = logging.getLogger("mock")


SINGLE_MUT_RE = re.compile(r"^([A-Z])([A-Za-z0-9])(\d+)([A-Z])$")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="生成假数据验证 pipeline 连通性")
    p.add_argument("--tasks", required=True,
                   help="prepare_inputs.py 生成的 tasks.csv")
    p.add_argument("--output", "-o", default="submitted_jobs.jsonl",
                   help="假 jsonl 路径（默认 submitted_jobs.jsonl）")
    p.add_argument("--result-dir", default="wemol_results",
                   help="假结果目录的根路径（与 config 的 result_dir 一致）")
    p.add_argument("--seed", type=int, default=42,
                   help="随机种子，方便重现")
    p.add_argument("--foldx-chains", default="M,N,S",
                   help="模拟生成 FoldX Interface_X_Y.csv 时用的链字母（逗号分隔，默认 M,N,S）")
    return p.parse_args()


def load_tasks(path: Path) -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def parse_params(s: str) -> dict[str, str]:
    out = {}
    for chunk in (s or "").split(";"):
        chunk = chunk.strip()
        if not chunk or "=" not in chunk:
            continue
        k, v = chunk.split("=", 1)
        out[k.strip()] = v.strip()
    return out


def write_flexddg_result(rdir: Path, mutation: str) -> None:
    """ddG.csv: mut,ddG 两列一行。"""
    (rdir).mkdir(parents=True, exist_ok=True)
    ddg = round(random.uniform(-5, 5), 4)
    with (rdir / "ddG.csv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["mut", "ddG"])
        w.writerow([mutation, ddg])


def write_foldx_result(rdir: Path, mut_file: Path | None, chains: list[str]) -> None:
    """
    给每对链生成一个 Interface_X_Y.csv，包含 WT 行 + 所有突变行。
    突变列表从 mut_file（mutlist_foldx.txt）读取（每行 `TM31A;`）。
    """
    rdir.mkdir(parents=True, exist_ok=True)
    # 读 FoldX mut 列表
    muts: list[str] = []
    if mut_file and mut_file.exists():
        for line in mut_file.read_text(encoding="utf-8").splitlines():
            line = line.strip().rstrip(";").strip()
            if line:
                muts.append(line)
    if not muts:
        log.warning(f"[FoldX mock] {mut_file} 为空或不存在，跳过")
        return

    pairs = [(chains[i], chains[j]) for i in range(len(chains)) for j in range(i + 1, len(chains))]
    fname_idx = 0
    for c1, c2 in pairs:
        wt_energy = round(random.uniform(-50, -10), 4)
        path = rdir / f"Interface_{c1}_{c2}.csv"
        with path.open("w", encoding="utf-8", newline="") as f:
            w = csv.writer(f)
            w.writerow(["Mutation", "File Name", "Chain1 Name", "Chain2 Name",
                        "Interaction Energy", "deltaEnergy"])
            w.writerow(["WT", "mock_Repair", c1, c2, wt_energy, 0.0])
            for mut in muts:
                fname_idx += 1
                delta = round(random.uniform(-3, 3), 4)
                w.writerow([mut, f"mutation_{fname_idx:03d}", c1, c2,
                            round(wt_energy + delta, 4), delta])


def write_fep_result(rdir: Path) -> None:
    """result.txt: 3 行 'X dG: <num>'。"""
    rdir.mkdir(parents=True, exist_ok=True)
    ligand = round(random.uniform(-100, -50), 6)
    complex_ = round(random.uniform(-100, -50), 6)
    final_ddg = round(random.uniform(-10, 10), 6)
    (rdir / "result.txt").write_text(
        f"ligand dG: {ligand}\n"
        f"complex dG: {complex_}\n"
        f"final ddG: {final_ddg}\n",
        encoding="utf-8",
    )


def main() -> int:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s",
                        datefmt="%H:%M:%S")
    args = parse_args()
    random.seed(args.seed)

    tasks_path = Path(args.tasks)
    tasks = load_tasks(tasks_path)
    if not tasks:
        log.error(f"任务清单为空: {tasks_path}")
        return 2

    base = Path(args.result_dir)
    base.mkdir(parents=True, exist_ok=True)

    out_path = Path(args.output)
    if out_path.exists():
        out_path.unlink()

    chains = [c.strip() for c in args.foldx_chains.split(",") if c.strip()]
    log.info(f"mock 假 FoldX 链对: {chains}")

    fake_job_id = 90000
    n = {"Flex DDG": 0, "Mutation Energy of Binding": 0, "Protein FEP": 0, "skip": 0}

    with out_path.open("w", encoding="utf-8") as fout:
        for task in tasks:
            module = task.get("module_name") or task.get("flow_name") or ""
            if module not in ("Flex DDG", "Mutation Energy of Binding", "Protein FEP"):
                n["skip"] += 1
                continue

            fake_job_id += 1
            params = parse_params(task.get("params", ""))
            rdir = base / f"{fake_job_id}-{module}" / f"{fake_job_id}-mock"

            if module == "Flex DDG":
                mut = task.get("mutation") or ""
                write_flexddg_result(rdir, mut)
            elif module == "Mutation Energy of Binding":
                mut_file = Path(params.get("Mutant File", "")) if params.get("Mutant File") else None
                write_foldx_result(rdir, mut_file, chains)
            elif module == "Protein FEP":
                write_fep_result(rdir)

            n[module] += 1

            rec = {
                "task_name":   task.get("task_name"),
                "mutation":    task.get("mutation"),
                "module_name": module,
                "module_id":   task.get("module_id") or None,
                "flow_name":   task.get("flow_name") or None,
                "flow_id":     task.get("flow_id") or None,
                "method_name": task.get("method_name") or None,
                "params":      params,
                "job_id":      fake_job_id,
                "submit_time": datetime.now().strftime("%Y-%m-%dT%H:%M:%S"),
                "status":      "Done",
                "fetch_time":  datetime.now().strftime("%Y-%m-%dT%H:%M:%S"),
                "result_dir":  f"{fake_job_id}-{module}",
                "error":       None,
                "_mock":       True,
            }
            fout.write(json.dumps(rec, ensure_ascii=False) + "\n")

    print()
    print("=" * 50)
    print(f"Mock 数据已生成")
    print("=" * 50)
    print(f"  假 jsonl       : {out_path}")
    print(f"  假结果目录根   : {base}")
    print(f"  FlexDDG 作业    : {n['Flex DDG']}")
    print(f"  FoldX   作业    : {n['Mutation Energy of Binding']}")
    print(f"  FEP     作业    : {n['Protein FEP']}")
    print(f"  跳过（非支持） : {n['skip']}")
    print()
    print("下一步：")
    print(f"  python aggregate_results.py -i {out_path} --result-dir {base} -o final_scores.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
