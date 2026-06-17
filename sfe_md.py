"""
SFE 流程 (1)(2)(3) 的 MD 编排：WT 与每个突变体各跑一次 MD，抽帧成单帧 PDB。

链路（每个目标一条，WT 无 PDB Mutation 这步）：
    [PDB Mutation 527] → 突变体 PDB → [flow 295 MD] → 轨迹 → [MD Trajectory 抽帧多帧PDB]
                                                                → 本地拆帧 frames_<tag>/<tag>_NN.pdb

设计为「可重复运行的断点续跑编排器」：状态记录在 --state jsonl，每跑一次就把每个目标
往前推进一个就绪的阶段（提交 / 查状态 / 下载 / 拆帧）。MD 在单卡上慢，跑一次推进不动就
等会儿再跑本脚本，直到所有目标 extracted。

flow 295 只覆盖必要参数（输入 PDB、链、模拟时长、抽帧 Skip/Type），其余保持平台默认。

⚠ 首次真实运行后需核对下载结果里的文件名：突变体 PDB 与多帧轨迹 PDB 的实际命名可能与
默认 glob 不符，用 --mutant-pdb-glob / --traj-glob 调整（脚本会把找到的候选打印出来）。

用法（测试小跑）：
    python sfe_md.py --config config.json \\
        --pdb test_files/7l7e-ab-ag-complex.pdb \\
        --mut-list prod/single_mut_list.txt \\
        --receptor-chain S --ligand-chain M,N \\
        --sim-ns 1 --skip-ps 500 \\
        --out-dir prod --frames-root prod/frames_mut --wt-frames-dir prod/frames_WT
    # 反复运行，直到输出 “全部目标已 extracted”
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
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
from sfe_split_frames import split_frames
from sfe_rechain import build_ref_order, rechain_pdb

# SDK run_job/run_job_use_flow 用的是「名字」不是 id
FLOW_MD_NAME = "Protein MD Simulation v2.1"   # flow 295
MOD_PDB_MUTATION_NAME = "PDB Mutation"        # module 527
TASK_PDB_PREPARE = "MD PDB Prepare"
TASK_MD_MDP = "MD MDP Generation"
TASK_TRAJECTORY = "MD Trajectory v2"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="SFE 的 MD 编排（WT+突变体各跑 MD 抽帧），断点续跑",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--config", required=True, help="登录配置 JSON")
    p.add_argument("--pdb", required=True, help="WT 复合物 PDB")
    p.add_argument("--mut-list", required=True, help="单点突变列表（每行 DM31A）")
    p.add_argument("--receptor-chain", default="S", help="MD PDB Prepare 的 Receptor Chain（默认 S=抗原）")
    p.add_argument("--ligand-chain", default="M,N", help="MD PDB Prepare 的 Ligand Chain（默认 M,N=抗体）")
    p.add_argument("--sim-ns", default="", help="MD 模拟时长 ns（空=平台默认；文档要 10）")
    p.add_argument("--skip-ps", default="500", help="MD Trajectory 抽帧间隔 ps（默认 500）")
    p.add_argument("--n-frames", type=int, default=21, help="期望帧数（仅用于拆帧告警，默认 21）")
    p.add_argument("--out-dir", default="prod", help="中间文件目录")
    p.add_argument("--wt-frames-dir", default="prod/frames_WT", help="WT 帧输出目录")
    p.add_argument("--frames-root", default="prod/frames_mut", help="突变体帧根目录（<root>/<mut>/）")
    p.add_argument("--state", default="sfe_md_state.jsonl", help="编排状态 jsonl")
    p.add_argument("--limit", type=int, default=0, help="只处理前 N 个突变（0=全部）")
    p.add_argument("--mutant-pdb-glob", default="*.pdb",
                   help="在 PDB Mutation 结果里找突变体 PDB 的 glob（默认 *.pdb）")
    p.add_argument("--traj-glob", default="md_center.pdb",
                   help="在 flow 295 结果里找多帧轨迹 PDB 的 glob"
                        "（默认 md_center.pdb，即 MD Trajectory v2 的多帧输出）")
    return p.parse_args()


def read_mut_list(path: Path, limit: int) -> list[str]:
    out: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            out.append(line)
    return out[:limit] if limit > 0 else out


def _flow_md_params(pdb_path: str, rec_chain: str, lig_chain: str,
                    sim_ns: str, skip_ps: str) -> dict:
    """flow 295 的 task-keyed 参数，只设必要项，其余留空走默认。"""
    params: dict[str, dict] = {
        TASK_PDB_PREPARE: {
            "PDB File": pdb_path,
            "Receptor Chain": rec_chain,
            "Ligand Chain": lig_chain,
        },
        TASK_TRAJECTORY: {
            "Type": "PDB",
            "Skip Time (ps)": skip_ps,
        },
    }
    if sim_ns:
        params[TASK_MD_MDP] = {"Simulation Time (ns)": sim_ns}
    return params


def _pick_one(result_dir: Path, glob: str, exclude: set[str], label: str) -> Path | None:
    """在结果目录里按 glob 找候选文件，排除已知输入名，返回唯一/最大候选并打印列表。"""
    cands = [p for p in result_dir.rglob(glob) if p.name not in exclude]
    if not cands:
        log.warning(f"[{label}] 在 {result_dir} 里按 {glob} 没找到候选")
        return None
    cands.sort(key=lambda p: p.stat().st_size, reverse=True)
    if len(cands) > 1:
        log.info(f"[{label}] 多个候选（取最大那个 {cands[0].name}）：{[p.name for p in cands[:6]]}")
    return cands[0]


def _state_index(records: list[dict]) -> dict[str, dict]:
    return {r["target"]: r for r in records}


def main() -> int:
    setup_console_logging(level=logging.INFO)
    args = parse_args()

    wt_pdb = Path(args.pdb)
    if not wt_pdb.exists():
        log.error(f"找不到 WT PDB: {wt_pdb}")
        return 2
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    mutfile_dir = out_dir / "sfe_md_mutfiles"
    mutfile_dir.mkdir(exist_ok=True)
    state_path = Path(args.state)
    result_base = Path(json.loads(Path(args.config).read_text(encoding="utf-8"))
                       .get("result_dir", "wemol_results"))

    muts = read_mut_list(Path(args.mut_list), args.limit)

    # 初始化/加载状态：WT 一个目标 + 每个突变一个目标
    records = read_jsonl(state_path)
    idx = _state_index(records)

    def ensure(target: str, kind: str, mut: str | None):
        if target not in idx:
            rec = {"target": target, "kind": kind, "mut": mut,
                   "stage": "init", "pdbmut_job": None, "md_job": None,
                   "mutant_pdb": None, "traj_pdb": None, "frames_dir": None,
                   "updated": now_iso()}
            records.append(rec)
            idx[target] = rec

    ensure("WT", "wt", None)
    for m in muts:
        ensure(m, "mut", m)

    session = build_session(args.config)
    changed = False

    def save():
        nonlocal changed
        if changed:
            rewrite_jsonl(state_path, records)

    for rec in records:
        target = rec["target"]
        kind = rec["kind"]
        stage = rec["stage"]

        # 已完成
        if stage == "extracted":
            continue

        # ── 突变体：先做 PDB Mutation 527 拿突变体结构 ────────────
        if kind == "mut" and rec["mutant_pdb"] is None:
            mut = rec["mut"]
            if rec["pdbmut_job"] is None:
                mutfile = mutfile_dir / f"{mut}.txt"
                mutfile.write_text(mut + "\n", encoding="utf-8", newline="\n")
                params = {"PDB File": wt_pdb.as_posix(), "Index Type": "POS",
                          "Mutation File": mutfile.as_posix()}
                jid = session.run_job(module_name=MOD_PDB_MUTATION_NAME,
                                      params=params, method_name="PDB Mutation",
                                      use_delay=False, result_dir_use_task_name=True)
                rec["pdbmut_job"] = jid
                rec["stage"] = "pdbmut_submitted"
                rec["updated"] = now_iso()
                changed = True
                log.info(f"[{target}] 提交 PDB Mutation -> job {jid}")
                continue
            # 已提交，查状态
            job = session.query_job_by_id(rec["pdbmut_job"])
            st = getattr(job, "status", None) if job else None
            if st == "Done":
                session.fetch_and_write_job_result(rec["pdbmut_job"], use_task_name=True)
                rdir = find_result_dir(result_base, rec["pdbmut_job"])
                pdb = _pick_one(rdir, args.mutant_pdb_glob, {wt_pdb.name}, f"{target}:mutantPDB") if rdir else None
                if pdb:
                    # 必须用 posix 正斜杠：服务器是 Linux，反斜杠路径会被搞乱
                    rec["mutant_pdb"] = pdb.as_posix()
                    rec["stage"] = "mutant_ready"
                    rec["updated"] = now_iso()
                    changed = True
                    log.info(f"[{target}] 突变体 PDB: {pdb}")
                else:
                    log.warning(f"[{target}] PDB Mutation Done 但没定位到突变体 PDB（调 --mutant-pdb-glob）")
                continue
            elif st in ("Abort", "Cancel"):
                rec["stage"] = f"pdbmut_{st}"
                changed = True
                log.error(f"[{target}] PDB Mutation {st}")
                continue
            else:
                log.info(f"[{target}] PDB Mutation 仍在跑 (status={st})")
                continue

        # ── 提交 MD flow 295 ──────────────────────────────────────
        md_input = wt_pdb.as_posix() if kind == "wt" else rec["mutant_pdb"]
        if rec["md_job"] is None:
            params = _flow_md_params(md_input, args.receptor_chain, args.ligand_chain,
                                     args.sim_ns, args.skip_ps)
            jid = session.run_job_use_flow(flow_name=FLOW_MD_NAME, params=params,
                                           use_delay=False, result_dir_use_task_name=True)
            rec["md_job"] = jid
            rec["stage"] = "md_submitted"
            rec["updated"] = now_iso()
            changed = True
            log.info(f"[{target}] 提交 MD flow 295 -> job {jid}")
            continue

        # ── MD 查状态 → Done 则下载+定位轨迹+拆帧 ─────────────────
        job = session.query_job_by_id(rec["md_job"])
        st = getattr(job, "status", None) if job else None
        if st == "Done":
            session.fetch_and_write_job_result(rec["md_job"], use_task_name=True)
            rdir = find_result_dir(result_base, rec["md_job"])
            traj = _pick_one(rdir, args.traj_glob, {wt_pdb.name}, f"{target}:traj") if rdir else None
            if not traj:
                log.warning(f"[{target}] MD Done 但没定位到多帧轨迹 PDB（调 --traj-glob）")
                continue
            rec["traj_pdb"] = str(traj)
            tag = "WT" if kind == "wt" else rec["mut"]
            frames_dir = (Path(args.wt_frames_dir) if kind == "wt"
                          else Path(args.frames_root) / rec["mut"])
            split_frames(Path(traj), tag, frames_dir, expect=args.n_frames)
            # MD 帧需还原成 Rosetta 兼容（链/编号/去氢/AMBER 命名归一），否则 FlexDDG 跑挂
            ref_order = build_ref_order(wt_pdb)
            n_fixed = 0
            for fp in sorted(frames_dir.glob("*.pdb")):
                rechain_pdb(fp, ref_order, strip_h=True)
                n_fixed += 1
            log.info(f"[{target}] 已 rechain {n_fixed} 帧（链/编号/去氢/命名归一）")
            rec["frames_dir"] = str(frames_dir)
            rec["stage"] = "extracted"
            rec["updated"] = now_iso()
            changed = True
            log.info(f"[{target}] 已抽帧 → {frames_dir}")
        elif st in ("Abort", "Cancel"):
            rec["stage"] = f"md_{st}"
            changed = True
            log.error(f"[{target}] MD flow {st}")
        else:
            log.info(f"[{target}] MD 仍在跑 (status={st})")

    save()

    # ── 汇总 ──────────────────────────────────────────────────
    done = sum(1 for r in records if r["stage"] == "extracted")
    print()
    print("=" * 50)
    print("MD 编排进度")
    print("=" * 50)
    for r in records:
        print(f"  {r['target']:<10} stage={r['stage']}")
    print(f"\n  已 extracted: {done}/{len(records)}")
    if done == len(records):
        print("\n✓ 全部目标已 extracted。下一步：")
        print(f"  python sfe_prepare.py --mut-list {args.mut_list} "
              f"--wt-frames-dir {args.wt_frames_dir} --wt-tag WT "
              f"--mut-frames-root {args.frames_root} --n-frames {args.n_frames} --out-dir {args.out_dir}")
    else:
        print("\n… 还有目标未完成（MD 慢）；等一会儿重跑本脚本继续推进。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
