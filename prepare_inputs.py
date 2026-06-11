"""
生成各模块需要的 mut 文件 + tasks.csv。突变来源支持两种模式：

【模式 A：自动生成】给定 PDB + 突变位点，直接枚举突变（零依赖，纯标准库）
    python prepare_inputs.py \\
        --pdb 7l7e-ab-ag-complex.pdb \\
        --muts "M31,M53,M55-56,M104-113,M116,N31-36,N38,N55-56,N59" \\
        --omit-aas CP \\
        --modules flexddg,foldx,fep \\
        --out-dir prod

【模式 B：使用现成文件】用已经生成好的突变列表（无任何额外依赖）
    python prepare_inputs.py \\
        --mut-list single_mut_list.txt \\
        --mut-list-foldx mut_list_for_foldx.txt \\
        --pdb 7l7e-ab-ag-complex.pdb \\
        --modules flexddg,foldx,fep \\
        --out-dir prod \\
        --limit 5

突变列表格式：每行一个突变 `DM31L`（原AA+链+UID+新AA）。

本脚本输出（写到 --out-dir）：
  - single_mut_list.txt / mut_list_for_foldx.txt （模式 A 自动生成）
  - Mut_Hchains.fasta / Mut_Lchains.fasta        （模式 A 自动输出，AbNatiV 用）
  - mutlists_flexddg/<mut>.txt   每个突变一个文件（FlexDDG 单元型）
  - mutlist_foldx.txt            （FoldX 用，每行带 ;）
  - tasks.csv                    submit_batch.py 直接吃

AbNatiV 的轻重链 FASTA：模式 A 会自动判定每条被突变链是重链还是轻链
  （基于抗体 J 区保守基序，无需 anarci），并输出对应 FASTA。
  如自动判定有误，可用 --heavy-chains / --light-chains 手动覆盖。
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

# ─────────────────────────────────────────────────────────────────────────────
# 零依赖突变生成器（逻辑对齐实验室 mutation_list.py，仅省略 H/L FASTA 输出）
# ─────────────────────────────────────────────────────────────────────────────

# 与实验室脚本保持相同顺序与集合
_NATURAL_AA = ["L","A","G","V","S","E","R","T","I","D","P","K","Q","N","F","Y","M","H","W","C"]

# PDB 三字母 → 单字母（标准 20 种）
_THREE_TO_ONE = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C",
    "GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I",
    "LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P",
    "SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V",
}


def _read_pdb_residues(pdb_file: Path) -> dict[str, dict[int, str]]:
    """
    纯标准库解析 PDB ATOM 行，返回 {chain: {resSeq(int): one_letter_aa}}。
    对应实验室脚本 get_res_from_pdb() 的残基字典部分（不含 chain_uid_pos）。
    每个 (chain, resSeq) 取首次出现的残基名（与 biopython 行为一致）。
    """
    result: dict[str, dict[int, str]] = {}
    for line in pdb_file.read_text(encoding="utf-8", errors="ignore").splitlines():
        if not line.startswith("ATOM"):
            continue
        resname = line[17:20].strip()
        one = _THREE_TO_ONE.get(resname)
        if one is None:
            continue                          # 非标准残基 / 水 / 配体跳过
        chain = line[21:22].strip()
        try:
            resseq = int(line[22:26].strip())
        except ValueError:
            continue
        chain_map = result.setdefault(chain, {})
        if resseq not in chain_map:           # 首次出现为准（与 biopython 一致）
            chain_map[resseq] = one
    return result


def _parse_muts_spec(spec: str) -> dict[str, list[int]]:
    """
    'M31,M53,M55-56,M104-113' → {'M': [31,53,55,56,104,...,113]}
    对应实验室脚本 single_mut() 里的 chain_mut_res_dict 构建逻辑。
    """
    chain_positions: dict[str, list[int]] = {}
    for raw in spec.strip().split(","):
        t = raw.strip()
        if not t:
            continue
        chain = t[0]
        rest  = t[1:]
        if "-" in rest:
            lo, hi = rest.split("-", 1)
            positions = list(range(int(lo), int(hi) + 1))
        else:
            positions = [int(rest)]
        chain_positions.setdefault(chain, []).extend(positions)
    return chain_positions


# 抗体 J 区保守基序：重链 W-G-x-G-[TS]，轻链 F-G-x-G-[TS]
_JREGION_RE = re.compile(r"[WF]G[A-Z]G[TS]")


def _guess_chain_type(seq: str) -> str | None:
    """
    基于抗体 J 区保守基序自动判定链类型。
    重链 J 区以 W 开头（WGxGT），轻链以 F 开头（FGxGT）。
    返回 'H' / 'L'；非抗体链（找不到该基序）返回 None。
    """
    m = _JREGION_RE.search(seq)
    if not m:
        return None
    return "H" if m.group(0)[0] == "W" else "L"


def _chain_seq_and_index(chain_map: dict[int, str]) -> tuple[str, dict[int, str]]:
    """
    由 {resSeq: aa}（文件顺序）得到 (整链序列, {resSeq: 0基位置索引})。
    对应实验室脚本里的 wt_seq 与 chain_uid_pos。
    """
    seq = "".join(chain_map.values())
    idx = {resseq: i for i, resseq in enumerate(chain_map.keys())}
    return seq, idx


def generate_mut_lists(
    pdb_file: Path,
    muts_spec: str,
    omit_aas: str,
    out_dir: Path,
    chain_types: dict[str, str] | None = None,
    limit: int = 0,
) -> tuple[Path, Path]:
    """
    等价于实验室 mutation_list.py 的核心逻辑（single_mut 函数），
    生成 single_mut_list.txt 和 mut_list_for_foldx.txt。

    自动判定每条被突变链的轻/重类型（基于 J 区保守基序），额外输出突变后整链
    序列的 Mut_Hchains.fasta / Mut_Lchains.fasta（AbNatiV 用）。
    chain_types 可显式覆盖自动判定，如 {'M':'H','N':'L'}。
    返回 (single_mut_list, mut_list_for_foldx) 的 Path。
    """
    residues       = _read_pdb_residues(pdb_file)
    chain_pos      = _parse_muts_spec(muts_spec)
    omit           = set(omit_aas.upper())
    target_aas     = [aa for aa in _NATURAL_AA if aa not in omit]
    overrides      = chain_types or {}

    print(f"  允许的目标氨基酸: {target_aas}")

    # 为所有被突变链预算整链序列与位置索引，并判定链类型
    seq_cache: dict[str, tuple[str, dict[int, str]]] = {}
    resolved_type: dict[str, str | None] = {}
    for chain in chain_pos:
        if chain not in residues:
            continue
        seq_cache[chain] = _chain_seq_and_index(residues[chain])
        if chain in overrides:
            resolved_type[chain] = overrides[chain]
            print(f"  链 {chain}: 指定为 {'重链(H)' if overrides[chain]=='H' else '轻链(L)'}")
        else:
            guess = _guess_chain_type(seq_cache[chain][0])
            resolved_type[chain] = guess
            label = {"H": "重链(H)", "L": "轻链(L)"}.get(guess, "非抗体链/无法判定(跳过FASTA)")
            print(f"  链 {chain}: 自动判定为 {label}")

    mutations: list[str] = []
    fasta_records: dict[str, list[tuple[str, str]]] = {"H": [], "L": []}

    for chain, positions in chain_pos.items():
        if chain not in residues:
            raise ValueError(
                f"PDB 中找不到链 {chain!r}，可用链: {sorted(residues.keys())}"
            )
        chain_map = residues[chain]
        ctype = resolved_type.get(chain)
        for uid in positions:
            wt = chain_map.get(uid)
            if wt is None:
                raise ValueError(
                    f"链 {chain} 位置 {uid} 在 PDB 中找不到标准氨基酸残基"
                )
            for aa in target_aas:
                if aa == wt:
                    continue
                mutation = f"{wt}{chain}{uid}{aa}"
                mutations.append(mutation)
                # 生成突变后整链序列写入对应 FASTA
                if ctype in ("H", "L"):
                    wt_seq, idx = seq_cache[chain]
                    pos = idx[uid]
                    mut_seq = wt_seq[:pos] + aa + wt_seq[pos + 1:]
                    fasta_records[ctype].append((mutation, mut_seq))

    # --limit：一致地截断突变列表与 FASTA（取前 N 个突变）
    if limit and limit > 0:
        kept = set(mutations[:limit])
        mutations = mutations[:limit]
        for ct in fasta_records:
            fasta_records[ct] = [(n, s) for n, s in fasta_records[ct] if n in kept]
        print(f"  ⓿ --limit {limit}：截断为前 {len(mutations)} 个突变（含 FASTA）")

    single_path = out_dir / "single_mut_list.txt"
    foldx_path  = out_dir / "mut_list_for_foldx.txt"
    single_path.write_text("\n".join(mutations) + "\n", encoding="utf-8", newline="\n")
    foldx_path.write_text("\n".join(f"{m};" for m in mutations) + "\n", encoding="utf-8", newline="\n")

    print(f"  生成突变数: {len(mutations)}")
    print(f"  → {single_path}")
    print(f"  → {foldx_path}")

    # 写 FASTA（标准格式：header 与序列分行）
    for ctype, fname in (("H", "Mut_Hchains.fasta"), ("L", "Mut_Lchains.fasta")):
        recs = fasta_records[ctype]
        if not recs:
            continue
        fa_path = out_dir / fname
        fa_path.write_text(
            "".join(f">{name}\n{seq}\n" for name, seq in recs),
            encoding="utf-8", newline="\n",
        )
        print(f"  → {fa_path} ({len(recs)} 条序列)")

    return single_path, foldx_path


# ─────────────────────────────────────────────────────────────────────────────
# 公共工具
# ─────────────────────────────────────────────────────────────────────────────

SINGLE_MUT_RE = re.compile(r"^([A-Z])([A-Za-z0-9])(\d+)([A-Z])$")


def parse_single_mut(name: str) -> tuple[str, str, int, str] | None:
    m = SINGLE_MUT_RE.match(name.strip())
    if not m:
        return None
    return m.group(1), m.group(2), int(m.group(3)), m.group(4)


def read_mut_list(path: Path) -> list[str]:
    if not path.exists():
        raise FileNotFoundError(f"找不到突变列表: {path}")
    out = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        out.append(line)
    return out


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="生成各模块的 mut 文件和 tasks.csv",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    # 模式 A
    p.add_argument("--muts", default=None,
                   help="【模式A】突变位点，如 'M31,M53,M55-56'（给了即自动生成突变列表）")
    p.add_argument("--omit-aas", default="CP",
                   help="【模式A】排除的目标氨基酸（默认 CP）")
    p.add_argument("--heavy-chains", default=None,
                   help="【模式A·可选】手动指定重链链名(逗号分隔)，覆盖自动判定")
    p.add_argument("--light-chains", default=None,
                   help="【模式A·可选】手动指定轻链链名(逗号分隔)，覆盖自动判定")
    # 模式 B
    p.add_argument("--mut-list", default=None,
                   help="【模式B】现成 single_mut_list.txt 路径")
    p.add_argument("--mut-list-foldx", default=None,
                   help="【模式B】现成 mut_list_for_foldx.txt（不给则从 single 自动合成）")
    # 公共
    p.add_argument("--pdb", required=True,
                   help="PDB 文件路径")
    p.add_argument("--modules", default="flexddg,foldx,fep",
                   help="逗号分隔的模块缩写：flexddg / foldx / fep / abnativ（默认前三个）")
    p.add_argument("--out-dir", default="prep",
                   help="输出目录（默认 ./prep）")
    p.add_argument("--limit", type=int, default=0,
                   help="只取前 N 个突变（0=全部；测试时用 --limit 5）")
    p.add_argument("--fep-type", default="B", choices=["S", "B"],
                   help="FEP Type：B=Binding affinity / S=Stability（默认 B）")
    p.add_argument("--heavy-nat", default="VH", choices=["VH", "VL"],
                   help="AbnatiV 重链 Nat 参数（默认 VH）")
    p.add_argument("--light-nat", default="VL", choices=["VH", "VL"],
                   help="AbnatiV 轻链 Nat 参数（默认 VL）")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    # ── 模块校验 ──────────────────────────────────────────────
    modules = {m.strip().lower() for m in args.modules.split(",") if m.strip()}
    unknown = modules - {"flexddg", "foldx", "fep", "abnativ"}
    if unknown:
        print(f"❌ 未识别的模块: {unknown}; 支持: flexddg / foldx / fep / abnativ", file=sys.stderr)
        return 2

    # ── 模式校验 ──────────────────────────────────────────────
    if not args.muts and not args.mut_list:
        print("❌ 必须给 --muts（模式A）或 --mut-list（模式B）之一", file=sys.stderr)
        return 2
    if args.muts and args.mut_list:
        print("⚠ 同时给了 --muts 和 --mut-list，以 --muts（模式A）为准")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # ── 取得突变列表 ──────────────────────────────────────────
    foldx_src: Path | None = None

    if args.muts:
        print("⚙ 模式A：从 PDB + 位点描述自动生成突变列表（零依赖）")
        pdb = Path(args.pdb)
        if not pdb.exists():
            print(f"❌ 找不到 PDB 文件: {pdb}", file=sys.stderr)
            return 2
        # 组装链类型映射（重/轻链 → H/L），用于输出 FASTA
        chain_types: dict[str, str] = {}
        if args.heavy_chains:
            for c in args.heavy_chains.split(","):
                if c.strip():
                    chain_types[c.strip()] = "H"
        if args.light_chains:
            for c in args.light_chains.split(","):
                if c.strip():
                    chain_types[c.strip()] = "L"

        # 模式A：生成时即按 --limit 一致截断（含 FASTA）
        single_path, foldx_path = generate_mut_lists(
            pdb, args.muts, args.omit_aas, out_dir, chain_types or None,
            limit=args.limit,
        )
        mut_list_path = single_path
        foldx_src = foldx_path
    else:
        print("⚙ 模式B：使用现成突变列表文件")
        mut_list_path = Path(args.mut_list)
        if args.mut_list_foldx:
            foldx_src = Path(args.mut_list_foldx)

    # ── 读突变列表 ────────────────────────────────────────────
    muts = read_mut_list(mut_list_path)
    if args.limit > 0:
        muts = muts[: args.limit]
        print(f"⓿ --limit {args.limit}：仅处理前 {len(muts)} 个突变")
    if not muts:
        print("❌ 突变列表为空", file=sys.stderr)
        return 2

    # 拆分单点 / 多点
    singles: list[tuple[str, tuple[str, str, int, str]]] = []
    multi:   list[str] = []
    for m in muts:
        parsed = parse_single_mut(m)
        if parsed:
            singles.append((m, parsed))
        else:
            multi.append(m)
    if multi:
        print(f"⚠ 检测到 {len(multi)} 个多点突变（FlexDDG / FEP 跳过，FoldX 保留）")

    pdb_path = args.pdb

    # ── FlexDDG：每个单点一个 mut 文件 ───────────────────────
    flex_files: dict[str, Path] = {}
    if "flexddg" in modules:
        flex_dir = out_dir / "mutlists_flexddg"
        flex_dir.mkdir(exist_ok=True)
        for mut_name, _ in singles:
            f = flex_dir / f"{mut_name}.txt"
            f.write_text(f"{mut_name}\n", encoding="utf-8", newline="\n")
            flex_files[mut_name] = f
        print(f"✓ FlexDDG: 生成 {len(flex_files)} 个 mut 文件 → {flex_dir}")

    # ── FoldX：所有突变一个文件 ───────────────────────────────
    foldx_file: Path | None = None
    if "foldx" in modules:
        if foldx_src and foldx_src.exists():
            foldx_file = out_dir / "mutlist_foldx.txt"
            if args.limit > 0:
                lines = [ln for ln in foldx_src.read_text(encoding="utf-8")
                         .splitlines() if ln.strip()][: args.limit]
                foldx_file.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
            else:
                # 读时 universal newline 归一为 \n，写时强制 LF（去掉源文件可能的 CRLF）
                normalized = "\n".join(foldx_src.read_text(encoding="utf-8").splitlines()) + "\n"
                foldx_file.write_text(normalized, encoding="utf-8", newline="\n")
            print(f"✓ FoldX: mut 文件 → {foldx_file}")
        else:
            # 没有 foldx 专用文件，从 single 合成（每行加 ;）
            foldx_file = out_dir / "mutlist_foldx.txt"
            foldx_file.write_text(
                "\n".join(f"{m};" for m in muts) + "\n", encoding="utf-8", newline="\n"
            )
            print(f"✓ FoldX: 从 single_mut_list 合成 → {foldx_file}")

    # ── tasks.csv ─────────────────────────────────────────────
    tasks_csv = out_dir / "tasks.csv"
    rows: list[dict[str, str]] = []

    if "flexddg" in modules:
        for mut_name, _ in singles:
            rows.append({
                "task_name":   f"flexddg_{mut_name}",
                "mutation":    mut_name,
                "module_name": "Flex DDG",
                "module_id":   "",
                "method_name": "",
                "flow_name":   "",
                "flow_id":     "",
                "params":      f"PDB File={pdb_path};Mut File={flex_files[mut_name].as_posix()}",
            })

    if "foldx" in modules and foldx_file:
        rows.append({
            "task_name":   "foldx_all",
            "mutation":    "*",
            "module_name": "Mutation Energy of Binding",
            "module_id":   "",
            "method_name": "Mutation Energy of Binding",
            "flow_name":   "",
            "flow_id":     "",
            "params":      f"PDB File={pdb_path};Mutant File={foldx_file.as_posix()};Reorder Mutation Number=no",
        })

    if "fep" in modules:
        for mut_name, (orig, chain, pos, new) in singles:
            rows.append({
                "task_name":   f"fep_{mut_name}",
                "mutation":    mut_name,
                "module_name": "Protein FEP",
                "module_id":   "",
                "method_name": "Single-point Mutation",
                "flow_name":   "",
                "flow_id":     "",
                "params":      f"PDB File={pdb_path};Mutation={orig}{pos}{new};Type={args.fep_type};Chain={chain}",
            })

    # AbnatiV：批量型，一份 FASTA 一个作业（重链 / 轻链各一个）
    if "abnativ" in modules:
        for fname, nat, tag in (("Mut_Hchains.fasta", args.heavy_nat, "H"),
                                ("Mut_Lchains.fasta", args.light_nat, "L")):
            fa = out_dir / fname
            if not fa.exists():
                print(f"⚠ AbnatiV: 找不到 {fa}（需模式A自动生成 FASTA），跳过 {tag} 链")
                continue
            rows.append({
                "task_name":   f"abnativ_{tag}",
                "mutation":    "*",
                "module_name": "AbnatiV",
                "module_id":   "",
                "method_name": "",
                "flow_name":   "",
                "flow_id":     "",
                "params":      f"Fasta File={fa.as_posix()};Nat={nat}",
            })

    fieldnames = ["task_name", "mutation", "module_name", "module_id",
                  "method_name", "flow_name", "flow_id", "params"]
    with tasks_csv.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)

    # ── 汇总 ──────────────────────────────────────────────────
    print()
    print("=" * 50)
    print(f"任务清单已生成: {tasks_csv}")
    print("=" * 50)
    n_flex  = sum(1 for r in rows if r["module_name"] == "Flex DDG")
    n_foldx = sum(1 for r in rows if r["module_name"] == "Mutation Energy of Binding")
    n_fep   = sum(1 for r in rows if r["module_name"] == "Protein FEP")
    n_abn   = sum(1 for r in rows if r["module_name"] == "AbnatiV")
    print(f"  FlexDDG : {n_flex} 个作业  (单元型，1 突变 1 作业)")
    print(f"  FoldX   : {n_foldx} 个作业  (批量型，1 作业管全部)")
    print(f"  FEP     : {n_fep} 个作业  (单元型，多点已跳过)")
    print(f"  AbnatiV : {n_abn} 个作业  (批量型，重/轻链各一个)")
    print(f"  合计    : {len(rows)} 个作业")
    print()
    print("下一步：")
    print(f"  python submit_batch.py --config config.json --tasks {tasks_csv} -o submitted_jobs.jsonl")
    return 0


if __name__ == "__main__":
    sys.exit(main())
