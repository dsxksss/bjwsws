"""
SFE 流程 (3.5)：把 GROMACS MD 帧的链/编号还原成原始 PDB 的，供 Rosetta FlexDDG 使用。

问题：MD Trajectory v2 出来的帧把链 ID 抹掉了（变空白），残基连续重编号成 1..N，
导致 FlexDDG 的 mut 文件（指定 链M+resSeq31）在帧里找不到对应残基 → Rosetta 跑挂。

修法：MD 帧的残基**数量与顺序**和参考 WT PDB 一致（前处理没增删残基），
所以按「第 k 个残基」位置，把参考 PDB 的 (chain, resSeq, iCode) 逐残基写回帧里。

零依赖。就地覆盖帧文件（强制 LF）。

用法：
    python sfe_rechain.py --reference test_files/7l7e-ab-ag-complex.pdb \\
        --frames-dir prod/frames_WT prod/frames_mut/DM31L prod/frames_mut/DM31A
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

_COORD = ("ATOM  ", "HETATM", "TER", "ANISOU")


def build_ref_order(ref_pdb: Path) -> list[tuple[str, str, str]]:
    """参考 PDB 里每个残基的 (chain, resSeq4, iCode)，按出现顺序去重。"""
    order: list[tuple[str, str, str]] = []
    last_key = None
    for ln in ref_pdb.read_text(encoding="utf-8", errors="ignore").splitlines():
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        chain = ln[21]
        resseq = ln[22:26]
        icode = ln[26] if len(ln) > 26 else " "
        key = (chain, resseq, icode)
        if key != last_key:
            order.append(key)
            last_key = key
    return order


# AMBER/GROMACS 质子化态残基名 → Rosetta 标准名
_RES_MAP = {
    "ASH": "ASP", "GLH": "GLU", "CYX": "CYS", "CYM": "CYS", "LYN": "LYS",
    "HID": "HIS", "HIE": "HIS", "HIP": "HIS", "HSD": "HIS", "HSE": "HIS", "HSP": "HIS",
}


def _fmt_atom4(name: str) -> str:
    """把原子名格式化回 4 字符 PDB 字段（单字母元素，前导空格）。"""
    if len(name) >= 4:
        return name[:4]
    return (" " + name).ljust(4)


def _normalize_names(ln: str) -> str:
    """归一 GROMACS 特有的残基名/原子名为 Rosetta 标准。"""
    resn = ln[17:20].strip()
    atom = ln[12:16].strip()
    new_atom = atom
    if atom in ("OC1", "OT1"):
        new_atom = "O"
    elif atom in ("OC2", "OT2"):
        new_atom = "OXT"
    elif resn == "ILE" and atom == "CD":
        new_atom = "CD1"
    new_resn = _RES_MAP.get(resn, resn)
    if new_atom != atom:
        ln = ln[:12] + _fmt_atom4(new_atom) + ln[16:]
    if new_resn != resn:
        ln = ln[:17] + f"{new_resn:>3}" + ln[20:]
    return ln


def _is_hydrogen(ln: str) -> bool:
    """PDB 行是否是氢原子：优先看元素列(77-78)，缺失则用原子名首字母启发式。"""
    elem = ln[76:78].strip() if len(ln) >= 78 else ""
    if elem:
        return elem == "H"
    name = ln[12:16].strip()
    # 原子名如 "1HB"/"HD21"/"H"：去掉前导数字后首字母为 H
    stripped = name.lstrip("0123456789")
    return stripped[:1] == "H"


def rechain_pdb(frame: Path, ref_order: list[tuple[str, str, str]],
                strip_h: bool = True) -> tuple[int, int]:
    """把 frame 的每个残基按位置改写成 ref_order 的链/编号（可选去氢）。返回 (改写残基数, 残基总数)。"""
    out: list[str] = []
    idx = -1
    last_frame_resid = None
    n_residues = 0
    for ln in frame.read_text(encoding="utf-8", errors="ignore").splitlines():
        if strip_h and ln.startswith(("ATOM  ", "HETATM", "ANISOU")) and _is_hydrogen(ln):
            continue  # 丢弃氢原子（Rosetta 要重原子、自己加氢）
        if ln.startswith(("ATOM  ", "HETATM", "ANISOU")):
            frame_resid = (ln[21], ln[22:27])  # 帧里的残基键（链+resSeq+icode）
            if frame_resid != last_frame_resid:
                last_frame_resid = frame_resid
                idx += 1
                n_residues += 1
            ln = _normalize_names(ln)
            if 0 <= idx < len(ref_order):
                chain, resseq, icode = ref_order[idx]
                ln = ln[:21] + chain + resseq + icode + ln[27:]
            out.append(ln)
        elif ln.startswith("TER"):
            # TER 行也带残基号，跟着上一个残基走（简单处理：用当前 idx）
            if 0 <= idx < len(ref_order):
                chain, resseq, icode = ref_order[idx]
                # TER 至少补到 27 列
                base = ln.ljust(27)
                ln = base[:21] + chain + resseq + icode + base[27:]
            out.append(ln)
        else:
            out.append(ln)
    frame.write_text("\n".join(out) + "\n", encoding="utf-8", newline="\n")
    return min(idx + 1, len(ref_order)), n_residues


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="把 MD 帧的链/编号还原成参考 PDB 的（供 FlexDDG）",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--reference", required=True, help="参考 WT PDB（提供链/编号顺序）")
    p.add_argument("--frames-dir", nargs="+", required=True,
                   help="一个或多个帧目录，处理其中所有 *.pdb（就地覆盖）")
    p.add_argument("--keep-h", action="store_true",
                   help="保留氢原子（默认去氢，Rosetta 要重原子）")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    ref = Path(args.reference)
    if not ref.exists():
        print(f"❌ 找不到参考 PDB: {ref}", file=sys.stderr)
        return 2
    ref_order = build_ref_order(ref)
    print(f"参考 {ref.name}: {len(ref_order)} 个残基")

    total = 0
    for d in args.frames_dir:
        dpath = Path(d)
        if not dpath.exists():
            print(f"⚠ 跳过不存在的目录: {dpath}")
            continue
        for pdb in sorted(dpath.glob("*.pdb")):
            n_re, n_res = rechain_pdb(pdb, ref_order, strip_h=not args.keep_h)
            total += 1
            flag = "" if n_res == len(ref_order) else f"  ⚠ 残基数 {n_res}≠参考 {len(ref_order)}"
            print(f"  ✓ {pdb}  ({n_re} 残基重链){flag}")
    print(f"共处理 {total} 个帧文件")
    return 0


if __name__ == "__main__":
    sys.exit(main())
