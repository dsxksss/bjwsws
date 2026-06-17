"""
SFE 流程 (3)：把 MD Trajectory v2 输出的多帧 PDB 拆成单帧 PDB。

MD Trajectory v2 (module 400, Type=PDB) 会把轨迹导成一个含多个 MODEL/ENDMDL 的
PDB 文件。本脚本按 MODEL 把它拆成 N 个单帧 PDB，命名 `<tag>_NN.pdb`（NN 从 01 起，
两位补零），对应文档「需要准备脚本来 model 进行拆分」。

文档协议：每条轨迹取 21 个构象（1 个初始构象 + 20 个 MD 帧，Skip Time=500ps）。
- 若多帧 PDB 已含 21 个 MODEL，直接拆成 21 个。
- 若 MD Trajectory 只给了 20 个 MD 帧，可用 --initial 把初始结构作为第 01 帧插到最前，
  MD 帧顺延为 02..21。

零依赖（纯标准库），输出强制 LF 换行（与项目其余部分一致，避免 Rosetta/FoldX 的 CRLF 坑）。

用法：
    # 直接按 MODEL 拆
    python sfe_split_frames.py --traj WT_traj.pdb --tag WT --out-dir prod/frames_WT

    # 初始结构单独给，作为第 01 帧
    python sfe_split_frames.py --traj DM31A_md.pdb --tag DM31A \\
        --initial DM31A_start.pdb --out-dir prod/frames_DM31A --expect 21
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

# 一帧里需要保留的记录类型（坐标 + 链接/终止），其余（MODEL/ENDMDL/CONECT 等）按需处理
_COORD_PREFIXES = ("ATOM  ", "HETATM", "TER", "ANISOU")
_PREAMBLE_PREFIXES = ("CRYST1", "SCALE", "ORIGX", "REMARK", "SEQRES", "HEADER", "TITLE")


def _split_models(text: str) -> tuple[list[str], list[list[str]]]:
    """
    返回 (preamble_lines, [model1_lines, model2_lines, ...])。
    - preamble：第一个 MODEL 之前的非坐标头部（CRYST1 等），拆分后会拼到每帧前面。
    - 每个 model_lines 是该 MODEL/ENDMDL 之间的原始行（不含 MODEL/ENDMDL 本身）。
    若文件无 MODEL 记录，则视作单帧：整份坐标作为一个 model。
    """
    lines = text.splitlines()
    preamble: list[str] = []
    models: list[list[str]] = []
    cur: list[str] | None = None
    seen_model = False

    for ln in lines:
        if ln.startswith("MODEL "):
            seen_model = True
            cur = []
            continue
        if ln.startswith("ENDMDL"):
            if cur is not None:
                models.append(cur)
                cur = None
            continue
        if cur is not None:
            cur.append(ln)
        elif not seen_model:
            # 第一个 MODEL 之前：收集头部（仅保留有意义的头，丢弃空行）
            if ln.startswith(_PREAMBLE_PREFIXES):
                preamble.append(ln)

    # 无 MODEL：整份坐标当作单帧
    if not seen_model:
        coord = [ln for ln in lines if ln.startswith(_COORD_PREFIXES)]
        head = [ln for ln in lines if ln.startswith(_PREAMBLE_PREFIXES)]
        return head, [coord] if coord else []

    return preamble, models


def _read_single_frame(pdb_path: Path) -> list[str]:
    """读一个单结构 PDB 的坐标记录（含 TER），用作初始帧。"""
    out: list[str] = []
    for ln in pdb_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if ln.startswith("MODEL ") or ln.startswith("ENDMDL"):
            continue
        if ln.startswith(_COORD_PREFIXES):
            out.append(ln)
    return out


def _write_frame(out_path: Path, preamble: list[str], coords: list[str]) -> None:
    body: list[str] = []
    body.extend(preamble)
    body.extend(ln for ln in coords if ln.startswith(_COORD_PREFIXES))
    body.append("END")
    out_path.write_text("\n".join(body) + "\n", encoding="utf-8", newline="\n")


def split_frames(
    traj: Path,
    tag: str,
    out_dir: Path,
    initial: Path | None = None,
    expect: int = 0,
) -> list[Path]:
    text = traj.read_text(encoding="utf-8", errors="ignore")
    preamble, models = _split_models(text)
    if not models:
        raise ValueError(f"{traj} 里没解析到任何坐标帧（MODEL/ATOM）")

    frames: list[list[str]] = []
    if initial is not None:
        init_coords = _read_single_frame(initial)
        if not init_coords:
            raise ValueError(f"--initial {initial} 没解析到坐标")
        frames.append(init_coords)
    frames.extend(models)

    out_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for i, coords in enumerate(frames, start=1):
        out_path = out_dir / f"{tag}_{i:02d}.pdb"
        _write_frame(out_path, preamble, coords)
        written.append(out_path)

    print(f"✓ 拆出 {len(written)} 帧 → {out_dir}/  ({tag}_01.pdb .. {tag}_{len(written):02d}.pdb)")
    if initial is not None:
        print(f"  (第 01 帧来自初始结构 {initial.name}，其余为 MD 帧)")
    if expect and len(written) != expect:
        print(f"  ⚠ 期望 {expect} 帧，实际 {len(written)} 帧——请核对 MD Trajectory 的 "
              f"Start/End/Skip Time 设置", file=sys.stderr)
    return written


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="把 MD Trajectory v2 的多帧 PDB 拆成单帧 <tag>_NN.pdb",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--traj", required=True, help="多帧 PDB（MD Trajectory v2 Type=PDB 输出）")
    p.add_argument("--tag", required=True, help="帧文件名前缀，如 WT 或 DM31A")
    p.add_argument("--out-dir", required=True, help="输出目录")
    p.add_argument("--initial", default=None,
                   help="（可选）初始构象 PDB，作为第 01 帧插到最前")
    p.add_argument("--expect", type=int, default=0,
                   help="（可选）期望帧数，用于校验（如 21），不符只告警不报错")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    traj = Path(args.traj)
    if not traj.exists():
        print(f"❌ 找不到轨迹 PDB: {traj}", file=sys.stderr)
        return 2
    initial = Path(args.initial) if args.initial else None
    if initial and not initial.exists():
        print(f"❌ 找不到初始结构: {initial}", file=sys.stderr)
        return 2
    split_frames(traj, args.tag, Path(args.out_dir), initial=initial, expect=args.expect)
    return 0


if __name__ == "__main__":
    sys.exit(main())
