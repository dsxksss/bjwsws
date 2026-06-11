#!/data/PRG/tools/miniconda3/bin/python
# Copyright (c) 2018-2023, Lawrence Livermore National Security, LLC
# SPDX-License-Identifier: MIT

"""
一键运行本地抗体设计流水线。

必需输入:
  (1) 2.2 单点突变五工具打分 CSV
  (2) PDB 结构 + 抗体链名（抗原链可省略，自动推断）

GP 仅训练一次；Step 03（生成）↔ Step 06（MEI 选择）可循环 bo_rounds 次（不重训 GP）。

用法:
  python run_pipeline.py \\
    --single-point-scores data/single_point_scores.csv \\
    --pdb structure/7l7e.pdb \\
    --antibody-chains M,N \\
    --antigen-chains S \\
    --output-dir work/run1 \\
    --bo-rounds 3 \\
    --num-sequences-per-round 2000 \\
    --mei-batch-size 50 \\
    --final-quota 20
"""

from __future__ import division, print_function

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

_REPO = Path(__file__).resolve().parents[1]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

import abag_ml.pareto_selection as ps  # noqa: E402

from local_pipeline.common.scoring import (  # noqa: E402
    combine_tool_scores, logistic_transform, normalize_sampling_weights, SCORE_COLUMNS,
    validate_single_point_scores_df,
)
from local_pipeline.common.menu_utils import (  # noqa: E402
    build_allowed_mutations_from_scores, save_allowed_mutations,
)
from local_pipeline.common.pdb_structure import (  # noqa: E402
    extract_master_from_pdb, extract_interface_pairs, resolve_antigen_chains,
    write_interface_pairs_json,
)
from local_pipeline.common.pipeline_io import (  # noqa: E402
    write_fasta, write_recommendations_fasta,
)
from local_pipeline.common.mutation_generator import generate_mutant_sequences  # noqa: E402
from local_pipeline.common.featurization import (  # noqa: E402
    load_interface_pairs, featurize_sequences_df,
)
from local_pipeline.common.gp_training import train_gp_model  # noqa: E402
from local_pipeline.common.mei_selection import mei_select_sequences  # noqa: E402
from local_pipeline.common.objective_aggregate import aggregate_objectives  # noqa: E402
from local_pipeline.common.sequence_utils import mutate_seq  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(description='Run full local antibody design pipeline')
    p.add_argument('--single-point-scores', required=True, help='2.2 单点打分 CSV')
    p.add_argument('--pdb', required=True, help='初始 PDB 结构')
    p.add_argument('--antibody-chains', required=True, help='重/轻链，逗号分隔，如 M,N')
    p.add_argument('--antigen-chains', default=None,
                   help='抗原链，逗号分隔，如 S；省略则取 PDB 中除抗体链外的所有链')
    p.add_argument('--output-dir', default=None, help='输出目录（默认 work/run_时间戳）')

    p.add_argument('--min-mutations', type=int, default=1, help='Step A 最少突变数')
    p.add_argument('--max-mutations', type=int, default=8, help='Step A 最多突变数')
    p.add_argument('--interface-cutoff', type=float, default=10.0)

    p.add_argument('--bo-rounds', type=int, default=3,
                   help='Step03↔Step06 循环次数（GP 不重训）')
    p.add_argument('--num-sequences-per-round', type=int, default=2000,
                   help='每轮生成候选序列数')
    p.add_argument('--mei-batch-size', type=int, default=50,
                   help='每轮 MEI 保留条数')

    p.add_argument('--gp-num-iters', type=int, default=500)
    p.add_argument('--final-quota', type=int, default=20,
                   help='Pareto 加权降维后最终推荐条数')
    p.add_argument('--seed', type=int, default=0)
    return p.parse_args()


def compute_sampling_weights(scores_df):
    df = scores_df.copy()
    df['location'] = df['location'].astype(int)
    df['mutationHumanReadable'] = df.apply(
        lambda r: '{}{}{}'.format(r['original_aa'], r['location'], r['mutant_aa']),
        axis=1,
    )
    for col in SCORE_COLUMNS:
        df['l_' + col] = df[col].apply(
            lambda v: logistic_transform(v) if pd.notna(v) else 0.0)
    df['samplingWeight'] = df.apply(combine_tool_scores, axis=1)
    return normalize_sampling_weights(df)


def build_single_point_seq_df(scores_df, master_seq):
    rows = []
    for _, r in scores_df.iterrows():
        mut = ('', str(int(r['location'])), r['original_aa'], r['mutant_aa'])
        seq = mutate_seq(master_seq, [mut])
        hr = '{}{}{}'.format(r['original_aa'], int(r['location']), r['mutant_aa'])
        rows.append({
            'sequence': seq,
            'mutationHumanReadable': hr,
            'num_mutations': 1,
        })
    return pd.DataFrame(rows)


def run_pareto_downselect(objectives_df, quota):
    # 目标列 → 打分方向
    obj_scorers = [
        ('sum_rosetta_flex', ps.simple_row_scorer),
        ('sum_foldx',        ps.simple_row_scorer),
        ('sum_sfe',          ps.simple_row_scorer),
        ('sum_abnativ',      ps.simple_negative_row_scorer),
        ('num_mutations',    ps.simple_row_scorer),
    ]
    for c, _ in obj_scorers:
        if c not in objectives_df.columns:
            objectives_df[c] = float('nan')

    # --- PATCH: 跳过整列全为 NaN 的目标（如 sfe 留空时的 sum_sfe），
    # 否则 NaN 会污染 Pareto 支配比较，导致选不出推荐。
    active = [(c, fn) for c, fn in obj_scorers if objectives_df[c].notna().any()]
    dom_funcs = [
        (lambda r=None, c=c, fn=fn: fn(r, column=c)) for c, fn in active
    ]
    pareto_idx = ps.get_pareto_rows(
        objectives_df, dominance_functions=dom_funcs,
        scalar_epsilons=[0.0] * len(dom_funcs), returnints=True,
    )
    df = objectives_df.copy()
    df['ParetoSet'] = [i in pareto_idx for i in range(len(df))]

    sel_cols = ['sum_rosetta_flex', 'sum_foldx', 'sum_abnativ', 'num_mutations']
    weights = [-1.0, -0.5, 2.0, -0.5]
    df['downselection_score'] = 0.0
    for col, w in zip(sel_cols, weights):
        # --- PATCH: NaN 列以 0 计入，避免整体打分变 NaN
        df['downselection_score'] += w * df[col].fillna(0.0)
    tmp = df['downselection_score'].copy()
    tmp[~df['ParetoSet']] = float('nan')
    best = []
    for _ in range(quota):
        if not tmp.notna().any():
            break
        i = int(np.nanargmax(tmp.values))
        best.append(i)
        tmp.iloc[i] = float('nan')
    df['Recommended'] = [i in best for i in range(len(df))]
    return df


def main():
    args = parse_args()
    out_dir = Path(args.output_dir) if args.output_dir else Path(
        'work/run_{}'.format(datetime.now().strftime('%Y%m%d_%H%M%S')))
    out_dir.mkdir(parents=True, exist_ok=True)

    ab_chains = [c.strip() for c in args.antibody_chains.split(',') if c.strip()]
    ag_chains = resolve_antigen_chains(args.pdb, ab_chains, args.antigen_chains)

    print('=== [1/7] 读取单点打分 & 从 PDB 提取 master ===')
    scores_df = pd.read_csv(args.single_point_scores)
    validate_single_point_scores_df(scores_df, antibody_chains=ab_chains)
    # --- BEGIN PATCH ---
    # 派生 mutationHumanReadable（original_aa + location + mutant_aa），与 featurization 产出的
    # 格式一致，供 train_gp_model 按此列 join 特征与标签。输入契约 CSV 不含此列，故此处补齐。
    if 'mutationHumanReadable' not in scores_df.columns:
        scores_df['mutationHumanReadable'] = scores_df.apply(
            lambda r: '{}{}{}'.format(
                r['original_aa'], int(r['location']), r['mutant_aa']),
            axis=1,
        )
    # --- END PATCH ---
    master_seq, mapping = extract_master_from_pdb(args.pdb, ab_chains)
    write_fasta(out_dir / 'master.fasta', 'master_from_pdb', master_seq)
    pd.DataFrame(mapping).to_csv(out_dir / 'master_residue_mapping.csv', index=False)

    allowed = build_allowed_mutations_from_scores(scores_df, master_seq)
    save_allowed_mutations(allowed, out_dir / 'allowed_mutations.json')

    print('=== [2/7] 采样权重 & 界面 pairs ===')
    weights_df = compute_sampling_weights(scores_df)
    weights_df.to_csv(out_dir / 'sampling_weights.csv', index=False)

    pairs_raw = extract_interface_pairs(
        args.pdb, ab_chains, ag_chains, cutoff=args.interface_cutoff)
    pairs_path = out_dir / 'interface_pairs.json'
    write_interface_pairs_json(
        args.pdb, ab_chains, ag_chains, pairs_raw, pairs_path,
        cutoff=args.interface_cutoff,
    )
    interface_pairs, _ = load_interface_pairs(str(pairs_path))

    print('=== [3/7] 单点特征化 & 训练 GP（仅一次）===')
    sp_seq_df = build_single_point_seq_df(scores_df, master_seq)
    sp_seq_df.to_csv(out_dir / 'single_point_sequences.csv', index=False)
    sp_feat_df = featurize_sequences_df(sp_seq_df, master_seq, interface_pairs)
    sp_feat_df.to_csv(out_dir / 'single_point_features.csv', index=False)

    gp_dir = out_dir / 'gp_model'
    train_gp_model(
        sp_feat_df, scores_df, gp_dir,
        join_on='mutationHumanReadable',
        num_iters=args.gp_num_iters,
    )
    print('GP model saved to {}'.format(gp_dir))

    print('=== [4/7] BO 循环: 生成 → 特征化 → MEI（{} 轮，不重训 GP）==='.format(
        args.bo_rounds))
    collected = []
    seen_sequences = set(sp_seq_df['sequence'].tolist())
    for rnd in range(args.bo_rounds):
        seed_r = args.seed + rnd
        print('  --- BO round {}/{} (seed={}) ---'.format(
            rnd + 1, args.bo_rounds, seed_r))
        gen_rows = generate_mutant_sequences(
            master_seq, weights_df, allowed,
            number_to_generate=args.num_sequences_per_round,
            min_locations=args.min_mutations,
            max_locations=args.max_mutations,
            exclude_sequences=seen_sequences,
            seed=seed_r,
        )
        if not gen_rows:
            print('  WARNING: round {} 未生成新序列'.format(rnd))
            continue
        gen_df = pd.DataFrame(gen_rows)
        gen_df['bo_round'] = rnd
        for s in gen_df['sequence']:
            seen_sequences.add(s)
        gen_path = out_dir / 'generated_round_{:02d}.csv'.format(rnd)
        gen_df.to_csv(gen_path, index=False)

        feat_df = featurize_sequences_df(gen_df, master_seq, interface_pairs)
        feat_path = out_dir / 'features_round_{:02d}.csv'.format(rnd)
        feat_df.to_csv(feat_path, index=False)

        mei_df, best_sf = mei_select_sequences(
            feat_df, gp_dir, args.mei_batch_size, scores_df)
        mei_df['bo_round'] = rnd
        mei_path = out_dir / 'mei_selected_round_{:02d}.csv'.format(rnd)
        mei_df.to_csv(mei_path, index=False)
        print('  selected {} sequences; best_so_far={:.4f}'.format(
            len(mei_df), best_sf))
        collected.append(mei_df)

    if not collected:
        raise RuntimeError('BO 循环未产生任何 MEI 候选')

    print('=== [5/7] 合并 MEI 候选 ===')
    pool_df = pd.concat(collected, ignore_index=True)
    pool_df = pool_df.drop_duplicates(subset=['sequence'], keep='first')
    pool_df.to_csv(out_dir / 'mei_pool_all_rounds.csv', index=False)

    print('=== [6/7] 聚合五工具 objective & Pareto 降维 ===')
    objectives_df = aggregate_objectives(pool_df, scores_df)
    objectives_df.to_csv(out_dir / 'multipoint_objectives.csv', index=False)

    final_df = run_pareto_downselect(objectives_df, args.final_quota)
    recommendations = final_df[final_df['Recommended']].copy()
    recommendations = recommendations.sort_values('downselection_score', ascending=False)
    rec_path = out_dir / 'final_recommendations.csv'
    recommendations.to_csv(rec_path, index=False)
    write_recommendations_fasta(
        recommendations, out_dir / 'final_recommendations.fasta')

    summary = {
        'single_point_scores': str(Path(args.single_point_scores).resolve()),
        'pdb': str(Path(args.pdb).resolve()),
        'antibody_chains': ab_chains,
        'antigen_chains': ag_chains,
        'master_length': len(master_seq),
        'n_single_point_mutations': len(scores_df),
        'bo_rounds': args.bo_rounds,
        'num_sequences_per_round': args.num_sequences_per_round,
        'mei_batch_size': args.mei_batch_size,
        'mei_pool_size': len(pool_df),
        'final_recommendations': len(recommendations),
        'output_dir': str(out_dir.resolve()),
    }
    with open(out_dir / 'pipeline_summary.json', 'w') as f:
        json.dump(summary, f, indent=2)

    print('=== [7/7] 完成 ===')
    print('推荐 {} 条突变组合 -> {}'.format(len(recommendations), rec_path))
    print('FASTA -> {}'.format(out_dir / 'final_recommendations.fasta'))
    if len(recommendations):
        print('\nTop 5 突变组合:')
        for _, row in recommendations.head(5).iterrows():
            print('  {}  (pred_mean={:.3f}, mei={:.3f})'.format(
                row['mutationHumanReadable'],
                row.get('pred_mean', float('nan')),
                row.get('mei_score', float('nan')),
            ))


if __name__ == '__main__':
    main()
