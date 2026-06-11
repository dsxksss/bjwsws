# 本地抗体设计流水线

无 HPC、无湿实验。假设 **2.2 单点五工具打分** 已完成（CSV）。

## 一键运行（推荐）

```bash
cd /data/PRG/tools/apps/CRPCA
export PYTHONPATH=$PWD

python local_pipeline/run_pipeline.py \
  --single-point-scores path/to/single_point_scores.csv \
  --pdb structure/7l7e.pdb \
  --antibody-chains M,N \
  --antigen-chains S \
  --output-dir work/my_run \
  --bo-rounds 3 \
  --num-sequences-per-round 2000 \
  --mei-batch-size 50 \
  --min-mutations 1 \
  --max-mutations 8 \
  --final-quota 20 \
  --seed 0
```

### 必需输入


| 参数                      | 说明                                   |
| ----------------------- | ------------------------------------ |
| `--single-point-scores` | 2.2 单点 CSV（见下表列定义）                   |
| `--pdb`                 | 初始共晶/ Fab PDB                        |
| `--antibody-chains`     | 重链,轻链… 按拼接顺序，如 `M,N`                 |
| `--antigen-chains`      | 抗原链，如 `S`；**可省略**，自动取 PDB 中除抗体链外的所有链 |


### 主要可调参数


| 参数                                    | 默认    | 说明                 |
| ------------------------------------- | ----- | ------------------ |
| `--bo-rounds`                         | 3     | Step03↔Step06 循环次数 |
| `--num-sequences-per-round`           | 2000  | 每轮生成候选数            |
| `--mei-batch-size`                    | 50    | 每轮 MEI 保留数         |
| `--min-mutations` / `--max-mutations` | 1 / 8 | 论文 Step A 突变数范围    |
| `--final-quota`                       | 20    | 最终推荐突变组合数          |
| `--gp-num-iters`                      | 500   | GP 训练迭代（**只训一次**）  |


### 最终输出


| 文件                            | 内容                                                        |
| ----------------------------- | --------------------------------------------------------- |
| `final_recommendations.csv`   | 推荐突变组合、`mutationHumanReadable`、全长 `sequence`、GP 预测、五工具聚合分 |
| `final_recommendations.fasta` | 同上，FASTA 格式                                               |
| `pipeline_summary.json`       | 运行摘要                                                      |
| `mei_pool_all_rounds.csv`     | 所有 BO 轮 MEI 候选合并                                          |
| `gp_model/`                   | 一次性训练的 GP（循环中不重训）                                         |


---

## 流程说明

```
单点 CSV + PDB
  → master.fasta（自动从 PDB 提取）
  → sampling_weights
  → interface_pairs.json
  → 单点特征 → 训练 GP（一次）
  → [循环 bo_rounds 次]
       Step A: k ~ Uniform[min, max]
       Step B: 加权采样 k 个突变 → 多点序列
       → 特征化 → MEI 选 Top-N
  → 合并 MEI 池 → 五工具求和 → Pareto 降维 → 最终推荐
```

**GP 不重训**；仅 **生成 + MEI** 按 `--bo-rounds` 循环，每轮用不同 seed 并排除已见序列。

---

## 2.2 单点 CSV 列定义


| 列                                               | 说明                                              |
| ----------------------------------------------- | ----------------------------------------------- |
| `mutation`                                      | 突变标签（如 `G112E`），便于人工查阅与下游追溯                      |
| `chain`                                         | PDB 链 ID，须属于 `--antibody-chains`（如 `M`、`N`）      |
| `location`                                      | master 序列 1-based 位置（`--antibody-chains` 按顺序拼接后的线性编号） |
| `original_aa`                                   | WT 氨基酸（须与 PDB 提取的 master 一致）                    |
| `mutant_aa`                                     | 突变氨基酸                                           |
| `sfe`, `fep`, `rosetta_flex`, `foldx`, `abnativ` | 五工具原始分数                                       |

`allowed_mutations.json` 会从 CSV **自动生成**。

### 从 42 个 Flex ddG 计算 SFE（step02c）

若已有 21 个 forward + 21 个 reverse 的 Rosetta Flex ddG，用宽表 CSV 输入：

```bash
python local_pipeline/step02c_compute_sfe_from_flex_ddg.py \
  --input local_pipeline/examples/sfe_flex_ddg.example.csv \
  --output work/sfe_scores.csv
```

宽表列定义（一行一个突变）：

```csv
mutation,chain,forward_01,forward_02,...,forward_21,reverse_01,...,reverse_21
G112E,M,-1.68,-2.22,...,reverse_21
```

`chain` 为 PDB 链 ID，标明突变所在抗体链（须与后续 2.2 CSV 中 `chain` 一致）。

公式：`SFE = (IQR_mean(forward) - IQR_mean(reverse)) / 2`。原始 `sfe_ddg` 越负表示结合越强；写入 2.2 CSV 时若需「越大越好」请加 `--invert-for-sampling`。

---

## 分步脚本（可选手动执行）


| 步骤  | 脚本                                       |
| --- | ---------------------------------------- |
| 02  | `step02_compute_sampling_weights.py`     |
| 02c | `step02c_compute_sfe_from_flex_ddg.py`   |
| 02a | `step02a_extract_master_fasta.py`        |
| 02b | `step02b_single_point_sequences.py`      |
| 03b | `step03b_extract_interface_pairs.py`     |
| 04  | `step04_featurize_sequences.py`          |
| 05  | `step05_train_gp_model.py`               |
| 03  | `step03_generate_sequences.py`（Step A+B） |
| 06  | `step06_mei_select_sequences.py`         |
| 07  | `step07_aggregate_objectives.py`         |
| 08  | `step08_pareto_select.py`                |


Step 03 参数：`--min-locations` / `--max-locations`（恢复论文 Step A）。

---

## 依赖

Python 3.8+，`numpy`, `pandas`, `torch`, `gpytorch`, `biopython`