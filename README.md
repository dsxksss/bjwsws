# Wemol 批量作业脚本 — 北京微生物所抗体设计

针对**抗体单点突变多模块打分**的批量提交 / 回收 / 汇总管线。

## 适用场景

- 输入：一份突变列表（如 ~n 个单点突变）+ PDB 复合物结构
- 跑 Wemol 平台多个模块（当前 MVP 含 **Flex DDG / Mutation Energy of Binding / Protein FEP**）
- 输出：一张汇总 CSV——第一列突变名，其余列各模块各打分

由于单条 Flex DDG 作业约 2 小时，全量作业需要数天，**提交和回收必须分两步、能断点续跑**。

---

## 文件结构

```
beijing_weishengwusuo/
├── README.md                  本文件
├── config_example.json        登录配置模板（复制为 config.json）
├── prepare_inputs.py          Step 0：突变位点/列表 → tasks.csv（内置零依赖突变生成）
├── submit_batch.py            Step 1：批量提交
├── fetch_batch.py             Step 2：批量回收 + 汇总成 scores.csv
├── build_pipeline_input.py    Step 3 桥接：scores.csv → run_pipeline 就绪的 single_point_scores.csv
├── _common.py                 共享工具（session / 解析器 / IO）
├── mutation_list.py           实验室原版突变生成脚本（备用，需 anarci+biopython，产 H/L FASTA）
├── _crpca_download/CRPCA/     CRPCA 序列优化流水线（step02c / run_pipeline 等，Step 3 用）
├── test_files/                示例 PDB + mut_list（一个真实跑过的小样）
└── _dev/                      开发期工具（mock_pipeline.py 等，非交付件）
```

---

## 准备工作（一次性）

### 1. 配置登录

复制 `config_example.json` 为 `config.json`，填入你的 Wemol 地址与账号：

```json
{
  "base_url": "your_wemol_base_url",
  "user": {
    "Name": "your_wemol_username",
    "Passwd": "your_wemol_password"
  },
  "result_dir": "wemol_results"
}
```

### 2. 准备突变列表（可选）

突变列表有两种来源，**默认推荐模式 A，无需任何额外准备**：

- **模式 A（推荐）**：不用提前准备。`prepare_inputs.py` 内置零依赖突变生成器，
  Step 0 直接给 `--muts` 位点描述即可，自动从 PDB 枚举突变。
- **模式 B**：如果你已经有现成的 `single_mut_list.txt`（每行一个 `DM31L`），
  Step 0 用 `--mut-list` 指定即可。

> AbNatiV 用的 `Mut_Hchains.fasta` / `Mut_Lchains.fasta` 也由模式 A 自动生成
> （自动判定轻/重链，无需 anarci，详见 Step 0）。
> `mutation_list.py` 是实验室原版脚本（需 anarci+biopython），仅作备用保留。

---

## 工作流

### Step 0：生成 tasks.csv

**模式 A（推荐，零依赖直接从位点生成）**：

```powershell
python prepare_inputs.py `
    --pdb test_files/7l7e-ab-ag-complex.pdb `
    --muts "M31,M53,M55-56,M104-113,M116,N31-36,N38,N55-56,N59" `
    --omit-aas CP `
    --out-dir prod `
    --modules flexddg,foldx,fep,abnativ
```

- `--muts` 位点描述：`链+位置`，支持范围，如 `M31,M53,M55-56`
- `--omit-aas CP` 排除的目标氨基酸（默认 `CP`，即不突变成 C / P）
- 位置编号 = PDB 文件 ATOM 行的残基序号（resSeq），与实验室脚本一致
- 会自动判定每条被突变链是重链还是轻链（基于抗体 J 区保守基序），
  输出 AbnatiV 用的 `Mut_Hchains.fasta` / `Mut_Lchains.fasta`；
  如判定有误可用 `--heavy-chains M --light-chains N` 手动覆盖
- `--modules` 可选 `flexddg / foldx / fep / abnativ`；选 `abnativ` 会自动用上面的
  FASTA 生成 VH/VL 两个作业（Nat 默认 VH/VL，可用 `--heavy-nat`/`--light-nat` 改）

**模式 B（用现成的突变列表文件）**：

```powershell
python prepare_inputs.py `
    --mut-list single_mut_list.txt `
    --mut-list-foldx mut_list_for_foldx.txt `
    --pdb test_files/7l7e-ab-ag-complex.pdb `
    --out-dir prod `
    --modules flexddg,foldx,fep
```

通用可选参数：

- `--limit 5` 仅取前 N 个突变（**强烈建议首次跑用 `--limit 5` 先验证**）
- `--modules flexddg` 只准备某个模块的任务
- `--fep-type B` FEP 类型，`B`=Binding affinity / `S`=Stability（默认 B）

输出：

```
prod/
├── single_mut_list.txt        （模式 A 自动生成）
├── mut_list_for_foldx.txt     （模式 A 自动生成）
├── Mut_Hchains.fasta          （模式 A 自动生成，重链突变序列，AbnatiV 用）
├── Mut_Lchains.fasta          （模式 A 自动生成，轻链突变序列，AbnatiV 用）
├── mutlists_flexddg/          每个突变一个 mut 文件（FlexDDG）
│   ├── DM31L.txt
│   ├── DM31V.txt
│   └── ...
├── mutlist_foldx.txt          所有突变一个文件（FoldX）
└── tasks.csv                  喂给 submit_batch 的任务清单
```

跑完会显示作业总数，**确认数量符合预期再进 Step 1**。

### Step 1：批量提交

```powershell
python submit_batch.py `
    --config config.json `
    --tasks prod/tasks.csv `
    -o prod/submitted.jsonl
```

- 每提交成功一条立即记录到 `submitted.jsonl`（断点续跑友好）
- 遇到平台并发上限（`JobRunMaxNumLimit`）自动 sleep 120s 重试该作业
- 中途 Ctrl+C 没事，**重跑同一条命令**：已提交的 task_name 会自动跳过

可选参数：

- `--delay 10` 每条提交后等待秒数（默认 10）
- `--retry-delay 120` 触限流后重试间隔（默认 120）

### Step 1.5：手动查进度

打开 Wemol 网页 → 作业列表，等所有作业跑完。

### Step 2：批量回收 + 汇总

```powershell
python fetch_batch.py `
    --config config.json `
    --input prod/submitted.jsonl `
    --output prod/scores.csv
```

依次做 3 件事：

1. 查每条 Submitted/Doing 作业状态
2. Done/Abort/Cancel 的下载结果，更新 `submitted.jsonl` 状态
3. 扫所有 Done 结果目录 → 各模块解析 → 输出 `scores.csv`

**可重复运行**：还有作业 Doing 时本次只汇总已完成部分，剩余下次再跑。

可选参数：

- `--no-aggregate` 仅 fetch 不汇总
- `--include-multipoint` 保留多点突变（默认仅单点）
- `--foldx-interface M_N` 指定 FoldX 哪个 interface 写入 `foldx` 列；不指定且存在多个 interface 时，`foldx` 暂时留空
- `--legacy-output` 输出旧版模块原始列（如 `flexddg_ddg` / `foldx_interface_*_ddg`）

默认输出会直接对齐 `run_pipeline.py --single-point-scores` 需要的格式：

```csv
mutation,chain,location,original_aa,mutant_aa,sfe,fep,rosetta_flex,foldx,abnativ
DM31L,M,31,D,L,,-8.123,-11.038,,0.9000
DM31V,M,31,D,V,,-7.456,-10.221,,0.8765
```

> `sfe` 暂时留空；`abnativ` 来自 AbnatiV 模块（取 `AbNatiV V* Score` 总分，4 位小数）。FoldX interface 未确认前默认留空，确认后加 `--foldx-interface M_S` 等填入。

⚠ 注意：`fetch_batch` 输出的 `location` 是 **PDB resSeq**，而下游 `run_pipeline.py` 需要的是
**master 序列线性位置**（重链+轻链拼接后的 1-based）。两者对重链一致、对轻链相差一个重链长度。
进入 Step 3 前**必须**经 `build_pipeline_input.py` 修正（见下）。

---

## 三、序列优化设计（Step 3）

把 Step 2 的 `scores.csv` 喂给 CRPCA 的 `run_pipeline.py`，做多点突变枚举 + GP/MEI + Pareto 筛选，
得到最终推荐突变组合。CRPCA 代码在 `_crpca_download/CRPCA/`。

### 3.1（可选）计算 SFE

SFE 需要每个突变 21 forward + 21 reverse 共 42 个 Rosetta Flex ddG（来自 MD 多构象采样，
目前未自动化）。**有了 42 值宽表后**，用 CRPCA 的 `step02c` 聚合（只需 numpy/pandas，可在本 venv 跑）：

```bash
PYTHONPATH=_crpca_download/CRPCA python _crpca_download/CRPCA/local_pipeline/step02c_compute_sfe_from_flex_ddg.py \
    --input prod/sfe_flex_ddg_wide.csv \
    --output prod/sfe_scores.csv
```

宽表列：`mutation,chain,forward_01..forward_21,reverse_01..reverse_21`（`mutation` 用与
`scores.csv` 一致的命名，如 `DM31L`）。公式 `SFE=(IQR_mean(forward)-IQR_mean(reverse))/2`。

### 3.2 生成 run_pipeline 就绪输入（必做）

`build_pipeline_input.py` 把 `scores.csv` 转成 `single_point_scores.csv`：修正 `location`（用 CRPCA
自己的 `extract_master_from_pdb`，轻链自动加偏移）、校验 `original_aa`、可选合并 SFE。需 `biopython`。

```bash
python build_pipeline_input.py \
    --scores prod/scores.csv \
    --pdb test_files/7l7e-ab-ag-complex.pdb \
    --antibody-chains M,N \
    --sfe prod/sfe_scores.csv \
    --output prod/single_point_scores.csv
```

- `--antibody-chains M,N`：重链在前、轻链在后，决定 master 拼接顺序与线性位置
- `--sfe` 可省略（没算 SFE 时 sfe 列留空，run_pipeline 会自动跳过该工具）
- `--strict`：`original_aa` 与 master 不一致时直接报错（默认仅警告并跳过该行）

### 3.3 运行优化流水线

`run_pipeline.py` 依赖 `torch` / `gpytorch` / `biopython` / `scipy` / `scikit-learn`，
CPU 即可（GP 数据量小）。本 venv 装好这些依赖后**直接跑**，无需容器：

```bash
# 一次性装依赖（biopython 固定 1.79；torch 用 CPU 版即可）
pip install torch --index-url https://download.pytorch.org/whl/cpu
pip install "biopython==1.79" gpytorch scipy scikit-learn

# 运行（真实规模 426 突变 / 3 轮 × 2000 序列 / GP 500 迭代 ≈ 70 秒）
PYTHONPATH=_crpca_download/CRPCA python _crpca_download/CRPCA/local_pipeline/run_pipeline.py \
    --single-point-scores prod/single_point_scores.csv \
    --pdb 7l7e-ab-ag-complex.pdb \
    --antibody-chains M,N \
    --antigen-chains S \
    --output-dir prod/optimization_run \
    --bo-rounds 3 --num-sequences-per-round 2000 \
    --mei-batch-size 50 --min-mutations 1 --max-mutations 8 \
    --final-quota 20 --seed 0
```

输出 `final_recommendations.csv`（推荐突变组合 + 全长序列 + 五工具聚合分）等，详见
`_crpca_download/CRPCA/local_pipeline/README.md`。

> CPU 实测约 70 秒，无需 GPU/容器。`Dockerfile.crpca` 保留作将来 GPU 容器化的参考（如需）。
>
> ⚠ CRPCA 代码含 3 处 AI 修复补丁（均标注 `BEGIN/END PATCH`），使 run_pipeline 可运行：
> - `abag_ml/pareto_selection.py`：缺失模块改可选导入
> - `local_pipeline/run_pipeline.py`：派生 `mutationHumanReadable` + Pareto 跳过全 NaN 目标
> - `local_pipeline/common/objective_aggregate.py`：按工具独立处理 NaN（支持 sfe 留空）

---

## 5 分钟最小可行验证

用项目自带的小样 + 你的 config.json 真跑一遍（1 个位点 Flex DDG 单模块，`--limit 1`）：

```powershell
python prepare_inputs.py --pdb test_files/7l7e-ab-ag-complex.pdb --muts "M31" --out-dir _try --modules flexddg --limit 1

python submit_batch.py --config config.json --tasks _try/tasks.csv -o _try/submitted.jsonl

# 等 ~2 小时（Flex DDG 单条耗时）
python fetch_batch.py --config config.json --input _try/submitted.jsonl --output _try/scores.csv

type _try\scores.csv
```

期望最后看到：

```
mutation,chain,location,original_aa,mutant_aa,sfe,fep,rosetta_flex,foldx,abnativ
DM31L,M,31,D,L,,,<-11 左右的数字>,,
```

---

## 常见问题

### Q1：怎么知道某个模块在 Wemol 上的真实名字？

跑：

```python
from _common import build_session
s = build_session("config.json")
mod = s.query_module_by_name("Flex DDG")    # 模糊查
for m in mod.methods:
    print(m.name)
    for a in (m.input_args or []):
        print(" ", a.name, "|", a.type, "| required:", a.required)
```

把打印的 `name` 一字不差填到 CSV 的 `params` 里。

### Q2：提交一半挂了 / Ctrl+C 了？

直接重跑 `submit_batch.py` 同一条命令——已成功提交的 task_name 会自动跳过。

### Q3：fetch 之后 scores.csv 里某列空着？

可能原因：

- 对应模块作业还没 Done（再等再跑 fetch）
- 模块名字符串没对上（看 fetch 日志里 `未注册解析器` warning）
- 结果文件路径变了（看 `找不到 ddG.csv` / `找不到 result.txt` warning）

### Q4：scores.csv 里突变数比预期少？

默认过滤了多点突变（如 `DM31L,SM55K`）。要保留就加 `--include-multipoint`。

### Q5：怎么加新模块（如 AbNatiV / PDB Mutation / MD）？

两步：

1. 在 `_common.py` 的 `PARSERS` dict 加一项：`"模块真名": parse_xxx`
2. 写对应的 `parse_xxx(result_dir, task) -> {mutation: {col_name: value}}` 函数

`prepare_inputs.py` 也要相应加分支处理新模块的参数格式。

> AbNatiV 用的 `Mut_Hchains.fasta` / `Mut_Lchains.fasta` 已由模式 A 自动生成，无需额外准备。

---

## 已知限制 / 待办

- 打分模块支持 **Flex DDG / Mutation Energy of Binding / Protein FEP / AbnatiV**，全部已打通
  - AbnatiV：提交侧自动用模式 A 产出的 FASTA 生成 VH/VL 两个作业；汇总侧解析
    `*_abnativ_seq_scores.csv` 的 `AbNatiV V* Score` 总分 → `abnativ` 列
  - ⚠ AbnatiV 模块需在 WeMol 把 `*_abnativ_seq_scores.csv` 绑定为输出文件，
    否则 SDK 下载不到（按残基的 `*_res_scores.csv` 没有序列总分，无法替代）
- **序列优化（Step 3）**：已打通。`build_pipeline_input.py`（location 修正 + SFE 合并）+
  `run_pipeline.py`（CPU ≈ 70 秒，本 venv 直接跑）。CRPCA 代码含 3 处 AI 修复补丁（见 3.3）
- **PDB Mutation / MD** 待模块上线后补 parser
- **SFE 流程**：`step02c` 已能算（输入 42 个 flex ddG 宽表）；但宽表的来源
  （MD → 21 构象 → forward/reverse flex ddG）尚未自动化，故 `sfe` 列暂空
- **多点突变**：FoldX 原生支持；Flex DDG 和 FEP 当前按单点逻辑处理，跳过多点
- **依赖**：`build_pipeline_input.py` 需 `biopython==1.79`（新版 biopython 移除了
  CRPCA 用到的 `three_to_one`，故固定 1.79）
- Flex DDG 内部有 stochastic sampling，同一突变多次运行 ddG 值有 ~1 单位浮动
