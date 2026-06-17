# Wemol 批量作业脚本 — 北京微生物所抗体设计

针对**抗体单点突变多模块打分**的批量提交 / 回收 / 汇总管线 + CRPCA 序列优化。
输入一份突变位点 + PDB 复合物结构，跑 WeMol 多模块（Flex DDG / FoldX / Protein FEP / AbNatiV / SFE），
输出一张汇总 CSV，再经 CRPCA 做组合突变枚举与筛选得到优势突变。

> 📖 **完整操作流程**（从初始位点选择 → 单点打分 → 序列优化的每一步命令）见 **[使用手册.md](使用手册.md)**。
> 本 README 是面向开发/维护的补充说明：文件结构、脚本参数速查、常见问题、已知限制。

---

## 文件结构

```
beijing_weishengwusuo/
├── 使用手册.md                 完整操作流程文档（交付/上手看这份）
├── README.md                  本文件（开发/维护补充）
├── config_example.json        登录配置模板（复制为 config.json）
├── prepare_inputs.py          生成突变列表 + tasks.csv（内置零依赖突变生成）
├── submit_batch.py            批量提交
├── fetch_batch.py             批量回收 + 汇总成 scores.csv
├── build_pipeline_input.py    scores.csv → run_pipeline 就绪的 single_point_scores.csv
├── sfe_md.py                  SFE：PDB Mutation + MD + 抽帧 + rechain 编排（断点续跑）
├── sfe_prepare.py             SFE：由构象帧生成 forward/reverse 的 Flex DDG 任务
├── sfe_fetch.py               SFE：回收 ddG → 42 列宽表 → step02c 出 sfe_scores.csv
├── sfe_split_frames.py        SFE 内部：多帧 PDB 拆单帧（被 sfe_md 调用）
├── sfe_rechain.py             SFE 内部：GROMACS 帧 → Rosetta 兼容（被 sfe_md 调用）
├── _common.py                 共享工具（session / 结果解析 / IO）
├── _crpca_download/CRPCA/     CRPCA 序列优化流水线（step02c / run_pipeline 等）
└── test_files/                示例 PDB（7l7e）+ mut_list
```

---

## 脚本与参数速查

操作顺序与完整示例见 [使用手册.md](使用手册.md)，这里只列各脚本的关键/可选参数。

### `prepare_inputs.py`（生成突变列表 + tasks.csv）

- 两种突变来源：**模式 A**（`--muts "M31,M53,M55-56"` 从 PDB 零依赖枚举，推荐）／**模式 B**（`--mut-list single_mut_list.txt` 用现成列表）
- `--omit-aas CP`：排除的目标氨基酸（默认 CP）
- `--modules flexddg,foldx,fep,abnativ`：任选组合
- `--limit N`：只取前 N 个突变（**首次/测试强烈建议**）
- `--fep-type B`：FEP 类型 B=Binding / S=Stability（默认 B）
- 自动按抗体 J 区基序判定重/轻链并输出 AbNatiV 用 FASTA；判定有误可用 `--heavy-chains` / `--light-chains` 覆盖
- 位置编号 = PDB ATOM 行的 resSeq

### `submit_batch.py`（批量提交）

- `--delay 10` 每条提交后等待秒数；`--retry-delay 120` 触并发上限后重试间隔
- 断点续跑：中断后重跑同一命令，已提交的 task_name 自动跳过

### `fetch_batch.py`（回收 + 汇总）

- 可重复运行：未完成的下次再收
- `--foldx-interface`：默认 `M_S+N_S`（抗体-抗原结合界面之和）；单界面用 `M_S`，多界面用 `+` 求和
- `--include-multipoint` 保留多点突变（默认仅单点）；`--no-aggregate` 仅 fetch 不汇总；`--legacy-output` 输出旧版原始列
- ⚠ 输出的 `location` 是 **PDB resSeq**，进入序列优化前**必须**经 `build_pipeline_input.py` 修正为 master 线性位置

### SFE 子流程（`sfe_md.py` → `sfe_prepare.py` → `sfe_fetch.py`）

- `sfe_md.py`：`--receptor-chain` / `--ligand-chain`、`--sim-ns`、`--skip-ps`、`--n-frames`；首次真跑后若定位不到结果文件，调 `--mutant-pdb-glob` / `--traj-glob`
- 测试小跑：`--sim-ns 1 --skip-ps 500 --n-frames 3 --limit 2`
- `sfe_fetch.py`：`--run-step02c` 直接出分；`--step02c-python` 指向装了 numpy/pandas 的解释器

### `build_pipeline_input.py`（序列优化就绪输入）

- `--antibody-chains M,N`（重链在前、轻链在后，决定 master 拼接顺序）；`--sfe` 可省略；`--strict` 校验 `original_aa` 不一致即报错

---

## CRPCA 本地适配补丁

`_crpca_download/CRPCA/` 含 3 处本地适配补丁（均标注 `BEGIN/END PATCH`），使 `run_pipeline` 可运行：

- `abag_ml/pareto_selection.py`：缺失模块改可选导入
- `local_pipeline/run_pipeline.py`：派生 `mutationHumanReadable` + Pareto 跳过全 NaN 目标
- `local_pipeline/common/objective_aggregate.py`：按工具独立处理 NaN（支持 sfe 留空）

`run_pipeline.py` 依赖 `torch` / `gpytorch` / `biopython==1.79` / `scipy` / `scikit-learn`，CPU 即可（实测 426 突变 ≈ 70 秒，无需 GPU/容器）。

---

## 常见问题

**Q：怎么知道某模块在 WeMol 上的真实名字？**
用 `_common.build_session` 登录后 `query_module_by_name("Flex DDG")` 查 methods/参数名，一字不差填到任务的 `params`。

**Q：提交一半挂了 / Ctrl+C 了？**
重跑同一条 `submit_batch.py`——已成功提交的 task_name 自动跳过。

**Q：fetch 后 scores.csv 某列空着？**
可能：对应作业还没 Done（再等再跑）；模块名没对上（看日志 `未注册解析器`）；结果文件路径变了（看 `找不到 ddG.csv` 等 warning）。

**Q：scores.csv 突变数比预期少？**
默认过滤多点突变，要保留加 `--include-multipoint`。

**Q：怎么加新模块的解析？**
在 `_common.py` 的 `PARSERS` 加 `"模块真名": parse_xxx`，并写 `parse_xxx(result_dir, task) -> {mutation: {col: value}}`；`prepare_inputs.py` 相应加该模块的参数分支。

---

## 已知限制 / 技术说明

- 打分模块 **Flex DDG / FoldX / Protein FEP / AbnatiV / SFE** 全部已打通。
  - AbnatiV 需在 WeMol 把 `*_seq_scores.csv` 绑定为输出文件，否则 SDK 下载不到（按残基的 `*_res_scores.csv` 无序列总分，无法替代）。
- **多点突变**：FoldX 原生支持；Flex DDG / FEP 当前按单点逻辑处理，跳过多点。
- **SFE 帧兼容**：GROMACS MD 帧与 Rosetta 不兼容（链/编号被改、带氢、AMBER 命名），`sfe_rechain.py` 负责还原链/编号、去氢、残基/原子名归一，否则 Flex DDG 会失败。
- **依赖分两套**：提交/回收用 `wemol-python-sdk`；聚合/优化用 `numpy pandas`（step02c）、`biopython==1.79`（build_pipeline_input，新版 biopython 移除了所需的 `three_to_one`）、`torch gpytorch scipy scikit-learn`（run_pipeline）。
- **算力约束**：Flex DDG 单作业按默认 nstruct 约半天，FEP/SFE 更重且吃 GPU。全量饱和扫描需相应规划 GPU 算力与时间。
- Flex DDG 内部有 stochastic sampling，同一突变多次运行 ddG 值有 ~1 单位浮动。
