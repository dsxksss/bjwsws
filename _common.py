"""共享工具：session 构建、任务清单加载、JSONL 读写、单作业提交分发。"""
from __future__ import annotations

import csv
import json
import logging
import re
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any, Callable

log = logging.getLogger("batch_jobs")


# ---------------------------------------------------------------------------
# 时间
# ---------------------------------------------------------------------------
def now_iso() -> str:
    return datetime.now().strftime("%Y-%m-%dT%H:%M:%S")


# ---------------------------------------------------------------------------
# Session
# ---------------------------------------------------------------------------
def build_session(config_path: str | Path):
    """
    读 JSON 配置，构建 UserSession。

    config.json 示例（账号密码登录）：
    {
      "base_url": "http://192.168.1.145:8200",
      "user":     {"Name": "your_name", "Passwd": "your_password"},
      "result_dir":     "wemol_results",
      "log_file":       "batch_jobs.log",
      "cache_file":     "batch_jobs_cache.json"
    }

    config.json 示例（ant_uid 登录）：
    {
      "base_url": "http://192.168.1.145:8200",
      "ant_uid":  "xxxx-session-token-xxxx",
      "user_agent": "Mozilla/5.0 ...",
      "result_dir":     "wemol_results",
      "log_file":       "batch_jobs.log",
      "cache_file":     "batch_jobs_cache.json"
    }
    """
    # SDK 懒加载：纯本地脚本（aggregate / mock）不需要装 SDK 也能跑
    from wemol_sdk import wemol, DEBUG_LOG_LEVEL
    from wemol_sdk.models.config import config_field
    from wemol_sdk.models.request_models import UserReq

    cfg_path = Path(config_path)
    if not cfg_path.exists():
        raise FileNotFoundError(f"找不到配置文件: {cfg_path}")
    cfg = json.loads(cfg_path.read_text(encoding="utf-8"))

    if "base_url" not in cfg:
        raise ValueError("config.json 缺少 base_url 字段")

    sdk_config: dict[str, Any] = {
        config_field.BASE_URL:             cfg["base_url"],
        config_field.RESULT_FILE_SAVE_DIR: cfg.get("result_dir", "wemol_results"),
        config_field.LOG_LEVEL:            DEBUG_LOG_LEVEL,
        config_field.LOG_FILE_SAVE_PATH:   cfg.get("log_file",   "batch_jobs.log"),
        config_field.CACHE_FILE_SAVE_PATH: cfg.get("cache_file", "batch_jobs_cache.json"),
    }
    if "user_agent" in cfg:
        sdk_config[config_field.USER_AGENT] = cfg["user_agent"]

    if "ant_uid" in cfg and cfg["ant_uid"]:
        log.info("使用 ant_uid 登录")
        return wemol.createUserSessionUseAntUid(ant_uid=cfg["ant_uid"], config=sdk_config)

    if "user" not in cfg:
        raise ValueError("config.json 必须提供 user (Name/Passwd) 或 ant_uid 二选一")
    log.info(f"使用账号密码登录: {cfg['user'].get('Name')}")
    return wemol.createUserSession(user=UserReq(**cfg["user"]), config=sdk_config)


# ---------------------------------------------------------------------------
# 任务清单加载
# ---------------------------------------------------------------------------
def _parse_csv_params(raw: str) -> dict[str, Any]:
    """CSV 单元格里 `k1=v1;k2=v2` 解析成 dict。空字符串返回 {}。"""
    out: dict[str, Any] = {}
    if not raw or not raw.strip():
        return out
    for chunk in raw.split(";"):
        chunk = chunk.strip()
        if not chunk:
            continue
        if "=" not in chunk:
            raise ValueError(f"params 段缺少 '=': {chunk!r}")
        k, v = chunk.split("=", 1)
        out[k.strip()] = v.strip()
    return out


def load_tasks(path: str | Path) -> list[dict[str, Any]]:
    """
    自动按后缀加载 .csv / .json 任务清单。

    每条任务统一规范成字段：
    - task_name   (必填)
    - module_name / module_id / flow_name / flow_id (四选一)
    - method_name (可选, 仅 run_job)
    - params      (dict)
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"找不到任务清单: {p}")

    suffix = p.suffix.lower()
    tasks: list[dict[str, Any]] = []

    if suffix == ".json":
        data = json.loads(p.read_text(encoding="utf-8"))
        if not isinstance(data, list):
            raise ValueError("JSON 任务清单顶层必须是数组")
        tasks = data

    elif suffix == ".csv":
        with p.open("r", encoding="utf-8-sig", newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                row = {k: (v.strip() if isinstance(v, str) else v) for k, v in row.items()}
                task: dict[str, Any] = {
                    "task_name":   row.get("task_name") or "",
                    "method_name": row.get("method_name") or None,
                    "params":      _parse_csv_params(row.get("params", "")),
                }
                # 把四个互斥字段按非空规则带入
                for key in ("module_name", "module_id", "flow_name", "flow_id"):
                    val = row.get(key)
                    if val:
                        task[key] = int(val) if key.endswith("_id") else val
                tasks.append(task)
    else:
        raise ValueError(f"不支持的任务清单格式: {suffix} (仅支持 .csv / .json)")

    # 基础校验
    for i, task in enumerate(tasks, start=1):
        if not task.get("task_name"):
            raise ValueError(f"第 {i} 条任务缺少 task_name")
        has_target = any(task.get(k) for k in ("module_name", "module_id", "flow_name", "flow_id"))
        if not has_target:
            raise ValueError(
                f"task_name={task['task_name']} 必须提供 module_name / module_id / flow_name / flow_id 之一"
            )
        task.setdefault("params", {})

    return tasks


# ---------------------------------------------------------------------------
# 单作业提交：根据字段自动选 run_job / run_job_use_flow
# ---------------------------------------------------------------------------
def submit_one(session, task: dict[str, Any]) -> int:
    """
    提交单条任务。返回 job_id（>=1 表示成功，-1 表示失败）。
    可能抛出 Exception（如 JobRunMaxNumLimit），调用方负责捕获并决定重试。
    """
    params = task.get("params") or {}

    # 当前 SDK(v2.0.0) 的 run_job / run_job_use_flow 只认「名字」，
    # 没有 module_id / flow_id 形参（传了会 TypeError）。任务清单统一用名字。
    if task.get("flow_id") or task.get("flow_name"):
        return session.run_job_use_flow(
            flow_name=task.get("flow_name"),
            params=params,
            use_delay=False,
            result_dir_use_task_name=True,
        )

    return session.run_job(
        module_name=task.get("module_name"),
        params=params,
        method_name=task.get("method_name"),
        use_delay=False,
        result_dir_use_task_name=True,
    )


# ---------------------------------------------------------------------------
# JSONL 读写
# ---------------------------------------------------------------------------
def append_jsonl(path: str | Path, record: dict[str, Any]) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("a", encoding="utf-8") as f:
        f.write(json.dumps(record, ensure_ascii=False) + "\n")


def read_jsonl(path: str | Path) -> list[dict[str, Any]]:
    p = Path(path)
    if not p.exists():
        return []
    out: list[dict[str, Any]] = []
    with p.open("r", encoding="utf-8") as f:
        for ln, line in enumerate(f, start=1):
            line = line.strip()
            if not line:
                continue
            try:
                out.append(json.loads(line))
            except json.JSONDecodeError as e:
                log.warning(f"{p}:{ln} 行 JSON 解析失败，已跳过: {e}")
    return out


def rewrite_jsonl(path: str | Path, records: list[dict[str, Any]]) -> None:
    """原地重写整个 jsonl（fetch 更新状态时用）。先写临时文件再替换以保证原子性。"""
    p = Path(path)
    tmp = p.with_suffix(p.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8") as f:
        for rec in records:
            f.write(json.dumps(rec, ensure_ascii=False) + "\n")
    tmp.replace(p)


def already_submitted(jsonl_path: str | Path, task_name: str) -> bool:
    """同 task_name 已存在且非 Failed 时视作已提交，直接跳过。"""
    for rec in read_jsonl(jsonl_path):
        if rec.get("task_name") == task_name and rec.get("status") != "Failed":
            return True
    return False


# ---------------------------------------------------------------------------
# 结果汇总：解析各模块的结果文件 → 拼成最终 scores.csv
# ---------------------------------------------------------------------------
def find_result_dir(base: Path, job_id: int) -> Path | None:
    """SDK 把作业结果落到 base/<job_id>-<module_safe>/ 下；用 glob 通配 module 名差异。"""
    matches = list(base.glob(f"{job_id}-*"))
    if not matches:
        matches = list(base.glob(f"{job_id}"))
    return matches[0] if matches else None


def parse_flexddg(result_dir: Path, task: dict) -> dict[str, dict[str, float]]:
    """Flex DDG: 找 ddG.csv，单行格式 mut,ddG。"""
    files = list(result_dir.rglob("ddG.csv"))
    if not files:
        log.warning(f"[Flex DDG] {result_dir} 找不到 ddG.csv")
        return {}
    out: dict[str, dict[str, float]] = {}
    with files[0].open(encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f):
            mut = row.get("mut", "").strip()
            ddg = row.get("ddG", "").strip()
            if not mut or not ddg:
                continue
            try:
                out[mut] = {"flexddg_ddg": float(ddg)}
            except ValueError:
                log.warning(f"[Flex DDG] {files[0]} ddG 非数值: {ddg!r}")
    return out


FOLDX_FNAME_RE = re.compile(r"^Interface_([A-Za-z0-9]+)_([A-Za-z0-9]+)\.csv$")
SINGLE_MUT_RE = re.compile(r"^([A-Z])([A-Za-z0-9])(\d+)([A-Z])$")
PIPELINE_FIELDNAMES = [
    "mutation",
    "chain",
    "location",
    "original_aa",
    "mutant_aa",
    "sfe",
    "fep",
    "rosetta_flex",
    "foldx",
    "abnativ",
]


def parse_foldx(result_dir: Path, task: dict) -> dict[str, dict[str, float]]:
    """
    Mutation Energy of Binding (FoldX): 每对链一个 Interface_X_Y.csv，
    取每行的 Mutation (canonical 名) 和 deltaEnergy。
    """
    out: dict[str, dict[str, float]] = defaultdict(dict)
    found = 0
    for csv_path in result_dir.rglob("Interface_*.csv"):
        m = FOLDX_FNAME_RE.match(csv_path.name)
        if not m:
            continue
        c1, c2 = m.group(1), m.group(2)
        col = f"foldx_interface_{c1}_{c2}_ddg"
        found += 1
        with csv_path.open(encoding="utf-8", newline="") as f:
            for row in csv.DictReader(f):
                mut = row.get("Mutation", "").strip().strip('"')
                if not mut or mut == "WT":
                    continue
                val = row.get("deltaEnergy", "").strip()
                try:
                    out[mut][col] = float(val)
                except ValueError:
                    log.warning(f"[FoldX] {csv_path} 行 {mut} deltaEnergy 非数值: {val!r}")
    if found == 0:
        log.warning(f"[FoldX] {result_dir} 找不到任何 Interface_*.csv")
    return dict(out)


def parse_fep(result_dir: Path, task: dict) -> dict[str, dict[str, float]]:
    """Protein FEP: result.txt 三行 `<name>: <num>`。文件本身无突变名，靠 task['mutation'] 关联。"""
    files = list(result_dir.rglob("result.txt"))
    if not files:
        log.warning(f"[FEP] {result_dir} 找不到 result.txt")
        return {}
    mut = task.get("mutation")
    if not mut or mut == "*":
        log.warning(f"[FEP] task_name={task.get('task_name')} 缺少有效 mutation 字段")
        return {}
    parsed: dict[str, float] = {}
    for line in files[0].read_text(encoding="utf-8").splitlines():
        if ":" not in line:
            continue
        k, v = line.split(":", 1)
        try:
            parsed[k.strip()] = float(v.strip())
        except ValueError:
            continue
    cols: dict[str, float] = {}
    if "ligand dG" in parsed:
        cols["fep_ligand_dg"] = parsed["ligand dG"]
    if "complex dG" in parsed:
        cols["fep_complex_dg"] = parsed["complex dG"]
    if "final ddG" in parsed:
        cols["fep_ddg"] = parsed["final ddG"]
    return {mut: cols} if cols else {}


# ---------------------------------------------------------------------------
# AbNatiV 人源性打分（WeMol 模块，批量型：一个作业一份 FASTA → 多条序列打分）
# ---------------------------------------------------------------------------
# 总分列形如 "AbNatiV VH Score" / "AbNatiV VKappa Score"（中间是单 token，无连字符）
# CDR/FR 子分列如 "AbNatiV CDR1-VH Score" / "AbNatiV FR-VH Score" 含连字符，排除
_ABNATIV_OVERALL_RE = re.compile(r"^AbNatiV V\w+ Score$")


def _parse_abnativ_seq_csv(csv_path: Path) -> dict[str, float]:
    """读一个 *_abnativ_seq_scores.csv，返回 {seq_id: 总分(4位小数)}，跳过 WT。"""
    out: dict[str, float] = {}
    with csv_path.open(encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        score_col = None
        for c in (reader.fieldnames or []):
            if _ABNATIV_OVERALL_RE.match((c or "").strip()):
                score_col = c
                break
        if not score_col:
            log.warning(f"[AbNatiV] {csv_path} 找不到总分列(AbNatiV V* Score)")
            return out
        for row in reader:
            seq_id = (row.get("seq_id") or "").strip()
            if not seq_id or seq_id == "WT":
                continue
            val = (row.get(score_col) or "").strip()
            try:
                out[seq_id] = round(float(val), 4)
            except ValueError:
                log.warning(f"[AbNatiV] {csv_path} 行 {seq_id} 分值非数值: {val!r}")
    return out


def parse_abnativ(result_dir: Path, task: dict) -> dict[str, dict[str, float]]:
    """
    AbnatiV: 结果目录里有按序列总分文件（一个作业一份 FASTA 的多条序列打分，批量型）。
    取每行的总分（AbNatiV V* Score），4 位小数，作为 Humanness 列。seq_id 列即突变名。

    文件名兼容两种：WeMol 模块输出 `H_seq_scores.csv` / `L_seq_scores.csv`，
    命令行 -oid 方式输出 `*_abnativ_seq_scores.csv`。统一用 `*_seq_scores.csv` 匹配
    （不会误匹配按残基的 `*_res_scores.csv`）。
    """
    out: dict[str, dict[str, float]] = defaultdict(dict)
    found = 0
    for csv_path in result_dir.rglob("*_seq_scores.csv"):
        found += 1
        for mut, score in _parse_abnativ_seq_csv(csv_path).items():
            out[mut]["Humanness"] = score
    if found == 0:
        log.warning(f"[AbnatiV] {result_dir} 找不到 *_seq_scores.csv")
    return dict(out)


# 模块名 → 解析器。键须与 WeMol 平台模块真名一致。
PARSERS: dict[str, Callable[[Path, dict], dict[str, dict[str, float]]]] = {
    "Flex DDG":                    parse_flexddg,
    "Mutation Energy of Binding":  parse_foldx,
    "Protein FEP":                 parse_fep,
    "AbnatiV":                     parse_abnativ,
}


def aggregate_to_csv(
    records: list[dict],
    result_dir_base: Path,
    output_csv: Path,
    include_multipoint: bool = False,
    pipeline_format: bool = True,
    foldx_interface: str | None = None,
) -> dict[str, int]:
    """
    扫所有 status=Done 的记录，调用各模块 parser，拼成总表写到 output_csv。
    返回统计 dict: {muts_parsed, columns, parsers_called, parser_missing}.
    """
    aggregated: dict[str, dict[str, float]] = defaultdict(dict)
    n_called = n_missing_parser = n_no_dir = 0

    for rec in records:
        if rec.get("status") != "Done":
            continue
        module = rec.get("module_name") or rec.get("flow_name")
        parser = PARSERS.get(module)
        if not parser:
            n_missing_parser += 1
            log.warning(f"⤳ 未注册解析器: {module} (task={rec.get('task_name')})")
            continue

        job_id = rec.get("job_id")
        rdir = find_result_dir(result_dir_base, job_id)
        if not rdir:
            n_no_dir += 1
            log.warning(f"✗ job_id={job_id} 找不到结果目录在 {result_dir_base}")
            continue

        try:
            partial = parser(rdir, rec)
        except Exception as e:
            log.error(f"✗ 解析 {module} (job {job_id}) 出错: {e}")
            continue
        n_called += 1
        for mut, cols in partial.items():
            aggregated[mut].update(cols)

    # 过滤多点（含逗号 / 下划线分隔）
    if not include_multipoint:
        aggregated = {m: v for m, v in aggregated.items() if "," not in m and "_" not in m}

    # 列名按模块前缀排序：flexddg → foldx → fep → 其他
    priority = ("flexddg_", "foldx_", "fep_")
    all_cols: set[str] = set()
    for v in aggregated.values():
        all_cols.update(v.keys())

    def _key(c: str) -> tuple[int, str]:
        for i, pre in enumerate(priority):
            if c.startswith(pre):
                return (i, c)
        return (len(priority), c)

    sorted_muts = sorted(aggregated.keys())

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    if pipeline_format:
        rows, n_foldx_ambiguous = _to_pipeline_rows(
            sorted_muts,
            aggregated,
            foldx_interface=foldx_interface,
        )
        with output_csv.open("w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=PIPELINE_FIELDNAMES)
            w.writeheader()
            w.writerows(rows)
        n_columns = len(PIPELINE_FIELDNAMES) - 1
    else:
        sorted_cols = sorted(all_cols, key=_key)
        with output_csv.open("w", encoding="utf-8", newline="") as f:
            w = csv.writer(f)
            w.writerow(["mutation"] + sorted_cols)
            for mut in sorted_muts:
                w.writerow([mut] + [aggregated[mut].get(c, "") for c in sorted_cols])
        n_columns = len(sorted_cols)
        n_foldx_ambiguous = 0

    return {
        "muts_parsed":    len(aggregated),
        "columns":        n_columns,
        "parsers_called": n_called,
        "parser_missing": n_missing_parser,
        "no_result_dir":  n_no_dir,
        "foldx_ambiguous": n_foldx_ambiguous,
    }


def _to_pipeline_rows(
    sorted_muts: list[str],
    aggregated: dict[str, dict[str, float]],
    foldx_interface: str | None = None,
) -> tuple[list[dict[str, object]], int]:
    rows: list[dict[str, object]] = []
    n_foldx_ambiguous = 0
    requested_foldx_cols = _foldx_interface_to_cols(foldx_interface)

    for mut in sorted_muts:
        cols = aggregated[mut]
        parsed = _parse_single_mutation(mut)
        foldx_value, ambiguous = _select_foldx_value(cols, requested_foldx_cols)
        if ambiguous:
            n_foldx_ambiguous += 1
        rows.append({
            "mutation": mut,
            "chain": parsed["chain"] if parsed else "",
            "location": parsed["location"] if parsed else "",
            "original_aa": parsed["original_aa"] if parsed else "",
            "mutant_aa": parsed["mutant_aa"] if parsed else "",
            "sfe": "",
            "fep": cols.get("fep_ddg", ""),
            "rosetta_flex": cols.get("flexddg_ddg", ""),
            "foldx": foldx_value,
            "abnativ": cols.get("Humanness", ""),
        })

    if n_foldx_ambiguous:
        log.warning(
            "FoldX has multiple interface columns for %s mutations; "
            "leaving foldx blank. Re-run with --foldx-interface such as M_N.",
            n_foldx_ambiguous,
        )
    return rows, n_foldx_ambiguous


def _parse_single_mutation(mut: str) -> dict[str, object] | None:
    m = SINGLE_MUT_RE.match(mut.strip())
    if not m:
        return None
    return {
        "original_aa": m.group(1),
        "chain": m.group(2),
        "location": int(m.group(3)),
        "mutant_aa": m.group(4),
    }


def _foldx_interface_to_cols(foldx_interface: str | None) -> list[str] | None:
    """
    解析 --foldx-interface。支持单个界面（如 M_S）或 '+' 连接的多个界面求和
    （如 M_S+N_S = 抗体两条链各自与抗原的界面之和）。返回列名列表，None 表示未指定。
    """
    if not foldx_interface:
        return None
    out: list[str] = []
    for part in foldx_interface.split("+"):
        normalized = part.strip().replace("-", "_")
        if not normalized:
            continue
        if normalized.startswith("foldx_interface_"):
            col = normalized if normalized.endswith("_ddg") else f"{normalized}_ddg"
        else:
            col = f"foldx_interface_{normalized}_ddg"
        out.append(col)
    return out or None


def _select_foldx_value(
    cols: dict[str, float],
    requested_cols: list[str] | None,
) -> tuple[float | str, bool]:
    # 指定了界面：取这些界面的和（缺失的列按缺省跳过；全缺则留空）
    if requested_cols:
        vals = [cols[c] for c in requested_cols if c in cols]
        if not vals:
            return "", False
        return sum(vals), False

    foldx_cols = sorted(c for c in cols if c.startswith("foldx_interface_"))
    if not foldx_cols:
        return "", False
    if len(foldx_cols) == 1:
        return cols.get(foldx_cols[0], ""), False
    return "", True


# ---------------------------------------------------------------------------
# 日志
# ---------------------------------------------------------------------------
def setup_console_logging(level: int = logging.INFO) -> None:
    """让脚本的 INFO 日志能打到控制台 (SDK 自己的 file handler 会另写文件)。"""
    root = logging.getLogger()
    if any(isinstance(h, logging.StreamHandler) for h in root.handlers):
        return
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s",
                                            datefmt="%H:%M:%S"))
    root.addHandler(handler)
    root.setLevel(level)
