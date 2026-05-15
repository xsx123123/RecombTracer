# RecombTracer

**版本：** 0.0.1  
[English](../README.md) | 中文

RecombTracer 是一个用于识别 MAGIC（多亲本高级世代互交）群体中的重组断点和局部祖先的工具包。它结合 PBWT 风格的染色体涂染与隐马尔可夫模型（HMM），以推断亲本单倍型片段和重组事件。

## 概述

该流程由两个核心模块组成：

- **`recombtracer.core.recombiner`** — PBWT 风格单倍型涂染，基于最大匹配片段将每个子代的 SNP 分配到最可能的亲本来源。
- **`recombtracer.core.hmm`** — HMM 细化，将嘈杂的 PBWT 分配平滑为连续的祖先区块，并以置信度评分调用重组断点。

此外还提供了一个实用模块，用于将标准 VCF 文件转换为分析流程所需的基于 NumPy 的输入格式。

---

## 依赖

- Python ^3.10
- numpy 2.2.6
- pandas 2.3.3
- scipy 1.15.3
- cyvcf2 0.33.0（用于 VCF 输入/输出）
- pyyaml >=6.0
- rich 15.0.0
- rich-gradient 0.3.12
- rich-argparse 1.8.0

通过 pip 安装（从源码可编辑安装）：

```bash
pip install -e .
```

或使用 Poetry：

```bash
poetry install
```

---

## VCF 转换（`recombtracer convert-vcf`）

如果你的数据从 VCF 文件开始（例如变异检测和过滤后），请使用 `recombtracer convert-vcf` 提取单倍型并将其转换为 `MagicRecombiner` 使用的矩阵格式。

### 输入要求

- VCF 需要建立索引（`.tbi` 或 `.csi`）。
- 所有位点应为纯合子（例如 `0/0` 或 `1/1`）。杂合位点应事先移除（例如 `bcftools view -e GT="het"`）。
- 缺失基因型（`./.`）默认会被丢弃，以避免分析错误。

### 命令行用法

一次转换一条染色体为压缩的 `.npz` 归档：

```bash
recombtracer convert-vcf test/subset_remove_het.vcf.gz \
    --parents 1,2,3,4,5,6,7,8,9,10,11 \
    --progeny 98,99,100,144,145,146 \
    --chrom LG1 \
    --out-dir ./magic_input
```

| 参数 | 描述 |
|------|------|
| `--parents` | **必需。** 逗号分隔的 founders/亲本样本名称。 |
| `--progeny` | 逗号分隔的子代样本名称。如果省略，自动使用 VCF 中所有非亲本样本。 |
| `--chrom` | **必需。** 要提取的染色体/ contig（例如 `LG1`）。 |
| `--out-dir` | `.npz` 文件的输出目录（默认：当前目录）。 |
| `--keep-missing` | 保留在某些样本中有缺失基因型的 SNP（默认：丢弃）。 |

输出文件（`{chrom}_magic.npz`）包含：

| 键 | 形状 / 类型 | 描述 |
|-----|------------|------|
| `parent_haps` | `(n_parents, n_snps)` int8 | 亲本单倍型矩阵。 |
| `parent_names` | str 列表 | 亲本样本名称。 |
| `progeny_{name}` | `(1, n_snps)` int8 | 每个子代个体的单倍型（纯合位点仅存储一个等位基因）。 |
| `progeny_names` | str 列表 | 子代样本名称。 |
| `positions` | `(n_snps,)` int32 | 基于1的基因组坐标。 |
| `chrom` | str | 染色体名称。 |

### Python API

```python
from recombtracer import vcf_to_magic_inputs, load_chromosome_npz

# 直接从 VCF 转换
data = vcf_to_magic_inputs(
    vcf_path="test/subset_remove_het.vcf.gz",
    parent_samples=["1", "2", "3", "4"],
    progeny_samples=["5", "6", "7", "8"],
    chrom="LG0",
)

# 或加载之前保存的 .npz 文件
# data = load_chromosome_npz("magic_input/LG0_magic.npz")

print(data["parent_haps"].shape)      # (n_parents, n_snps)
print(data["progeny_haps"]["5"].shape) # (1, n_snps)
```

---

## 重组分析 CLI（`recombtracer run`）

将 VCF 转换为 `.npz` 后，从命令行运行完整的 PBWT + HMM 流程：

```bash
recombtracer run magic_input/LG1_magic.npz \
    --out-dir ./results \
    --smooth-window 5 \
    --save-raw
```

| 参数 | 描述 |
|------|------|
| `npz` | **必需。** `convert-vcf` 生成的输入 `.npz` 文件。 |
| `--out-dir` | 结果 CSV 的输出目录（默认：当前目录）。 |
| `--smooth-window` | PBWT 中值滤波窗口大小（默认：5）。 |
| `--min-segment-snps` | 每个祖先片段的最小 SNP 数（默认：5）。 |
| `--min-segment-bp` | 每个祖先片段的最小碱基对长度（默认：1000）。 |
| `--min-posterior` | 保留断点的最小 HMM 后验概率（默认：0.8）。 |
| `--save-raw` | 同时保存原始 PBWT 结果与 HMM 细化结果。 |
| `--progeny` | 要分析的逗号分隔子代列表（默认：`.npz` 中所有子代）。 |
| `--haplotype` | 仅分析特定单倍型索引（默认：全部）。 |

### 输出文件

对每个子代 × 单倍型，你将获得：

| 文件 | 描述 |
|------|------|
| `{name}_hap{h}_{chrom}_hmm_viterbi.csv` | 每个 SNP 的 Viterbi 亲本分配和后验概率。 |
| `{name}_hap{h}_{chrom}_hmm_segments.csv` | HMM 平滑后的连续祖先片段。 |
| `{name}_hap{h}_{chrom}_hmm_recombinations.csv` | 带置信度评分的 HMM 过滤重组断点。 |
| `{name}_hap{h}_{chrom}_paint.csv` | *（带 `--save-raw`）* 原始 PBWT 每个 SNP 的亲本分配。 |
| `{name}_hap{h}_{chrom}_segments_raw.csv` | *（带 `--save-raw`）* 原始 PBWT 祖先片段。 |
| `{name}_hap{h}_{chrom}_recombinations_raw.csv` | *（带 `--save-raw`）* 原始 PBWT 断点调用。 |
| `summary_{chrom}.csv` | 汇总表：每个个体的断点和片段计数。 |

---

## 完整流水线 CLI（`recombtracer pipeline`）

通过一条命令运行完整的工作流——VCF 转换后紧跟 PBWT + HMM 分析：

```bash
recombtracer pipeline test/subset_remove_het.vcf.gz \
    --parents 1,2,3,4,5,6,7,8,9,10,11 \
    --chrom LG1 \
    --out-dir ./results \
    --save-raw
```

这等价于使用相同参数依次运行 `convert-vcf` 和 `run`。`run` 中的所有算法和过滤参数在此同样可用。

---

## 重组分析（`recombtracer.core.recombiner`）

获得单倍型矩阵后，运行 PBWT 涂染：

```python
from recombtracer import MagicRecombiner

recombiner = MagicRecombiner(
    parent_haps=data["parent_haps"],
    parent_names=data["parent_names"],
    positions=data["positions"],
    chrom=data["chrom"],
)

# 分析一个子代
prog_name = data["progeny_names"][0]
prog_hap = data["progeny_haps"][prog_name]  # 形状 (1, n_snps)

paint_df = recombiner.paint_progeny(
    progeny_haps=prog_hap,
    progeny_name=prog_name,
    smooth_window=5,
)

# 提取连续祖先片段
segments = recombiner.extract_segments(paint_df, min_segment_snps=5)

# 调用重组断点
rec_df = recombiner.call_recombinations(segments)
```

### 输出列

`paint_df`（每个 SNP 的亲本分配）：
- `chrom`, `position`, `progeny`, `haplotype`
- `parent` — 推断的亲本来源
- `confidence` — PBWT 匹配置信度

`segments`（连续区块）：
- `chrom`, `start_pos`, `end_pos`, `parent`, `haplotype`
- `score` — 片段平均置信度

`rec_df`（断点调用）：
- `chrom`, `position`, `haplotype`
- `left_parent`, `right_parent`
- `confidence` — 两侧片段的平均评分

---

## HMM 细化（`recombtracer.core.hmm`）

为了平滑嘈杂的 PBWT 调用并获得后验概率，运行 HMM：

```python
from recombtracer import run_hmm_refinement

viterbi_df, segments_df, rec_df = run_hmm_refinement(
    paint_df=paint_df[paint_df["haplotype"] == 0],
    parent_haps=data["parent_haps"],
    parent_names=data["parent_names"],
    progeny_hap=prog_hap[0],  # 展平为 (n_snps,)
    progeny_name=prog_name,
    chrom=data["chrom"],
)
```

- `viterbi_df` — 每个 SNP 最可能的亲本和后验概率。
- `segments_df` — 平滑后的祖先片段。
- `rec_df` — HMM 过滤的重组断点。

---

## 完整工作流示例

### 快速 CLI 流水线（一键运行）

```bash
# 完整流水线：VCF → .npz → PBWT + HMM
recombtracer pipeline test/subset_remove_het.vcf.gz \
    --parents 1,2,3,4,5,6,7,8,9,10,11 \
    --chrom LG1 \
    --out-dir ./results \
    --save-raw
```

### 分步 CLI

```bash
# 1. 将 VCF 转换为 NumPy 输入
recombtracer convert-vcf test/subset_remove_het.vcf.gz \
    --parents 1,2,3,4,5,6,7,8,9,10,11 \
    --chrom LG1 \
    --out-dir ./magic_input

# 2. 运行 PBWT + HMM 分析
recombtracer run magic_input/LG1_magic.npz \
    --out-dir ./results \
    --save-raw
```

### Python API（用于自定义分析）

如果你需要对参数或下游分析进行程序化控制：

```python
from recombtracer import (
    load_chromosome_npz,
    MagicRecombiner,
    run_hmm_refinement,
)

data = load_chromosome_npz("magic_input/LG1_magic.npz")

recombiner = MagicRecombiner(
    parent_haps=data["parent_haps"],
    parent_names=data["parent_names"],
    positions=data["positions"],
    chrom=data["chrom"],
)

for prog_name in data["progeny_names"]:
    prog_hap = data["progeny_haps"][prog_name]

    paint_df = recombiner.paint_progeny(prog_hap, progeny_name=prog_name)
    segments = recombiner.extract_segments(paint_df)
    rec_raw = recombiner.call_recombinations(segments)

    viterbi_df, seg_hmm, rec_hmm = run_hmm_refinement(
        paint_df[paint_df["haplotype"] == 0],
        data["parent_haps"],
        data["parent_names"],
        prog_hap[0],
        progeny_name=prog_name,
        chrom=data["chrom"],
    )

    print(f"{prog_name}: {len(rec_raw)} 原始断点, {len(rec_hmm)} HMM 断点")
```

---

## 项目结构

```
.
├── recombtracer/              # 主程序包
│   ├── __init__.py            # 公共 API 导出
│   ├── cli.py                 # 统一命令行接口
│   ├── core/                  # 核心分析模块
│   │   ├── recombiner.py      # PBWT 染色体涂染
│   │   ├── hmm.py             # HMM 平滑与断点调用
│   │   ├── vcf.py             # VCF → NumPy 转换工具
│   │   ├── convert.py         # convert-vcf 的 CLI 处理
│   │   ├── run.py             # run 的 CLI 处理
│   │   └── pipeline.py        # pipeline 的 CLI 处理
│   ├── config/                # 程序包配置
│   │   ├── software.yaml      # 软件元数据
│   │   └── default.yaml       # 默认分析参数
│   └── utils/                 # 工具模块（日志、Logo 等）
├── test/
│   └── subset_remove_het.vcf.gz   # 示例 VCF（已移除杂合位点）
├── pyproject.toml             # 基于 Poetry 的构建配置
└── README.md
```

---

## 许可证

MIT
