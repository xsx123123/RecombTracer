# RecombTracer 优化 TODO

> 生成时间：2026-05-20

---

## 🔴 高优先级（建议尽快修）

### 1. `__all__` 与实际导入名不匹配
- **文件**：`recombtracer/__init__.py`
- **问题**：`__all__` 第 56 行写的是 `"extract_chromosome"`，但实际导入的是 `extract_homozygous_chromosome`
- **影响**：`from recombtracer import *` 行为异常
- **修复**：统一名称，或修改 `__all__` 中的字符串

### 2. 测试架构几乎是空的
- **文件/目录**：`test_pbwt_simple.py`、`test/`
- **问题**：
  - `test_pbwt_simple.py` 放在根目录，不是 pytest 发现的格式
  - `test/` 目录里有数据生成脚本和大量 CSV，但没有真正的 `test_*.py` 测试文件
  - `test/` 应该改名叫 `tests/`（pytest 推荐）
- **修复**：
  - [ ] 将 `test/` 重命名为 `tests/`
  - [ ] CSV 测试数据移到 `tests/fixtures/` 或 `tests/data/`
  - [ ] 把 `test_pbwt_simple.py` 改为 pytest 风格（`test_` 前缀 + `assert`）
  - [ ] 为 `generate_report_test_data.py` 补充对应的断言测试

### 3. `config/` 配置重复
- **文件**：`recombtracer/config/default.yaml`、`recombtracer/config/software.yaml`
- **问题**：`default.yaml` 里已经包含了 `software:` 段，但同时又单独有一个 `software.yaml`，两者内容高度重叠
- **修复**：合并为一个配置文件，或明确分工（如 `software.yaml` 只存元数据，`default.yaml` 只存运行参数）

### 4. `pyproject.toml` 格式问题
- **文件**：`pyproject.toml`
- **问题**：`authors = ["JZHANG zhangjian199567@outlook.com"]` 不是 Poetry 推荐格式
- **修复**：改为标准格式：
  ```toml
  authors = ["JZHANG <zhangjian199567@outlook.com>"]
  ```

---

## 🟡 中优先级（代码整洁 & 可维护性）

### 5. `utils/logo.py` 代码清理
- **文件**：`recombtracer/utils/logo.py`
- **问题**：
  - 大量被注释掉的 `ascii_type_2/3/5/8/10-16`，git 历史里都有，可以删掉
  - 第 158 行 `ascii_type_10` 里有 `_123123` 残留字符
  - `rice_color` 参数名疑似 typo，应该是 `rich_color`？
  - `GradientText` 导入后从未使用
- **修复**：清理无用代码和残留字符，修正参数名

### 6. `core/run.py` 里的 hardcode 逻辑
- **文件**：`recombtracer/core/run.py`
- **问题**：
  ```python
  paint_hap = paint_df[paint_df["haplotype"] == 0].copy()
  paint_hap["haplotype"] = h
  ```
  这里依赖了 `paint_progeny` 传入单条 haplotype 时内部编号为 0 的假设
- **修复**：加断言保护，或让 `paint_progeny` 支持传入单 haplotype 时保留原索引

### 7. `report/generator.py` 的重复解析逻辑
- **文件**：`recombtracer/report/generator.py`
- **问题**：`rec_files` 和 `seg_files` 的文件名解析代码几乎一模一样
- **修复**：提取成 `_parse_filename(basename, suffix)` 辅助函数

### 8. 空的 `config/__init__.py`
- **文件**：`recombtracer/config/__init__.py`
- **问题**：完全为空文件
- **修复**：添加 docstring，或暴露便捷导入

---

## 🟢 低优先级（锦上添花）

### 9. `recombtracer/__init__.py` 的导入策略
- **文件**：`recombtracer/__init__.py`
- **问题**：`import recombtracer` 会立即加载 `core.recombiner`、`core.hmm`、`core.vcf` 等所有重模块，即使只想看个版本号也会拖慢启动
- **修复**：考虑对重模块使用延迟导入（lazy import），或把 `__version__` 的加载和重模块解耦

### 10. `.gitignore` 完善
- **问题**：可能遗漏了一些常见忽略项
- **建议补充**：
  ```gitignore
  __pycache__/
  *.egg-info/
  dist/
  build/
  *.html        # 如果 test_report.html 是生成的
  .pytest_cache/
  ```

### 11. `README.md` 与 `docs/README_zh.md` 同步
- **文件**：`README.md`、`docs/README_zh.md`
- **问题**：英文 README 里的 `loguru 0.7.2` 也需要同步改成 `>=0.7.2`
- **修复**：检查并同步中英文 README 的依赖版本、安装说明等内容
