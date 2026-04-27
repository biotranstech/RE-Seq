# RE-Seq

RE-Seq is a repository for sequencing analysis workflow organization, historical pipeline preservation, and Python package refactoring.

当前仓库已经完成第一阶段结构整理，分为两个核心部分：

- `workflow/`：保存历史流程、原始脚本和旧版分析逻辑
- `python_package/`：保存重构和标准化的 Python 包版本

本仓库的目标不是简单保存脚本，而是将原有流程逐步整理为一个**可维护、可安装、可扩展、可复用**的标准化分析工具仓库。

---

## 1. Repository Overview

很多生物信息分析项目在早期开发时，通常采用“脚本 + 手工目录管理 + 流程拼接”的方式完成。  
这种方式虽然能快速实现分析目标，但随着项目规模变大，往往会出现以下问题：

- 目录结构不清晰
- 配置和路径耦合严重
- 不便迁移和复用
- 不利于版本管理
- 难以标准化部署
- 新老流程混杂，维护成本高

因此，本仓库进行了结构重组，将内容拆分为两个层级：

1. **历史流程层（workflow）**
2. **标准化开发层（python_package）**

这种组织方式既能保留原始流程，又能为后续 Python 包开发和正式发布打下基础。

---

## 2. Current Repository Structure

当前仓库顶层结构如下：

```text
RE-Seq/
├── workflow/
├── python_package/
└── README.md
