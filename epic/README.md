# Epic backlog（规划任务清单）

本目录用于维护 SHUD-up 的 **Epic backlog**：把 `docs/spec.md` 中的 To-Be 目标拆分为可执行的 issue，并持续追踪进度。

- 计划基线（Spec）：[`docs/spec.md`](../docs/spec.md)
- 当前实现（As-Is）：[`docs/model_upgrade.md`](../docs/model_upgrade.md)

## Backlog 维护方式

目前采用“以 Markdown 文件为单元”的方式维护 backlog，并可用脚本批量创建 GitHub issues：

- Backlog 文件命名：`epic/backlog-*.md`
  - 文件首行 `# ...` 作为 issue 标题
  - 其余内容作为 issue body
- 创建 issues 脚本：`epic/create_issues_from_backlog.sh`
  - 依赖：`gh`（GitHub CLI）已登录、`git`
  - 会为每个 `backlog-*.md` 创建一个 issue
  - 会在 `epic/issue-map.tsv` 记录 “文件路径 -> issue URL -> issue number” 的映射（自动生成）

运行示例（在仓库根目录）：

```bash
bash epic/create_issues_from_backlog.sh
```

## 注意事项

- 本目录下的 `issue-map.tsv` 为脚本生成产物；是否提交到仓库由团队约定决定。
- 若仓库中暂时没有 `backlog-*.md` 文件：表示 backlog 可能仍在 GitHub issues/看板中维护，后续再逐步回填到本目录以实现“文档即 backlog”。
