# CLAUDE.md

> **NOTE:** This is version 2 — a reconstructed stand-in for the original .MD file,
> which is temporarily unavailable (it lives on a local machine in Calgary).
> Use this file as the active instructions until the original is restored.
> When the original is recovered, compare the two and merge anything missing.

## 1. Skills — always use them (non-negotiable)

Before starting **any** task, check the available skills and invoke every one that
matches the work. Do not do a task inline when a skill exists for it.

Big notes:

- **Always scan the skills list first.** If a skill matches the request — even
  partially — invoke it before writing any other response about the task.
- **Code changes:** run `verify` before committing nontrivial changes, and
  `code-review` on the diff. Use `simplify` to clean up changed code.
- **Reviews:** use `review` for GitHub pull requests, `security-review` for
  security-sensitive changes on the current branch.
- **New projects/repos:** use `init` to generate codebase documentation.
- **Running the app:** use `run` when launching or screenshotting the project.
- **Anything LLM/Claude-API related:** read `claude-api` before answering — never
  answer model/pricing/API questions from memory.
- **Charts and visualization:** read `dataviz` before writing any chart code.
- **Recurring tasks:** use `loop`; **harness/settings changes:** use `update-config`.

If a skill exists and wasn't used, that is a mistake — go back and use it.

## 2. Agents — use them no matter what (non-negotiable)

Delegating to agents is mandatory, not optional. Do not keep work inline that an
agent can carry.

- **Use agents for every task where they are available** — research, codebase
  search, planning, multi-step implementation. This is non-negotiable.
- **Explore agent:** any broad codebase search or fan-out across many files.
- **Plan agent:** designing the implementation strategy before nontrivial changes.
- **general-purpose agent:** complex research and multi-step task execution.
- Continue an existing agent with its context (SendMessage) rather than
  re-spawning cold when following up on the same thread of work.
- If in doubt whether a task warrants an agent: it does. Spawn it.

## 3. General working notes

- Commit and push work to the designated feature branch; never push to a
  different branch without explicit permission.
- Keep this file up to date: when new standing instructions come up in
  conversation, add them here so they survive across sessions.
