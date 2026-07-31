# CLAUDE.md

## Purpose

This repo holds personal Claude Code configuration used from mobile/web
sessions (Claude Code on the web). It is intentionally separate from the
local Claude config on the main machine — nothing here should assume or
depend on that machine's setup. There's no app code here, just config.

## Working notes

- Develop on a `claude/*` feature branch and push there; never push to a
  different branch without explicit permission.
- Prefer small, direct changes — this repo is config, not a product
  codebase, so keep additions (settings, skills, hooks) minimal and easy
  to reason about from a phone.
- `.claude/settings.json` carries a small read-only permission allowlist
  (git status/diff/log/show/branch/remote, ls, pwd, Read/Glob/Grep) to
  cut down on permission prompts during mobile sessions. Expand it
  deliberately, not by default — nothing destructive or network-facing
  should be auto-approved.
- When a new standing preference or instruction comes up in conversation,
  add it here so it survives across sessions.
