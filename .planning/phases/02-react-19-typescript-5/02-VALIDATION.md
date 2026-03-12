---
phase: 2
slug: react-19-typescript-5
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-03-12
---

# Phase 2 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | Vitest 4.1.0 |
| **Config file** | `vite.config.ts` (test section) |
| **Quick run command** | `cd frontend && npm test` |
| **Full suite command** | `cd frontend && npx tsc --noEmit && npm test && npm run build` |
| **Estimated runtime** | ~30 seconds |

---

## Sampling Rate

- **After every task commit:** Run `cd frontend && npx tsc --noEmit && npm test`
- **After every plan wave:** Run `cd frontend && npx tsc --noEmit && npm test && npm run build && npm run lint`
- **Before `/gsd:verify-work`:** Full suite must be green + visual verification of all 9 pages + DevTools console clean
- **Max feedback latency:** 30 seconds

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|-----------|-------------------|-------------|--------|
| 02-01-01 | 01 | 0 | FRAME-03 | build | `cd frontend && npx tsc --noEmit` | No (W0 creates tsconfig) | ⬜ pending |
| 02-01-02 | 01 | 1 | FRAME-01 | smoke | `cd frontend && npm test` | Yes (App.test.tsx) | ⬜ pending |
| 02-01-03 | 01 | 1 | FRAME-01 | build | `cd frontend && npm run build` | N/A | ⬜ pending |
| 02-02-01 | 02 | 1 | FRAME-03 | lint | `cd frontend && npm run lint` | N/A | ⬜ pending |
| 02-02-02 | 02 | 1 | FRAME-05 | build | `cd frontend && npm run build` | N/A | ⬜ pending |
| 02-02-03 | 02 | 1 | FRAME-05 | manual | Visual verification of page transitions, layout animations | N/A | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

- [ ] `tsconfig.json` + `tsconfig.app.json` + `tsconfig.node.json` — must be created before any TS compilation
- [ ] `src/types/global.d.ts` — must exist for 3Dmol, Vite env vars, CSS module imports
- [ ] `@types/react@^19` + `@types/react-dom@^19` — must be installed before TS compilation
- [ ] Update `package.json` scripts: lint/format globs from `.js/.jsx` to `.ts/.tsx`

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| No console errors/warnings on all 9 pages | FRAME-01 | Runtime console output cannot be captured in unit tests | Open each page in browser, check DevTools console for errors/warnings |
| Motion animations work correctly | FRAME-05 | Visual behavior cannot be asserted programmatically | Navigate between pages, verify transitions are smooth; check layout animations on tab-based pages |

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 30s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
