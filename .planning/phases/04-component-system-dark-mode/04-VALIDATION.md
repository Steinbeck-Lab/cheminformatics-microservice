---
phase: 4
slug: component-system-dark-mode
status: draft
nyquist_compliant: true
wave_0_complete: false
created: 2026-03-12
---

# Phase 4 -- Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | Vitest 4.1.0 + Testing Library |
| **Config file** | vite.config.mts (test section) |
| **Quick run command** | `cd frontend && npm test` |
| **Full suite command** | `cd frontend && npm run test:coverage` |
| **Estimated runtime** | ~15 seconds |

---

## Sampling Rate

- **After every task commit:** Run `cd frontend && npm test`
- **After every plan wave:** Run `cd frontend && npm run test:coverage`
- **Before `/gsd:verify-work`:** Full suite must be green
- **Max feedback latency:** 15 seconds

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|-----------|-------------------|-------------|--------|
| 04-01-00 | 01 | 1 | ALL | scaffold | `cd frontend && npx vitest run src/__tests__/ --reporter=verbose` | Created by Task 0 | ⬜ pending |
| 04-01-01 | 01 | 1 | COMP-01 | unit | `cd frontend && npx vitest run src/__tests__/utils.test.ts -x` | Created by 04-01 Task 0, filled by Task 1 | ⬜ pending |
| 04-01-02 | 01 | 1 | COMP-02 | unit | `cd frontend && npx vitest run src/__tests__/components/button.test.tsx -x` | Created by 04-01 Task 0, filled by Task 1 | ⬜ pending |
| 04-01-03 | 01 | 1 | COMP-03 | unit | `cd frontend && npx vitest run src/__tests__/components/card.test.tsx -x` | Created by 04-01 Task 0, filled by Task 1 | ⬜ pending |
| 04-01-04 | 01 | 1 | COMP-04 | unit | `cd frontend && npx vitest run src/__tests__/components/dialog.test.tsx -x` | Created by 04-01 Task 0 | ⬜ pending |
| 04-01-05 | 01 | 1 | COMP-05 | unit | `cd frontend && npx vitest run src/__tests__/components/input.test.tsx -x` | Created by 04-01 Task 0, filled by Task 1 | ⬜ pending |
| 04-02-01 | 02 | 1 | COMP-06 | smoke | `cd frontend && npx vitest run src/__tests__/icon-migration.test.ts -x` | Created by 04-01 Task 0, filled by 04-02 Task 2 | ⬜ pending |
| 04-02-02 | 02 | 1 | COMP-07 | smoke | `cd frontend && npx vitest run src/__tests__/App.test.tsx -x` | ✅ | ⬜ pending |
| 04-03-01 | 03 | 2 | THEME-01 | unit | `cd frontend && npx vitest run src/__tests__/context/theme.test.tsx -x` | Created by 04-01 Task 0, filled by 04-01 Task 2 | ⬜ pending |
| 04-03-02 | 03 | 2 | THEME-02 | unit | `cd frontend && npx vitest run src/__tests__/context/theme.test.tsx -x` | Created by 04-01 Task 0, filled by 04-01 Task 2 | ⬜ pending |
| 04-03-03 | 03 | 2 | THEME-03 | unit | `cd frontend && npx vitest run src/__tests__/context/theme.test.tsx -x` | Created by 04-01 Task 0, filled by 04-01 Task 2 | ⬜ pending |
| 04-04-01 | 04 | 2 | THEME-04 | manual | Visual inspection | N/A | ⬜ pending |
| 04-04-02 | 04 | 2 | THEME-05 | manual | Visual inspection | N/A | ⬜ pending |

*Status: ⬜ pending / ✅ green / ❌ red / ⚠️ flaky*

---

## Wave 0 Requirements

All Wave 0 test stubs are created by **Plan 04-01, Task 0** (first task executed in the phase):

- [ ] `src/__tests__/utils.test.ts` -- stubs for COMP-01 (cn() utility), filled in by 04-01 Task 1
- [ ] `src/__tests__/components/button.test.tsx` -- stubs for COMP-02, filled in by 04-01 Task 1
- [ ] `src/__tests__/components/card.test.tsx` -- stubs for COMP-03, filled in by 04-01 Task 1
- [ ] `src/__tests__/components/dialog.test.tsx` -- stubs for COMP-04, remains as todo until 04-04
- [ ] `src/__tests__/components/input.test.tsx` -- stubs for COMP-05, filled in by 04-01 Task 1
- [ ] `src/__tests__/icon-migration.test.ts` -- stubs for COMP-06, filled in by 04-02 Task 2
- [ ] `src/__tests__/context/theme.test.tsx` -- stubs for THEME-01/02/03, filled in by 04-01 Task 2

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| CSS variables change between light/dark | THEME-04 | CSS variables applied by browser runtime | Toggle dark mode, inspect computed styles on body/root |
| Pages render correctly in both modes | THEME-05 | Visual layout correctness requires human judgment | Visit all 9 pages in light and dark mode, check for contrast/readability |

---

## Validation Sign-Off

- [x] All tasks have `<automated>` verify or Wave 0 dependencies
- [x] Sampling continuity: no 3 consecutive tasks without automated verify
- [x] Wave 0 covers all MISSING references (created by Plan 04-01 Task 0)
- [x] No watch-mode flags
- [x] Feedback latency < 15s
- [x] `nyquist_compliant: true` set in frontmatter

**Approval:** pending execution
