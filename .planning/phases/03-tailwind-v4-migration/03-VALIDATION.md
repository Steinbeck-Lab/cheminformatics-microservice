---
phase: 3
slug: tailwind-v4-migration
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-03-12
---

# Phase 3 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | Vitest 4.1.0 |
| **Config file** | `frontend/vite.config.mts` (test section) |
| **Quick run command** | `cd frontend && npx vite build` |
| **Full suite command** | `cd frontend && npx vitest run` |
| **Estimated runtime** | ~15 seconds (build), ~10 seconds (tests) |

---

## Sampling Rate

- **After every task commit:** Run `cd frontend && npx vite build` (verifies CSS pipeline intact)
- **After every plan wave:** Run `cd frontend && npx vitest run` + visual spot-check
- **Before `/gsd:verify-work`:** Full 9-page manual visual verification + build succeeds + tests pass
- **Max feedback latency:** 15 seconds

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|-----------|-------------------|-------------|--------|
| 03-01-01 | 01 | 1 | FRAME-02 | smoke | `cd frontend && npx vite build` | N/A (build) | pending |
| 03-01-02 | 01 | 1 | FRAME-02 | smoke | `cd frontend && npx vitest run` | setup.ts, App.test.tsx | pending |
| 03-02-01 | 02 | 1 | FRAME-04 | manual | Visual verification of all 9 pages | N/A (manual) | pending |
| 03-02-02 | 02 | 1 | FRAME-04 | smoke | `cd frontend && npx vite build 2>&1` (no warnings) | N/A (build output) | pending |

*Status: pending / green / red / flaky*

---

## Wave 0 Requirements

Existing infrastructure covers all phase requirements.

- Vitest already configured in `vite.config.mts`
- Smoke tests exist (`App.test.tsx`)
- Vite build serves as CSS pipeline validation

*No Wave 0 test stubs needed.*

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| All 9 pages render with correct styling | FRAME-04 | No visual regression testing infrastructure (no screenshot comparison) | Open each page in browser, verify: colors match, layouts intact, dark mode works, cards/buttons/forms styled correctly |
| Custom color shades preserved | FRAME-02 | Visual accuracy of custom hex values | Check dark mode backgrounds (gray-950), blue accent shades in buttons and links |
| Component class specificity maintained | FRAME-04 | CSS specificity changes are visual | Verify btn-*, card-*, glass classes render identically to v3 |
| Border default color preserved | FRAME-04 | Default border color changed from gray-200 to currentColor in v4 | Check all borders on cards, inputs, dividers for unexpected color changes |

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 15s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
