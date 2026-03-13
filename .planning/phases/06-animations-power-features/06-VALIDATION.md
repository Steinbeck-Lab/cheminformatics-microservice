---
phase: 6
slug: animations-power-features
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-03-13
---

# Phase 6 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | Vitest 4.1.0 + @testing-library/react 16.3.2 |
| **Config file** | `frontend/vite.config.mts` (test section) |
| **Quick run command** | `cd frontend && npx vitest run --reporter=verbose` |
| **Full suite command** | `cd frontend && npx vitest run` |
| **Estimated runtime** | ~15 seconds |

---

## Sampling Rate

- **After every task commit:** Run `cd frontend && npx vitest run --reporter=verbose`
- **After every plan wave:** Run `cd frontend && npx vitest run`
- **Before `/gsd:verify-work`:** Full suite must be green
- **Max feedback latency:** 15 seconds

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|-----------|-------------------|-------------|--------|
| 06-01-01 | 01 | 1 | ANIM-01 | integration | `cd frontend && npx vitest run src/__tests__/components/page-transitions.test.tsx` | ❌ W0 | ⬜ pending |
| 06-01-02 | 01 | 1 | ANIM-02 | unit | `cd frontend && npx vitest run src/__tests__/components/button-animations.test.tsx` | ❌ W0 | ⬜ pending |
| 06-01-03 | 01 | 1 | ANIM-03 | unit | `cd frontend && npx vitest run src/__tests__/components/focus-ring.test.tsx` | ❌ W0 | ⬜ pending |
| 06-01-04 | 01 | 1 | ANIM-04 | manual-only | N/A (animation timing not testable in jsdom) | N/A | ⬜ pending |
| 06-01-05 | 01 | 1 | ANIM-05 | manual-only | N/A (stagger timing not testable in jsdom) | N/A | ⬜ pending |
| 06-01-06 | 01 | 1 | ANIM-06 | manual-only | N/A (layout animation not testable in jsdom) | N/A | ⬜ pending |
| 06-02-01 | 02 | 2 | POWER-01 | integration | `cd frontend && npx vitest run src/__tests__/components/command-palette.test.tsx` | ❌ W0 | ⬜ pending |
| 06-02-02 | 02 | 2 | POWER-02 | integration | `cd frontend && npx vitest run src/__tests__/components/smiles-preview.test.tsx` | ❌ W0 | ⬜ pending |
| 06-03-01 | 03 | 2 | POWER-03 | unit | `cd frontend && npx vitest run src/__tests__/components/breadcrumbs.test.tsx` | ❌ W0 | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

- [ ] `src/__tests__/components/page-transitions.test.tsx` — stubs for ANIM-01 (AnimatedOutlet renders with key)
- [ ] `src/__tests__/components/button-animations.test.tsx` — stubs for ANIM-02 (CSS class presence)
- [ ] `src/__tests__/components/focus-ring.test.tsx` — stubs for ANIM-03 (focus ring classes)
- [ ] `src/__tests__/components/command-palette.test.tsx` — stubs for POWER-01 (open/close, search, navigation)
- [ ] `src/__tests__/components/smiles-preview.test.tsx` — stubs for POWER-02 (debounce, show/hide, error handling)
- [ ] `src/__tests__/components/breadcrumbs.test.tsx` — stubs for POWER-03 (section/tool display, visibility rules)

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| Tab content cross-fades on switch | ANIM-04 | Animation timing not testable in jsdom | Switch between tabs on ChemPage, ConvertPage, DepictPage, ToolsPage — verify ~150ms cross-fade, no abrupt content swap |
| Results lists stagger in | ANIM-05 | Stagger timing not testable in jsdom | Navigate to a page that displays results — verify items animate in with staggered timing |
| Layout changes animate | ANIM-06 | Layout animation not testable in jsdom | Trigger expand/collapse or show/hide interactions — verify smooth layout transitions |

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 15s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
