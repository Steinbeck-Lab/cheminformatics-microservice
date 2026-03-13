---
phase: 5
slug: loading-states-ux
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-03-13
---

# Phase 5 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | Vitest 4.1.0 |
| **Config file** | `frontend/vite.config.mts` (test section) |
| **Quick run command** | `cd frontend && npx vitest run` |
| **Full suite command** | `cd frontend && npx vitest run --coverage` |
| **Estimated runtime** | ~15 seconds |

---

## Sampling Rate

- **After every task commit:** Run `cd frontend && npx vitest run`
- **After every plan wave:** Run `cd frontend && npx vitest run --coverage`
- **Before `/gsd:verify-work`:** Full suite must be green
- **Max feedback latency:** 15 seconds

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|-----------|-------------------|-------------|--------|
| 05-01-01 | 01 | 1 | LOAD-01 | unit | `cd frontend && npx vitest run src/__tests__/components/glass-skeleton.test.tsx -t "shimmer"` | ❌ W0 | ⬜ pending |
| 05-01-02 | 01 | 1 | LOAD-02 | unit | `cd frontend && npx vitest run src/__tests__/components/sonner.test.tsx` | ❌ W0 | ⬜ pending |
| 05-01-03 | 01 | 1 | LOAD-03 | unit | `cd frontend && npx vitest run src/__tests__/components/glass-error-card.test.tsx` | ❌ W0 | ⬜ pending |
| 05-01-04 | 01 | 1 | LOAD-04 | unit | `cd frontend && npx vitest run src/__tests__/components/tool-loading.test.tsx` | ❌ W0 | ⬜ pending |
| 05-02-01 | 02 | 1 | UX-01 | manual-only | Manual viewport testing at 375px/768px | N/A | ⬜ pending |
| 05-02-02 | 02 | 1 | UX-02 | unit | `cd frontend && npx vitest run src/__tests__/App.test.tsx -t "lazy"` | ❌ W0 | ⬜ pending |
| 05-02-03 | 02 | 1 | UX-03 | unit | `cd frontend && npx vitest run src/__tests__/components/navigation.test.tsx` | ❌ W0 | ⬜ pending |
| 05-02-04 | 02 | 1 | UX-04 | smoke | `cd frontend && npm run build 2>&1 \| grep -c "dist/assets/"` | N/A | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

- [ ] `src/__tests__/components/glass-skeleton.test.tsx` — stubs for LOAD-01
- [ ] `src/__tests__/components/sonner.test.tsx` — stubs for LOAD-02
- [ ] `src/__tests__/components/glass-error-card.test.tsx` — stubs for LOAD-03
- [ ] `src/__tests__/components/tool-loading.test.tsx` — stubs for LOAD-04
- [ ] `src/__tests__/components/navigation.test.tsx` — stubs for UX-03

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| Responsive layout no horizontal overflow at 375px/768px | UX-01 | Viewport testing requires browser rendering; JSDOM can't measure overflow | 1. Open app in Chrome DevTools 2. Toggle device toolbar to 375px (iPhone SE) 3. Navigate all pages, verify no horizontal scrollbar 4. Repeat at 768px (iPad) |
| Build produces multiple chunks | UX-04 | Build output verification, not runtime behavior | 1. Run `cd frontend && npm run build` 2. Verify `dist/assets/` contains multiple .js files (not single monolith) 3. Verify main entry chunk < 400KB gzipped |

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 15s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
