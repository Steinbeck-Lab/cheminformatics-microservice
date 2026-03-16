---
phase: 1
slug: vite-migration
status: draft
nyquist_compliant: false
wave_0_complete: false
created: 2026-03-12
---

# Phase 1 — Validation Strategy

> Per-phase validation contract for feedback sampling during execution.

---

## Test Infrastructure

| Property | Value |
|----------|-------|
| **Framework** | Vitest 4.1.0 + @testing-library/react 16.3.2 |
| **Config file** | vite.config.js (test section) -- Wave 0 creates this |
| **Quick run command** | `cd frontend && npx vitest run` |
| **Full suite command** | `cd frontend && npx vitest run --coverage` |
| **Estimated runtime** | ~15 seconds |

---

## Sampling Rate

- **After every task commit:** Run `cd frontend && npx vite build && npm audit --audit-level=low`
- **After every plan wave:** Run `cd frontend && npx vitest run --coverage && npx vite build`
- **Before `/gsd:verify-work`:** Full suite must be green + `npm audit` clean + Docker build succeeds
- **Max feedback latency:** 30 seconds

---

## Per-Task Verification Map

| Task ID | Plan | Wave | Requirement | Test Type | Automated Command | File Exists | Status |
|---------|------|------|-------------|-----------|-------------------|-------------|--------|
| 01-01-01 | 01 | 0 | BUILD-01 | smoke | `cd frontend && npx vite build` | Wave 0 creates | ⬜ pending |
| 01-01-02 | 01 | 0 | BUILD-06 | unit | `cd frontend && npx vitest run src/__tests__/env.test.js` | Wave 0 creates | ⬜ pending |
| 01-02-01 | 02 | 1 | BUILD-02 | security | `cd frontend && npm audit --audit-level=low` | N/A -- npm command | ⬜ pending |
| 01-02-02 | 02 | 1 | BUILD-03 | manual-only | Inspect package.json for `"overrides"` key | N/A -- file inspection | ⬜ pending |
| 01-02-03 | 02 | 1 | BUILD-04 | manual-only | `cd frontend && npm outdated` | N/A -- npm command | ⬜ pending |
| 01-03-01 | 03 | 2 | BUILD-05 | integration | `docker build -t test-frontend frontend/` | N/A -- Docker command | ⬜ pending |

*Status: ⬜ pending · ✅ green · ❌ red · ⚠️ flaky*

---

## Wave 0 Requirements

- [ ] `frontend/vite.config.js` -- Vite + Vitest configuration (test section with jsdom environment)
- [ ] `frontend/src/__tests__/setup.js` -- jest-dom matcher setup for Vitest
- [ ] `frontend/src/__tests__/App.test.jsx` -- smoke test: App renders without errors
- [ ] `frontend/src/__tests__/env.test.js` -- verify env vars use import.meta.env pattern
- [ ] `frontend/eslint.config.mjs` -- ESLint v9 flat config
- [ ] Framework install: `npm install -D vitest @vitest/coverage-v8 jsdom @testing-library/react @testing-library/jest-dom`

*Wave 0 creates the test infrastructure before any feature work begins.*

---

## Manual-Only Verifications

| Behavior | Requirement | Why Manual | Test Instructions |
|----------|-------------|------------|-------------------|
| No overrides in package.json | BUILD-03 | Simple file inspection, no runtime behavior | Open `package.json`, confirm no `"overrides"` key exists |
| Dependencies at latest stable | BUILD-04 | Version currency is a point-in-time check | Run `npm outdated`, confirm no outdated packages |
| All 9 pages load correctly | BUILD-01 | Visual UI verification across routes | Start dev server, navigate to each of the 9 routes, confirm rendering |
| Docker frontend serves app | BUILD-05 | End-to-end container verification | Build and run Docker image, confirm app loads in browser |

---

## Validation Sign-Off

- [ ] All tasks have `<automated>` verify or Wave 0 dependencies
- [ ] Sampling continuity: no 3 consecutive tasks without automated verify
- [ ] Wave 0 covers all MISSING references
- [ ] No watch-mode flags
- [ ] Feedback latency < 30s
- [ ] `nyquist_compliant: true` set in frontmatter

**Approval:** pending
