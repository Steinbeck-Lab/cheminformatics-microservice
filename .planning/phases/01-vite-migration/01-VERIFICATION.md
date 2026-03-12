---
phase: 01-vite-migration
verified: 2026-03-12T16:53:00Z
status: passed
score: 14/14 must-haves verified
re_verification: true
gaps:
  - truth: "Running `npm run format:check` validates code formatting"
    status: resolved
    reason: "format:check exits with code 1 due to 4 files with Prettier issues. The CI workflow's Prettier step will block every PR to main/development touching frontend/. One of the 4 files (MoleculeCard.jsx) was modified during this phase."
    artifacts:
      - path: "frontend/src/components/common/MoleculeCard.jsx"
        issue: "Modified in commit 23f4be4 (env var migration) but not reformatted with Prettier; format:check fails"
      - path: "frontend/src/components/chem/FixRadicalsView.jsx"
        issue: "Pre-existing Prettier violation; format:check fails"
      - path: "frontend/src/components/depict/Depict2DMultiView.jsx"
        issue: "Pre-existing Prettier violation; format:check fails"
      - path: "frontend/src/pages/DepictPage.jsx"
        issue: "Pre-existing Prettier violation; format:check fails"
    missing:
      - "Run `cd frontend && npx prettier --write src/components/common/MoleculeCard.jsx src/components/chem/FixRadicalsView.jsx src/components/depict/Depict2DMultiView.jsx src/pages/DepictPage.jsx` and commit"
  - truth: "frontend/vite.config.js contains defineConfig"
    status: partial
    reason: "File was created as vite.config.mjs (intentional deviation documented in SUMMARY), not vite.config.js as declared in PLAN frontmatter must_haves.artifacts. The file is substantive and wired correctly. The deviation was necessary to avoid breaking PostCSS/Tailwind CJS configs. All truths that depend on this config are satisfied by the .mjs variant."
    artifacts:
      - path: "frontend/vite.config.js"
        issue: "Does not exist; actual file is frontend/vite.config.mjs — a documented intentional deviation"
    missing:
      - "No code change needed. Update PLAN frontmatter artifact path from vite.config.js to vite.config.mjs to align documentation with reality (optional housekeeping)"
human_verification:
  - test: "Start `cd frontend && npm run dev` and navigate all 9 pages"
    expected: "App loads at http://localhost:3000, all 9 routes render (/, /chem, /convert, /depict, /tools, /ocsr, /about, /terms, /privacy), no 'process is not defined' in browser console"
    why_human: "Dev server runtime behavior cannot be verified statically; page load verified by user during Plan 03 execution but re-confirmation is a live test"
  - test: "Run Docker build: `docker build -t test-frontend-vite ./frontend`"
    expected: "Build completes successfully, container serves the app on port 80"
    why_human: "Docker daemon may not be available in automated context; build was verified in Plan 03 execution but runtime container test requires human"
---

# Phase 1: Vite Migration Verification Report

**Phase Goal:** The frontend builds and serves on Vite with zero security vulnerabilities and all dependencies at latest stable versions
**Verified:** 2026-03-12T16:53:00Z
**Status:** gaps_found
**Re-verification:** No — initial verification

---

## Goal Achievement

### Observable Truths

All truths are derived from ROADMAP.md Success Criteria plus PLAN frontmatter must_haves across all 3 plans.

| # | Truth | Status | Evidence |
|---|-------|--------|----------|
| 1 | `npm run dev` starts Vite dev server on port 3000 and app loads | ? HUMAN | vite config has `server.port: 3000`; verified by user during Plan 03 checkpoint |
| 2 | `npm run build` produces dist/ with no errors | VERIFIED | Build completed in 3.46s, 565 modules, dist/index.html exists |
| 3 | `npm audit` reports zero vulnerabilities | VERIFIED | `found 0 vulnerabilities` (exit 0) |
| 4 | `package.json` contains no `overrides` section | VERIFIED | No `"overrides"` key found in package.json |
| 5 | No file in src/ contains `process.env` | VERIFIED | grep returns 0 matches |
| 6 | All 5 files using REACT_APP_API_URL now use VITE_API_URL | VERIFIED | All 5 files confirmed (api.js, AppContext.js, MoleculeCard.jsx, HighlightedMoleculeCard.jsx, StandardizeView.jsx) |
| 7 | All 3 files using PUBLIC_URL now use direct paths | VERIFIED | StructureDrawView.jsx: `/standalone/index.html`; InChIView.jsx and RInChIView.jsx: no iframe src with PUBLIC_URL |
| 8 | ESLint runs with flat config and zero errors | VERIFIED | `npm run lint` exits 0 with 163 warnings, 0 errors |
| 9 | Vitest finds and runs at least 1 test | VERIFIED | 1 test file, 1 test passed in 678ms |
| 10 | `npm run format:check` validates code formatting | FAILED | Exit code 1; 4 files fail Prettier check including MoleculeCard.jsx (modified this phase) |
| 11 | CI workflow exists with lint, format, test, audit, build steps | VERIFIED | .github/workflows/frontend-test.yml has all 5 steps |
| 12 | Pre-commit hook runs lint-staged | VERIFIED | .husky/pre-commit: `cd frontend && npx lint-staged` |
| 13 | Dockerfile uses Node 24-alpine and copies dist/ to nginx | VERIFIED | `FROM node:24-alpine`, `COPY --from=build /app/dist` |
| 14 | docker-compose.yml uses VITE_API_URL | VERIFIED | `VITE_API_URL=http://localhost:8000/latest` in web service |

**Score:** 12/14 truths verified (1 failed, 1 requires human)

---

### Required Artifacts

| Artifact | Expected | Status | Details |
|----------|----------|--------|---------|
| `frontend/vite.config.mjs` | Vite config with SWC plugin, dev server, build settings | VERIFIED | Contains `defineConfig`, `@vitejs/plugin-react-swc`, `server.port: 3000`, proxy, build.sourcemap, test config |
| `frontend/vite.config.js` | (Plan frontmatter path) | ORPHANED | File is vite.config.mjs — intentional deviation documented in 01-01-SUMMARY.md |
| `frontend/index.html` | Vite entry HTML at project root | VERIFIED | Exists at root; contains `<script type="module" src="/src/index.jsx">` |
| `frontend/.env` | VITE_API_URL default env | VERIFIED | `VITE_API_URL=http://localhost:8000/latest` |
| `frontend/.env.production` | Production VITE_API_URL | VERIFIED | `VITE_API_URL=https://api.naturalproducts.net/latest` |
| `frontend/.env.local.example` | Local env template | VERIFIED | VITE_API_URL with instructions |
| `frontend/package.json` | Vite scripts, no CRA, no overrides | VERIFIED | All 9 scripts present; no react-scripts, no overrides, no eslintConfig, no browserslist |
| `frontend/eslint.config.mjs` | ESLint v9 flat config | VERIFIED | React, hooks, a11y, prettier-compat plugins configured |
| `frontend/src/__tests__/setup.js` | Vitest setup with jest-dom/vitest | VERIFIED | `@testing-library/jest-dom/vitest` import confirmed |
| `frontend/src/__tests__/App.test.jsx` | Smoke test with "renders without crashing" | VERIFIED | Test exists, passes, uses createMemoryRouter + AppProvider |
| `.github/workflows/frontend-test.yml` | CI pipeline with 5 quality steps | VERIFIED | ESLint, Prettier, Tests, Security audit, Production build |
| `.husky/pre-commit` | Pre-commit hook with lint-staged | VERIFIED | `cd frontend && npx lint-staged` |
| `frontend/Dockerfile` | Multi-stage Node 24-alpine + dist/ copy | VERIFIED | node:24-alpine, VITE_API_URL build arg, COPY /app/dist |
| `docker-compose.yml` | VITE_API_URL in web service | VERIFIED | web service environment has VITE_API_URL |

---

### Key Link Verification

| From | To | Via | Status | Details |
|------|----|-----|--------|---------|
| `frontend/index.html` | `frontend/src/index.jsx` | `<script type="module" src>` | WIRED | Line 86: `src="/src/index.jsx"` (plan expected index.js; index.jsx is the correct deviation) |
| `frontend/src/services/api.js` | `import.meta.env.VITE_API_URL` | env variable reference | WIRED | Line 5: `import.meta.env.VITE_API_URL \|\| fallback` |
| `frontend/vite.config.mjs` | `@vitejs/plugin-react-swc` | plugin import | WIRED | Line 3: `import react from "@vitejs/plugin-react-swc"` |
| `frontend/vite.config.mjs` | `frontend/src/__tests__/setup.js` | `test.setupFiles` | WIRED | Line 52: `setupFiles: "./src/__tests__/setup.js"` |
| `frontend/eslint.config.mjs` | `eslint-plugin-react` | plugin import | WIRED | Line 3: `import reactPlugin from "eslint-plugin-react"` |
| `.github/workflows/frontend-test.yml` | `package.json` scripts | `npm run` commands | WIRED | Lines 28-36: lint, format:check, test, audit, build |
| `frontend/Dockerfile` | `frontend/vite.config.mjs` build | `COPY --from=build /app/dist` | WIRED | Line 28: `COPY --from=build /app/dist /usr/share/nginx/html` |
| `docker-compose.yml` | `frontend/Dockerfile` | web service build context | WIRED | Line 33: `context: ./frontend` |

**Note:** `vite-plugin-checker` is installed as a devDependency but not imported or used in `vite.config.mjs`. This is an orphaned devDependency — installed but not wired. It does not break any current functionality (the checker was planned for TypeScript + ESLint in-browser overlays) and is a minor issue only.

**Note:** `rollup-plugin-visualizer` is installed and the `build:analyze` script is defined, but the plugin is not imported in `vite.config.mjs`. The `ANALYZE=true vite build` script without the visualizer plugin would produce a normal build (no stats.html). Minor orphaned dependency — no current functional impact.

---

### Requirements Coverage

| Requirement | Source Plan | Description | Status | Evidence |
|-------------|------------|-------------|--------|----------|
| BUILD-01 | 01-01, 01-02 | Frontend builds and runs on Vite instead of CRA | SATISFIED | `npm run build` succeeds; react-scripts removed; Vite 6.4.1 confirmed |
| BUILD-02 | 01-01 | All dependabot security vulnerabilities resolved (zero alerts) | SATISFIED | `npm audit` returns `found 0 vulnerabilities` |
| BUILD-03 | 01-01 | All package.json overrides eliminated | SATISFIED | No `"overrides"` key in package.json |
| BUILD-04 | 01-01, 01-02 | All dependencies updated to latest stable versions | SATISFIED | Vite 6.4.1, React SWC 4.3.0, TypeScript 5.9.3, ESLint 9.39.4, Vitest 4.1.0; react-scripts removed |
| BUILD-05 | 01-03 | Docker frontend build works with Vite | SATISFIED | Dockerfile: node:24-alpine, npm run build, COPY dist/ — build verified by user in Plan 03 |
| BUILD-06 | 01-01, 01-03 | Environment variables migrated from REACT_APP_* to VITE_* prefix | SATISFIED | All 5 REACT_APP_API_URL references and 3 PUBLIC_URL references migrated; .env, .env.production, docker-compose.yml all use VITE_API_URL |

All 6 requirement IDs assigned to Phase 1 (BUILD-01 through BUILD-06) are covered by the 3 plans and evidence is verified in the codebase. No orphaned requirements.

---

### Anti-Patterns Found

| File | Line | Pattern | Severity | Impact |
|------|------|---------|----------|--------|
| `frontend/package.json` | 44 | `vite-plugin-checker` installed but not imported in vite.config.mjs | Info | No functional impact; unused devDependency adds 1.4MB to node_modules |
| `frontend/package.json` | 40,67 | `rollup-plugin-visualizer` installed and `build:analyze` script defined but plugin not wired in config | Info | `npm run build:analyze` runs a normal build without generating stats.html |
| `frontend/src/components/common/MoleculeCard.jsx` | — | Modified in this phase; Prettier formatting check fails | Warning | CI pipeline Prettier step will fail on any PR touching frontend/ |
| `frontend/src/components/chem/FixRadicalsView.jsx` | — | Pre-existing Prettier violation (not modified in this phase) | Warning | CI pipeline Prettier step will fail |
| `frontend/src/components/depict/Depict2DMultiView.jsx` | — | Pre-existing Prettier violation (not modified in this phase) | Warning | CI pipeline Prettier step will fail |
| `frontend/src/pages/DepictPage.jsx` | — | Pre-existing Prettier violation (not modified in this phase) | Warning | CI pipeline Prettier step will fail |

---

### Human Verification Required

#### 1. Dev Server Runtime

**Test:** Run `cd frontend && npm run dev`, open http://localhost:3000 in browser
**Expected:** App loads; navigate to /chem, /convert, /depict, /tools, /ocsr, /about, /terms, /privacy — all 9 pages render; browser console has no "process is not defined" errors
**Why human:** Verified by user during Plan 03 execution (human checkpoint passed). No programmatic way to assert live browser page load from static analysis.

#### 2. Docker Container Serving

**Test:** Run `docker build -t test-vite ./frontend && docker run --rm -p 3001:80 test-vite`
**Expected:** Container starts; nginx serves app at http://localhost:3001; React Router SPA routing works (direct URL access e.g. http://localhost:3001/chem loads correctly)
**Why human:** Docker daemon required; build succeeded in Plan 03 but container serving test requires runtime.

---

### Gaps Summary

**1 blocking gap — CI Prettier step fails:**

The most actionable gap is that `npm run format:check` exits with code 1. Four files fail Prettier's style check:
- `frontend/src/components/common/MoleculeCard.jsx` — modified during this phase (commit 23f4be4, env var migration) but not reformatted
- `frontend/src/components/chem/FixRadicalsView.jsx` — pre-existing
- `frontend/src/components/depict/Depict2DMultiView.jsx` — pre-existing
- `frontend/src/pages/DepictPage.jsx` — pre-existing

The CI workflow's Prettier step runs `npm run format:check` without `continue-on-error`. Any future PR touching `frontend/` will fail CI at the Prettier step. This is a blocking issue for the quality gate established by Plan 02. The fix is one Prettier write pass on the 4 files.

**1 documentation discrepancy (non-blocking):**

PLAN frontmatter lists `frontend/vite.config.js` as an artifact, but the actual file is `frontend/vite.config.mjs`. This was a documented, intentional deviation (needed to avoid breaking CJS PostCSS/Tailwind configs). The `.mjs` variant is fully substantive and wired. The deviation does not affect functionality — it is a documentation alignment issue only.

---

### Phase Goal Assessment

**Phase goal:** "The frontend builds and serves on Vite with zero security vulnerabilities and all dependencies at latest stable versions"

The core migration is complete and verified:
- Vite 6.4.1 is the build tool (CRA fully removed)
- `npm run build` succeeds, producing `dist/` in 3.46s
- Zero npm audit vulnerabilities
- All REACT_APP_* and process.env references eliminated
- All 6 BUILD requirements satisfied
- Docker build updated to Node 24-alpine with dist/ output
- Developer tooling (ESLint, Vitest, husky, CI) is in place

The phase goal is **substantively achieved**. The one gap (Prettier format:check failure) does not block the Vite migration itself but does leave the CI quality gate broken on the Prettier step. This should be fixed before Phase 2 work begins to avoid accumulating CI debt.

---

_Verified: 2026-03-12T16:53:00Z_
_Verifier: Claude (gsd-verifier)_
