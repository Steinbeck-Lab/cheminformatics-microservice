---
phase: 01-vite-migration
plan: 01
subsystem: ui
tags: [vite, react, swc, cra-migration, esbuild, npm-audit]

# Dependency graph
requires: []
provides:
  - Vite 6.x build toolchain replacing CRA
  - Zero npm audit vulnerabilities
  - All env vars migrated from REACT_APP_* to VITE_*
  - Working production build producing dist/ directory
affects: [01-02-PLAN, 01-03-PLAN, docker-frontend]

# Tech tracking
tech-stack:
  added: [vite@6.4.1, "@vitejs/plugin-react-swc@4.3.0", vite-plugin-checker@0.11.0, rollup-plugin-visualizer@7.0.1, typescript@5.9.3]
  removed: [react-scripts@5.0.1, web-vitals@2.1.4, "@testing-library/jest-dom", "@testing-library/react", "@testing-library/user-event"]
  patterns: [vite-esm-build, import-meta-env, jsx-in-js-custom-plugin]

key-files:
  created:
    - frontend/vite.config.mjs
    - frontend/index.html
    - frontend/src/index.jsx
  modified:
    - frontend/package.json
    - frontend/src/services/api.js
    - frontend/src/context/AppContext.js
    - frontend/src/components/common/MoleculeCard.jsx
    - frontend/src/components/common/HighlightedMoleculeCard.jsx
    - frontend/src/components/chem/StandardizeView.jsx
    - frontend/src/components/depict/StructureDrawView.jsx
    - frontend/src/components/tools/InChIView.jsx
    - frontend/src/components/tools/RInChIView.jsx
    - frontend/.env.production
    - frontend/.env.local.example
    - frontend/.gitignore
    - frontend/tailwind.config.js

key-decisions:
  - "Used vite.config.mjs (ESM) instead of vite.config.js to avoid needing type:module in package.json, which would break CJS PostCSS/Tailwind configs"
  - "Added custom Vite plugin (treat-js-as-jsx) to handle JSX in .js files rather than renaming all .js files to .jsx"
  - "Renamed entry point index.js to index.jsx because Vite HTML entry parsing requires explicit .jsx extension for JSX content"
  - "File locations differed from plan: MoleculeCard at common/ not chem/, HighlightedMoleculeCard at common/ not chem/, StructureDrawView at depict/ not tools/, InChIView at tools/ not convert/, RInChIView at tools/ not convert/"
  - "StandardizeView.jsx found as 5th REACT_APP_API_URL reference (not ConvertPage.jsx as plan suggested)"

patterns-established:
  - "JSX-in-JS handling: Custom Vite plugin using transformWithEsbuild for .js files containing JSX"
  - "Environment variables: All frontend env vars prefixed with VITE_ and accessed via import.meta.env"
  - "Dev server proxy: /v1 and /latest routes proxied to localhost:8000"

requirements-completed: [BUILD-01, BUILD-02, BUILD-03, BUILD-04, BUILD-06]

# Metrics
duration: 13min
completed: 2026-03-12
---

# Phase 1 Plan 01: CRA to Vite Migration Summary

**Replaced CRA with Vite 6.4.1 + SWC, migrated all env vars from REACT_APP_* to VITE_*, achieved zero npm audit vulnerabilities**

## Performance

- **Duration:** 13 min
- **Started:** 2026-03-12T14:59:54Z
- **Completed:** 2026-03-12T15:12:48Z
- **Tasks:** 2
- **Files modified:** 18

## Accomplishments
- Removed react-scripts and all CRA dependencies, eliminating all 16 npm audit vulnerabilities
- Installed Vite 6.4.1 with SWC-based React plugin for fast builds (3.8s production build)
- Migrated 9 process.env references across 8 source files to Vite-compatible import.meta.env
- Removed package.json overrides, eslintConfig, browserslist, and homepage (CRA-specific cruft)
- Created Vite config with dev server proxy, hidden sourcemaps, and bundle analysis support

## Task Commits

Each task was committed atomically:

1. **Task 1: Uninstall CRA, install Vite stack, restructure package.json** - `a0018d7` (feat)
2. **Task 2: Create Vite config, migrate index.html, migrate all env vars** - `23f4be4` (feat)

## Files Created/Modified
- `frontend/vite.config.mjs` - Vite build config with SWC plugin, proxy, JSX-in-JS custom plugin
- `frontend/index.html` - Vite entry HTML moved from public/ to project root
- `frontend/src/index.jsx` - Renamed from index.js for Vite HTML entry compatibility
- `frontend/package.json` - Cleaned: no overrides, no eslintConfig, no browserslist, Vite scripts
- `frontend/src/services/api.js` - REACT_APP_API_URL -> VITE_API_URL, NODE_ENV -> import.meta.env.DEV
- `frontend/src/context/AppContext.js` - REACT_APP_API_URL -> VITE_API_URL
- `frontend/src/components/common/MoleculeCard.jsx` - REACT_APP_API_URL -> VITE_API_URL
- `frontend/src/components/common/HighlightedMoleculeCard.jsx` - REACT_APP_API_URL -> VITE_API_URL
- `frontend/src/components/chem/StandardizeView.jsx` - REACT_APP_API_URL -> VITE_API_URL
- `frontend/src/components/depict/StructureDrawView.jsx` - PUBLIC_URL -> direct path
- `frontend/src/components/tools/InChIView.jsx` - PUBLIC_URL -> direct path
- `frontend/src/components/tools/RInChIView.jsx` - PUBLIC_URL -> direct path
- `frontend/.env.production` - REACT_APP_API_URL -> VITE_API_URL
- `frontend/.env.local.example` - Updated to VITE_ prefix
- `frontend/.gitignore` - /build -> /dist, added stats.html
- `frontend/tailwind.config.js` - Updated content path from ./public/index.html to ./index.html

## Decisions Made
- Used `vite.config.mjs` instead of `vite.config.js` to avoid adding `"type": "module"` to package.json, which would have required renaming PostCSS and Tailwind configs to `.cjs`
- Created a custom Vite plugin (`treat-js-as-jsx`) that uses `transformWithEsbuild` to handle JSX in `.js` files during build, rather than renaming all `.js` files to `.jsx` (massive rename would have been disruptive)
- Renamed only the entry point `index.js` to `index.jsx` because Vite's HTML entry parsing requires explicit `.jsx` extension before any plugins run

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Entry point index.js required rename to index.jsx**
- **Found during:** Task 2 (Vite config creation)
- **Issue:** Vite's HTML parser runs before any transform plugins, so it cannot parse JSX in `.js` entry files
- **Fix:** Renamed `src/index.js` to `src/index.jsx` and updated HTML script src
- **Files modified:** frontend/src/index.jsx, frontend/index.html
- **Verification:** Build succeeds with 534 modules transformed
- **Committed in:** 23f4be4

**2. [Rule 3 - Blocking] JSX in .js files not transformed during Vite build**
- **Found during:** Task 2 (Vite config creation)
- **Issue:** esbuild only applies JSX transform to `.jsx`/`.tsx` files by default; the codebase uses JSX in `.js` files (CRA convention)
- **Fix:** Added custom Vite plugin that pre-transforms `.js` files with `transformWithEsbuild` using JSX loader
- **Files modified:** frontend/vite.config.mjs
- **Verification:** Build succeeds, all 534 modules transformed correctly
- **Committed in:** 23f4be4

**3. [Rule 3 - Blocking] Plan file paths did not match actual codebase locations**
- **Found during:** Task 2 (env var migration)
- **Issue:** Plan referenced files at `components/chem/MoleculeCard.jsx`, `components/chem/HighlightedMoleculeCard.jsx`, `components/tools/StructureDrawView.jsx`, `components/convert/RInChIView.jsx`, `components/convert/InChIView.jsx` -- actual locations differ (common/, depict/, tools/)
- **Fix:** Used grep to find all actual `process.env` references and migrated at correct locations. Also found `StandardizeView.jsx` as the 5th REACT_APP_API_URL reference (not ConvertPage.jsx)
- **Files modified:** All 8 source files at their correct locations
- **Verification:** `grep -rn "process\.env" src/` returns zero results
- **Committed in:** 23f4be4

**4. [Rule 1 - Bug] vite.config.js with type:module broke PostCSS config**
- **Found during:** Task 2 (build testing)
- **Issue:** Adding `"type": "module"` to package.json caused PostCSS/Tailwind CJS configs to fail with `module is not defined in ES module scope`
- **Fix:** Used `.mjs` extension for vite config instead, avoiding the need for `"type": "module"` entirely
- **Files modified:** frontend/vite.config.mjs (renamed from .js)
- **Verification:** Build succeeds, PostCSS/Tailwind configs load correctly
- **Committed in:** 23f4be4

---

**Total deviations:** 4 auto-fixed (1 bug, 3 blocking)
**Impact on plan:** All auto-fixes were necessary to achieve a working build. No scope creep.

## Issues Encountered
- The `@vitejs/plugin-react-swc` build-mode plugin has no `transform` hook -- it only configures esbuild via the `config` hook. This means JSX transformation during production builds relies entirely on esbuild, which only handles `.jsx`/`.tsx` by default. Required adding a custom plugin for `.js` JSX handling.
- Vite's internal `vite:build-import-analysis` plugin parses files with Rollup's parser before esbuild transforms are applied, creating a chicken-and-egg problem for the HTML entry point file specifically.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Vite build toolchain fully operational, ready for Plan 02 (ESLint flat config + Vitest setup)
- The `test`, `lint`, and `prepare` scripts in package.json reference tools not yet installed (vitest, eslint, husky) -- Plan 02 will install these
- PostCSS and Tailwind configs remain as CJS (`.js` with `module.exports`) -- this will be addressed in Phase 3 (Tailwind v4 migration)

---
## Self-Check: PASSED

- frontend/vite.config.mjs: FOUND
- frontend/index.html: FOUND
- frontend/src/index.jsx: FOUND
- SUMMARY.md: FOUND
- Commit a0018d7 (Task 1): FOUND
- Commit 23f4be4 (Task 2): FOUND

---
*Phase: 01-vite-migration*
*Completed: 2026-03-12*
