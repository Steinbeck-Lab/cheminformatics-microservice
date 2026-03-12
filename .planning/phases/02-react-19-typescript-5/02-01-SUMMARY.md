---
phase: 02-react-19-typescript-5
plan: 01
subsystem: ui
tags: [react-19, typescript-5, vite, eslint, strict-mode]

# Dependency graph
requires:
  - phase: 01-vite-migration
    provides: Vite build system, ESM module structure, React 18 baseline
provides:
  - React 19.2.4 runtime with updated @types/react@19
  - Full TypeScript strict-mode codebase (61 .ts/.tsx files, 0 errors)
  - Shared API type definitions (frontend/src/types/api.ts, molecule.ts)
  - TypeScript-aware ESLint config with typescript-eslint
  - Three-file tsconfig reference pattern (root, app, node)
affects: [02-02-motion-migration, frontend-components, frontend-services]

# Tech tracking
tech-stack:
  added: [react@19.2.4, react-dom@19.2.4, "@types/react@19", "@types/react-dom@19", "@types/dompurify", typescript-eslint@8]
  removed: ["@headlessui/react@1.7.19"]
  patterns: [typescript-strict-mode, three-file-tsconfig-reference, typed-service-layer, typed-context-pattern]

key-files:
  created:
    - frontend/tsconfig.json
    - frontend/tsconfig.app.json
    - frontend/src/types/api.ts
    - frontend/src/types/molecule.ts
    - frontend/src/types/global.d.ts
  modified:
    - frontend/package.json
    - frontend/vite.config.mts
    - frontend/eslint.config.mjs
    - frontend/tsconfig.node.json
    - frontend/index.html
    - frontend/src/services/api.ts
    - frontend/src/services/chemService.ts
    - frontend/src/services/convertService.ts
    - frontend/src/services/depictService.ts
    - frontend/src/services/ocsrService.ts
    - frontend/src/services/toolsService.ts
    - frontend/src/context/AppContext.tsx

key-decisions:
  - "Used vite.config.mts extension (not .ts) to force ESM resolution without requiring type:module in package.json"
  - "Removed dead @headlessui/react to unblock React 19 peer dependency resolution"
  - "TypeScript strict mode with 0 errors achieved on first pass — no relaxation needed"
  - "Kept framer-motion imports as-is (migration to motion/react deferred to Plan 02-02)"

patterns-established:
  - "Three-file tsconfig: root references, tsconfig.app.json for src/, tsconfig.node.json for vite config"
  - "Typed service layer: explicit param/return types with shared API type interfaces"
  - "Typed context pattern: createContext<T | undefined>(undefined) with runtime guard hook"
  - "Error handling pattern: error instanceof Error with typed catch blocks"

requirements-completed: [FRAME-01, FRAME-03]

# Metrics
duration: 10min
completed: 2026-03-12
---

# Phase 02 Plan 01: React 19 + TypeScript Strict Mode Summary

**React 19.2.4 upgrade with full TypeScript strict-mode conversion of 61 files, shared API type system, and typed service layer**

## Performance

- **Duration:** 10 min
- **Started:** 2026-03-12T18:10:00Z
- **Completed:** 2026-03-12T18:20:00Z
- **Tasks:** 2
- **Files modified:** 75

## Accomplishments
- Upgraded React 18 to 19.2.4 with matching @types packages, resolving CVE-2025-55182
- Converted all 61 frontend files from .js/.jsx to .ts/.tsx with git mv (preserving history)
- Created comprehensive shared type system: 20+ API response interfaces in src/types/api.ts, molecule data types in src/types/molecule.ts
- Typed all 6 service files and AppContext with explicit parameter/return types
- Achieved 0 tsc errors in strict mode on first attempt, 0 ESLint errors
- Production build passes (567 modules, 3.5s), test suite passes (1/1)

## Task Commits

Each task was committed atomically:

1. **Task 1: Install React 19, create TypeScript infrastructure, rename all files** - `ca39fb2` (feat)
2. **Task 2: Add TypeScript types to services, context, and config** - `43285af` (feat)

## Files Created/Modified

**Type definitions (created):**
- `frontend/src/types/api.ts` - 20+ shared API response interfaces (DescriptorResult, DepictionOptions, SugarRemovalOptions, etc.)
- `frontend/src/types/molecule.ts` - RecentMolecule, MoleculeData, ApiConfig, AppContextValue interfaces
- `frontend/src/types/global.d.ts` - $3Dmol CDN ambient declaration, Vite env vars

**TypeScript config (created):**
- `frontend/tsconfig.json` - Root references file
- `frontend/tsconfig.app.json` - App config: strict, ES2020, jsx react-jsx, noUncheckedIndexedAccess

**Config updates (modified):**
- `frontend/package.json` - React 19.2.4, TypeScript type packages, script globs .ts/.tsx
- `frontend/vite.config.mts` - Renamed from .mjs, removed treat-js-as-jsx plugin, .ts setupFiles
- `frontend/eslint.config.mjs` - typescript-eslint integration, @typescript-eslint/no-unused-vars
- `frontend/tsconfig.node.json` - Updated include to vite.config.mts
- `frontend/index.html` - Script src /src/index.tsx

**Typed services (modified):**
- `frontend/src/services/api.ts` - AxiosInstance, AxiosError typed interceptors
- `frontend/src/services/chemService.ts` - 14 functions with typed params/returns
- `frontend/src/services/convertService.ts` - 12 functions with typed params/returns
- `frontend/src/services/depictService.ts` - All functions using DepictionOptions/EnhancedDepictionOptions
- `frontend/src/services/ocsrService.ts` - File|Blob param, OCSRResult return
- `frontend/src/services/toolsService.ts` - SugarRemovalOptions, ExtractionOptions, ChemicalFilterOptions

**Typed context (modified):**
- `frontend/src/context/AppContext.tsx` - createContext<AppContextValue | undefined> with runtime guard

**Renamed (61 files):**
- All .js -> .ts, .jsx -> .tsx via git mv across src/ (pages, components, hooks, utils, services, context, tests)

## Decisions Made
- **vite.config.mts extension:** Vite bundles .ts config using CJS require() when package.json lacks "type":"module". The old .mjs forced ESM; .mts does the same for TypeScript. This avoids ERR_REQUIRE_ESM without changing package.json.
- **Removed @headlessui/react:** Dead dependency (not imported anywhere) that blocked React 19 peer resolution. Verified no imports exist before removal.
- **No tsconfig relaxation needed:** strict mode with noUncheckedIndexedAccess passed with 0 errors, so noUnusedLocals/noUnusedParameters were kept enabled.
- **framer-motion imports preserved:** Plan 02-02 handles the migration to motion/react; changing imports here would create merge conflicts.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Removed @headlessui/react to unblock npm install**
- **Found during:** Task 1 (npm install)
- **Issue:** @headlessui/react@1.7.19 has peer react@"^16 || ^17 || ^18", blocking @types/react@19 installation with ERESOLVE
- **Fix:** Verified no imports exist in codebase, ran `npm uninstall @headlessui/react`
- **Files modified:** frontend/package.json, frontend/package-lock.json
- **Verification:** npm install succeeds, no runtime errors
- **Committed in:** ca39fb2 (Task 1 commit)

**2. [Rule 3 - Blocking] Used .mts extension for vite config to fix ESM resolution**
- **Found during:** Task 2 (vite.config conversion)
- **Issue:** Renaming vite.config.mjs to vite.config.ts caused ERR_REQUIRE_ESM because Vite uses CJS require() for .ts config files
- **Fix:** Used .mts extension instead of .ts to force ESM resolution (equivalent to old .mjs behavior)
- **Files modified:** frontend/vite.config.mts, frontend/tsconfig.node.json
- **Verification:** npm test, npm run build both pass
- **Committed in:** 43285af (Task 2 commit)

---

**Total deviations:** 2 auto-fixed (2 blocking issues, Rule 3)
**Impact on plan:** Both auto-fixes were necessary to unblock execution. No scope creep.

## Issues Encountered
- The large chunk warning from Vite build (1,271 KB index bundle) is pre-existing and unrelated to this migration. Code-splitting is a future optimization.
- 118 ESLint warnings (mostly react/no-unescaped-entities in JSX strings) are pre-existing; 0 errors.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- TypeScript infrastructure complete, ready for Plan 02-02 (framer-motion to motion/react migration)
- All 61 files are .ts/.tsx with strict mode; future plans can add types incrementally
- Shared type system in src/types/ provides foundation for component prop typing

## Self-Check: PASSED

All 13 key files verified present. Both task commits (ca39fb2, 43285af) verified in git log.

---
*Phase: 02-react-19-typescript-5*
*Completed: 2026-03-12*
