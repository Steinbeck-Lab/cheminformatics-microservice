---
phase: 03-tailwind-v4-migration
plan: 02
subsystem: ui
tags: [tailwindcss, v4-migration, edge-cases, visual-verification, css-variables]

# Dependency graph
requires:
  - phase: 03-tailwind-v4-migration/plan-01
    provides: Tailwind v4 CSS-first config, automated class renames across 50 files
provides:
  - All v4 migration edge cases resolved (transform-none, ring defaults, bg-opacity merging)
  - Zero visual regressions confirmed across all 9 pages
  - Complete Tailwind v3 to v4 migration verified and production-ready
affects: [04-component-system]

# Tech tracking
tech-stack:
  added: []
  patterns: [slash-modifier opacity syntax, individual transform resets for v4, CSS custom properties for ring styling]

key-files:
  created: []
  modified:
    - frontend/src/styles/tailwind.css
    - frontend/src/pages/ChemPage.tsx
    - frontend/src/components/chem/AllFiltersView.tsx

key-decisions:
  - "Replaced md:transform-none with md:translate-x-0 md:translate-y-0 md:scale-100 md:rotate-0 for v4 individual transform resets"
  - "Converted all bg-opacity-* two-class patterns to slash modifier syntax across 21 TSX files"
  - "Replaced --tw-ring-color/opacity CSS variables with standard ring-* utilities in tailwind.css"

patterns-established:
  - "Opacity: use slash modifier (bg-red-900/30) not separate bg-opacity-* class"
  - "Transforms: use individual properties (translate-x-0) not transform-none shorthand"
  - "Ring styling: use ring-* utilities not --tw-* CSS variables"

requirements-completed: [FRAME-04]

# Metrics
duration: 5min
completed: 2026-03-12
---

# Phase 3 Plan 02: Edge Case Audit + Visual Verification Summary

**All Tailwind v4 edge cases resolved (transform-none, ring defaults, bg-opacity merging across 23 files) with zero visual regressions confirmed across all 9 pages**

## Performance

- **Duration:** 5 min (Task 1 execution) + human verification pause
- **Started:** 2026-03-12T19:18:29Z
- **Completed:** 2026-03-12T19:30:36Z
- **Tasks:** 2
- **Files modified:** 23

## Accomplishments
- All `transform-none` usages replaced with individual transform resets (translate-x-0, translate-y-0, scale-100, rotate-0) in ChemPage.tsx
- All 21 `bg-opacity-*` two-class patterns converted to v4 slash modifier syntax across TSX files
- `--tw-ring-color` and `--tw-ring-opacity` CSS variables replaced with standard ring utilities in tailwind.css
- Human verified all 9 pages with zero visual regressions in both light and dark modes

## Task Commits

Each task was committed atomically:

1. **Task 1: Manual edge case audit and fixes** - `169ea0c` (fix)
2. **Task 2: Visual verification of all 9 pages** - human-verify checkpoint (approved, no commit needed)

## Files Created/Modified
- `frontend/src/styles/tailwind.css` - Replaced --tw-ring-color/opacity with standard ring utilities
- `frontend/src/pages/ChemPage.tsx` - Replaced md:transform-none with individual transform resets
- `frontend/src/components/chem/AllFiltersView.tsx` - bg-opacity to slash modifier conversion
- `frontend/src/components/chem/ClassyfireView.tsx` - bg-opacity to slash modifier conversion
- `frontend/src/components/chem/CoconutPreProcessingView.tsx` - bg-opacity to slash modifier conversion
- `frontend/src/components/chem/DescriptorsView.tsx` - bg-opacity to slash modifier conversion
- `frontend/src/components/chem/ErtlFunctionalGroupView.tsx` - bg-opacity to slash modifier conversion
- `frontend/src/components/chem/FixRadicalsView.tsx` - bg-opacity to slash modifier conversion
- `frontend/src/components/chem/HOSECodeView.tsx` - bg-opacity to slash modifier conversion
- `frontend/src/components/chem/NPlikenessView.tsx` - bg-opacity to slash modifier conversion
- `frontend/src/components/chem/PubChemLookupView.tsx` - bg-opacity to slash modifier conversion
- `frontend/src/components/chem/StandardizedTautomerView.tsx` - bg-opacity to slash modifier conversion
- `frontend/src/components/chem/StandardizeView.tsx` - bg-opacity to slash modifier conversion
- `frontend/src/components/chem/StereoisomersView.tsx` - bg-opacity to slash modifier conversion
- `frontend/src/components/chem/StructureErrorView.tsx` - bg-opacity to slash modifier conversion
- `frontend/src/components/chem/TanimotoView.tsx` - bg-opacity to slash modifier conversion
- `frontend/src/components/convert/FormatConversionView.tsx` - bg-opacity to slash modifier conversion
- `frontend/src/components/convert/Mol2DView.tsx` - bg-opacity to slash modifier conversion
- `frontend/src/components/convert/Mol3DView.tsx` - bg-opacity to slash modifier conversion
- `frontend/src/components/depict/Depict2DMultiView.tsx` - bg-opacity to slash modifier conversion
- `frontend/src/components/depict/Depict3DView.tsx` - bg-opacity to slash modifier conversion
- `frontend/src/components/depict/MoleculeDepiction3D.tsx` - bg-opacity to slash modifier conversion
- `frontend/src/components/tools/SugarRemovalView.tsx` - bg-opacity to slash modifier conversion

## Decisions Made
- Replaced `md:transform-none` with `md:translate-x-0 md:translate-y-0 md:scale-100 md:rotate-0` to cover all transform axes -- the sidebar in ChemPage uses translateX for mobile slide-in, so all axes are reset for clean desktop layout
- Converted all `bg-opacity-*` two-class patterns to slash modifier syntax (e.g., `dark:bg-red-900 dark:bg-opacity-30` became `dark:bg-red-900/30`) -- this is the v4-idiomatic approach and eliminates deprecated utilities
- Replaced `--tw-ring-color` and `--tw-ring-opacity` CSS variables with standard `ring-*` utilities -- v4 may not expose internal `--tw-*` variables the same way

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
None - all edge cases identified in the plan were found and fixed as expected.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Tailwind v4 migration fully complete -- all edge cases resolved, zero visual regressions
- CSS-first configuration stable and production-ready
- Ready for Phase 4: Component System + Dark Mode (shadcn/ui integration)
- Note: Local development requires Node 20+ for Tailwind v4's Oxide engine

## Self-Check: PASSED

- FOUND: 03-02-SUMMARY.md
- FOUND: commit 169ea0c (Task 1)
- FOUND: frontend/src/styles/tailwind.css
- FOUND: frontend/src/pages/ChemPage.tsx
- FOUND: frontend/src/components/chem/AllFiltersView.tsx

---
*Phase: 03-tailwind-v4-migration*
*Completed: 2026-03-12*
