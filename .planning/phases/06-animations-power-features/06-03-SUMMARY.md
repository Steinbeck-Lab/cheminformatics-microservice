---
phase: 06-animations-power-features
plan: 03
subsystem: ui
tags: [cmdk, command-palette, breadcrumbs, navigation, react-router, motion]

# Dependency graph
requires:
  - phase: 06-01
    provides: AnimatedOutlet, page transitions, tab cross-fades
provides:
  - Cmd+K command palette with fuzzy search across pages, tools, and molecules
  - Breadcrumb navigation for tool pages (section > tool hierarchy)
  - Navigation registry as single source of truth for all navigable items
affects: []

# Tech tracking
tech-stack:
  added: [cmdk]
  patterns: [navigation-registry-pattern, global-keyboard-shortcut-pattern]

key-files:
  created:
    - frontend/src/data/navigationRegistry.ts
    - frontend/src/components/ui/command.tsx
    - frontend/src/components/common/CommandPalette.tsx
    - frontend/src/components/common/Breadcrumbs.tsx
    - frontend/src/__tests__/components/command-palette.test.tsx
    - frontend/src/__tests__/components/breadcrumbs.test.tsx
  modified:
    - frontend/src/App.tsx
    - frontend/src/__tests__/setup.ts
    - frontend/package.json

key-decisions:
  - "Navigation registry uses exact tool IDs from each page's tabs array -- verified all 25 tools match"
  - "Added ResizeObserver and scrollIntoView mocks to test setup for cmdk jsdom compatibility"
  - "Breadcrumbs placed inside <main> before <Suspense> to leverage existing header clearance padding"
  - "shadcn CLI still writes to literal @/ directory -- manually moved command.tsx (same pattern as Phase 04)"

patterns-established:
  - "Navigation registry pattern: centralized page/tool/molecule data exported from navigationRegistry.ts"
  - "Global keyboard shortcut pattern: useEffect on document keydown for Cmd+K toggle"

requirements-completed: [POWER-01, POWER-03]

# Metrics
duration: 7min
completed: 2026-03-13
---

# Phase 06 Plan 03: Command Palette + Breadcrumbs Summary

**Cmd+K command palette with cmdk-powered fuzzy search across 7 pages, 25 tools, and 6 example molecules, plus breadcrumb navigation with animated tool name cross-fade on tool pages**

## Performance

- **Duration:** 7 min
- **Started:** 2026-03-13T22:30:13Z
- **Completed:** 2026-03-13T22:37:44Z
- **Tasks:** 2
- **Files modified:** 10

## Accomplishments
- Command palette opens with Cmd+K (Mac) / Ctrl+K (Windows), showing categorized results (Pages, Tools, Example Molecules, Recent)
- Navigation registry provides single source of truth with 7 pages, 25 tools, 6 molecules and keyword data for fuzzy search
- Breadcrumbs appear on /chem, /convert, /depict, /tools routes showing section > tool hierarchy with 150ms cross-fade
- Breadcrumbs hidden on non-tool pages and mobile (<768px), CommandPalette uses glass-bold dialog styling
- 14 new tests (6 command-palette + 8 breadcrumbs) all passing, production build succeeds

## Task Commits

Each task was committed atomically:

1. **Task 1: Navigation registry + shadcn Command component + CommandPalette** - `048b356` (feat)
2. **Task 2: Breadcrumbs component + App.tsx wiring** - `5c95d44` (feat)

## Files Created/Modified
- `frontend/src/data/navigationRegistry.ts` - Centralized pages, tools, molecules data with NavEntry/MoleculeEntry types and helper functions
- `frontend/src/components/ui/command.tsx` - shadcn/ui Command component wrapping cmdk with glass-bold dialog styling
- `frontend/src/components/common/CommandPalette.tsx` - Global Cmd+K palette with categorized search and navigation
- `frontend/src/components/common/Breadcrumbs.tsx` - Section > Tool breadcrumbs with animated cross-fade
- `frontend/src/App.tsx` - Wired CommandPalette and Breadcrumbs into Layout
- `frontend/src/__tests__/setup.ts` - Added ResizeObserver and scrollIntoView mocks for cmdk
- `frontend/src/__tests__/components/command-palette.test.tsx` - 6 tests for keyboard shortcuts and content rendering
- `frontend/src/__tests__/components/breadcrumbs.test.tsx` - 8 tests for route filtering and display logic
- `frontend/package.json` - Added cmdk dependency

## Decisions Made
- Navigation registry uses exact tool IDs verified from each page's tabs array (hosecodes not hose-codes, structureexplorer not structure-explorer, etc.)
- Added ResizeObserver and Element.scrollIntoView mocks to test setup -- cmdk requires both, jsdom provides neither
- shadcn CLI still outputs to literal @/ directory -- moved command.tsx manually (same approach as Phase 04)
- Breadcrumbs placed inside `<main>` before `<Suspense>` to leverage existing `main > *` padding-top for header clearance

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] shadcn CLI output to literal @/ directory**
- **Found during:** Task 1 (shadcn command component installation)
- **Issue:** shadcn CLI created files in literal `@/components/ui/` instead of `src/components/ui/`
- **Fix:** Copied command.tsx to correct location, removed literal @/ directory
- **Files modified:** frontend/src/components/ui/command.tsx
- **Verification:** Import resolves correctly, build passes
- **Committed in:** 048b356 (Task 1 commit)

**2. [Rule 3 - Blocking] jsdom missing ResizeObserver and scrollIntoView**
- **Found during:** Task 1 (command palette tests)
- **Issue:** cmdk requires ResizeObserver and Element.scrollIntoView, jsdom implements neither
- **Fix:** Added ResizeObserver class mock and scrollIntoView noop to test setup.ts
- **Files modified:** frontend/src/__tests__/setup.ts
- **Verification:** All 6 command palette tests pass
- **Committed in:** 048b356 (Task 1 commit)

---

**Total deviations:** 2 auto-fixed (2 blocking)
**Impact on plan:** Both auto-fixes necessary for test environment compatibility. No scope creep.

## Issues Encountered
- Pre-existing test failures in molecule-card-enhanced.test.tsx (missing ComparisonProvider) -- not caused by this plan, out of scope

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Phase 06 is now complete (3/3 plans done)
- All animation polish, inline SMILES preview, command palette, and breadcrumb features are implemented
- Ready for Phase 07 or final project completion

---
*Phase: 06-animations-power-features*
*Completed: 2026-03-13*
