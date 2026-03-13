---
phase: 06-animations-power-features
plan: 02
subsystem: ui
tags: [react, debounce, smiles, 2d-depiction, glass-skeleton, preview]

# Dependency graph
requires:
  - phase: 04.1-visual-design-system-glassmorphism-liquid-glass-bento-grids
    provides: GlassSkeleton component, glass-bold CSS utility
  - phase: 05-loading-states-ux
    provides: GlassSkeleton feedback component
provides:
  - useDebounce generic hook for delayed state updates
  - SMILESPreview inline 2D structure preview component
  - Integrated preview in all 23 tool views via SMILESInput
affects: [06-animations-power-features]

# Tech tracking
tech-stack:
  added: []
  patterns: [useDebounce hook, inline image preview with glass skeleton loading, withFakeTimers test helper]

key-files:
  created:
    - frontend/src/hooks/useDebounce.ts
    - frontend/src/components/common/SMILESPreview.tsx
    - frontend/src/__tests__/components/smiles-preview.test.tsx
  modified:
    - frontend/src/components/common/SMILESInput.tsx

key-decisions:
  - "withFakeTimers helper pattern for vitest -- restores real timers before test returns to prevent global afterEach cleanup timeout with React effect teardowns"
  - "SMILESPreview uses img src={url} pattern with onLoad/onError for zero-fetch overhead -- browser handles the image request natively"

patterns-established:
  - "withFakeTimers: wrap test body to ensure vi.useRealTimers() runs before global cleanup, preventing hook timeouts"
  - "Inline preview: self-contained component that manages its own debounce/loading/error state, requires zero parent logic"

requirements-completed: [POWER-02]

# Metrics
duration: 6min
completed: 2026-03-13
---

# Phase 06 Plan 02: Inline SMILES Preview Summary

**Debounced 2D structure preview below SMILES input with glass skeleton loading, silent error hiding, and 500ms keystroke debounce across all 23 tool views**

## Performance

- **Duration:** 6 min
- **Started:** 2026-03-13T22:19:50Z
- **Completed:** 2026-03-13T22:26:49Z
- **Tasks:** 2
- **Files modified:** 4

## Accomplishments
- Generic useDebounce hook reusable across the frontend
- SMILESPreview component with glass skeleton loading state, image display, and silent error hiding
- Integrated into SMILESInput so all 23 tool views automatically get inline 2D structure previews
- 7 tests covering debounce behavior, empty/whitespace input, loading skeleton, successful image load, error hiding, and keystroke debounce

## Task Commits

Each task was committed atomically:

1. **Task 1: useDebounce hook + SMILESPreview component** - `de4813d` (feat - TDD)
2. **Task 2: Integrate SMILESPreview into SMILESInput** - `4b896dc` (feat)

_Note: Task 1 used TDD flow -- tests written first (RED), then implementation (GREEN)_

## Files Created/Modified
- `frontend/src/hooks/useDebounce.ts` - Generic debounce hook using useState + useEffect
- `frontend/src/components/common/SMILESPreview.tsx` - Inline 2D structure preview with glass skeleton loading
- `frontend/src/__tests__/components/smiles-preview.test.tsx` - 7 tests for debounce + preview behavior
- `frontend/src/components/common/SMILESInput.tsx` - Added SMILESPreview import and JSX insertion

## Decisions Made
- Used withFakeTimers helper pattern to wrap test bodies instead of beforeEach/afterEach for fake timer management. This avoids a vitest/React interaction where global afterEach cleanup() hangs when fake timers are active during React effect teardowns.
- SMILESPreview uses native img src={url} with onLoad/onError callbacks rather than fetch() + blob URL. This leverages the browser's image loading pipeline and requires no additional state management for the HTTP request lifecycle.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed vitest fake timer interaction with React cleanup**
- **Found during:** Task 1 (TDD GREEN phase)
- **Issue:** Using beforeEach(vi.useFakeTimers()) + afterEach(cleanup/useRealTimers) caused "Hook timed out in 10000ms" errors because React's effect teardown callbacks require real setTimeout to flush, but fake timers intercept them
- **Fix:** Created withFakeTimers() helper that wraps test body with try/finally ensuring vi.useRealTimers() runs before the test function returns, so the global afterEach cleanup() executes with real timers
- **Files modified:** frontend/src/__tests__/components/smiles-preview.test.tsx
- **Verification:** All 7 tests pass consistently with no timeouts
- **Committed in:** de4813d (Task 1 commit)

---

**Total deviations:** 1 auto-fixed (1 bug fix)
**Impact on plan:** Test infrastructure fix needed for correct fake timer + React interaction. No scope creep.

## Issues Encountered
- Pre-existing test failure in molecule-card-enhanced.test.tsx (needs ComparisonProvider context) -- unrelated to this plan's changes, out of scope.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- useDebounce hook available for reuse in other components
- SMILESPreview integrated into all tool views via shared SMILESInput
- Ready for Plan 03 (next features in Phase 06)

---
*Phase: 06-animations-power-features*
*Completed: 2026-03-13*

## Self-Check: PASSED
- All 4 created files exist on disk
- Both task commits (de4813d, 4b896dc) found in git log
