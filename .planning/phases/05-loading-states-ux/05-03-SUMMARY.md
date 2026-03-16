---
phase: 05-loading-states-ux
plan: 03
subsystem: ui
tags: [react-lazy, suspense, code-splitting, vite-chunks, responsive, scrollbar-hide, navigation, tap-targets]

# Dependency graph
requires:
  - phase: 05-loading-states-ux
    plan: 01
    provides: RouteLoadingFallback component, Sonner Toaster wrapper
provides:
  - Route-level code splitting with React.lazy for all 9 pages
  - Vite vendor chunk splitting (vendor-react, vendor-motion, vendor-ui)
  - Toaster wired into Layout for global toast notifications
  - Suspense boundary with RouteLoadingFallback in Layout
  - CaffeineMolecule3D hidden on mobile for GPU performance
  - scrollbar-hide utility for horizontal tab strips
  - 44px minimum tap targets on mobile nav and hamburger
  - Active nav link visual indicator (desktop pill + mobile border)
  - Navigation wayfinding tests and lazy-loading Suspense test
affects: [05-04, 05-05]

# Tech tracking
tech-stack:
  added: []
  patterns: [route-level-code-splitting, vendor-chunk-splitting, scrollbar-hide-utility, mobile-responsive-guards]

key-files:
  created:
    - frontend/src/__tests__/components/navigation.test.tsx
  modified:
    - frontend/src/App.tsx
    - frontend/vite.config.mts
    - frontend/src/pages/HomePage.tsx
    - frontend/src/pages/ConvertPage.tsx
    - frontend/src/pages/DepictPage.tsx
    - frontend/src/pages/ToolsPage.tsx
    - frontend/src/components/common/Header.tsx
    - frontend/src/components/common/Navigation.tsx
    - frontend/src/styles/tailwind.css
    - frontend/src/__tests__/App.test.tsx

key-decisions:
  - "Main entry chunk reduced from 1268KB to 321KB via React.lazy code splitting"
  - "Tab pages use mobile dropdown (md:hidden) and desktop horizontal tabs -- scrollbar-hide only needed on desktop overflow-x-auto container"
  - "Hamburger button bumped from h-8 w-8 (32px) to h-11 w-11 min-h-[44px] for WCAG tap target compliance"
  - "Lazy-loading test uses simulated lazy component to avoid jsdom Canvas/IntersectionObserver limitations"

patterns-established:
  - "React.lazy for all page-level routes with Suspense fallback={RouteLoadingFallback} in Layout"
  - "scrollbar-hide @utility in tailwind.css for hidden scrollbar on overflow-x-auto containers"
  - "min-h-[44px] on all mobile interactive elements for WCAG touch target compliance"

requirements-completed: [UX-01, UX-02, UX-03, UX-04]

# Metrics
duration: 6min
completed: 2026-03-13
---

# Phase 5 Plan 03: Code Splitting, Responsive Fixes & Navigation Summary

**Route-level React.lazy code splitting reducing main chunk from 1268KB to 321KB, with vendor chunking, mobile responsive fixes, and navigation wayfinding tests**

## Performance

- **Duration:** 6 min
- **Started:** 2026-03-13T17:16:28Z
- **Completed:** 2026-03-13T17:22:28Z
- **Tasks:** 2
- **Files modified:** 15

## Accomplishments
- Replaced all 9 eager page imports with React.lazy + dynamic import(), reducing main entry chunk from 1268KB to 321KB
- Configured Vite manualChunks producing vendor-react (78KB), vendor-motion (96KB), vendor-ui (60KB) separate chunks
- Wired RouteLoadingFallback and Toaster into Layout (both from Plan 01)
- Hidden CaffeineMolecule3D on mobile (<768px) via `hidden md:block` for GPU performance
- Added scrollbar-hide utility and applied to tab strips on ConvertPage, DepictPage, ToolsPage
- Bumped mobile hamburger and nav links to 44px minimum tap targets
- Created 3 navigation wayfinding tests and 3 App/lazy-loading tests (all passing)

## Task Commits

Each task was committed atomically:

1. **Task 1: Code splitting (React.lazy + Suspense + Toaster) and vendor chunking** - `902700e` (feat)
2. **Task 2: Responsive fixes, navigation wayfinding, scrollbar-hide, and tests** - `ff9da0e` (feat)

## Files Created/Modified
- `frontend/src/App.tsx` - 9 React.lazy imports, Suspense wrapping Outlet, Toaster in Layout
- `frontend/vite.config.mts` - manualChunks for vendor-react, vendor-motion, vendor-ui
- `frontend/src/pages/HomePage.tsx` - CaffeineMolecule3D wrapped in hidden md:block
- `frontend/src/pages/ConvertPage.tsx` - flex-nowrap scrollbar-hide on desktop tab strip
- `frontend/src/pages/DepictPage.tsx` - flex-nowrap scrollbar-hide on desktop tab strip
- `frontend/src/pages/ToolsPage.tsx` - flex-nowrap scrollbar-hide on desktop tab strip
- `frontend/src/components/common/Header.tsx` - Hamburger button 44px tap target
- `frontend/src/components/common/Navigation.tsx` - Mobile nav 44px tap target + border-l-4 active indicator
- `frontend/src/styles/tailwind.css` - scrollbar-hide @utility
- `frontend/src/__tests__/components/navigation.test.tsx` - 3 tests for active nav highlighting and tap targets
- `frontend/src/__tests__/App.test.tsx` - 3 tests including lazy-loading Suspense verification

## Decisions Made
- Main entry chunk reduced from 1268KB to 321KB via React.lazy code splitting -- well under the 400KB target
- Tab pages already use mobile dropdown (md:hidden) and desktop horizontal tabs, so scrollbar-hide only applied to the desktop overflow-x-auto container
- Hamburger button bumped from h-8/w-8 (32px) to h-11/w-11 with min-h-[44px] for WCAG tap target compliance
- Lazy-loading test uses a simulated lazy component with Suspense to avoid jsdom Canvas/IntersectionObserver limitations when rendering full App with Three.js and motion/react

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Fixed extra closing brace in 4 chem view components**
- **Found during:** Task 2 (build verification)
- **Issue:** DescriptorsView, CoconutPreProcessingView, PubChemLookupView, ClassyfireView all had `/>}}` (extra `}`) at end of GlassErrorCard JSX, causing esbuild transform error
- **Fix:** Removed the extra `}` in each file
- **Files modified:** frontend/src/components/chem/DescriptorsView.tsx, CoconutPreProcessingView.tsx, PubChemLookupView.tsx, ClassyfireView.tsx
- **Verification:** Build succeeds after fix
- **Committed in:** ff9da0e (Task 2 commit)

**2. [Rule 3 - Blocking] Fixed malformed import in AllFiltersView.tsx**
- **Found during:** Task 2 (build verification)
- **Issue:** Line 18 had `, Loader2 }` starting with comma instead of proper import list continuation
- **Fix:** Reformatted import to `Loader2,` on its own line
- **Files modified:** frontend/src/components/chem/AllFiltersView.tsx
- **Verification:** Build succeeds after fix
- **Committed in:** ff9da0e (Task 2 commit)

---

**Total deviations:** 2 auto-fixed (both Rule 3 blocking)
**Impact on plan:** Both were pre-existing syntax errors preventing production build. Required fixing to verify Task 1 and Task 2 work. No scope creep.

## Issues Encountered
- Pre-existing test failures in `molecule-card-enhanced.test.tsx` (7 tests) -- same as documented in 05-01 summary, unrelated to this plan

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- All feedback components from Plan 01 are now wired into the app (Toaster in Layout, RouteLoadingFallback in Suspense)
- Code splitting active -- each page loads as a separate chunk on navigation
- Responsive mobile fixes complete -- ready for remaining UX plans (05-04, 05-05)

## Self-Check: PASSED

All 12 files verified present. Both commit hashes (902700e, ff9da0e) verified in git log.

---
*Phase: 05-loading-states-ux*
*Completed: 2026-03-13*
