---
phase: 05-loading-states-ux
plan: 01
subsystem: ui
tags: [sonner, react-error-boundary, skeleton, toast, glassmorphism, shimmer, error-handling, axios]

# Dependency graph
requires:
  - phase: 04.1-visual-design-system
    provides: glass-bold utility, clay-interactive, cn() util, motion/react, shadcn Button
provides:
  - GlassSkeleton component with 5 size variants and shimmer animation
  - GlassErrorCard with retry button and motion entrance animation
  - ToolSkeleton with 4 content-shaped skeleton variants
  - EmptyState for tools with no results
  - RouteLoadingFallback for route-level code splitting Suspense
  - Glass-styled Sonner Toaster wrapper
  - Contextual error-messages.ts with per-domain messages
  - Axios response interceptor with deduplicated toast for network/5xx/429 errors
affects: [05-02, 05-03, 05-04, 05-05]

# Tech tracking
tech-stack:
  added: [sonner, react-error-boundary]
  patterns: [glass-shimmer-skeleton, contextual-error-messages, axios-toast-interceptor]

key-files:
  created:
    - frontend/src/components/feedback/GlassSkeleton.tsx
    - frontend/src/components/feedback/GlassErrorCard.tsx
    - frontend/src/components/feedback/ToolSkeleton.tsx
    - frontend/src/components/feedback/EmptyState.tsx
    - frontend/src/components/feedback/RouteLoadingFallback.tsx
    - frontend/src/components/ui/sonner.tsx
    - frontend/src/lib/error-messages.ts
    - frontend/src/__tests__/components/glass-skeleton.test.tsx
    - frontend/src/__tests__/components/glass-error-card.test.tsx
    - frontend/src/__tests__/components/sonner.test.tsx
    - frontend/src/__tests__/components/tool-loading.test.tsx
  modified:
    - frontend/src/styles/tailwind.css
    - frontend/src/services/api.ts
    - frontend/package.json

key-decisions:
  - "Force-add lib/error-messages.ts due to root .gitignore lib/ pattern (same approach as Phase 04)"
  - "Used data-testid attributes on all feedback components for reliable test selectors"

patterns-established:
  - "GlassSkeleton as base building block: composable via variant prop, used by ToolSkeleton and RouteLoadingFallback"
  - "Axios interceptor deduplication via Sonner toast id parameter (network-error, rate-limit)"
  - "Domain-scoped error messages: getErrorMessage(domain, error) inspects error type for contextual messaging"

requirements-completed: [LOAD-01, LOAD-02, LOAD-03, LOAD-04]

# Metrics
duration: 4min
completed: 2026-03-13
---

# Phase 5 Plan 01: Feedback Component Library Summary

**Glass shimmer skeletons, error cards with retry, Sonner toast system, and Axios interceptor for global network/5xx error toasts**

## Performance

- **Duration:** 4 min
- **Started:** 2026-03-13T17:07:32Z
- **Completed:** 2026-03-13T17:12:26Z
- **Tasks:** 2
- **Files modified:** 14

## Accomplishments
- Created 5 feedback components (GlassSkeleton, GlassErrorCard, ToolSkeleton, EmptyState, RouteLoadingFallback) with glass aesthetic
- Installed Sonner and created glass-styled Toaster wrapper with dark mode support
- Created contextual error message mapping with 6 domains (chem, convert, depict, tools, ocsr, network)
- Wired Axios response interceptor with deduplicated toast notifications for network/5xx/429 errors
- Added 13 passing tests across 4 test files with real assertions

## Task Commits

Each task was committed atomically:

1. **Task 1: Install Sonner + shimmer CSS + feedback components + error messages** - `d2ed85d` (feat)
2. **Task 2: Wire Axios interceptor + test stubs for all feedback components** - `4e80c11` (feat)

## Files Created/Modified
- `frontend/src/components/feedback/GlassSkeleton.tsx` - Base glass shimmer skeleton with 5 size variants
- `frontend/src/components/feedback/GlassErrorCard.tsx` - Error card with AlertTriangle icon, message, and retry button
- `frontend/src/components/feedback/ToolSkeleton.tsx` - Content-shaped skeleton variants for tool output areas
- `frontend/src/components/feedback/EmptyState.tsx` - Empty state with Beaker icon for tools with no results
- `frontend/src/components/feedback/RouteLoadingFallback.tsx` - Inline Suspense fallback with skeleton blocks
- `frontend/src/components/ui/sonner.tsx` - Glass-styled Sonner Toaster with dark mode from AppContext
- `frontend/src/lib/error-messages.ts` - Domain-scoped contextual error messages with getErrorMessage()
- `frontend/src/styles/tailwind.css` - Added shimmer keyframes and animate-shimmer animation
- `frontend/src/services/api.ts` - Added toast import and network/5xx/429 error toast notifications
- `frontend/package.json` - Added sonner and react-error-boundary dependencies
- `frontend/src/__tests__/components/glass-skeleton.test.tsx` - 4 tests for shimmer, variants, className
- `frontend/src/__tests__/components/glass-error-card.test.tsx` - 4 tests for message, retry, click
- `frontend/src/__tests__/components/sonner.test.tsx` - 1 test for Toaster render
- `frontend/src/__tests__/components/tool-loading.test.tsx` - 4 tests for ToolSkeleton, EmptyState, RouteLoadingFallback

## Decisions Made
- Force-add `lib/error-messages.ts` due to root `.gitignore` `lib/` pattern (same established approach from Phase 04)
- Used `data-testid` attributes on all feedback components for reliable test selectors independent of CSS classes
- Kept pre-existing `molecule-card-enhanced.test.tsx` failures untouched (out of scope -- logged below)

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
- Pre-existing test failures in `molecule-card-enhanced.test.tsx` (7 tests failing) -- not related to this plan's changes, out of scope.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- All feedback components ready for consumption by downstream plans (05-02 through 05-05)
- GlassSkeleton, GlassErrorCard, ToolSkeleton, EmptyState, RouteLoadingFallback all importable from `@/components/feedback/`
- Sonner Toaster importable from `@/components/ui/sonner` -- needs to be mounted in App.tsx (Plan 05-02+ scope)
- getErrorMessage importable from `@/lib/error-messages`
- Axios interceptor active -- toast notifications fire on network/5xx/429 errors globally

## Self-Check: PASSED

All 12 files verified present. Both commit hashes (d2ed85d, 4e80c11) verified in git log.

---
*Phase: 05-loading-states-ux*
*Completed: 2026-03-13*
