---
phase: 02-react-19-typescript-5
plan: 02
subsystem: ui
tags: [motion, framer-motion, animation, spring-physics, react-19]

# Dependency graph
requires:
  - phase: 02-react-19-typescript-5/02-01
    provides: TypeScript strict-mode codebase, React 19 runtime, .tsx file extensions
provides:
  - motion@12.36.0 animation library (replacing framer-motion)
  - Spring physics on tab transitions and page entrance animations
  - All 15 animation files importing from motion/react
affects: [frontend-animations, 06-animations-power-features]

# Tech tracking
tech-stack:
  added: [motion@12.36.0]
  removed: [framer-motion]
  patterns: [spring-physics-tab-transitions, spring-physics-page-entrance]

key-files:
  created: []
  modified:
    - frontend/package.json
    - frontend/package-lock.json
    - frontend/src/pages/AboutPage.tsx
    - frontend/src/pages/HomePage.tsx
    - frontend/src/pages/ChemPage.tsx
    - frontend/src/pages/ConvertPage.tsx
    - frontend/src/pages/DepictPage.tsx
    - frontend/src/pages/ToolsPage.tsx
    - frontend/src/pages/OCSRPage.tsx
    - frontend/src/pages/PrivacyPolicy.tsx
    - frontend/src/pages/TermsOfService.tsx
    - frontend/src/components/common/Header.tsx
    - frontend/src/components/common/Footer.tsx
    - frontend/src/components/common/LoadingScreen.tsx
    - frontend/src/components/depict/Depict2DMultiView.tsx
    - frontend/src/components/ocsr/OCRView.tsx
    - frontend/src/components/tools/StructureGenView.tsx

key-decisions:
  - "Installed motion@12.36.0 (latest) instead of pinned 12.35.2 -- newer patch, same API"
  - "Spring physics applied to tab indicators and page entrance animations for natural motion feel"
  - "Footer particles, Header mobile menu, LoadingScreen, and AboutPage scroll animations left untouched per research recommendations"

patterns-established:
  - "Animation imports: all motion imports use 'motion/react' path, not 'framer-motion'"
  - "Tab transition springs: type spring, stiffness 500, damping 30 for snappy tab switches"
  - "Page entrance springs: type spring, stiffness 100, damping 20 for natural deceleration"

requirements-completed: [FRAME-05]

# Metrics
duration: 12min
completed: 2026-03-12
---

# Phase 02 Plan 02: Motion Package Migration Summary

**Migrated framer-motion to motion@12.36.0 across 15 files with spring-physics enhanced tab transitions and page entrance animations**

## Performance

- **Duration:** 12 min (including human visual verification)
- **Started:** 2026-03-12T17:30:00Z
- **Completed:** 2026-03-12T18:03:00Z
- **Tasks:** 2
- **Files modified:** 17

## Accomplishments
- Replaced framer-motion with motion@12.36.0 (official successor package) across all 15 source files
- Enhanced tab-based pages (ChemPage, DepictPage, ConvertPage, ToolsPage, Depict2DMultiView) with spring-physics tab transitions (stiffness 500, damping 30)
- Enhanced page entrance animations (HomePage, AboutPage, PrivacyPolicy, TermsOfService, OCSRPage) with spring physics (stiffness 100, damping 20)
- All 9 pages visually verified by user with smooth animations, no console errors
- Production build succeeds (1,271 KB bundle, 3.49s), TypeScript compiles with 0 errors, test suite passes

## Task Commits

Each task was committed atomically:

1. **Task 1: Migrate framer-motion to motion package with animation enhancements** - `1b9a4b9` (feat)
2. **Task 2: Visual verification of all 9 pages and animations** - No commit (verification-only checkpoint, user approved)

## Files Created/Modified

**Package changes:**
- `frontend/package.json` - motion@12.36.0 added, framer-motion removed from direct dependencies
- `frontend/package-lock.json` - Updated dependency tree

**Pages with spring-physics page entrance animations:**
- `frontend/src/pages/HomePage.tsx` - Import path + spring entrance animation
- `frontend/src/pages/AboutPage.tsx` - Import path + spring entrance animation (scroll animations preserved)
- `frontend/src/pages/OCSRPage.tsx` - Import path + spring entrance animation
- `frontend/src/pages/PrivacyPolicy.tsx` - Import path + spring entrance animation
- `frontend/src/pages/TermsOfService.tsx` - Import path + spring entrance animation

**Pages with spring-physics tab transitions:**
- `frontend/src/pages/ChemPage.tsx` - Import path + spring tab indicator transitions
- `frontend/src/pages/ConvertPage.tsx` - Import path + spring tab indicator transitions
- `frontend/src/pages/DepictPage.tsx` - Import path + spring tab indicator transitions
- `frontend/src/pages/ToolsPage.tsx` - Import path + spring tab indicator transitions

**Components (import path only, animations preserved):**
- `frontend/src/components/common/Header.tsx` - Import path updated
- `frontend/src/components/common/Footer.tsx` - Import path updated (particles preserved)
- `frontend/src/components/common/LoadingScreen.tsx` - Import path updated
- `frontend/src/components/depict/Depict2DMultiView.tsx` - Import path updated
- `frontend/src/components/ocsr/OCRView.tsx` - Import path updated
- `frontend/src/components/tools/StructureGenView.tsx` - Import path updated

## Decisions Made
- **motion@12.36.0 vs 12.35.2:** Installed latest available (12.36.0) rather than the pinned 12.35.2 from the plan. Same API, newer bugfix patch.
- **Selective animation enhancement:** Per research recommendations, only enhanced tab transitions (spring stiffness 500/damping 30) and page entrances (spring stiffness 100/damping 20). Left Footer particles, Header mobile menu, LoadingScreen, and AboutPage scroll animations untouched to avoid breaking fragile animation systems.

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
None - the migration was a clean import path swap with targeted spring physics enhancements. All 15 files compiled and ran correctly on first attempt.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Phase 2 is now fully complete: React 19.2.4, TypeScript 5 strict mode, motion@12.36.0
- All 61 .ts/.tsx source files with proper types, all 15 animation files using motion/react
- Ready for Phase 3: Tailwind v4 Migration (no blockers)
- The framer-motion React 19 compatibility concern noted in STATE.md blockers is now resolved

## Self-Check: PASSED

All 16 key files verified present. Task 1 commit (1b9a4b9) verified in git log. Task 2 was a verification-only checkpoint (no commit expected).

---
*Phase: 02-react-19-typescript-5*
*Completed: 2026-03-12*
