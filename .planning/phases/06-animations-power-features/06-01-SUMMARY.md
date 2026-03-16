---
phase: 06-animations-power-features
plan: 01
subsystem: ui
tags: [motion, css-animations, page-transitions, button-interactions, framer-motion, tailwind-v4]

# Dependency graph
requires:
  - phase: 04-component-system-dark-mode
    provides: Button component with cva variants, clay-btn/clay-pressed CSS utilities
  - phase: 04.1-visual-design-system
    provides: glass-bold utility, GradientMesh, tab page structure
provides:
  - btn-hover-lift CSS utility for hover lift + glow + press animations
  - focus-ring-animate CSS utility for animated focus ring transitions
  - AnimatedOutlet component for spring-physics page transitions
  - Simplified 150ms cross-fade contentVariants for all tab pages
  - staggerContainer/staggerItem animation variants for list animations
affects: [06-02, 06-03]

# Tech tracking
tech-stack:
  added: []
  patterns: [AnimatePresence mode=sync page transitions, CSS @utility for interactive animations, layout prop for automatic layout shift animations]

key-files:
  created:
    - frontend/src/components/common/AnimatedOutlet.tsx
    - frontend/src/__tests__/components/button-animations.test.tsx
    - frontend/src/__tests__/components/page-transitions.test.tsx
  modified:
    - frontend/src/styles/tailwind.css
    - frontend/src/components/ui/button.tsx
    - frontend/src/App.tsx
    - frontend/src/pages/ChemPage.tsx
    - frontend/src/pages/ConvertPage.tsx
    - frontend/src/pages/DepictPage.tsx
    - frontend/src/pages/ToolsPage.tsx

key-decisions:
  - "btn-hover-lift and focus-ring-animate implemented as @utility directives (not @layer components) for Tailwind v4 compatibility"
  - "AnimatedOutlet uses mode=sync per user decision with absolute positioning on exit to prevent vertical stacking"
  - "Tab content simplified from spring+scale/x-slide to pure 150ms opacity cross-fade per ANIM-04 user decision"

patterns-established:
  - "CSS @utility for interactive animation: btn-hover-lift combines transform + box-shadow transitions"
  - "focus-ring-animate: animated outline-based focus ring replacing static Tailwind ring classes"
  - "AnimatePresence mode=sync + absolute exit positioning: pattern for overlap page transitions"
  - "layout prop on content motion.divs: enables automatic layout shift animations"

requirements-completed: [ANIM-01, ANIM-02, ANIM-03, ANIM-04, ANIM-05, ANIM-06]

# Metrics
duration: 4min
completed: 2026-03-13
---

# Phase 06 Plan 01: Animation Polish Summary

**Page transitions with spring physics via AnimatedOutlet, button hover lift+glow+press via CSS utilities, animated focus rings, 150ms tab cross-fades, stagger variants, and layout animations across all 4 tab pages**

## Performance

- **Duration:** 4 min
- **Started:** 2026-03-13T22:19:35Z
- **Completed:** 2026-03-13T22:24:23Z
- **Tasks:** 3
- **Files modified:** 10

## Accomplishments
- Button hover shows translateY(-2px) lift with primary-colored glow shadow, press shows scale(0.97) -- all via pure CSS @utility (no JS)
- Focus rings animate in/out smoothly over 200ms using outline transitions instead of static Tailwind ring classes
- Page navigation triggers cross-fade with subtle vertical rise animation using spring physics (stiffness 100, damping 20) with sync overlap
- All 4 tab pages (Chem, Convert, Depict, Tools) simplified from spring+scale/x-slide to clean 150ms opacity cross-fade
- Stagger animation variants defined in all pages for future list animation use
- Layout prop added to content motion.divs for automatic expand/collapse animations
- 11 new tests (6 button + 5 page transition) all passing

## Task Commits

Each task was committed atomically:

1. **Task 1: CSS animation utilities + Button component update** - `ab2f1cd` (feat)
2. **Task 2: AnimatedOutlet page transitions + App.tsx integration** - `f2e783a` (feat)
3. **Task 3: Tab content cross-fade + stagger/layout animations in all 4 tab pages** - `6fc6a6c` (feat)

## Files Created/Modified
- `frontend/src/styles/tailwind.css` - Added btn-hover-lift and focus-ring-animate @utility blocks
- `frontend/src/components/ui/button.tsx` - Replaced hover:scale/active:scale with btn-hover-lift focus-ring-animate
- `frontend/src/components/common/AnimatedOutlet.tsx` - New: page transition wrapper with AnimatePresence mode=sync
- `frontend/src/App.tsx` - Replaced bare Outlet with AnimatedOutlet in Layout
- `frontend/src/pages/ChemPage.tsx` - Simplified contentVariants, added stagger variants and layout prop
- `frontend/src/pages/ConvertPage.tsx` - Simplified contentVariants, added stagger variants and layout prop
- `frontend/src/pages/DepictPage.tsx` - Simplified contentVariants, added stagger variants and layout prop
- `frontend/src/pages/ToolsPage.tsx` - Simplified contentVariants, added stagger variants and layout prop
- `frontend/src/__tests__/components/button-animations.test.tsx` - New: 6 tests for button animation classes
- `frontend/src/__tests__/components/page-transitions.test.tsx` - New: 5 tests for AnimatedOutlet

## Decisions Made
- btn-hover-lift and focus-ring-animate use @utility (not @layer components) since they need to work as single-class utilities in Tailwind v4
- AnimatedOutlet uses mode="sync" per user decision in CONTEXT.md for overlap transitions; exiting page gets absolute positioning to prevent two-page vertical stacking
- Tab content cross-fade simplified to pure opacity 150ms per ANIM-04 user decision (removed spring, scale, and x-axis movement)
- Ghost and link button variants preserve hover:scale-100/active:scale-100 to neutralize btn-hover-lift transforms

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness
- Animation foundation complete with CSS utilities, page transitions, and tab cross-fades
- staggerContainer/staggerItem variants ready for use in Plans 02 and 03
- layout prop enables automatic layout animations for content size changes

## Self-Check: PASSED

- All 10 files verified present on disk
- All 3 task commits verified in git history (ab2f1cd, f2e783a, 6fc6a6c)
- Production build succeeds
- All 11 new tests pass (6 button + 5 page transition)

---
*Phase: 06-animations-power-features*
*Completed: 2026-03-13*
