---
phase: 04-component-system-dark-mode
plan: 05
subsystem: ui
tags: [shadcn-ui, button, card, react, framer-motion, asChild, gap-closure]

# Dependency graph
requires:
  - phase: 04-component-system-dark-mode (plan 04)
    provides: "shadcn/ui base migration, Button/Input/Select/Sheet components across most pages"
provides:
  - "Complete shadcn Button coverage across all 9 pages (zero hand-rolled button-styled elements)"
  - "Card component usage for OCSRPage info panels, HomePage feature/molecule cards, AboutPage funder logos, Footer resource links"
  - "COMP-02, COMP-03, COMP-07 fully satisfied (previously PARTIAL)"
affects: [05-loading-states-ux, 06-polish-performance]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "motion.div wrapper + Button asChild + child anchor/Link for animated button links"
    - "Card with py-0 gap-0 overrides for custom padding via CardContent"

key-files:
  created: []
  modified:
    - frontend/src/pages/OCSRPage.tsx
    - frontend/src/pages/HomePage.tsx
    - frontend/src/pages/AboutPage.tsx
    - frontend/src/components/common/Footer.tsx

key-decisions:
  - "Used motion.div wrapper pattern (not motion.create(Button)) for animated button links per research findings"
  - "Kept feature-card-enhanced class on Card element to preserve ::before pseudo-element border reveal CSS effect"
  - "Funder logos wrapped in Card (not Button) since they contain images, not button text"

patterns-established:
  - "motion.div + Button asChild: All animated external links use motion.div for animation wrapping and Button asChild for semantic button rendering"
  - "Card py-0 gap-0: Card components with custom internal padding use py-0 gap-0 overrides with CardContent for padding control"

requirements-completed: [COMP-01, COMP-02, COMP-03, COMP-04, COMP-05, COMP-06, COMP-07, THEME-01, THEME-02, THEME-03, THEME-04, THEME-05]

# Metrics
duration: 7min
completed: 2026-03-13
---

# Phase 4 Plan 5: Gap Closure Summary

**Migrated all remaining hand-rolled button-styled anchors and card-styled divs to shadcn/ui Button and Card across OCSRPage, HomePage, AboutPage, and Footer -- zero motion.a elements remain in codebase**

## Performance

- **Duration:** 7 min
- **Started:** 2026-03-13T09:16:41Z
- **Completed:** 2026-03-13T09:23:46Z
- **Tasks:** 2
- **Files modified:** 4

## Accomplishments
- OCSRPage: 2 motion.a buttons migrated to Button asChild, 3 card containers migrated to Card + CardContent
- HomePage: 3 buttons (CTA + 2 secondary) migrated to Button asChild, 5 feature cards and 3 molecule cards wrapped in Card
- AboutPage: 2 button-styled motion.a links converted to Button asChild, 3 funder logo links wrapped in Card
- Footer: 4 resource link motion.a elements converted to motion.div + Card + a
- Zero motion.a elements remain in any page or Footer component
- All existing motion animations (variants, whileHover, stagger, spring physics) preserved
- Build succeeds, all 16 tests pass

## Task Commits

Each task was committed atomically:

1. **Task 1: Migrate OCSRPage and HomePage to shadcn Button and Card** - `5c97df4` (feat)
2. **Task 2: Migrate AboutPage button-styled links and Footer resource cards** - `2f8a69c` (feat)

## Files Created/Modified
- `frontend/src/pages/OCSRPage.tsx` - Added Button + Card imports; 2 buttons + 3 cards migrated to shadcn
- `frontend/src/pages/HomePage.tsx` - Added Button + Card imports; 3 buttons + 8 cards migrated to shadcn
- `frontend/src/pages/AboutPage.tsx` - Added Card import; 2 buttons + 3 funder logos migrated to shadcn
- `frontend/src/components/common/Footer.tsx` - Added Card import; 4 resource links migrated to shadcn

## Decisions Made
- Used motion.div wrapper pattern (not motion.create(Button)) for animated button links, consistent with research recommendation to avoid Radix Slot + motion component conflicts
- Kept feature-card-enhanced CSS class on Card element to preserve the ::before pseudo-element border reveal effect on hover
- Funder logos wrapped in Card (not Button) since they are image-only links, not button-styled text links -- semantically correct

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered
None

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- All COMP-* and THEME-* requirements fully satisfied
- Phase 4 (Component System + Dark Mode) is complete with all 5 plans executed
- Ready for Phase 5 (Loading States + UX)

---
*Phase: 04-component-system-dark-mode*
*Completed: 2026-03-13*
