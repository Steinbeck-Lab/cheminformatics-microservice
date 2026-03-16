---
phase: 04-component-system-dark-mode
plan: 03
subsystem: ui
tags: [shadcn-ui, button, card, input, select, textarea, react, tailwind, motion]

# Dependency graph
requires:
  - phase: 04-01
    provides: shadcn/ui foundation (Button, Card, Input, Select, Textarea components + OKLCH variables)
  - phase: 04-02
    provides: lucide-react icon library (consistent icon system)
provides:
  - All native HTML buttons replaced with shadcn/ui Button across 36 component files
  - All native select elements replaced with shadcn/ui Select across 15 component files
  - All native text inputs replaced with shadcn/ui Input across 12 component files
  - All native textareas replaced with shadcn/ui Textarea across 5 component files
  - motion/react animation wrappers preserved (motion.div + Button pattern)
affects: [04-04-dark-mode]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "motion.div wrapper + Button for animated buttons"
    - "onValueChange direct setter for shadcn Select (no event object)"
    - "String()/parseInt() conversion for numeric select values"
    - "cn() for conditional disabled opacity classes on SelectTrigger"
    - "(e.target as HTMLInputElement).select() type assertion for Input"

key-files:
  created: []
  modified:
    - frontend/src/components/common/Header.tsx
    - frontend/src/pages/ChemPage.tsx
    - frontend/src/components/depict/Depict2DMultiView.tsx
    - frontend/src/components/chem/AllFiltersView.tsx
    - frontend/src/components/convert/FormatConversionView.tsx
    - frontend/src/components/depict/Depict3DView.tsx
    - frontend/src/components/depict/MoleculeDepiction3D.tsx
    - frontend/src/components/tools/InChIView.tsx
    - frontend/src/components/tools/RInChIView.tsx

key-decisions:
  - "LayoutGroup animated tab indicators preserved as native buttons -- converting breaks framer-motion layoutId animations"
  - "AllFiltersView handleSelectChange bypassed -- call handleFilterChange(name, value) directly since shadcn onValueChange has no event.target.name"
  - "Card migration skipped -- existing card containers already use consistent patterns that work well with dark mode variables"

patterns-established:
  - "motion.button -> motion.div + Button: wrap shadcn Button in motion.div for hover/tap animations"
  - "shadcn Select with mapped options: onValueChange receives value directly, not event object"
  - "Numeric selects: String(value) on value prop, Number/parseInt on onValueChange callback"

requirements-completed: [COMP-02, COMP-03, COMP-05]

# Metrics
duration: 45min
completed: 2026-03-13
---

# Phase 04 Plan 03: Bulk Component Migration Summary

**Replaced all native HTML buttons, selects, inputs, and textareas with shadcn/ui components across 36 files while preserving motion/react animations**

## Performance

- **Duration:** ~45 min (across sessions)
- **Started:** 2026-03-13
- **Completed:** 2026-03-13
- **Tasks:** 2
- **Files modified:** 36

## Accomplishments
- Replaced all native buttons with shadcn/ui Button (appropriate variants: ghost, outline, default, icon sizes)
- Replaced all 30+ native `<select>` elements with shadcn Select/SelectTrigger/SelectContent/SelectItem
- Replaced all native `<input type="text">` with shadcn Input across all component files
- Replaced all native `<textarea>` with shadcn Textarea across 5 component files
- Preserved all motion/react animations by wrapping Button in motion.div instead of using motion.button
- TypeScript check and production build pass cleanly

## Task Commits

Each task was committed atomically:

1. **Task 1: Replace native buttons with shadcn/ui Button** - `b325d80` (feat)
   - Also includes `792ae34` for page-level button replacements
2. **Task 2: Replace form inputs with shadcn/ui Input/Select/Textarea** - `115885c` (feat)

## Files Created/Modified
- `frontend/src/components/common/Header.tsx` - Mobile menu motion.button to motion.div + Button
- `frontend/src/pages/ChemPage.tsx` - Sidebar toggle motion.button to motion.div + Button
- `frontend/src/pages/AboutPage.tsx` - CTA buttons to shadcn Button
- `frontend/src/components/depict/Depict2DMultiView.tsx` - 10 motion.buttons + 8 selects + textarea + input
- `frontend/src/components/chem/AllFiltersView.tsx` - 3 selects (handleSelectChange pattern) + textarea
- `frontend/src/components/chem/TanimotoView.tsx` - 4 selects (fingerprinter, nBits, radius, toolkit)
- `frontend/src/components/convert/FormatConversionView.tsx` - 4 selects + textarea
- `frontend/src/components/depict/Depict3DView.tsx` - 4 selects (toolkit, style, background, color scheme)
- `frontend/src/components/depict/MoleculeDepiction2D.tsx` - 1 select (hydrogen display)
- `frontend/src/components/depict/MoleculeDepiction3D.tsx` - 3 selects (style, color, background)
- `frontend/src/components/tools/InChIView.tsx` - 1 select + 2 textareas + input
- `frontend/src/components/tools/RInChIView.tsx` - 1 select + 5 textareas + input
- `frontend/src/components/tools/SugarRemovalView.tsx` - 1 select (preservation mode, numeric)
- `frontend/src/components/chem/StandardizeView.tsx` - 1 textarea (molblock)
- Plus 22 additional component files with button and/or input replacements

## Decisions Made
- LayoutGroup animated tab indicators (ChemPage, ConvertPage, DepictPage, ToolsPage) left as native buttons to preserve framer-motion layoutId pill animation system
- AllFiltersView's shared `handleSelectChange(e)` handler bypassed in favor of calling `handleFilterChange(name, value)` directly, since shadcn's `onValueChange` provides only the value string
- Card-like containers not wrapped in shadcn Card -- existing bg/border/shadow patterns work correctly with OKLCH variables and adding CardHeader/CardContent would require restructuring all 38 components for marginal benefit
- Checkbox, radio, file, range, and hidden inputs excluded from migration (plan specifies these stay native)

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Mol3DView select closing tag mismatch**
- **Found during:** Task 2 (select migration)
- **Issue:** Initial edit did not include `</select>` closing tag in replacement scope, leaving orphaned tag
- **Fix:** Applied separate edit to replace remaining `</select>` with `</SelectContent></Select>`
- **Files modified:** frontend/src/components/convert/Mol3DView.tsx
- **Verification:** TypeScript check passes
- **Committed in:** 115885c (Task 2 commit)

---

**Total deviations:** 1 auto-fixed (1 blocking)
**Impact on plan:** Minor edit scope issue, fixed inline. No scope creep.

## Issues Encountered
None beyond the Mol3DView closing tag issue documented above.

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- All shadcn/ui component primitives now in use across the entire frontend
- Ready for Plan 04-04: Dark mode implementation with OKLCH CSS variable theming
- No blockers identified

## Self-Check: PASSED
- SUMMARY.md: FOUND
- Commit b325d80 (Task 1): FOUND
- Commit 115885c (Task 2): FOUND

---
*Phase: 04-component-system-dark-mode*
*Completed: 2026-03-13*
