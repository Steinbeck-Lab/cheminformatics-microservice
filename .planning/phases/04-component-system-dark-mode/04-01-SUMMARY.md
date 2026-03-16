---
phase: 04-component-system-dark-mode
plan: 01
subsystem: ui
tags: [shadcn-ui, tailwind-css, oklch, dark-mode, fouc-prevention, lucide-react, class-variance-authority]

# Dependency graph
requires:
  - phase: 03-tailwind-v4-migration
    provides: "Tailwind v4 CSS-first architecture with custom variants"
provides:
  - "shadcn/ui component infrastructure (components.json, cn(), 8 UI components)"
  - "OKLCH CSS variables for light/dark themes mapped to Tailwind utilities"
  - "FOUC prevention blocking script in index.html"
  - "System preference detection in AppContext"
  - "Wave 0 test stubs for all COMP-* and THEME-* requirements"
affects: [04-component-system-dark-mode, 05-ux-polish]

# Tech tracking
tech-stack:
  added: [lucide-react, class-variance-authority, clsx, tailwind-merge, tw-animate-css, shadcn, radix-ui]
  patterns: [shadcn-ui-component-system, oklch-css-variables, fouc-prevention-script, system-preference-detection]

key-files:
  created:
    - frontend/components.json
    - frontend/src/lib/utils.ts
    - frontend/src/components/ui/button.tsx
    - frontend/src/components/ui/card.tsx
    - frontend/src/components/ui/dialog.tsx
    - frontend/src/components/ui/sheet.tsx
    - frontend/src/components/ui/input.tsx
    - frontend/src/components/ui/select.tsx
    - frontend/src/components/ui/textarea.tsx
    - frontend/src/components/ui/tabs.tsx
    - frontend/src/__tests__/utils.test.ts
    - frontend/src/__tests__/components/button.test.tsx
    - frontend/src/__tests__/components/card.test.tsx
    - frontend/src/__tests__/components/dialog.test.tsx
    - frontend/src/__tests__/components/input.test.tsx
    - frontend/src/__tests__/icon-migration.test.ts
    - frontend/src/__tests__/context/theme.test.tsx
  modified:
    - frontend/src/styles/tailwind.css
    - frontend/index.html
    - frontend/src/context/AppContext.tsx
    - frontend/tsconfig.app.json
    - frontend/vite.config.mts
    - frontend/package.json
    - frontend/src/pages/ChemPage.tsx
    - frontend/src/pages/OCSRPage.tsx
    - frontend/src/pages/ConvertPage.tsx
    - frontend/src/pages/DepictPage.tsx
    - frontend/src/pages/ToolsPage.tsx
    - frontend/src/components/tools/StructureGenView.tsx
    - frontend/src/__tests__/setup.ts

key-decisions:
  - "Used force-add for src/lib/utils.ts since root .gitignore has lib/ pattern"
  - "shadcn CLI created files in literal @/ directory -- manually moved to src/components/ui/"
  - "Restored react-icons and @fortawesome after shadcn CLI unexpectedly removed them"
  - "Added matchMedia mock to test setup.ts for jsdom compatibility"
  - "Removed inline CSS variable definitions from ConvertPage, DepictPage, ToolsPage style blocks"
  - "Replaced btn/btn-primary/btn-secondary classes in StructureGenView with inline Tailwind utilities"

patterns-established:
  - "shadcn/ui components: use @/components/ui/ path with cn() for class merging"
  - "CSS variables: all theme colors use shadcn standard OKLCH naming (--background, --foreground, --primary, etc.)"
  - "FOUC prevention: blocking ES5 script in index.html reads localStorage darkMode key"
  - "System preference: AppContext detects prefers-color-scheme when no localStorage value exists"

requirements-completed: [COMP-01, THEME-02, THEME-03, THEME-04]

# Metrics
duration: 18min
completed: 2026-03-12
---

# Phase 04 Plan 01: shadcn/ui Foundation Summary

**shadcn/ui component system with OKLCH CSS variables, FOUC-free dark mode, system preference detection, and 7 Wave 0 test stubs**

## Performance

- **Duration:** 18 min
- **Started:** 2026-03-12T21:12:18Z
- **Completed:** 2026-03-12T21:30:22Z
- **Tasks:** 3 (+ 1 deviation fix)
- **Files modified:** 29 (17 created, 12 modified)

## Accomplishments
- shadcn/ui infrastructure fully configured: components.json, cn() utility, @/ path alias in both Vite and Vitest
- 8 shadcn/ui components generated: button, card, dialog, sheet, input, select, textarea, tabs
- CSS variables rewritten from custom names to shadcn/ui standard OKLCH format with light/dark palettes
- FOUC prevention script added to index.html -- prevents flash of wrong theme on page load
- AppContext enhanced with system preference detection via matchMedia for first-time visitors
- 7 Wave 0 test stub files created covering all COMP-* and THEME-* requirements
- All old CSS variable references (--bg-primary, --text-primary, etc.) removed from 6 TSX files
- All @utility blocks for btn, card, form replaced -- shadcn/ui components handle these now

## Task Commits

Each task was committed atomically:

1. **Task 0: Wave 0 test stubs** - `9200ae2` (test)
2. **Task 1: shadcn/ui dependencies + path aliases + cn()** - `7724b13` (feat)
3. **Fix: restore react-icons/fortawesome** - `24ffe0b` (fix)
4. **Task 2: CSS rewrite + FOUC + AppContext** - `1d0ba77` (feat)

## Files Created/Modified

### Created
- `frontend/components.json` - shadcn/ui configuration (new-york style, rsc: false, slate base)
- `frontend/src/lib/utils.ts` - cn() utility (clsx + tailwind-merge)
- `frontend/src/components/ui/button.tsx` - shadcn Button with variant system
- `frontend/src/components/ui/card.tsx` - shadcn Card with header/content/footer
- `frontend/src/components/ui/dialog.tsx` - shadcn Dialog (Radix-based)
- `frontend/src/components/ui/sheet.tsx` - shadcn Sheet (side panel)
- `frontend/src/components/ui/input.tsx` - shadcn Input
- `frontend/src/components/ui/select.tsx` - shadcn Select
- `frontend/src/components/ui/textarea.tsx` - shadcn Textarea
- `frontend/src/components/ui/tabs.tsx` - shadcn Tabs
- `frontend/src/__tests__/utils.test.ts` - cn() tests (3 real assertions)
- `frontend/src/__tests__/components/button.test.tsx` - Button render test
- `frontend/src/__tests__/components/card.test.tsx` - Card render test
- `frontend/src/__tests__/components/dialog.test.tsx` - Dialog stubs
- `frontend/src/__tests__/components/input.test.tsx` - Input render test
- `frontend/src/__tests__/icon-migration.test.ts` - Icon migration stubs
- `frontend/src/__tests__/context/theme.test.tsx` - Theme tests (THEME-02, THEME-03 real assertions)

### Modified
- `frontend/src/styles/tailwind.css` - Complete rewrite to OKLCH variables
- `frontend/index.html` - FOUC prevention script added
- `frontend/src/context/AppContext.tsx` - System preference detection added
- `frontend/tsconfig.app.json` - @/ path alias added
- `frontend/vite.config.mts` - @/ resolve alias for Vite and Vitest
- `frontend/package.json` - New dependencies installed
- `frontend/src/__tests__/setup.ts` - matchMedia mock for jsdom
- `frontend/src/pages/ChemPage.tsx` - Old CSS vars replaced
- `frontend/src/pages/OCSRPage.tsx` - Old CSS vars replaced
- `frontend/src/pages/ConvertPage.tsx` - Old CSS vars + inline styles removed
- `frontend/src/pages/DepictPage.tsx` - Old CSS vars + inline styles removed
- `frontend/src/pages/ToolsPage.tsx` - Old CSS vars + inline styles removed
- `frontend/src/components/tools/StructureGenView.tsx` - Old CSS vars + btn classes replaced

## Decisions Made
- Used `git add -f` for `src/lib/utils.ts` because root `.gitignore` contains `lib/` pattern
- shadcn CLI generated components into a literal `@/` directory instead of resolving the path alias -- manually moved files to `src/components/ui/`
- Restored react-icons and @fortawesome after shadcn CLI unexpectedly removed them from package.json (still needed by source files until Plan 04-02)
- Added global `window.matchMedia` mock to test setup.ts since jsdom doesn't implement it
- Removed redundant inline `<style>` CSS variable definitions from ConvertPage, DepictPage, ToolsPage (now centralized in tailwind.css)
- Changed `isDarkMode` default from `true` to `false` since the FOUC prevention script handles the initial dark class

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] shadcn CLI generated files in wrong directory**
- **Found during:** Task 1 (component generation)
- **Issue:** `npx shadcn@latest add` created files in a literal `@/` directory instead of resolving the path alias
- **Fix:** Manually moved all 8 component files from `frontend/@/components/ui/` to `frontend/src/components/ui/`
- **Files modified:** 8 component files relocated
- **Committed in:** 7724b13

**2. [Rule 3 - Blocking] shadcn CLI removed react-icons and @fortawesome from package.json**
- **Found during:** Task 1 (after component generation)
- **Issue:** shadcn unexpectedly removed react-icons and @fortawesome from dependencies even though source files still import from them
- **Fix:** Re-installed react-icons@^4.12.0, @fortawesome/free-solid-svg-icons@^6.7.2, @fortawesome/react-fontawesome@^0.2.2
- **Files modified:** package.json, package-lock.json
- **Committed in:** 24ffe0b

**3. [Rule 3 - Blocking] jsdom lacks window.matchMedia implementation**
- **Found during:** Task 2 (theme tests)
- **Issue:** After adding matchMedia call to AppContext, all tests using AppProvider failed with "window.matchMedia is not a function"
- **Fix:** Added global matchMedia mock to test setup.ts
- **Files modified:** frontend/src/__tests__/setup.ts
- **Committed in:** 1d0ba77

---

**Total deviations:** 3 auto-fixed (all Rule 3 - blocking issues)
**Impact on plan:** All fixes necessary for correct operation. No scope creep.

## Issues Encountered
- Lint-staged prettier reverted tailwind.css changes during one commit attempt (stash/restore cycle) -- rewrote the file and committed successfully on retry
- ESLint flagged constant boolean expressions in test assertions (`false && "hidden"`) -- fixed by using variables

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- shadcn/ui component infrastructure complete -- Plan 04-02 can proceed with icon migration
- All 8 UI components generated and type-checked
- CSS variables use shadcn/ui standard naming -- components will pick up theme automatically
- Wave 0 test stubs in place for all requirements
- Build, typecheck, and tests all pass (7 test files, 16 passing, 13 todo)

## Self-Check: PASSED

- All 17 created files: FOUND
- All 13 modified files: FOUND
- All 4 commits (9200ae2, 7724b13, 24ffe0b, 1d0ba77): FOUND

---
*Phase: 04-component-system-dark-mode*
*Completed: 2026-03-12*
