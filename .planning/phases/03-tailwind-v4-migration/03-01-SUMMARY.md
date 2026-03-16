---
phase: 03-tailwind-v4-migration
plan: 01
subsystem: ui
tags: [tailwindcss, vite, css-variables, css-first-config, postcss-removal]

# Dependency graph
requires:
  - phase: 02-react-19-typescript-5
    provides: TypeScript config, Vite build pipeline, React 19 components
provides:
  - Tailwind CSS v4 with @tailwindcss/vite plugin
  - CSS-first configuration (@theme, @plugin, @custom-variant, @utility)
  - All v3-to-v4 class renames applied across 50 TSX files
  - Custom theme values preserved (colors, fonts, animations)
  - Legacy JS config and PostCSS pipeline removed
affects: [03-02-PLAN, 04-component-system]

# Tech tracking
tech-stack:
  added: [tailwindcss@4.2.1, @tailwindcss/vite@4.2.1, @tailwindcss/forms@0.5.11, @tailwindcss/typography@0.5.19]
  patterns: [CSS-first config via @theme/@plugin/@custom-variant, @utility for custom utilities, var(--color-*) instead of theme()]

key-files:
  created: []
  modified:
    - frontend/vite.config.mts
    - frontend/src/styles/tailwind.css
    - frontend/src/styles/animations.css
    - frontend/package.json

key-decisions:
  - "Used @utility directive for component classes (btn, card, form) because v4 @layer components cannot use @apply to reference same-layer classes"
  - "Kept component classes as @utility rather than inlining -- Phase 4 replaces them with shadcn/ui anyway"
  - "Scoped global * transition-colors to interactive elements only (a, button, input, select, textarea, [role=button])"
  - "Removed all unused pattern definitions and classes (dots-bg, grid-bg, mesh-bg, noise-bg) -- zero TSX references"
  - "Dropped all custom shadows from @theme (inner-lg, inner-xl, soft-xl, blue-glow) -- zero TSX usage"
  - "Kept custom gray-950 (#0f1521) to preserve visual continuity -- differs from v4 default #030712"
  - "glass and text-gradient converted to @utility with dark mode via :where() selector"

patterns-established:
  - "@theme block: custom colors, fonts, animations defined in CSS not JS"
  - "@plugin directive: plugins loaded in CSS via @plugin instead of require() in JS"
  - "@custom-variant: dark mode via class-based &:where(.dark, .dark *)"
  - "@utility: custom single-class utilities with dark mode via :where() nesting"
  - "var(--color-*): all theme references use CSS variables, not theme() function"

requirements-completed: [FRAME-02, FRAME-04]

# Metrics
duration: 19min
completed: 2026-03-12
---

# Phase 3 Plan 01: Tailwind v4 Migration Summary

**Tailwind CSS v4 with CSS-first config (@theme, @plugin, @custom-variant), @tailwindcss/vite plugin, automated class renames across 50 files, and all 32 theme() calls replaced with var()**

## Performance

- **Duration:** 19 min
- **Started:** 2026-03-12T18:59:13Z
- **Completed:** 2026-03-12T19:18:29Z
- **Tasks:** 2
- **Files modified:** 54

## Accomplishments
- Tailwind CSS v4 installed with @tailwindcss/vite plugin replacing entire PostCSS pipeline
- All 50 TSX files class-renamed by upgrade tool (bg-gradient-to -> bg-linear-to, outline-none -> outline-hidden, shadow shifts, bg-opacity -> slash modifier, flex-shrink -> shrink, etc.)
- CSS-first configuration with @theme (custom colors, fonts, animations), @plugin (forms, typography), @custom-variant (class-based dark mode)
- All 32 theme() calls replaced with var(--color-*) CSS variables
- Legacy tailwind.config.js and postcss.config.js deleted
- glass and text-gradient utilities converted to @utility with :where() dark mode
- Unused pattern definitions and custom shadows removed
- Build passes, all tests pass

## Task Commits

Each task was committed atomically:

1. **Task 1: Run upgrade tool and swap packages** - `98c26fe` (feat)
2. **Task 2: Refine CSS configuration and consolidate** - `e531d83` (feat)

## Files Created/Modified
- `frontend/vite.config.mts` - Added @tailwindcss/vite plugin import and registration
- `frontend/src/styles/tailwind.css` - Complete CSS-first Tailwind v4 configuration
- `frontend/src/styles/animations.css` - Removed duplicate bounce-slow (now in @theme)
- `frontend/package.json` - Updated tailwindcss v4, added @tailwindcss/vite, removed postcss/autoprefixer
- `frontend/package-lock.json` - Updated dependency tree
- `frontend/tailwind.config.js` - DELETED (replaced by @theme in CSS)
- `frontend/postcss.config.js` - DELETED (replaced by @tailwindcss/vite)
- `frontend/src/components/chem/AllFiltersView.tsx` - Fixed theme(space.12) -> var(--spacing-12)
- 47 additional TSX files - Automated class renames by upgrade tool

## Decisions Made
- Used @utility directive for component classes (btn, card, form) because v4's @layer components cannot use @apply to reference same-layer classes -- this is a v4 requirement, not a design choice
- Kept custom gray-950 (#0f1521) in @theme to preserve visual continuity -- v4's built-in gray-950 is #030712 which would be a visible regression on the dark mode background
- Scoped global `* { transition-colors }` to interactive elements only -- was causing unnecessary transitions on all DOM elements
- Removed all pattern definitions (dots-bg, grid-bg, mesh-bg, noise-bg) and custom shadows -- confirmed zero usage in TSX files
- Converted glass and text-gradient to @utility with dark mode via `&:where(.dark, .dark *)` -- survives past Phase 4

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Node 20 required for native bindings**
- **Found during:** Task 1 (build verification)
- **Issue:** npm install ran on Node 18 but @tailwindcss/oxide requires Node 20+ native bindings. Build failed with "Cannot find native binding" error
- **Fix:** Clean reinstall with `rm -rf node_modules && npm install` under Node 20 via nvm
- **Files modified:** frontend/node_modules/ (reinstall)
- **Verification:** Build succeeds with Node 20
- **Committed in:** 98c26fe (Task 1 commit)

**2. [Rule 1 - Bug] Upgrade tool created broken @utility dark block**
- **Found during:** Task 2 (CSS refinement)
- **Issue:** The upgrade tool converted `.dark .theme-toggle-button:hover` into a separate `@utility dark` block which is semantically broken -- "dark" is not a valid custom utility name and it conflicted with the dark mode variant
- **Fix:** Moved theme-toggle-button styles to @layer components where complex selectors are supported
- **Files modified:** frontend/src/styles/tailwind.css
- **Verification:** Build succeeds, theme-toggle-button styles preserved
- **Committed in:** e531d83 (Task 2 commit)

**3. [Rule 1 - Bug] Upgrade tool couldn't auto-convert JS config**
- **Found during:** Task 1 (upgrade tool execution)
- **Issue:** The upgrade tool added `@config '../../tailwind.config.js'` fallback instead of converting the JS config to CSS, because the config was too complex to auto-convert
- **Fix:** Manually deleted the @config reference and created the @theme block with all custom values from the old tailwind.config.js
- **Files modified:** frontend/src/styles/tailwind.css
- **Verification:** Build succeeds with CSS-first config only
- **Committed in:** 98c26fe (removed @config), e531d83 (added @theme block)

---

**Total deviations:** 3 auto-fixed (2 bugs, 1 blocking)
**Impact on plan:** All auto-fixes necessary for correct v4 operation. No scope creep.

## Issues Encountered
- Upgrade tool's template migration step took extremely long (10+ minutes CPU) on the codebase -- possibly scanning node_modules. The tool was killed after the template migration had already completed.
- The upgrade tool converted component classes to @utility which is actually the correct v4 approach (not @layer components as the plan suggested) -- @layer components can't use @apply to reference same-layer classes in v4

## User Setup Required
None - no external service configuration required.

## Next Phase Readiness
- Tailwind v4 CSS pipeline fully operational
- CSS-first configuration in place -- ready for Phase 4 shadcn/ui integration
- Plan 03-02 can now perform edge case audit and 9-page visual verification
- Note: Local development requires Node 20+ (nvm use 20) for Tailwind v4's Oxide engine

---
*Phase: 03-tailwind-v4-migration*
*Completed: 2026-03-12*
