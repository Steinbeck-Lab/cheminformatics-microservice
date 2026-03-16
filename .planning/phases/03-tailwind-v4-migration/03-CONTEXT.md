# Phase 3: Tailwind v4 Migration - Context

**Gathered:** 2026-03-12
**Status:** Ready for planning

<domain>
## Phase Boundary

Upgrade Tailwind CSS v3 to v4 with CSS-first configuration, replace JS config and PostCSS setup with @tailwindcss/vite plugin, resolve all theme() function usages and class renames, and verify zero visual regressions across all 9 pages. No new visual features — this is a pure infrastructure migration.

</domain>

<decisions>
## Implementation Decisions

### Custom Color Palette
- Keep ALL custom color shades (gray.750/850/950, blue.350/450/550/650/750/850) in v4's @theme block
- Remove any that v4 now includes by default (e.g., gray.950 is built-in in v4)
- Claude's discretion on color format: evaluate whether to keep exact hex values or convert to OKLCH for consistency with v4's native palette

### CSS Variable Theming
- Replace all 32 `theme()` function calls in tailwind.css with v4-native CSS variables (e.g., `var(--color-slate-100)`) — Claude decides exact approach
- Inline `theme(space.12)` in AllFiltersView.tsx: replace with `var(--spacing-12)` (v4 CSS variables)
- Dark mode: Claude decides between class-based (`@custom-variant`) or media-query approach, prioritizing Phase 4 compatibility (dark mode toggle planned there)
- Keep current CSS variable names (--bg-primary, --text-primary, etc.) — do NOT pre-align with shadcn/ui naming; Phase 4 handles that transition
- Preserve exact specificity behavior of component classes (btn, card, glass) — test and fix any v4 @layer/@apply specificity changes
- Claude's discretion on: global `* { transition-colors }` rule (scope to interactive elements vs keep global), pattern CSS variables (audit usage, keep only what's used)

### Vite Plugin Setup
- Use @tailwindcss/vite plugin (standard approach) — add to vite.config.mts plugins array
- Remove postcss.config.js entirely — uninstall postcss, autoprefixer, and old tailwindcss package
- @tailwindcss/vite handles vendor prefixing internally — no separate autoprefixer needed

### @apply & Component Classes
- Claude's discretion on whether to migrate @apply component classes as-is or strip them — Phase 4 replaces btn/card/form classes with shadcn/ui
- Glass and text-gradient utility classes: Claude decides whether to keep as @apply or convert to @utility/@layer — these survive past Phase 4
- Claude's discretion on @utility vs @layer utilities syntax for custom utilities
- Claude's discretion on converting @apply in @layer base to plain CSS where v4 recommends it

### Typography Plugin
- Migrate typography plugin config (dark prose overrides) to v4's @plugin + CSS-based configuration
- Phase 4 will simplify or replace with shadcn/ui typography styles — full config preserved for now
- Keep @tailwindcss/forms plugin — migrate to @plugin directive loading
- Keep @tailwindcss/typography plugin — migrate to @plugin directive loading

### Class Rename Audit
- Full audit across ALL TSX files for v3-to-v4 breaking class changes (shadow, ring, spacing, etc.)
- Use official Tailwind upgrade tool (`npx @tailwindcss/upgrade`) as starting point, then manual review and fixes
- Zero visual regressions is the hard success criteria

### Animation Consolidation
- Consolidate and deduplicate animations: tailwind.config.js animations (pulse-slow, pulse-fast, bounce-slow) + animations.css keyframes — one source of truth
- Claude's discretion on whether to keep animations.css separate or merge into tailwind.css

### Visual Verification
- Manual 9-page verification: systematically check all pages in browser for correct colors, layouts, dark mode, form inputs, cards, typography
- Same approach as Phase 1 verification (proven effective)

### Fonts
- Claude's discretion: evaluate whether to keep Inter (sans) and JetBrains Mono (mono) font stacks in @theme or use v4 defaults, prioritizing visual continuity

### Shadows
- Claude's discretion: audit actual usage of custom shadows (inner-lg, inner-xl, soft-xl, blue-glow), keep only referenced ones

### IDE Support
- Claude's discretion: update .vscode/settings.json for Tailwind v4 IntelliSense if settings file already exists

### Claude's Discretion
- Color format (hex vs OKLCH) for custom shades
- Dark mode variant approach (class-based vs media-query) — must support Phase 4 toggle
- Global transition-colors rule scoping
- Pattern CSS variable retention (audit usage)
- @apply component class migration strategy (migrate vs strip)
- Custom utility syntax (@utility vs @layer)
- Base style @apply-to-plain-CSS conversion
- Animation file organization (separate vs merged)
- Font stack retention
- Shadow retention (usage-based audit)
- VS Code IntelliSense settings

</decisions>

<specifics>
## Specific Ideas

- "Avoid band-aids, shortcuts — this is a production system, accuracy first" (carried from Phase 1)
- "Make everything future proof" (carried from Phase 2)
- Use official Tailwind upgrade tool as starting point, then manual review — not fully manual migration
- Phase 4 dependency: dark mode variant choice must support the toggle UI planned there
- Phase 4 dependency: CSS variable names stay as current (--bg-primary etc.) — shadcn/ui naming happens in Phase 4
- Animation consolidation: one source of truth for all animation keyframes and classes

</specifics>

<code_context>
## Existing Code Insights

### Reusable Assets
- `tailwind.config.js`: Extensive custom theme (colors, spacing, shadows, fonts, animations, typography prose overrides) — all must migrate to @theme or CSS
- `postcss.config.js`: Simple CJS config (tailwindcss + autoprefixer) — will be deleted
- `tailwind.css`: 32 theme() calls, 20+ @apply component classes, CSS variables for light/dark theming, pattern definitions
- `animations.css`: 10 keyframe animations, 12 animation classes, staggered list animations — pure CSS, minimal v4 impact
- `vite.config.mts`: Already has @vitejs/plugin-react-swc — @tailwindcss/vite adds alongside

### Established Patterns
- Dark mode via `darkMode: 'class'` with `.dark` class toggling
- CSS variables (--bg-primary, --text-primary, etc.) for theme-aware component styling
- @layer base/components/utilities structure in tailwind.css
- Component classes via @apply (btn-*, card-*, form-*, glass, text-gradient)
- Plugins: @tailwindcss/forms (global form resets) + @tailwindcss/typography (prose styling)

### Integration Points
- `AllFiltersView.tsx`: 2 inline `theme(space.12)` usages that will break in v4
- `tailwind.css` imported in `index.tsx` — entry point stays the same
- `animations.css` imported in `index.tsx` — separate file, minimal migration
- `vite.config.mts` plugins array — add @tailwindcss/vite here
- `package.json` — swap tailwindcss v3 for v4, remove postcss + autoprefixer

</code_context>

<deferred>
## Deferred Ideas

- Bento grid layouts for pages — Phase 5 (UX) or Phase 6 (Power Features)
- Liquid glass visual effect — Phase 6 (Animations + Power Features)
- Claymorphism/glassmorphism combination styles — Phase 6 (Animations + Power Features)
- New animation features beyond current set — Phase 6 (Animations + Power Features)

</deferred>

---

*Phase: 03-tailwind-v4-migration*
*Context gathered: 2026-03-12*
