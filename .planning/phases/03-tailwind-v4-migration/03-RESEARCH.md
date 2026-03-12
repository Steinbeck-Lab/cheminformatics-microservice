# Phase 3: Tailwind v4 Migration - Research

**Researched:** 2026-03-12
**Domain:** Tailwind CSS v3-to-v4 migration, CSS-first configuration, Vite plugin integration
**Confidence:** HIGH

## Summary

Tailwind CSS v4 is a major rewrite with a CSS-first configuration model that replaces the JavaScript `tailwind.config.js` with CSS directives (`@theme`, `@plugin`, `@custom-variant`, `@utility`). The `@tailwind base/components/utilities` directives are replaced with a single `@import "tailwindcss"`, and the `theme()` function is replaced with native CSS variables (e.g., `var(--color-slate-500)`). The Vite plugin `@tailwindcss/vite` replaces PostCSS-based processing entirely.

This project has **significant migration surface**: 32 `theme()` calls in tailwind.css, 115 `bg-gradient-to-*` usages across 24 files (renamed to `bg-linear-to-*`), 62 `bg-opacity-*` usages across 21 files (removed, replaced with slash modifier), ~95 `flex-shrink-0`/`flex-grow` usages across 40+ files (renamed to `shrink-0`/`grow`), 344 shadow utility usages across 44 files (naming shifted down one size), and 115 `outline-none` usages across 37 files (renamed to `outline-hidden`). The `@tailwindcss/upgrade` CLI tool handles ~90% of class renames automatically across TSX files, but the tool requires Node.js 20+ and the project currently runs Node.js 18.20.8.

**Primary recommendation:** Run `npx @tailwindcss/upgrade` (after ensuring Node 20+ is available) as the first step to handle mechanical renames, then manually migrate the CSS file (`tailwind.css`) theme configuration, plugin loading, and custom component/utility classes. Follow up with thorough manual verification of all 9 pages.

<user_constraints>

## User Constraints (from CONTEXT.md)

### Locked Decisions
- Keep ALL custom color shades (gray.750/850/950, blue.350/450/550/650/750/850) in v4's @theme block
- Remove any that v4 now includes by default (e.g., gray.950 is built-in in v4)
- Replace all 32 `theme()` function calls in tailwind.css with v4-native CSS variables
- Inline `theme(space.12)` in AllFiltersView.tsx: replace with `var(--spacing-12)`
- Keep current CSS variable names (--bg-primary, --text-primary, etc.) -- do NOT pre-align with shadcn/ui naming
- Preserve exact specificity behavior of component classes (btn, card, glass)
- Use @tailwindcss/vite plugin -- add to vite.config.mts plugins array
- Remove postcss.config.js entirely -- uninstall postcss, autoprefixer, and old tailwindcss package
- Migrate typography plugin config (dark prose overrides) to v4's @plugin + CSS-based configuration
- Keep @tailwindcss/forms plugin -- migrate to @plugin directive loading
- Keep @tailwindcss/typography plugin -- migrate to @plugin directive loading
- Full audit across ALL TSX files for v3-to-v4 breaking class changes
- Use official Tailwind upgrade tool (`npx @tailwindcss/upgrade`) as starting point, then manual review
- Zero visual regressions is the hard success criteria
- Consolidate and deduplicate animations from tailwind.config.js + animations.css
- Manual 9-page verification

### Claude's Discretion
- Color format (hex vs OKLCH) for custom shades
- Dark mode variant approach (class-based vs media-query) -- must support Phase 4 toggle
- Global transition-colors rule scoping
- Pattern CSS variable retention (audit usage)
- @apply component class migration strategy (migrate vs strip)
- Custom utility syntax (@utility vs @layer)
- Base style @apply-to-plain-CSS conversion
- Animation file organization (separate vs merged)
- Font stack retention
- Shadow retention (usage-based audit)
- VS Code IntelliSense settings
- Glass and text-gradient utility classes: keep as @apply or convert to @utility/@layer

### Deferred Ideas (OUT OF SCOPE)
- Bento grid layouts for pages -- Phase 5 (UX) or Phase 6 (Power Features)
- Liquid glass visual effect -- Phase 6 (Animations + Power Features)
- Claymorphism/glassmorphism combination styles -- Phase 6 (Animations + Power Features)
- New animation features beyond current set -- Phase 6 (Animations + Power Features)

</user_constraints>

<phase_requirements>

## Phase Requirements

| ID | Description | Research Support |
|----|-------------|-----------------|
| FRAME-02 | Tailwind CSS upgraded from v3.3 to v4.x with CSS-first config | Full migration path documented: @tailwindcss/vite plugin, @import directive, @theme block, @plugin directive, removal of JS config and PostCSS. Standard stack and architecture patterns sections cover this comprehensively. |
| FRAME-04 | All Tailwind v3 class renames and breaking changes resolved | Complete class rename mapping documented with exact counts per category (shadow, gradient, opacity, flex-shrink, outline). Upgrade tool handles ~90% automatically; remaining manual fixes itemized in pitfalls section. |

</phase_requirements>

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| tailwindcss | ^4.2.1 | CSS framework (v4, CSS-first config) | Latest stable; Oxide engine (Rust), 2-5x faster builds |
| @tailwindcss/vite | ^4.2.1 | Vite integration plugin | Official Vite plugin, replaces PostCSS pipeline entirely |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| @tailwindcss/forms | ^0.5.11 | Form element base resets | Already used; load via `@plugin` directive in CSS |
| @tailwindcss/typography | ^0.5.19 | Prose/typography styling | Already used; load via `@plugin` directive in CSS |

### Packages to Remove
| Package | Reason |
|---------|--------|
| tailwindcss (v3) | Replaced by tailwindcss v4 |
| postcss | No longer needed -- @tailwindcss/vite handles everything |
| autoprefixer | No longer needed -- @tailwindcss/vite handles vendor prefixing |

### Files to Delete
| File | Reason |
|------|--------|
| `tailwind.config.js` | Replaced by CSS-first @theme/@plugin directives |
| `postcss.config.js` | No longer needed with @tailwindcss/vite |

**Installation:**
```bash
cd frontend
npm uninstall tailwindcss postcss autoprefixer
npm install -D tailwindcss@latest @tailwindcss/vite@latest @tailwindcss/forms@latest @tailwindcss/typography@latest
```

## Architecture Patterns

### CSS Entry Point Migration (tailwind.css)

**v3 (current):**
```css
@tailwind base;
@tailwind components;
@tailwind utilities;

@layer base { ... }
@layer components { ... }
@layer utilities { ... }
```

**v4 (target):**
```css
@import "tailwindcss";

@plugin "@tailwindcss/forms";
@plugin "@tailwindcss/typography";

@custom-variant dark (&:where(.dark, .dark *));

@theme {
  --color-gray-750: #2b3544;
  --color-gray-850: #1a2231;
  --color-blue-350: #7dabf8;
  --color-blue-450: #4f85e6;
  --color-blue-550: #3a72d6;
  --color-blue-650: #2a5db8;
  --color-blue-750: #214a98;
  --color-blue-850: #183878;

  --font-sans: 'Inter', ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, 'Noto Sans', sans-serif;
  --font-mono: 'JetBrains Mono', ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, 'Liberation Mono', 'Courier New', monospace;

  --animate-pulse-slow: pulse 3s cubic-bezier(0.4, 0, 0.6, 1) infinite;
  --animate-pulse-fast: pulse 1s cubic-bezier(0.4, 0, 0.6, 1) infinite;
  --animate-bounce-slow: bounce 2s infinite;
}

@layer base {
  :root { /* CSS variables for theming */ }
  .dark { /* Dark mode CSS variables */ }
  body { ... }
}

/* Component classes via @layer components or @utility */
/* Custom utilities via @utility */
```

### Vite Configuration

```typescript
// vite.config.mts
import { defineConfig } from "vite";
import react from "@vitejs/plugin-react-swc";
import tailwindcss from "@tailwindcss/vite";

export default defineConfig({
  plugins: [
    tailwindcss(),
    react({
      jsxRuntime: "automatic",
    }),
  ],
  // ... rest of config
});
```

### theme() Replacement Pattern

Every `theme("colors.slate.100")` becomes `var(--color-slate-100)`:

```css
/* v3 */
--bg-primary: theme("colors.slate.100");
--text-accent: theme("colors.sky.700");

/* v4 */
--bg-primary: var(--color-slate-100);
--text-accent: var(--color-sky-700);
```

For `theme("colors.white")` and `theme("colors.gray.950")` (custom color):
```css
/* v4 */
--input-bg: var(--color-white);
--bg-primary: var(--color-gray-950); /* custom color in @theme */
```

### Dark Mode Configuration

**Recommendation: Use class-based dark mode via @custom-variant.**

Phase 4 requires a dark mode toggle UI, which means class-based switching (not media-query). Use:

```css
@custom-variant dark (&:where(.dark, .dark *));
```

This replaces `darkMode: 'class'` from the JS config. The `.dark` class on `<html>` continues to work exactly as before. The `dark:` variant in all TSX classes remains unchanged.

### Plugin Loading

```css
/* v4 CSS-first plugin loading */
@plugin "@tailwindcss/forms";
@plugin "@tailwindcss/typography";
```

The typography plugin's prose dark mode overrides from the JS config (custom prose colors for headings, links, code, etc.) cannot be passed via `@plugin`. Two approaches:

1. **Use `prose-invert` class** for dark mode (built-in to typography plugin) -- this provides a complete dark mode for prose content. Apply `dark:prose-invert` alongside `prose`.
2. **Override with CSS custom properties** -- the typography plugin exposes CSS variables like `--tw-prose-body`, `--tw-prose-headings`, `--tw-prose-links` that can be overridden in a `.dark .prose` rule.

**Recommendation:** Use `prose dark:prose-invert` as the standard approach. The existing custom prose overrides in tailwind.config.js are dark-theme-specific colors that `prose-invert` already handles. If exact color matching is needed, override with CSS variables.

### bg-opacity Migration Pattern

The codebase has a consistent pattern of `dark:bg-red-900 dark:bg-opacity-30` (62 occurrences across 21 files). The upgrade tool converts these to the slash modifier:

```html
<!-- v3 -->
<div class="dark:bg-red-900 dark:bg-opacity-30">

<!-- v4 -->
<div class="dark:bg-red-900/30">
```

**Important:** The `bg-opacity-*` pattern is used ONLY in conjunction with `dark:` variants in this codebase, always following a `dark:bg-{color}-900` class. The upgrade tool should handle this automatically.

### Class Rename Reference

| v3 Class | v4 Class | Count in Codebase | Files |
|----------|----------|-------------------|-------|
| `bg-gradient-to-*` | `bg-linear-to-*` | 115 | 24 |
| `flex-shrink-0` | `shrink-0` | ~80 | 40+ |
| `flex-grow` | `grow` | ~15 | 10+ |
| `outline-none` | `outline-hidden` | 115 | 37 |
| `shadow-sm` | `shadow-xs` | varies | 44 |
| `shadow` (bare) | `shadow-sm` | varies | 44 |
| `shadow-md` | (unchanged) | varies | 44 |
| `shadow-lg` | (unchanged) | varies | 44 |
| `bg-opacity-*` | slash modifier | 62 | 21 |
| `ring` (bare) | `ring-3` | varies | - |
| `transform-none` | individual resets | 1 | ChemPage.tsx |
| `rounded-sm` | `rounded-xs` | varies | - |
| `rounded` (bare) | `rounded-sm` | varies | - |
| `blur-sm` | `blur-xs` | varies | - |
| `blur` (bare) | `blur-sm` | varies | - |
| `backdrop-blur-sm` | `backdrop-blur-xs` | varies | 14 |
| `backdrop-blur` (bare) | `backdrop-blur-sm` | varies | 14 |

### Custom Utilities Migration

**glass and text-gradient utilities (survive past Phase 4):**

**Recommendation: Convert to `@utility` directive** since v4 uses native CSS cascade layers and the old `@layer utilities` behavior has changed.

```css
/* v4 */
@utility glass {
  background-color: rgb(255 255 255 / 0.7);
  backdrop-filter: blur(16px);
  border: 1px solid rgb(203 213 225 / 0.5);

  &:where(.dark, .dark *) {
    background-color: rgb(30 41 59 / 0.7);
    backdrop-filter: blur(24px);
    border-color: rgb(51 65 85 / 0.5);
  }
}

@utility text-gradient {
  color: transparent;
  background-clip: text;
  background-image: linear-gradient(to right, var(--color-sky-600), var(--color-indigo-600));

  &:where(.dark, .dark *) {
    background-image: linear-gradient(to right, var(--color-sky-400), var(--color-cyan-400));
  }
}
```

**Component classes (btn-*, card-*, form-*) -- will be replaced in Phase 4:**

**Recommendation: Keep as `@apply` within `@layer components`** for now. These are replaced by shadcn/ui in Phase 4, so minimal effort is appropriate. Fix any `@apply` breakage (e.g., v4 class renames within `@apply` rules) but don't restructure.

### @apply Behavior in v4

In v4, `@apply` still works within `@layer` blocks, but class names used in `@apply` must be valid v4 class names. This means:
- `@apply bg-gradient-to-r` must become `@apply bg-linear-to-r`
- `@apply flex-shrink-0` must become `@apply shrink-0`
- `@apply outline-none` must become `@apply outline-hidden`
- `@apply shadow-sm` must become `@apply shadow-xs`

The upgrade tool does NOT always catch `@apply` renames in CSS files -- manual verification required.

### Animation Consolidation

Currently animations are split between:
1. `tailwind.config.js` -- 3 animation definitions (pulse-slow, pulse-fast, bounce-slow)
2. `tailwind.css` -- 2 keyframe animations (subtle-pulse, mesh-gradient-move)
3. `animations.css` -- 10 keyframe animations, 12 animation classes

**Recommendation: Keep animations.css as a separate file.** Rationale:
- animations.css is pure CSS keyframes and classes -- no Tailwind-specific syntax, minimal v4 impact
- Only `bounce-slow` from animations.css overlaps with tailwind.config.js (same name, same animation)
- The 3 tailwind.config.js animations should move to the `@theme` block
- The 2 tailwind.css keyframes (subtle-pulse, mesh-gradient-move) stay in tailwind.css
- No duplication to consolidate beyond the bounce-slow overlap

### Custom Theme Values -- Audit Results

**Custom shadows: NONE used in TSX files.** `inner-lg`, `inner-xl`, `soft-xl`, `blue-glow` have zero references in TSX. **Recommendation: Drop all custom shadows from @theme.**

**Pattern CSS variables: NOT used in TSX files directly.** The pattern variables (--dots-pattern, --grid-pattern, --mesh-pattern, --noise-pattern) and their component classes (dots-bg, grid-bg, mesh-bg, noise-bg) are defined in tailwind.css but have zero usage in any TSX file. **Recommendation: Remove pattern definitions and their component classes.**

**Font stacks: Used extensively.** `font-mono` appears in 85 references across 39 files, `Inter` is set as body font. **Recommendation: Keep custom font stacks in @theme.**

**Custom spacing (72, 84, 96, 128): Check needed.** Many of these values may now be built-in to v4 (v4 has more default spacing values). Keep only truly custom values.

**Custom opacity (15, 35, 85, 95): Check needed.** v4 supports arbitrary values more easily. May not need explicit theme entries.

### gray.950 Note

Tailwind v4 includes `gray-950` as a built-in shade. The project's custom `gray-950: '#0f1521'` has a DIFFERENT hex value from v4's built-in `gray-950` (#030712). Decision needed:
- If visual continuity matters: keep custom value in @theme (it will override the built-in)
- If aligning with v4 palette: remove custom value and accept the color shift

**Recommendation: Keep the custom value (#0f1521) to preserve visual continuity** -- this is the `--bg-primary` dark mode background, and changing it would be a visible regression.

### Anti-Patterns to Avoid
- **Running upgrade tool on Node 18:** The `@tailwindcss/upgrade` tool requires Node.js 20+. The project currently runs Node 18.20.8. Must use `nvm use 20` or similar before running the tool.
- **Leaving postcss.config.js:** Will cause conflicts with @tailwindcss/vite. Must be deleted.
- **Using `@layer utilities` for custom utilities in v4:** v4 uses native CSS cascade layers; use `@utility` directive instead.
- **Assuming upgrade tool catches everything:** It handles ~90% of class renames in TSX files but may miss: @apply references in CSS, inline styles using theme(), dynamic class construction via template literals.
- **Using `theme()` in CSS:** The `theme()` function is deprecated in v4. Use `var(--color-*)` etc. Only exception: `@media` queries where CSS variables don't work -- there use `theme(--breakpoint-xl)` (note: different syntax, uses `--` prefix).

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Class renames | Manual find-replace across 47 files | `npx @tailwindcss/upgrade` | Handles 90%+ of mechanical renames across TSX, finds edge cases |
| Vendor prefixing | autoprefixer config | @tailwindcss/vite | Built-in, zero config |
| PostCSS pipeline | postcss.config.js + plugins | @tailwindcss/vite | Single plugin replaces entire pipeline |
| Dark mode prose | Custom typography theme config | `prose dark:prose-invert` | Built-in dark mode for typography plugin |
| CSS variable generation | Manual --color-* definitions | @theme block | Tailwind v4 auto-generates CSS variables for all theme values |

**Key insight:** The `@tailwindcss/upgrade` tool is the single most important time-saver. It automatically updates dependencies, rewrites CSS directives, converts the JS config to @theme, and renames classes in template files. Manual migration without it would require touching 47+ files for class renames alone.

## Common Pitfalls

### Pitfall 1: Node.js Version Incompatibility
**What goes wrong:** `npx @tailwindcss/upgrade` fails silently or with errors on Node 18
**Why it happens:** The upgrade tool requires Node.js 20+; project uses Node 18.20.8
**How to avoid:** Use `nvm use 20` (or `nvm use 22`) before running the upgrade tool. The Dockerfile uses Node 24-alpine which is fine, but local development needs attention.
**Warning signs:** Cryptic errors from the upgrade tool, or the tool seemingly running but not modifying files

### Pitfall 2: Border Color Default Change
**What goes wrong:** All `border` utilities without explicit colors suddenly show `currentColor` instead of gray-200
**Why it happens:** v4 changed the default border color from gray-200 to currentColor
**How to avoid:** Add a base style to preserve v3 behavior:
```css
@layer base {
  *, ::after, ::before, ::backdrop, ::file-selector-button {
    border-color: var(--color-gray-200, currentColor);
  }
}
```
**Warning signs:** Borders appearing darker/different color than expected, especially on cards and containers

### Pitfall 3: bg-opacity + dark:bg-color Ordering
**What goes wrong:** After upgrade, `dark:bg-red-900/30` doesn't produce the expected semi-transparent background
**Why it happens:** The upgrade tool may not correctly merge separate `dark:bg-red-900 dark:bg-opacity-30` into a single `dark:bg-red-900/30` class when they appear as separate tokens in a className string
**How to avoid:** Verify the upgrade tool output for each `bg-opacity` instance; the pattern `dark:bg-{color}-900 dark:bg-opacity-{n}` must become `dark:bg-{color}-900/{n}` (a single class, not two)
**Warning signs:** Error messages or dark mode backgrounds being fully opaque or fully transparent

### Pitfall 4: @apply Using Renamed Classes
**What goes wrong:** CSS @apply rules silently fail because they reference v3 class names
**Why it happens:** The upgrade tool focuses on template files (TSX) but may not fully update @apply references in CSS files
**How to avoid:** Manually audit all @apply rules in tailwind.css after running the upgrade tool
**Warning signs:** Component classes (btn-primary, card, etc.) losing their styling

### Pitfall 5: Shadow Size Shift Confusion
**What goes wrong:** Shadows are smaller than expected after migration
**Why it happens:** v4 shifted all shadow sizes down: `shadow-sm` (v3) = `shadow-xs` (v4), `shadow` (v3) = `shadow-sm` (v4). With 344 shadow usages, even one missed rename is visible.
**How to avoid:** The upgrade tool handles this, but verify shadow appearance on cards, buttons, and modals
**Warning signs:** Cards looking "flatter" than before, subtle visual regression in depth

### Pitfall 6: ring Utility Default Change
**What goes wrong:** `ring` no longer produces a 3px blue ring, now produces 1px currentColor
**Why it happens:** v4 changed ring defaults: width from 3px to 1px, color from blue-500 to currentColor
**How to avoid:** Replace `ring` with `ring-3 ring-blue-500` where the old behavior is desired; `ring-2` remains unchanged
**Warning signs:** Focus rings being invisible or wrong color

### Pitfall 7: space-between and divide Selector Change
**What goes wrong:** `space-y-*` and `divide-y-*` may affect elements differently (applies margin-bottom to all-but-last instead of margin-top to all-but-first)
**Why it happens:** v4 changed the CSS selector from `:not([hidden]) ~ :not([hidden])` to `:not(:last-child)`
**How to avoid:** This is mostly backward-compatible but may cause subtle spacing differences with dynamically shown/hidden elements. Visual verification catches this.
**Warning signs:** Spacing before/after specific list items looking different

### Pitfall 8: transform-none No Longer Works
**What goes wrong:** `md:transform-none` in ChemPage.tsx (sidebar reset) stops working
**Why it happens:** v4 removed `transform-none`; individual transform properties must be reset independently
**How to avoid:** Replace `md:transform-none` with `md:translate-x-0 md:translate-y-0` or `md:scale-100` depending on which transform is being reset
**Warning signs:** Sidebar stuck in wrong position on medium+ screens

### Pitfall 9: CSS Variable Syntax in Arbitrary Values
**What goes wrong:** `bg-[--brand-color]` bracket syntax stops working
**Why it happens:** v4 changed arbitrary value syntax from square brackets `[--var]` to parentheses `(--var)` for CSS variables
**How to avoid:** The codebase uses `bg-[var(--bg-primary)]` format (with explicit `var()`) which should still work. But check for any `[--var-name]` patterns without `var()`.
**Warning signs:** Arbitrary values not applying

### Pitfall 10: Tailwind Internal Variables (--tw-*)
**What goes wrong:** Direct references to `--tw-ring-color`, `--tw-ring-opacity`, `--tw-gradient-stops` may break
**Why it happens:** v4 may generate different internal variable names
**How to avoid:** The tailwind.css file directly sets `--tw-ring-color` and `--tw-ring-opacity` in several places. These need testing and potentially replacing with the official v4 approach.
**Warning signs:** Focus rings, gradients not working as expected

## Code Examples

### Complete vite.config.mts (target state)
```typescript
// Source: Official Tailwind v4 docs + existing project config
/// <reference types="vitest/config" />
import { defineConfig } from "vite";
import react from "@vitejs/plugin-react-swc";
import tailwindcss from "@tailwindcss/vite";

export default defineConfig({
  plugins: [
    tailwindcss(),
    react({
      jsxRuntime: "automatic",
    }),
  ],
  server: {
    port: 3000,
    proxy: {
      "/v1": { target: "http://localhost:8000", changeOrigin: true },
      "/latest": { target: "http://localhost:8000", changeOrigin: true },
    },
  },
  build: {
    sourcemap: "hidden",
    outDir: "dist",
  },
  test: {
    globals: true,
    environment: "jsdom",
    setupFiles: "./src/__tests__/setup.ts",
    css: true,
    coverage: {
      provider: "v8",
      reporter: ["text", "lcov"],
    },
  },
});
```

### theme() to var() Conversion (complete mapping)
```css
/* Source: Project tailwind.css analysis + Tailwind v4 docs */

/* :root (light mode) */
--bg-primary: var(--color-slate-100);
--bg-secondary: var(--color-white);
--text-primary: var(--color-slate-900);
--text-secondary: var(--color-slate-600);
--text-accent: var(--color-sky-700);
--text-accent-hover: var(--color-sky-800);
--border-primary: var(--color-slate-300);
--border-secondary: var(--color-slate-200);
--link-primary: var(--color-sky-600);
--link-hover: var(--color-sky-700);
--input-bg: var(--color-white);
--input-border: var(--color-slate-300);
--input-text: var(--color-slate-900);
--input-focus-ring: var(--color-sky-500);
--placeholder-color: var(--color-slate-400);

/* .dark */
--bg-primary: var(--color-gray-950);
--bg-secondary: var(--color-slate-800);
--text-primary: var(--color-slate-100);
--text-secondary: var(--color-slate-300);
--text-accent: var(--color-sky-400);
--text-accent-hover: var(--color-sky-300);
--border-primary: var(--color-slate-700);
--border-secondary: var(--color-slate-800);
--link-primary: var(--color-sky-400);
--link-hover: var(--color-sky-300);
--input-bg: var(--color-slate-800);
--input-border: var(--color-slate-700);
--input-text: var(--color-white);
--input-focus-ring: var(--color-sky-400);
--placeholder-color: var(--color-slate-500);

/* In tailwind.css utilities section */
--tw-ring-color: var(--color-sky-500);  /* was theme("colors.sky.500") */
--tw-ring-color: var(--color-sky-400);  /* was theme("colors.sky.400") */
```

### AllFiltersView.tsx Fix
```tsx
// Source: Project AllFiltersView.tsx lines 756, 785
// v3 (current -- broken in v4)
className="sticky left-[calc(theme(space.12))] ..."

// v4 (fixed)
className="sticky left-[calc(var(--spacing-12))] ..."
```

### @theme Block for Custom Values
```css
/* Source: Tailwind v4 docs + project tailwind.config.js */
@theme {
  /* Custom colors -- only non-default shades */
  --color-gray-750: #2b3544;
  --color-gray-850: #1a2231;
  --color-gray-950: #0f1521;  /* Keep: different from v4 default #030712 */

  --color-blue-350: #7dabf8;
  --color-blue-450: #4f85e6;
  --color-blue-550: #3a72d6;
  --color-blue-650: #2a5db8;
  --color-blue-750: #214a98;
  --color-blue-850: #183878;

  /* Font stacks */
  --font-sans: 'Inter', ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, 'Noto Sans', sans-serif;
  --font-mono: 'JetBrains Mono', ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, 'Liberation Mono', 'Courier New', monospace;

  /* Custom animations */
  --animate-pulse-slow: pulse 3s cubic-bezier(0.4, 0, 0.6, 1) infinite;
  --animate-pulse-fast: pulse 1s cubic-bezier(0.4, 0, 0.6, 1) infinite;
  --animate-bounce-slow: bounce 2s infinite;
}
```

## State of the Art

| Old Approach (v3) | Current Approach (v4) | When Changed | Impact |
|--------------------|-----------------------|--------------|--------|
| `tailwind.config.js` for theme | `@theme {}` block in CSS | v4.0 (Jan 2025) | Config files eliminated |
| `@tailwind base/components/utilities` | `@import "tailwindcss"` | v4.0 | Single import |
| `theme("colors.x.y")` | `var(--color-x-y)` | v4.0 | All theme values are CSS variables |
| `require('@tailwindcss/forms')` in JS | `@plugin "@tailwindcss/forms"` in CSS | v4.0 | CSS-first plugin loading |
| `darkMode: 'class'` in JS | `@custom-variant dark (...)` | v4.0 | CSS-first variant config |
| `@layer utilities { .foo {} }` | `@utility foo {}` | v4.0 | Native cascade layers |
| PostCSS + autoprefixer pipeline | @tailwindcss/vite plugin | v4.0 | Simpler, faster build |
| `bg-opacity-50` utility | `/50` modifier on color | v4.0 | Opacity modifiers only |
| `flex-shrink-0` / `flex-grow` | `shrink-0` / `grow` | v4.0 | Shorter names |
| `bg-gradient-to-r` | `bg-linear-to-r` | v4.0 | CSS spec alignment |
| `shadow-sm` / `shadow` | `shadow-xs` / `shadow-sm` | v4.0 | Size naming shifted |
| `outline-none` | `outline-hidden` | v4.0 | `outline-none` now removes outline entirely |
| `ring` = 3px blue | `ring` = 1px currentColor | v4.0 | Must use `ring-3` for old behavior |

**Browser requirements:** Safari 16.4+, Chrome 111+, Firefox 128+ (uses @property, color-mix())

## Discretion Recommendations

### Color Format: Keep hex
Keep exact hex values for custom shades. OKLCH conversion would change the visual appearance slightly and add complexity for zero user-visible benefit in this phase.

### Dark Mode: Class-based via @custom-variant
```css
@custom-variant dark (&:where(.dark, .dark *));
```
This preserves the existing `.dark` class toggling behavior and is required for Phase 4's dark mode toggle UI.

### Global transition-colors: Scope to interactive elements
The current `* { @apply transition-colors duration-200 ease-in-out; }` rule applies transitions to ALL elements. This is expensive and causes jank on large DOM trees. Scope to interactive elements:
```css
a, button, input, select, textarea, [role="button"] {
  transition-property: color, background-color, border-color, outline-color, text-decoration-color, fill, stroke;
  transition-duration: 200ms;
  transition-timing-function: ease-in-out;
}
```

### Pattern CSS variables: Remove
Zero usage in TSX files. The dots-bg, grid-bg, mesh-bg, noise-bg classes are defined but never referenced. Remove them.

### @apply component classes: Migrate as-is
btn-*, card-*, form-* classes are replaced by shadcn/ui in Phase 4. Minimal effort: fix @apply class renames (e.g., `@apply shadow-sm` -> `@apply shadow-xs`), keep the structure. Don't convert to @utility.

### Custom utility syntax: Use @utility
For `glass` and `text-gradient` (which survive past Phase 4), use `@utility` directive. This is the v4-recommended approach and works correctly with the native cascade layers.

### Base style @apply: Convert where straightforward
Where v4 recommends plain CSS instead of @apply in @layer base, convert. For example:
```css
/* Keep @apply where it references Tailwind utilities */
body { @apply min-h-screen antialiased; }

/* Convert where it's just setting properties */
h1, h2, h3, h4, h5, h6 { font-weight: 600; color: var(--text-primary); }
```

### Animation file: Keep separate
animations.css is pure CSS, no Tailwind-specific syntax. Merging adds complexity for no benefit. Remove `animate-bounce-slow` from animations.css (redundant with @theme animation).

### Font stacks: Keep in @theme
85 `font-mono` references across 39 files. Inter/JetBrains Mono are integral to the design. Keep them.

### Shadows: Drop custom shadows
Zero usage in TSX. Remove inner-lg, inner-xl, soft-xl, blue-glow from @theme.

### VS Code IntelliSense: Skip
No .vscode directory exists in the frontend. Don't create one -- not needed for this migration.

## Open Questions

1. **@tailwindcss/upgrade tool + Node 18**
   - What we know: Tool requires Node 20+, project uses Node 18.20.8 locally
   - What's unclear: Whether the project's dev environment has Node 20+ available via nvm
   - Recommendation: The executor should check for `nvm` and use Node 20+ for the upgrade tool step only. If nvm is unavailable, perform manual class renames (slower but doable).

2. **Typography prose dark overrides -- exact color matching**
   - What we know: `prose-invert` provides built-in dark mode for typography content
   - What's unclear: Whether `prose-invert` colors match the current custom prose colors exactly
   - Recommendation: Accept `prose-invert` defaults. The current custom colors were hand-picked to be "nice dark mode prose" -- `prose-invert` was designed for the same purpose by the Tailwind team. Visual verification will catch any issues.

3. **Internal --tw-* variable compatibility**
   - What we know: The CSS file directly manipulates `--tw-ring-color`, `--tw-ring-opacity`, `--tw-gradient-stops`
   - What's unclear: Whether v4 generates the same internal variable names
   - Recommendation: Test these specific features after migration. The subtle-pulse animation and bg-gradient-radial utility reference these variables. May need to use v4's official ring/gradient syntax instead.

## Validation Architecture

### Test Framework
| Property | Value |
|----------|-------|
| Framework | Vitest 4.1.0 |
| Config file | `vite.config.mts` (test section) |
| Quick run command | `cd frontend && npx vitest run` |
| Full suite command | `cd frontend && npx vitest run --coverage` |

### Phase Requirements -> Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| FRAME-02 | Tailwind v4 installed, Vite plugin configured, JS config removed | smoke | `cd frontend && npx vite build` (build succeeds = CSS pipeline works) | N/A -- build is the test |
| FRAME-02 | CSS import and @theme block correctly processed | smoke | `cd frontend && npx vitest run` (existing smoke test renders without crash) | setup.ts, App.test.tsx exist |
| FRAME-04 | All class renames applied correctly | manual | Visual verification of all 9 pages | N/A -- manual visual check |
| FRAME-04 | No build warnings about unknown classes | smoke | `cd frontend && npx vite build 2>&1` (check for warnings) | N/A -- build output |

### Sampling Rate
- **Per task commit:** `cd frontend && npx vite build` (verifies CSS pipeline)
- **Per wave merge:** `cd frontend && npx vitest run` + visual spot-check
- **Phase gate:** Full 9-page manual visual verification, build succeeds, existing tests pass

### Wave 0 Gaps
None -- existing test infrastructure (Vitest smoke test + Vite build) covers the automated verification needs. The primary verification method for this phase is manual visual inspection, which is appropriate given the nature of CSS migration (no programmatic way to verify visual correctness without screenshot comparison, which is not set up).

## Sources

### Primary (HIGH confidence)
- [Tailwind CSS v4 Upgrade Guide](https://tailwindcss.com/docs/upgrade-guide) - Complete breaking changes list, class renames, migration steps
- [Tailwind CSS v4 Dark Mode docs](https://tailwindcss.com/docs/dark-mode) - @custom-variant configuration, class-based and media-query approaches
- [Tailwind CSS v4.0 Blog Post](https://tailwindcss.com/blog/tailwindcss-v4) - Architecture overview, Oxide engine, CSS-first philosophy

### Secondary (MEDIUM confidence)
- [@tailwindcss/vite npm](https://www.npmjs.com/package/@tailwindcss/vite) - Latest version 4.2.1
- [@tailwindcss/forms npm](https://www.npmjs.com/package/@tailwindcss/forms) - Latest version 0.5.11, @plugin directive usage
- [@tailwindcss/typography npm](https://www.npmjs.com/package/@tailwindcss/typography) - Latest version 0.5.19
- [Tailwind CSS forms plugin guide](https://benjamincrozat.com/tailwind-css-forms-plugin) - @plugin strategy option in v4
- [Tailwind CSS v4 Migration Best Practices](https://www.digitalapplied.com/blog/tailwind-css-v4-2026-migration-best-practices) - Community migration patterns

### Tertiary (LOW confidence)
- None -- all critical claims verified with official sources

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - official docs confirm versions, installation, and approach
- Architecture: HIGH - patterns derived from official upgrade guide and verified against project structure
- Pitfalls: HIGH - derived from official breaking changes list cross-referenced with actual codebase usage counts
- Discretion recommendations: MEDIUM - informed opinions based on project structure, but visual outcomes need verification

**Research date:** 2026-03-12
**Valid until:** 2026-04-12 (Tailwind v4 is stable, not fast-moving at this point)
