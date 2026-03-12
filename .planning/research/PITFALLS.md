# Domain Pitfalls

**Domain:** CRA-to-Vite + React 18-to-19 + Tailwind v3-to-v4 + shadcn/ui migration
**Researched:** 2026-03-12
**Confidence:** HIGH for Tailwind v4 and React 19 pitfalls (verified from official upgrade guides); MEDIUM for shadcn/ui and framer-motion pitfalls (training data)

## Critical Pitfalls

Mistakes that cause rewrites, broken builds, or major regressions.

### Pitfall 1: Tailwind v4 Class Renames Break Existing UI Silently

**What goes wrong:** Tailwind v4 renamed several utility classes (shadow-sm -> shadow-xs, rounded-sm -> rounded-xs, blur-sm -> blur-xs, ring -> ring-3, outline-none -> outline-hidden). These are NOT build errors -- the old class names simply produce no CSS output. Your UI looks broken (no shadows, no rounded corners, no rings) but there are no errors in the console or build log.
**Why it happens:** Tailwind generates CSS only for classes it recognizes. Unrecognized classes are silently ignored.
**Consequences:** Subtle visual regressions across the entire app. Borders, shadows, rounded corners, focus rings all disappear. Hard to catch without visual regression testing.
**Prevention:**
1. Run `npx @tailwindcss/upgrade` -- the automated tool handles most class renames in template files.
2. After the tool runs, grep for known renamed classes: `shadow-sm`, `rounded-sm`, `blur-sm`, `ring ` (bare), `outline-none`, `drop-shadow-sm`.
3. Visually inspect every page after migration. Screenshot before and after.
**Detection:** Components that look "flat" (missing shadows, borders, rounded corners) without any error messages.
**Confidence:** HIGH (verified from official Tailwind v3->v4 upgrade guide)

### Pitfall 2: CSS Variable Theme Functions Break in Tailwind v4

**What goes wrong:** The existing `tailwind.css` uses `theme("colors.slate.100")` inside CSS variable definitions:
```css
:root {
  --bg-primary: theme("colors.slate.100");
}
```
Tailwind v4's CSS-first approach uses `@theme` blocks and CSS custom properties. The old `theme()` function may not resolve the same way with the `@tailwindcss/vite` plugin.
**Why it happens:** `theme()` was a PostCSS plugin feature. The new Vite plugin processes CSS differently.
**Consequences:** All CSS variables for theming resolve to empty/invalid values. Dark mode completely breaks. Every component using `var(--bg-primary)` etc. gets no styling.
**Prevention:**
1. Replace `theme("colors.slate.100")` references with actual color values or Tailwind v4 CSS variable references.
2. In Tailwind v4, colors are exposed as CSS variables: `var(--color-slate-100)`. Use these.
3. Or define with raw values: `--bg-primary: #f1f5f9;`
**Detection:** All backgrounds are transparent/white. Dark mode toggle does nothing visually.
**Confidence:** HIGH (verified: Tailwind v4 upgrade guide confirms migration to CSS-first config)

### Pitfall 3: Environment Variable Prefix Change Causes Silent API Failures

**What goes wrong:** CRA uses `REACT_APP_*` prefix. Vite uses `VITE_*`. If you forget to rename env vars, `import.meta.env.VITE_API_URL` returns `undefined`, and the fallback URL is used. Works locally, breaks in production Docker.
**Why it happens:** Vite strips all environment variables that don't start with `VITE_` for security.
**Consequences:** API calls go to wrong URL. Network errors in production with no clear cause.
**Prevention:**
1. Create `.env` files with `VITE_` prefix
2. Update ALL `process.env.REACT_APP_*` to `import.meta.env.VITE_*`
3. Update Dockerfile ARG/ENV from `REACT_APP_API_URL` to `VITE_API_URL`
4. Files to update: `services/api.js`, `context/AppContext.js`, `Dockerfile`, `.env.production`
5. `grep -r "REACT_APP" frontend/` should return zero results after migration
**Detection:** API calls fail in production but work in development. Network tab shows requests to wrong URL.
**Confidence:** HIGH

### Pitfall 4: Vite Output Directory and index.html Location Change

**What goes wrong:** CRA puts `index.html` in `public/` and outputs to `build/`. Vite requires `index.html` in the project root and outputs to `dist/`. If not adjusted, Vite can't find the entry point and Docker copies from wrong directory.
**Consequences:** `vite build` fails (can't find index.html), or Docker container serves empty directory.
**Prevention:**
1. Move `public/index.html` to `frontend/index.html`
2. Update HTML: remove `%PUBLIC_URL%`, add `<script type="module" src="/src/main.jsx"></script>`
3. Update Dockerfile: `COPY --from=build /app/dist /usr/share/nginx/html`
**Detection:** Build fails immediately or Docker serves blank page.
**Confidence:** HIGH

### Pitfall 5: @layer utilities/components Syntax Change in Tailwind v4

**What goes wrong:** Existing `tailwind.css` uses `@layer utilities { .glass { ... } }` and `@layer components { .btn { ... } }`. Tailwind v4 changes `@layer utilities` to `@utility` for single-utility definitions.
**Consequences:** Custom utilities (`.glass`, `.text-gradient`, `.bg-gradient-radial`, `.transform-style-3d`) stop working. Components using these classes lose their styling.
**Prevention:**
1. Convert to `@utility glass { ... }` syntax
2. The automated upgrade tool handles some of this
3. `@layer components` may be kept or converted
**Detection:** Custom utility classes produce no styling.
**Confidence:** HIGH (verified from Tailwind v4 upgrade guide)

### Pitfall 6: shadcn/ui CSS Variables Conflict with Existing Dark Mode

**What goes wrong:** The existing app has three competing theming approaches: (1) `isDarkMode ? "bg-X" : "bg-Y"` ternaries in JSX, (2) Tailwind `dark:bg-X` classes, (3) CSS variables in `:root` and `.dark` scopes. shadcn/ui adds a fourth: its own CSS variable system. Without cleanup, colors become unpredictable.
**Consequences:** Text invisible against background, incorrect contrast, jarring flashes on mode switch.
**Prevention:** Establish a single strategy: use shadcn/ui CSS variables as source of truth. Tailwind's `dark:` variant triggers from `.dark` class on `<html>` (already how AppContext works). Remove `isDarkMode` ternaries gradually.
**Detection:** Test every page in both light and dark mode after each component migration.
**Confidence:** MEDIUM-HIGH

## Moderate Pitfalls

### Pitfall 7: react-router-dom Future Flags with React 19

**What goes wrong:** The current `App.js` sets `v7_startTransition: true` and `v7_relativeSplatPath: true`. These may interact with React 19's transition system differently.
**Prevention:** Test all routing after React 19 upgrade: navigation, back/forward, URL params (`:convertId`, `:depictId`, `:toolId`).
**Confidence:** MEDIUM

### Pitfall 8: framer-motion and React 19 Compatibility

**What goes wrong:** framer-motion uses `forwardRef` internally, which React 19 deprecates. If framer-motion's version doesn't support React 19's ref-as-prop, you get deprecation warnings or broken animations.
**Prevention:** framer-motion v12.x (current) likely supports React 19. Verify changelog. Upgrade to latest if warnings appear.
**Confidence:** MEDIUM

### Pitfall 9: Default Border Color Change in Tailwind v4

**What goes wrong:** Default border color changed from `gray-200` to `currentColor`. All `border` classes without explicit color will have borders matching text color.
**Prevention:** Add base layer override:
```css
@layer base {
  *, ::after, ::before {
    border-color: var(--color-gray-200, currentColor);
  }
}
```
**Detection:** Borders appear too dark/prominent, especially in dark mode.
**Confidence:** HIGH (verified from Tailwind v4 upgrade guide)

### Pitfall 10: shadcn/ui CLI Assumes TypeScript

**What goes wrong:** `npx shadcn@latest init` generates TypeScript files by default. This project uses JavaScript.
**Prevention:** During `shadcn init`, select JavaScript option (or configure `components.json`). Alternatively, generate as TypeScript and rename to .jsx (removing type annotations is trivial).
**Detection:** Import errors from generated components.
**Confidence:** MEDIUM

### Pitfall 11: Important Modifier Syntax Reversal

**What goes wrong:** `!flex` (v3) becomes `flex!` (v4). Important overrides stop working silently.
**Prevention:** Grep for `!` prefix in class names. Automated upgrade tool handles most cases.
**Confidence:** HIGH (verified from Tailwind v4 upgrade guide)

### Pitfall 12: Node.js 20+ Required for Tailwind Upgrade Tool

**What goes wrong:** `npx @tailwindcss/upgrade` requires Node.js 20+. Current Dockerfile uses `node:18-alpine`.
**Prevention:** Run upgrade tool locally with Node 20+. Update Dockerfile to `node:20-alpine`.
**Confidence:** HIGH (verified from Tailwind v4 upgrade guide)

### Pitfall 13: CSS Arbitrary Value Syntax Change

**What goes wrong:** `bg-[--brand-color]` becomes `bg-(--brand-color)` in v4. The current codebase uses `@apply placeholder-[var(--placeholder-color)]` and `focus:ring-offset-[var(--bg-primary)]`.
**Prevention:** Grep for `[--` and `[var(--` patterns. Update to parentheses syntax.
**Confidence:** HIGH (verified from Tailwind v4 upgrade guide)

### Pitfall 14: @testing-library/react Version Mismatch

**What goes wrong:** Current `@testing-library/react` v13 does not support React 19.
**Prevention:** Upgrade to `@testing-library/react@^16` alongside React 19 upgrade.
**Confidence:** MEDIUM-HIGH

## Minor Pitfalls

### Pitfall 15: CRA browserslist/eslintConfig Ignored by Vite

**Prevention:** Remove `browserslist` and `eslintConfig` from `package.json` during Vite migration.

### Pitfall 16: web-vitals No Longer Auto-Wired

**Prevention:** Remove `web-vitals` unless explicitly needed.

### Pitfall 17: Vite Dev Server Port Default

**Prevention:** Set `server.port: 3000` in `vite.config.js` to match CRA default.

### Pitfall 18: react-icons to lucide-react Mapping

**What goes wrong:** The codebase uses 20+ icons from `react-icons/hi` (Heroicons). Not all have direct lucide-react equivalents.
**Prevention:** Create mapping table before migration. Keep react-icons temporarily for unmapped icons.

## Phase-Specific Warnings

| Phase Topic | Likely Pitfall | Mitigation |
|-------------|---------------|------------|
| CRA to Vite | Pitfalls 3, 4, 15, 16, 17 | Env var grep, move index.html, update Dockerfile, clean package.json |
| React 18 to 19 | Pitfalls 7, 8, 14 | Run codemod, test routing, verify framer-motion, upgrade testing-library |
| Tailwind v3 to v4 | Pitfalls 1, 2, 5, 9, 11, 12, 13 | Run `@tailwindcss/upgrade`, visual regression test all pages |
| shadcn/ui | Pitfalls 6, 10, 18 | Establish single theming strategy, configure JS output, create icon mapping |
| Dark mode | Pitfall 6 | Remove competing approaches, use CSS variables only |
| Docker/CI | Pitfalls 3, 4, 12 | Update Dockerfile, test Docker build immediately after Vite migration |

## Sources

- Tailwind CSS v3->v4 upgrade guide: https://tailwindcss.com/docs/upgrade-guide (verified via WebFetch, HIGH confidence)
- React 19 upgrade guide: https://react.dev/blog/2024/04/25/react-19-upgrade-guide (verified via WebFetch, HIGH confidence)
- Vite migration patterns (training data, MEDIUM-HIGH confidence)
- Direct codebase analysis of all 55+ frontend source files (HIGH confidence)
- shadcn/ui patterns (training data, MEDIUM confidence)
