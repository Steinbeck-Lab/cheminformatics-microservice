# Domain Pitfalls

**Domain:** Frontend modernization -- CRA to Vite, React 18 to 19, Tailwind v3 to v4, shadcn/ui integration
**Researched:** 2026-03-12
**Sources:** React 19 official blog (react.dev), Tailwind CSS v4 upgrade guide (tailwindcss.com/docs/upgrade-guide), codebase analysis of 55+ JSX/JS source files

---

## Critical Pitfalls

Mistakes that cause rewrites, broken builds, or silent visual regressions across the entire application.

---

### Pitfall 1: Environment Variable Prefix -- Silent Runtime Failures After CRA-to-Vite Migration

**What goes wrong:** After migrating to Vite, all `process.env.REACT_APP_*` references silently resolve to `undefined` at runtime. The app builds without errors but the API calls fail because the base URL is `undefined`. CRA uses `process.env.REACT_APP_*` while Vite uses `import.meta.env.VITE_*`.

**Why it happens:** Vite does not polyfill `process.env`. The build succeeds because JavaScript does not throw on accessing properties of `undefined` objects -- `process.env` becomes `{}` or `undefined` depending on Vite config, and property access just returns `undefined`.

**Consequences:** API calls silently fail in production. The fallback URL (`https://dev.api.naturalproducts.net/latest`) may mask the issue in development but the production `.env.production` value is lost.

**Specific impact on this codebase (9 affected locations):**
- `services/api.js:5` -- `process.env.REACT_APP_API_URL` (API base URL)
- `context/AppContext.js:20` -- `process.env.REACT_APP_API_URL` (context API config)
- `components/common/MoleculeCard.jsx:61` -- `process.env.REACT_APP_API_URL`
- `components/common/HighlightedMoleculeCard.jsx:45` -- `process.env.REACT_APP_API_URL`
- `components/chem/StandardizeView.jsx:18` -- `process.env.REACT_APP_API_URL`
- `services/api.js:33` -- `process.env.NODE_ENV` (Vite uses `import.meta.env.MODE`)
- `components/tools/RInChIView.jsx:1437` -- `process.env.PUBLIC_URL`
- `components/tools/InChIView.jsx:1959` -- `process.env.PUBLIC_URL`
- `components/depict/StructureDrawView.jsx:669` -- `process.env.PUBLIC_URL`

**Prevention:**
1. Rename all env files: `REACT_APP_API_URL` becomes `VITE_API_URL`
2. Global find-and-replace `process.env.REACT_APP_` with `import.meta.env.VITE_`
3. Replace `process.env.NODE_ENV` with `import.meta.env.MODE`
4. Replace `process.env.PUBLIC_URL` with empty string (Vite serves from root by default) or `import.meta.env.BASE_URL`
5. Update `.env.local.example`, `.env.production`, and `.env.local`
6. Update Dockerfile ARG/ENV from `REACT_APP_API_URL` to `VITE_API_URL`
7. Update docker-compose.yml environment variables (line 37)
8. Update GitHub Actions workflows (dev-build.yml line 69, prod-build.yml line 91)
9. Run `grep -r "REACT_APP\|PUBLIC_URL\|process\.env" frontend/src/` -- should return zero results after migration

**Detection:** API calls return network errors or `undefined` URLs after migration. Test the app against the live API immediately after the Vite migration.

**Phase:** CRA-to-Vite migration (must be the very first step)

**Confidence:** HIGH -- verified against codebase grep and official Vite documentation

---

### Pitfall 2: Tailwind v4 `theme()` Function Removal Breaks Entire CSS Variable System

**What goes wrong:** The entire dark mode CSS variable system in `tailwind.css` uses `theme("colors.slate.100")` syntax (32 occurrences) to set CSS custom properties. Tailwind v4 replaces the `theme()` function with native CSS variables, so `theme("colors.slate.100")` no longer works. All CSS variables break, and the entire light/dark theming system silently produces invalid CSS.

**Why it happens:** Tailwind v4 moves from JavaScript-based configuration to CSS-first configuration. The `theme()` CSS function is replaced by direct access to CSS variables like `var(--color-slate-100)`. The upgrade tool may not correctly migrate `theme()` calls inside `:root` and `.dark` custom property declarations.

**Consequences:** Every component that uses `var(--bg-primary)`, `var(--text-primary)`, `var(--border-primary)`, etc. loses its styling. Since the CSS variable declarations are all in `@layer base`, this breaks the entire application appearance.

**Specific impact on this codebase:**
- `tailwind.css` lines 10-24: 15 light mode `theme()` calls (e.g., `--bg-primary: theme("colors.slate.100")`)
- `tailwind.css` lines 43-57: 15 dark mode `theme()` calls
- `tailwind.css` lines 245, 249: ring color `theme()` calls
- `components/chem/AllFiltersView.jsx` lines 756, 785: inline `theme(space.12)` in className attributes

**Prevention:**
1. Replace all `theme("colors.slate.100")` with `var(--color-slate-100)` (Tailwind v4 auto-generates these CSS variables)
2. Replace `theme(space.12)` in JSX className strings with the actual value (`3rem`)
3. Run the `@tailwindcss/upgrade` tool first, then manually verify the CSS variable declarations
4. Test both light and dark modes visually after migration

**Detection:** Visual regression -- colors revert to browser defaults or become transparent. Check the `:root` and `.dark` CSS variable declarations in browser DevTools.

**Phase:** Tailwind v3 to v4 migration

**Confidence:** HIGH -- verified against official Tailwind v4 upgrade guide and codebase grep

---

### Pitfall 3: Tailwind v4 `@tailwind` Directive Removal and `@layer` Behavior Change

**What goes wrong:** The `@tailwind base; @tailwind components; @tailwind utilities;` directives are removed in v4, replaced by `@import "tailwindcss";`. Additionally, custom CSS inside `@layer utilities` must use the new `@utility` directive for custom utility definitions. The existing `tailwind.css` file uses all three `@tailwind` directives and has extensive custom CSS in all three layers.

**Why it happens:** Tailwind v4 is now a full CSS-first framework. The old directives are no longer recognized. Custom utilities in `@layer utilities` that use `@apply` may not work correctly because the `@utility` directive is now required for custom utility definitions.

**Consequences:** Build fails with CSS parsing errors, or worse, custom component classes like `.btn`, `.btn-primary`, `.card`, `.glass`, `.text-gradient` silently stop working.

**Specific impact on this codebase:**
- `tailwind.css` lines 2-4: Three `@tailwind` directives that must be replaced
- `tailwind.css` line 7: `@layer base` block with ~120 lines of base styles
- `tailwind.css` line 134: `@layer components` with `.btn`, `.btn-primary`, `.btn-secondary`, `.btn-danger`, `.card`, `.card-header`, `.card-body`, `.card-footer`, `.form-group`, `.form-label`, pattern utility classes
- `tailwind.css` line 214: `@layer utilities` with `.glass`, `.text-gradient`, keyframes, `.bg-gradient-radial`, `.transform-style-3d`

**Prevention:**
1. Replace the three `@tailwind` directives with single `@import "tailwindcss";`
2. Keep `@layer base` styles as-is (this still works in v4)
3. Convert custom utilities in `@layer utilities` to use `@utility` directive:
   ```css
   @utility glass {
     @apply bg-white/70 backdrop-blur-lg border border-slate-300/50;
     @apply dark:bg-slate-800/70 dark:backdrop-blur-xl dark:border-slate-700/50;
   }
   ```
4. Test all custom component classes (`.btn-primary`, `.card`, `.glass`, `.text-gradient`) after migration
5. Run the Tailwind upgrade tool and carefully review its diff

**Detection:** Build errors referencing `@tailwind`, or missing styles on elements using custom component classes.

**Phase:** Tailwind v3 to v4 migration

**Confidence:** HIGH -- verified against official Tailwind v4 upgrade guide

---

### Pitfall 4: Tailwind v4 Deprecated Opacity Utilities -- 40+ Occurrences in This Codebase

**What goes wrong:** Tailwind v4 removes `bg-opacity-*`, `text-opacity-*`, `ring-opacity-*`, and `border-opacity-*` utilities. This codebase uses `dark:bg-opacity-30`, `dark:bg-opacity-20`, `dark:bg-opacity-50`, `bg-opacity-75`, and `ring-opacity-5` extensively across every component (40+ occurrences across 25+ files).

**Why it happens:** Tailwind v4 replaces these with the modern slash syntax: `bg-red-900/30` instead of `bg-red-900 bg-opacity-30`. The old utilities are completely removed, not just deprecated.

**Consequences:** All opacity modifiers silently fail -- elements that should be semi-transparent become fully opaque, causing a dramatically different (worse) visual appearance. Alert boxes, modals, overlays, and info panels all lose their subtle transparency.

**Specific patterns requiring find-and-replace:**
- `dark:bg-red-900 dark:bg-opacity-30` becomes `dark:bg-red-900/30` (error alerts across ~15 files)
- `dark:bg-blue-900 dark:bg-opacity-20` becomes `dark:bg-blue-900/20` (info panels across ~15 files)
- `dark:bg-yellow-900 dark:bg-opacity-20` and `dark:bg-opacity-30` (warning panels)
- `dark:bg-green-900 dark:bg-opacity-20` (success panels)
- `bg-gray-500 bg-opacity-75` becomes `bg-gray-500/75` (AllFiltersView.jsx:609 modal overlay)
- `dark:bg-gray-900 dark:bg-opacity-75` becomes `dark:bg-gray-900/75` (same modal)
- `ring-black ring-opacity-5` becomes `ring-black/5` (SMILESInput.jsx:346,385 dropdown shadow)
- `dark:bg-blue-900 dark:bg-opacity-50` becomes `dark:bg-blue-900/50` (ClassyfireView.jsx:521 badges)
- `--tw-ring-opacity: 0.5` in tailwind.css:115 must be replaced with appropriate slash syntax

**Prevention:**
1. Run the `@tailwindcss/upgrade` tool -- it should handle most of these automatically
2. Verify with a global regex search for `opacity-\d+` in className strings post-upgrade
3. Manually fix any the tool misses -- especially combined patterns where `bg-opacity` is on a separate utility from the color class (e.g., `dark:bg-red-900 dark:bg-opacity-30` where the tool must merge two separate classes)
4. The `--tw-ring-opacity` internal variable in tailwind.css needs manual attention

**Detection:** Visual regression -- semi-transparent backgrounds become fully opaque. Compare screenshots before and after migration.

**Phase:** Tailwind v3 to v4 migration

**Confidence:** HIGH -- verified against official Tailwind v4 upgrade guide and exhaustive codebase grep

---

### Pitfall 5: Tailwind v4 Renamed Utility Scale -- `shadow`, `rounded`, `blur`, `ring` All Shift

**What goes wrong:** Tailwind v4 renames the default-size utilities in a shifted scale. `shadow` becomes `shadow-sm`, `shadow-sm` becomes `shadow-xs`, `rounded` becomes `rounded-sm`, `rounded-sm` becomes `rounded-xs`, `blur` becomes `blur-sm`, `ring` (3px) becomes `ring-3` (default is now 1px), and `outline-none` becomes `outline-hidden`. These changes affect virtually every component in the codebase.

**Why it happens:** Tailwind v4 normalized its sizing scales. The "default" (no suffix) size was inconsistent across utilities. The shift makes the naming consistent but breaks every existing usage.

**Consequences:** Subtle but pervasive visual differences. Shadows become slightly smaller. Border radii shift. Ring widths change from 3px to 1px. Focus outlines disappear if using `outline-none` (now means `outline-style: none` instead of visually hidden). These changes are individually small but collectively make the UI look "off."

**Specific high-impact changes for this codebase:**
- `shadow` (used on many cards/alerts) now maps to a smaller shadow
- `rounded-md`, `rounded-lg` (heavily used) are unaffected, but any bare `rounded` shifts
- `ring-2` (used in focus states via tailwind.css) is unaffected, but bare `ring` shifts from 3px to 1px
- `focus:outline-none` (used in `.btn` class and base input styles) needs to become `focus:outline-hidden`
- `backdrop-blur-lg`, `backdrop-blur-xl` (used in `.glass` utility) are unaffected

**Prevention:**
1. The `@tailwindcss/upgrade` tool should handle most renames automatically
2. After running the tool, search for bare `shadow\b`, `rounded\b`, `blur\b`, `ring\b` without suffixes to verify they were updated
3. Test focus states specifically -- the `outline-none` to `outline-hidden` change affects keyboard navigation visibility
4. Compare shadow depths visually before/after

**Detection:** Subtle visual differences in shadow depth, border radius, and ring width. Focus states may look different or disappear. Use browser DevTools to compare computed styles.

**Phase:** Tailwind v3 to v4 migration

**Confidence:** HIGH -- verified against official Tailwind v4 upgrade guide

---

### Pitfall 6: `index.html` Location and `%PUBLIC_URL%` Template Syntax in Vite

**What goes wrong:** CRA expects `public/index.html` with `%PUBLIC_URL%` template variables. Vite expects `index.html` at the project root (not in `public/`) and does not process `%PUBLIC_URL%` templates. If you leave index.html in `public/`, Vite will not find it. If you move it but keep `%PUBLIC_URL%`, the templates render as literal text in the HTML.

**Why it happens:** Vite's dev server serves `index.html` from the project root as the entry point. It does not have a template variable system like CRA. Static assets in `public/` are still served, but `index.html` itself must be at root.

**Consequences:** Vite dev server shows a blank page or 404. Favicon, manifest, and apple-touch-icon links break. The 3Dmol.js script tag and Matomo analytics script survive (they use absolute URLs) but asset references with `%PUBLIC_URL%` break.

**Specific impact on this codebase:**
- `public/index.html` line 5: `%PUBLIC_URL%/favicon.ico` -- breaks
- `public/index.html` line 17: `%PUBLIC_URL%/img/logo_small.png` -- breaks
- `public/index.html` line 18: `%PUBLIC_URL%/manifest.json` -- breaks
- The inline CSP meta tag (line 10) and Matomo script (lines 57-73) are fine
- The 3Dmol.js CDN script (lines 77-81) is fine
- A `<script type="module" src="/src/index.jsx"></script>` entry must be added manually before `</body>`

**Prevention:**
1. Move `public/index.html` to `frontend/index.html` (project root)
2. Replace all `%PUBLIC_URL%` with empty string or `/` (Vite serves public assets from root)
3. Add `<script type="module" src="/src/index.jsx"></script>` before `</body>`
4. Remove the CRA-injected runtime chunk comments if present
5. Update the Dockerfile `COPY` path since the build output directory changes from `build/` to `dist/`

**Detection:** Blank page on `npm run dev`, or broken favicon/manifest in production.

**Phase:** CRA-to-Vite migration (first step)

**Confidence:** HIGH -- well-documented Vite behavior

---

### Pitfall 7: Dockerfile Build Output Directory Change (`build/` to `dist/`)

**What goes wrong:** CRA outputs to `build/`. Vite outputs to `dist/`. The Dockerfile copies `COPY --from=build /app/build /usr/share/nginx/html` which either fails or produces an empty nginx directory.

**Why it happens:** Different default output directories between CRA and Vite.

**Consequences:** Docker-deployed frontend shows a blank nginx page or 403 error. This breaks both dev and prod deployments, and also breaks CI/CD pipelines.

**Specific impact on this codebase:**
- `frontend/Dockerfile` line 27: `COPY --from=build /app/build /usr/share/nginx/html` must change to `/app/dist`
- `frontend/Dockerfile` line 17: `ARG REACT_APP_API_URL` must become `ARG VITE_API_URL`
- `frontend/Dockerfile` line 18: `ENV REACT_APP_API_URL=$REACT_APP_API_URL` must become `ENV VITE_API_URL=$VITE_API_URL`
- `docker-compose.yml` line 37: `REACT_APP_API_URL` environment variable must become `VITE_API_URL`
- GitHub Actions `dev-build.yml` line 69: `REACT_APP_API_URL` build arg must become `VITE_API_URL`
- GitHub Actions `prod-build.yml` line 91: same change

**Prevention:**
1. Update Dockerfile: change `COPY --from=build /app/build` to `COPY --from=build /app/dist`
2. Alternatively, configure Vite to output to `build/` via `build.outDir: 'build'` in vite.config.js
3. Update all `REACT_APP_*` ARGs/ENVs to `VITE_*` in Dockerfile, docker-compose.yml, and GitHub Actions
4. Update Dockerfile base image from `node:18-alpine` to `node:20-alpine` (Node 20 is LTS and recommended for Vite + Tailwind v4)
5. Test the Docker build locally before merging

**Detection:** Docker container serves blank page. Check nginx logs for 404/403 errors.

**Phase:** CRA-to-Vite migration

**Confidence:** HIGH -- standard Vite default behavior

---

### Pitfall 8: Tailwind v4 Config File Incompatibility -- JS Config No Longer Auto-Detected

**What goes wrong:** Tailwind v4 does not automatically detect `tailwind.config.js`. The extensive configuration in this project (custom colors, fonts, spacing, shadows, animations, plugins) is silently ignored. The app builds but all custom theme values disappear.

**Why it happens:** Tailwind v4 uses CSS-first configuration via `@theme` directive in CSS files. The JavaScript config file is supported via the `@config` directive, but it is not auto-detected. Additionally, plugins like `@tailwindcss/forms` and `@tailwindcss/typography` have new v4-compatible versions that work differently.

**Consequences:** Custom color shades (gray-750, gray-850, gray-950, blue-350 through blue-850), custom spacing, custom animations, custom shadows, custom font families, and plugin styles all stop working. The typography plugin's custom prose styling is lost.

**Specific theme values at risk (from `tailwind.config.js`):**
- 6 custom gray shades (750, 850, 950)
- 6 custom blue shades (350-850)
- Custom spacing values (72, 84, 96, 128)
- Custom maxWidth (8xl)
- Custom opacity values (15, 35, 85, 95)
- Custom borderRadius extensions (xl, 2xl, 3xl)
- Custom boxShadow (inner-lg, inner-xl, soft-xl, blue-glow)
- Custom font families (Inter, JetBrains Mono -- full stack definitions)
- Custom animations (pulse-slow, pulse-fast, bounce-slow)
- Custom z-index values (60-100)
- Custom backdrop blur (xs)
- `@tailwindcss/forms` plugin
- `@tailwindcss/typography` plugin with extensive dark prose customization

**Prevention:**
1. **Option A (recommended):** Migrate the JS config to CSS `@theme` directive in the main CSS file:
   ```css
   @import "tailwindcss";
   @theme {
     --color-gray-750: #2b3544;
     --color-gray-850: #1a2231;
     --color-gray-950: #0f1521;
     --color-blue-350: #7dabf8;
     /* ... all custom values ... */
     --font-sans: 'Inter', ui-sans-serif, system-ui, /* ... */;
     --font-mono: 'JetBrains Mono', ui-monospace, /* ... */;
   }
   @plugin "@tailwindcss/forms";
   @plugin "@tailwindcss/typography";
   ```
2. **Option B (quick migration):** Add `@config "./tailwind.config.js";` to the CSS file to keep the old config. This is a compatibility shim, not a long-term solution.
3. Install v4-compatible plugin versions: `@tailwindcss/forms` and `@tailwindcss/typography` both have v4 versions
4. The `darkMode: 'class'` config option is now handled separately (see Pitfall 16)
5. Remove `autoprefixer` from PostCSS config -- Tailwind v4 includes autoprefixing

**Detection:** Custom colors render as the Tailwind default palette (or not at all). Custom shadows and animations disappear. Typography plugin styles are missing.

**Phase:** Tailwind v3 to v4 migration

**Confidence:** HIGH -- verified against official Tailwind v4 upgrade guide

---

### Pitfall 9: `darkMode: 'class'` Strategy Must Be Explicitly Configured in Tailwind v4

**What goes wrong:** Tailwind v4 removes the `darkMode` config option from `tailwind.config.js`. By default, v4 uses `@media (prefers-color-scheme: dark)` for the `dark:` variant. This codebase uses class-based dark mode (`darkMode: 'class'` in config, toggled via `document.documentElement.classList` in AppContext). If not reconfigured, the manual dark mode toggle stops working entirely -- dark mode only activates based on OS preference.

**Why it happens:** Tailwind v4 defaults to media-query-based dark mode. Class-based dark mode must be explicitly configured via CSS.

**Consequences:** The dark mode toggle button in the header has no effect. The app only respects the OS-level dark mode preference. Every `dark:` variant in every component (hundreds of occurrences) responds to OS preference instead of the class toggle.

**Specific impact on this codebase:**
- `context/AppContext.js` lines 82-91: Toggling `document.documentElement.classList.add/remove("dark")` has no effect on Tailwind utilities
- Every `dark:` variant in every component file (50+ files) now responds to OS preference instead of the class toggle
- The CSS variables in `.dark {}` selector in `tailwind.css` lines 41-72 still work (plain CSS), but the Tailwind `dark:` variant classes don't match the class

**Prevention:**
1. Add this to the main CSS file (after the `@import "tailwindcss"` line):
   ```css
   @custom-variant dark (&:where(.dark, .dark *));
   ```
2. This makes the `dark:` variant match the `.dark` class on any ancestor element, restoring the v3 behavior
3. Keep the existing AppContext toggle logic unchanged -- it already manages the `dark` class on `documentElement`

**Detection:** Toggle the dark mode switch -- if nothing changes visually (Tailwind classes don't respond), the class-based dark mode is not configured.

**Phase:** Tailwind v3 to v4 migration

**Confidence:** HIGH -- verified against official Tailwind v4 upgrade guide

---

## Moderate Pitfalls

Mistakes that cause bugs or significant rework but don't require full rewrites.

---

### Pitfall 10: Tailwind v4 Default Border/Ring Color Change Causes Subtle Visual Breakage

**What goes wrong:** In Tailwind v3, bare `border` defaults to `gray-200`. In v4, it defaults to `currentColor`. This means every `border` utility without an explicit color now uses the text color instead of light gray, making borders much more prominent and visually jarring.

**Why it happens:** Tailwind v4 changed the defaults to be more CSS-native. `currentColor` is the CSS default for border-color, and v4 aligns with this.

**Consequences:** Cards, inputs, dividers, and panels that rely on bare `border` (without `border-gray-200`) suddenly have dark borders instead of subtle gray ones.

**Prevention:**
1. Search for `\bborder\b` in JSX files that lack an accompanying `border-{color}` class
2. Add a compatibility base style:
   ```css
   @layer base {
     *, ::after, ::before, ::backdrop, ::file-selector-button {
       border-color: var(--color-gray-200, currentColor);
     }
   }
   ```
3. Or explicitly add `border-gray-200 dark:border-gray-700` to every bare `border` usage

**Detection:** Borders become noticeably darker/thicker-looking because they match text color instead of subtle gray.

**Phase:** Tailwind v3 to v4 migration

**Confidence:** HIGH -- verified against official Tailwind v4 upgrade guide

---

### Pitfall 11: Headless UI v1 Incompatible with React 19

**What goes wrong:** The project lists `@headlessui/react` v1.7.17 as a dependency. Even though no imports were found in the current source files, it remains in the dependency tree and may cause peer dependency conflicts with React 19.

**Why it happens:** Headless UI v1 was designed for React 17/18. React 19 changes ref handling and deprecates several APIs that Headless UI v1 depends on internally.

**Prevention:**
1. Since no imports were found, simply remove `@headlessui/react` from package.json
2. shadcn/ui (which uses Radix UI primitives) replaces any need for Headless UI
3. If any component does use it, migrate to shadcn/ui equivalents

**Detection:** Peer dependency warnings during `npm install` after React 19 upgrade. Runtime errors if any component imports Headless UI.

**Phase:** React 18 to 19 migration (or during dependency cleanup)

**Confidence:** MEDIUM -- no imports found but dependency exists; verify before removing

---

### Pitfall 12: React 19 `Context.Provider` Deprecation in AppContext

**What goes wrong:** The codebase uses `<AppContext.Provider value={contextValue}>` in `context/AppContext.js:113`. React 19 deprecates `Context.Provider` in favor of rendering `<Context>` directly. While this still works with a deprecation warning, it produces console noise.

**Why it happens:** React 19 simplifies the Context API.

**Prevention:**
1. Change `<AppContext.Provider value={contextValue}>{children}</AppContext.Provider>` to `<AppContext value={contextValue}>{children}</AppContext>`
2. This is a simple one-line change

**Detection:** Console deprecation warning about `Context.Provider`.

**Phase:** React 18 to 19 migration

**Confidence:** HIGH -- verified against React 19 official release notes

---

### Pitfall 13: CRA Testing Setup Incompatible with Vite

**What goes wrong:** The project uses `@testing-library/react` v13 and `@testing-library/jest-dom` v5 with CRA's built-in Jest configuration. After removing `react-scripts`, the test runner (`npm test` which runs `react-scripts test`) breaks completely. The `eslintConfig` section referencing `react-app` also breaks.

**Why it happens:** CRA bundles Jest, Babel transforms, and test configuration. Removing `react-scripts` removes all of this. Vite projects typically use Vitest instead of Jest.

**Consequences:** `npm test` fails. CI pipeline may break if it runs frontend tests.

**Specific impact on this codebase:**
- `package.json` scripts: `"test": "react-scripts test"` must be replaced
- `package.json` eslintConfig: `"extends": ["react-app", "react-app/jest"]` must be replaced
- `package.json` browserslist: not used by Vite (can be removed or kept for other tools)
- Testing libraries need version updates for React 19 compatibility

**Prevention:**
1. Install Vitest: `npm install -D vitest @testing-library/react@latest @testing-library/jest-dom@latest jsdom`
2. Configure Vitest in `vite.config.js`:
   ```js
   test: {
     globals: true,
     environment: 'jsdom',
     setupFiles: './src/setupTests.js',
   }
   ```
3. Update `"test"` script to `"vitest"` or `"vitest run"` for CI
4. Replace ESLint config with a standalone ESLint configuration (eslint.config.js)
5. Update `@testing-library/react` to v16+ for React 19 support

**Detection:** `npm test` fails immediately after removing react-scripts.

**Phase:** CRA-to-Vite migration

**Confidence:** HIGH -- standard CRA migration requirement

---

### Pitfall 14: PostCSS Configuration Must Change for Tailwind v4

**What goes wrong:** The current PostCSS setup uses `tailwindcss` and `autoprefixer` as plugins. Tailwind v4 replaces both with a single `@tailwindcss/postcss` plugin (or preferably the `@tailwindcss/vite` Vite plugin). Using the old PostCSS config produces either build errors or no Tailwind processing at all.

**Why it happens:** Tailwind v4 includes autoprefixing built-in and has a completely new PostCSS plugin.

**Prevention:**
1. **Preferred for Vite:** Use `@tailwindcss/vite` as a Vite plugin instead of PostCSS:
   ```js
   import tailwindcss from "@tailwindcss/vite";
   export default defineConfig({
     plugins: [react(), tailwindcss()],
   });
   ```
2. Remove `autoprefixer` from devDependencies (Tailwind v4 includes it)
3. Remove or simplify `postcss.config.js` (not needed if using the Vite plugin)
4. Remove `postcss` from devDependencies (Vite includes it)

**Detection:** Tailwind classes not being processed. Raw `@import "tailwindcss"` appearing in browser styles.

**Phase:** Tailwind v3 to v4 migration (done alongside Vite migration)

**Confidence:** HIGH -- verified against official Tailwind v4 upgrade guide

---

### Pitfall 15: shadcn/ui Requires Path Aliases and Specific Setup

**What goes wrong:** shadcn/ui initialization (`npx shadcn@latest init`) expects path aliases (`@/` mapping to `src/`) and a specific Tailwind CSS setup. This codebase uses JavaScript (.js/.jsx) with no path aliases configured. Running `shadcn init` may fail or produce components with broken imports.

**Why it happens:** shadcn/ui generates components using `@/components/ui/...` import paths and depends on a `cn()` utility from `@/lib/utils`.

**Consequences:** Generated components have broken imports. The `cn()` utility function cannot be found. Component files are placed in wrong directories.

**Prevention:**
1. Configure path aliases in `vite.config.js`:
   ```js
   resolve: {
     alias: {
       "@": path.resolve(__dirname, "./src"),
     },
   }
   ```
2. Configure path aliases in `jsconfig.json` (or `tsconfig.json` if converting to TS):
   ```json
   {
     "compilerOptions": {
       "baseUrl": ".",
       "paths": { "@/*": ["./src/*"] }
     }
   }
   ```
3. Create `src/lib/utils.js` with the `cn()` helper (uses `clsx` + `tailwind-merge`)
4. Install required dependencies: `clsx`, `tailwind-merge`, `class-variance-authority`
5. During `shadcn init`, configure for your specific setup (JS vs TS, Tailwind v4)

**Detection:** Import errors referencing `@/components/ui/...` or missing `cn` utility.

**Phase:** shadcn/ui integration (must happen after Vite and Tailwind v4 migrations)

**Confidence:** MEDIUM -- verify current shadcn/ui Vite installation docs at time of migration

---

### Pitfall 16: `flex-shrink-0` and `flex-grow` Renamed in Tailwind v4

**What goes wrong:** Tailwind v4 renames `flex-shrink-*` to `shrink-*` and `flex-grow-*` to `grow-*`. The old names are removed. This codebase uses `flex-shrink-0` in 50+ locations (icons, layout containers, sidebars) and `flex-grow` in 10+ locations.

**Why it happens:** Tailwind v4 shortened these utility names for consistency. The old names no longer generate CSS.

**Consequences:** Icons in alert messages lose their fixed size and get squished. Layout containers flex unexpectedly. The sidebar on ChemPage may collapse. Footer dividers lose their grow behavior.

**Prevention:**
1. The `@tailwindcss/upgrade` tool should handle these automatically
2. Global find-and-replace: `flex-shrink-0` to `shrink-0`, `flex-grow` to `grow`
3. Verify with regex search after migration: `grep -r "flex-shrink\|flex-grow" frontend/src/`

**Detection:** Icons in alert/error messages appear squished or misaligned. Layout containers flex when they shouldn't.

**Phase:** Tailwind v3 to v4 migration

**Confidence:** HIGH -- verified against official Tailwind v4 upgrade guide and codebase grep (60+ occurrences)

---

### Pitfall 17: Vite Does Not Support `require()` -- Config File ESM Conflict

**What goes wrong:** The `tailwind.config.js` uses `require('@tailwindcss/forms')` and `require('@tailwindcss/typography')` with `module.exports`. If the project adds `"type": "module"` to package.json (which Vite scaffolding often does), `require()` and `module.exports` break.

**Why it happens:** Vite is ESM-native. Adding `"type": "module"` makes all `.js` files treated as ES modules where `require()` is not available.

**Prevention:**
1. If migrating to Tailwind v4's CSS-first config, this becomes moot -- plugins are imported via CSS: `@plugin "@tailwindcss/forms";`
2. If keeping the JS config temporarily, ensure package.json does NOT have `"type": "module"`, or rename config files to `.cjs` extensions
3. Best approach: migrate fully to Tailwind v4 CSS-first config and avoid the issue entirely

**Detection:** Build error: `require is not defined in ES module scope`.

**Phase:** CRA-to-Vite migration / Tailwind v3 to v4 migration (depends on order)

**Confidence:** MEDIUM -- depends on whether `"type": "module"` is added

---

### Pitfall 18: shadcn/ui CSS Variables Conflict with Existing Dark Mode System

**What goes wrong:** The existing app has CSS variables in `:root` and `.dark` scopes for theming. shadcn/ui brings its own CSS variable system (HSL-based colors like `--background`, `--foreground`, `--primary`, etc.). Without reconciliation, the two systems conflict -- components reference different variables, and color consistency breaks.

**Why it happens:** shadcn/ui was designed with its own theming approach. When added to an existing themed app, the variable namespaces can overlap or diverge.

**Consequences:** Inconsistent colors between existing components and new shadcn/ui components. Some elements may be invisible (same background and text color).

**Prevention:**
1. When running `shadcn init`, map its CSS variables to your existing color scheme
2. Consolidate into one system: either adopt shadcn/ui's variable naming or map shadcn/ui variables to your existing ones
3. Test both light and dark mode after each component migration
4. The existing `var(--bg-primary)` / `var(--text-primary)` variables can coexist with shadcn's if properly coordinated

**Detection:** Test every page in both light and dark mode after adding shadcn/ui components. Look for invisible text or jarring color mismatches.

**Phase:** shadcn/ui integration

**Confidence:** MEDIUM-HIGH -- based on architectural analysis of the existing theming system

---

### Pitfall 19: CSS Arbitrary Value Syntax Change for CSS Variables

**What goes wrong:** Tailwind v4 changes the syntax for CSS variables in arbitrary values from square brackets to parentheses: `bg-[--brand-color]` becomes `bg-(--brand-color)`. The current codebase uses `@apply placeholder-[var(--placeholder-color)]` and `focus:ring-offset-[var(--bg-primary)]` in the base CSS.

**Why it happens:** Tailwind v4 adopted a new syntax to distinguish CSS variables from other arbitrary values.

**Prevention:**
1. Search for `[var(--` and `[--` patterns in CSS and JSX files
2. Replace with parentheses syntax: `bg-(--brand-color)` or `bg-(var(--brand-color))`
3. The upgrade tool should handle most of these

**Detection:** Elements using CSS variable arbitrary values lose their styling.

**Phase:** Tailwind v3 to v4 migration

**Confidence:** HIGH -- verified against official Tailwind v4 upgrade guide

---

## Minor Pitfalls

Issues that cause confusion or minor bugs but are quickly fixable.

---

### Pitfall 20: React 19 `ref` Callback Cleanup Functions

**What goes wrong:** React 19 allows ref callbacks to return cleanup functions. If any ref callback has an implicit return (arrow function without braces), React 19 interprets the return value as a cleanup function, which may cause unexpected behavior.

**Prevention:**
1. This codebase uses no `forwardRef` (verified by grep), which is good
2. Search for ref callbacks with implicit returns: `ref={el => someVar = el}` should become `ref={el => { someVar = el }}`
3. Use the React codemod: `npx types-react-codemod@latest no-implicit-ref-callback-return`

**Detection:** Console warnings about unexpected ref cleanup return values.

**Phase:** React 18 to 19 migration

**Confidence:** HIGH -- verified against React 19 official release notes

---

### Pitfall 21: Vite's JSX Transform Requires Explicit Configuration for .js Files

**What goes wrong:** Vite with `@vitejs/plugin-react` expects JSX in `.jsx` or `.tsx` files by default. This codebase has two entry files with `.js` extensions that contain JSX: `index.js` and `App.js`.

**Prevention:**
1. Rename `src/index.js` to `src/index.jsx` and `src/App.js` to `src/App.jsx`
2. Also check `src/context/AppContext.js` and hook files -- rename if they contain JSX
3. Update the `<script>` tag in `index.html` to reference the new filename
4. Alternatively, configure Vite to treat `.js` as `.jsx`:
   ```js
   esbuild: { loader: { '.js': 'jsx' } }
   ```
   But renaming is cleaner.

**Detection:** Build error: `The JSX syntax is not enabled for .js files`.

**Phase:** CRA-to-Vite migration

**Confidence:** HIGH -- standard Vite behavior with React

---

### Pitfall 22: Global `* { transition-colors }` Rule Conflicts with shadcn/ui

**What goes wrong:** The codebase has `* { @apply transition-colors duration-200 ease-in-out; }` in the base layer (tailwind.css line 128-130). This applies color transitions to EVERY element. When shadcn/ui components are added, their own transition definitions conflict with this global rule, creating janky double-transitions.

**Prevention:**
1. Remove the global `* { transition-colors }` rule
2. Add transitions only where needed using Tailwind's `transition-colors` utility on specific elements
3. shadcn/ui components include their own transitions -- let them handle their own animation

**Detection:** Sluggish or double transitions on interactive elements after adding shadcn/ui.

**Phase:** Tailwind v4 migration (address preemptively before shadcn/ui)

**Confidence:** MEDIUM -- based on architectural analysis

---

### Pitfall 23: Tailwind v4 Preflight Changes -- Button Cursor and Placeholder Color

**What goes wrong:** Tailwind v4 Preflight no longer adds `cursor: pointer` to buttons (uses `cursor: default` instead). Placeholder text changes from gray-400 to current text color at 50% opacity. Buttons throughout the app lose their pointer cursor.

**Prevention:**
1. Add a base style to restore pointer cursor:
   ```css
   @layer base {
     button, [role="button"] { cursor: pointer; }
   }
   ```
2. Check placeholder appearance in both light and dark modes

**Detection:** Hovering over buttons shows default cursor instead of pointer.

**Phase:** Tailwind v3 to v4 migration

**Confidence:** HIGH -- verified against official Tailwind v4 upgrade guide

---

### Pitfall 24: Vite's Public Asset Handling for Standalone HTML Iframe

**What goes wrong:** Three components embed `public/standalone/index.html` in an iframe using `process.env.PUBLIC_URL`. After Vite migration, this path must change.

**Specific impact:**
- `StructureDrawView.jsx:669` -- iframe src for SMILES editor
- `RInChIView.jsx:1437` -- iframe src for SMILES editor
- `InChIView.jsx:1959` -- iframe src for SMILES editor

**Prevention:**
1. Replace `${process.env.PUBLIC_URL}/standalone/index.html` with `/standalone/index.html`
2. Ensure `public/standalone/` directory is not processed by Vite's build pipeline (files in `public/` are copied as-is by default)
3. Test that the standalone editor loads in both dev and production

**Detection:** Iframe shows blank page or 404 after migration.

**Phase:** CRA-to-Vite migration

**Confidence:** HIGH

---

### Pitfall 25: Node.js 20+ Required for Tailwind v4 Upgrade Tool

**What goes wrong:** `npx @tailwindcss/upgrade` requires Node.js 20+. The current Dockerfile uses `node:18-alpine`.

**Prevention:**
1. Run upgrade tool locally with Node 20+
2. Update Dockerfile to `node:20-alpine` (Node 20 is current LTS)
3. Verify CI/CD workflows use Node 20+

**Detection:** Upgrade tool fails with version requirement error.

**Phase:** Tailwind v3 to v4 migration (prerequisite)

**Confidence:** HIGH -- verified against official Tailwind v4 upgrade guide

---

### Pitfall 26: react-router-dom Future Flags May Interact with React 19 Transitions

**What goes wrong:** The current `App.js` sets `v7_startTransition: true` and `v7_relativeSplatPath: true`. These may interact with React 19's transition system differently than React 18.

**Prevention:** Test all routing after React 19 upgrade: navigation, back/forward, URL params (`:convertId`, `:depictId`, `:toolId`). Monitor for transition-related rendering glitches.

**Detection:** Pages flash or re-render unexpectedly during navigation.

**Phase:** React 18 to 19 migration

**Confidence:** MEDIUM

---

### Pitfall 27: web-vitals and CRA-Specific Packages No Longer Needed

**What goes wrong:** `web-vitals` was auto-wired by CRA's template. After removing react-scripts, the auto-wiring is gone. The package just sits unused.

**Prevention:** Remove `web-vitals` from dependencies during Vite migration. Remove `browserslist` config from package.json. Remove `eslintConfig` from package.json.

**Phase:** CRA-to-Vite migration (cleanup)

**Confidence:** HIGH

---

## Phase-Specific Warnings

| Phase Topic | Likely Pitfall | Mitigation |
|---|---|---|
| **CRA to Vite** | Env vars silently fail (#1) | Find-replace all 9 `process.env` locations before testing |
| **CRA to Vite** | index.html location and templates (#6) | Move file, strip `%PUBLIC_URL%`, add module script tag |
| **CRA to Vite** | Dockerfile output directory (#7) | Change `build` to `dist` in COPY command, update all env var names |
| **CRA to Vite** | Test runner dies (#13) | Set up Vitest before removing react-scripts |
| **CRA to Vite** | JSX in .js files (#21) | Rename index.js and App.js to .jsx |
| **CRA to Vite** | Standalone iframe paths (#24) | Replace PUBLIC_URL with root-relative paths |
| **CRA to Vite** | Cleanup (#27) | Remove web-vitals, browserslist, eslintConfig |
| **Tailwind v3 to v4** | theme() function breaks (#2) | Replace with var(--color-*) syntax |
| **Tailwind v3 to v4** | @tailwind directives removed (#3) | Replace with @import "tailwindcss" |
| **Tailwind v3 to v4** | Opacity utilities removed (#4) | Convert 40+ occurrences to slash syntax |
| **Tailwind v3 to v4** | Utility scale renamed (#5) | Run upgrade tool, verify shadow/rounded/ring |
| **Tailwind v3 to v4** | JS config not auto-detected (#8) | Migrate to @theme CSS directive |
| **Tailwind v3 to v4** | Class-based dark mode breaks (#9) | Add @custom-variant dark rule |
| **Tailwind v3 to v4** | Default border color change (#10) | Add compatibility base style |
| **Tailwind v3 to v4** | PostCSS config must change (#14) | Use @tailwindcss/vite plugin instead |
| **Tailwind v3 to v4** | flex-shrink/grow renamed (#16) | Find-replace 60+ occurrences |
| **Tailwind v3 to v4** | ESM require conflict (#17) | Use CSS-first config to avoid |
| **Tailwind v3 to v4** | CSS variable syntax change (#19) | Replace bracket with parentheses syntax |
| **Tailwind v3 to v4** | Global transition rule (#22) | Remove `*` transition-colors |
| **Tailwind v3 to v4** | Button cursor default change (#23) | Add base cursor:pointer rule |
| **Tailwind v3 to v4** | Node.js version requirement (#25) | Upgrade to Node 20+ |
| **React 18 to 19** | Context.Provider deprecated (#12) | One-line change in AppContext.js |
| **React 18 to 19** | Headless UI v1 incompatible (#11) | Remove (shadcn/ui replaces it) |
| **React 18 to 19** | Ref callback returns (#20) | Run codemod for implicit returns |
| **React 18 to 19** | Router future flags (#26) | Test all navigation paths |
| **shadcn/ui** | Requires path aliases (#15) | Configure @/ alias in Vite + jsconfig/tsconfig |
| **shadcn/ui** | CSS variable conflicts (#18) | Reconcile theming systems during init |
| **shadcn/ui** | Global transitions conflict (#22) | Remove `*` transition-colors first |

---

## Recommended Phase Ordering Based on Pitfall Dependencies

1. **CRA to Vite first** -- Pitfalls #1, #6, #7, #13, #17, #21, #24, #27 are all blocked on this. Must be done before anything else because the build pipeline changes. This is where the highest-risk silent failures live (env vars, Dockerfile).

2. **Tailwind v3 to v4 second** -- Pitfalls #2, #3, #4, #5, #8, #9, #10, #14, #16, #19, #22, #23, #25 form the largest and most dangerous cluster. Do this as a dedicated migration step with thorough visual testing. Use the `@tailwindcss/upgrade` tool to automate what it can, then manually fix the CSS variable system and dark mode configuration. This phase has the most pitfalls because Tailwind v4 is the most breaking change of the four migrations.

3. **React 18 to 19 third** -- Pitfalls #11, #12, #20, #26 are relatively contained. This can be done safely after the build pipeline is stable. React 19 is the least breaking migration since no `forwardRef` usage exists in the codebase.

4. **shadcn/ui last** -- Pitfalls #15, #18, #22 require that Vite, Tailwind v4, and path aliases are all in place first. This is an additive step that should not break existing functionality if the theming system is properly reconciled.

---

## Sources

- React 19 release blog: https://react.dev/blog/2024/12/05/react-19 (HIGH confidence -- official, fetched and verified)
- Tailwind CSS v4 upgrade guide: https://tailwindcss.com/docs/upgrade-guide (HIGH confidence -- official, fetched and verified)
- Codebase analysis: exhaustive grep of 55+ source files for affected patterns (HIGH confidence -- direct evidence with file:line references)
- Vite documentation and CRA migration patterns: based on training data (MEDIUM confidence -- verify current Vite docs at time of migration)
- shadcn/ui setup requirements: based on training data (MEDIUM confidence -- verify current shadcn/ui Vite installation docs at time of migration)
