# Technology Stack

**Project:** Cheminformatics Microservice Frontend Modernization
**Researched:** 2026-03-12

## Current Stack (Baseline)

| Technology | Version | Status |
|------------|---------|--------|
| React | 18.2.0 | Outdated (19.x current) |
| react-scripts (CRA) | 5.0.1 | **Deprecated, unmaintained** |
| Tailwind CSS | 3.3.5 | Outdated (v4.x current) |
| TypeScript | 4.9.5 | Installed but **not used** (all files are .js/.jsx) |
| react-router-dom | 6.30.3 | Current, keep |
| framer-motion | 12.6.2 | Current, keep |
| axios | 1.13.5 | Current, keep |
| react-dropzone | 14.2.3 | Current, keep |
| react-icons | 4.12.0 | Current, keep |
| @headlessui/react | 1.7.17 | Installed but **not imported anywhere** |
| @fortawesome/* | 6.7.2 / 0.2.2 | Used in 1 file only (StructureVisualizerView.jsx) |
| dompurify | 3.3.1 | Current, keep |
| jszip | 3.10.1 | Current, keep |
| @tailwindcss/forms | 0.5.7 | Needs v4 migration |
| @tailwindcss/typography | 0.5.10 | Needs v4 migration |
| autoprefixer | 10.4.16 | **Remove** (Tailwind v4 handles this) |
| postcss | 8.4.31 | **Remove** (Tailwind v4 Vite plugin replaces this) |
| prettier | 3.6.2 | Current, keep |

## Recommended Stack

### Build Tool

| Technology | Version | Purpose | Why |
|------------|---------|---------|-----|
| Vite | ^6.x | Build tool, dev server, HMR | Replaces deprecated CRA. 10-100x faster cold start than webpack. Native ESM. Official React recommendation. **Confidence: HIGH** (verified: Tailwind docs show `@tailwindcss/vite` as primary integration path) |
| @vitejs/plugin-react | ^4.x | React Fast Refresh for Vite | Official Vite plugin for React. Uses Babel for JSX transform. Use this over plugin-react-swc because shadcn/ui and the ecosystem assumes Babel compatibility. **Confidence: HIGH** |

### Core Framework

| Technology | Version | Purpose | Why |
|------------|---------|---------|-----|
| React | ^19.2.0 | UI library | Latest stable. Includes useActionState, useFormStatus, useOptimistic, `use()` API, ref-as-prop (no more forwardRef), `<Context>` as provider, document metadata hoisting. React Compiler v1.0 available. **Confidence: HIGH** (verified: React blog confirms 19.2 released Oct 2025) |
| react-dom | ^19.2.0 | React DOM renderer | Must match React version. **Confidence: HIGH** |

### Styling

| Technology | Version | Purpose | Why |
|------------|---------|---------|-----|
| tailwindcss | ^4.2 | Utility-first CSS | CSS-first configuration (no more tailwind.config.js). `@theme` blocks replace JS config. `@import "tailwindcss"` replaces `@tailwind` directives. **Confidence: HIGH** (verified: Tailwind installation docs show v4.2) |
| @tailwindcss/vite | ^4.2 | Vite integration | First-class Vite plugin, replaces PostCSS-based setup. Faster than PostCSS path. **Confidence: HIGH** (verified: Tailwind installation docs) |

### Component Library

| Technology | Version | Purpose | Why |
|------------|---------|---------|-----|
| shadcn/ui | latest (CLI-based) | Pre-built, customizable components | Not a dependency - components are copied into your project via CLI. Built on Radix UI primitives + Tailwind CSS. Supports React 19 + Tailwind v4. Gives you ownership of component code. Perfect for dark mode (uses CSS variables). **Confidence: MEDIUM** (based on training data + ecosystem knowledge; shadcn/ui updated for TW v4 support in early 2025) |

shadcn/ui dependencies that will be installed:

| Library | Purpose | Why |
|---------|---------|-----|
| @radix-ui/* | Accessible UI primitives | Underlying primitives for shadcn/ui components (Dialog, Dropdown, Tooltip, etc.) |
| class-variance-authority (cva) | Component variant management | Defines component variants declaratively |
| clsx | Conditional class joining | Standard utility for conditional Tailwind classes |
| tailwind-merge | Tailwind class deduplication | Prevents conflicting Tailwind classes (e.g., `p-2 p-4` resolves to `p-4`) |
| lucide-react | Icon library | shadcn/ui's default icon set. Replaces react-icons and FontAwesome |

### Routing (Keep)

| Technology | Version | Purpose | Why |
|------------|---------|---------|-----|
| react-router-dom | ^6.30.3 | Client-side routing | Already current. Using `createBrowserRouter` data API. No change needed. **Confidence: HIGH** |

### HTTP Client (Keep)

| Technology | Version | Purpose | Why |
|------------|---------|---------|-----|
| axios | ^1.13.5 | HTTP client | Already current. Well-configured with interceptors. Keep as-is. **Confidence: HIGH** |

### Animation (Keep)

| Technology | Version | Purpose | Why |
|------------|---------|---------|-----|
| framer-motion | ^12.x | Animations | Already current. Widely used. Works with React 19. **Confidence: HIGH** |

### Dark Mode

| Technology | Version | Purpose | Why |
|------------|---------|---------|-----|
| (built-in) | N/A | Dark mode toggle | Tailwind v4 dark mode via `@media (prefers-color-scheme: dark)` or `dark:` variant with class strategy. shadcn/ui uses CSS variables for theming. The existing `AppContext` already has dark mode state + localStorage persistence. **Confidence: HIGH** |

### Testing

| Technology | Version | Purpose | Why |
|------------|---------|---------|-----|
| vitest | ^3.x | Test runner | Replaces Jest (CRA bundled). Native Vite integration, same config, faster. **Confidence: MEDIUM** (version from training data) |
| @testing-library/react | ^16.x | React testing utilities | Upgrade to v16 for React 19 compatibility. **Confidence: MEDIUM** |
| @testing-library/jest-dom | ^6.x | DOM matchers | Keep using, works with Vitest. **Confidence: HIGH** |
| jsdom | latest | DOM environment for tests | Vitest needs this configured as test environment. **Confidence: HIGH** |

### Linting and Formatting

| Technology | Version | Purpose | Why |
|------------|---------|---------|-----|
| eslint | ^9.x | Linting | Flat config format. Replace CRA's eslint-config-react-app. **Confidence: MEDIUM** |
| prettier | ^3.6.2 | Code formatting | Already current. Keep. **Confidence: HIGH** |

## Dependencies to REMOVE

| Package | Why Remove |
|---------|-----------|
| react-scripts | Deprecated CRA. Replaced by Vite. |
| @headlessui/react | Installed but never imported. Dead dependency. |
| @fortawesome/free-solid-svg-icons | Used in 1 file. Replace with lucide-react (shadcn/ui default). |
| @fortawesome/react-fontawesome | Used in 1 file. Replace with lucide-react. |
| autoprefixer | Tailwind v4 handles vendor prefixing internally. |
| postcss | Tailwind v4 Vite plugin replaces PostCSS pipeline. |
| @tailwindcss/forms | Absorbed into Tailwind v4 core or replaced by shadcn/ui form components. |
| @tailwindcss/typography | Must check v4 compatibility. May need `@tailwindcss/typography` v4 version or replace with prose styles. |
| web-vitals | CRA-specific. Remove unless explicitly needed. |
| typescript | Not used (no .ts/.tsx files). Remove to avoid confusion. Add back only if migrating to TypeScript. |

## Dependencies to ADD

```bash
# Core (install in frontend/)
npm install react@^19.2.0 react-dom@^19.2.0

# Build tooling
npm install -D vite@^6 @vitejs/plugin-react@^4

# Tailwind CSS v4
npm install -D tailwindcss@^4.2 @tailwindcss/vite@^4.2

# shadcn/ui utilities (installed by shadcn CLI, listed for reference)
npm install class-variance-authority clsx tailwind-merge lucide-react

# Radix UI primitives (installed per-component by shadcn CLI)
# e.g., npm install @radix-ui/react-dialog @radix-ui/react-dropdown-menu etc.

# Testing
npm install -D vitest@^3 @testing-library/react@^16 @testing-library/jest-dom@^6 jsdom
```

## Alternatives Considered

| Category | Recommended | Alternative | Why Not |
|----------|-------------|-------------|---------|
| Build tool | Vite | Turbopack, Rspack | Vite has the largest ecosystem, best plugin support, and is the official Tailwind/shadcn recommendation. Turbopack is Next.js-specific. Rspack is newer with smaller ecosystem. |
| Component library | shadcn/ui | Radix Themes, Mantine, Chakra UI | shadcn/ui gives you source code ownership, is Tailwind-native (no runtime CSS-in-JS), and is the most popular Tailwind component system. Radix Themes is opinionated on styling. Mantine/Chakra use CSS-in-JS which conflicts with Tailwind approach. |
| Icons | lucide-react | react-icons, heroicons | lucide-react is shadcn/ui's default. Consistent design language. Tree-shakeable. react-icons is already in the project but lucide integrates better with shadcn. However, react-icons can be kept alongside if migration is gradual. |
| CSS framework | Tailwind v4 | Tailwind v3 (stay) | v4 is stable, has better performance, eliminates PostCSS dependency with Vite plugin, and is required for shadcn/ui's latest templates. Staying on v3 means maintaining a deprecated config approach. |
| Test runner | Vitest | Jest | Vitest is Vite-native, uses the same config, faster, and designed for Vite projects. Jest requires separate webpack/babel config after removing CRA. |
| Dark mode | Tailwind class strategy + CSS vars | next-themes | next-themes is Next.js-focused. The existing AppContext already handles dark mode state. shadcn/ui's CSS variable approach works perfectly with Tailwind's class-based dark mode. |
| React version | React 19 | Stay on React 18 | React 19 is stable for 15+ months. Breaking changes are minimal for this codebase (no class components, no string refs, no legacy context, already using createRoot). The automated codemod handles most migration. |

## Migration Order (Critical)

The upgrades have dependencies. This order prevents broken intermediate states:

```
Phase 1: CRA -> Vite (foundational, everything else depends on this)
   |
   v
Phase 2: React 18 -> React 19 (independent of Tailwind, but do after Vite)
   |
   v
Phase 3: Tailwind v3 -> v4 (uses @tailwindcss/vite plugin, requires Vite)
   |
   v
Phase 4: Add shadcn/ui (requires Tailwind v4 + React 19)
   |
   v
Phase 5: Dark mode polish (requires shadcn/ui CSS variable system)
```

**Why this order:**

1. **Vite first** because Tailwind v4's best integration is the `@tailwindcss/vite` plugin, and all subsequent tooling (Vitest, shadcn CLI) assumes a Vite project.

2. **React 19 second** because it is independent of Tailwind version but shadcn/ui's latest components assume React 19 APIs.

3. **Tailwind v4 third** because shadcn/ui's latest CLI and component templates target Tailwind v4's CSS variable system.

4. **shadcn/ui fourth** because it requires both React 19 and Tailwind v4 to be in place for its components to work correctly.

5. **Dark mode last** because it depends on shadcn/ui's CSS variable theming system being fully integrated.

## Env Variable Migration

CRA uses `REACT_APP_*` prefix. Vite uses `VITE_*` prefix.

| Current | New | Where Used |
|---------|-----|-----------|
| `REACT_APP_API_URL` | `VITE_API_URL` | `services/api.js`, `context/AppContext.js` |
| `process.env.REACT_APP_*` | `import.meta.env.VITE_*` | All env var references |
| `process.env.NODE_ENV` | `import.meta.env.MODE` | Error logging in api.js |

Files that need env var updates:
- `frontend/src/services/api.js` (2 references)
- `frontend/src/context/AppContext.js` (1 reference)
- `frontend/Dockerfile` (ARG/ENV names)
- `frontend/.env.production` (if exists)

## Docker Changes Required

Current Dockerfile uses:
- `node:18-alpine` -- upgrade to `node:20-alpine` (Tailwind v4 upgrade tool requires Node 20+)
- `ARG REACT_APP_API_URL` -- change to `ARG VITE_API_URL`
- `COPY --from=build /app/build` -- change to `COPY --from=build /app/dist` (Vite outputs to `dist/`, not `build/`)

## Vite Configuration (Reference)

```javascript
// vite.config.js
import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'
import tailwindcss from '@tailwindcss/vite'
import path from 'path'

export default defineConfig({
  plugins: [
    react(),
    tailwindcss(),
  ],
  resolve: {
    alias: {
      '@': path.resolve(__dirname, './src'),
    },
  },
  server: {
    port: 3000,       // Match current CRA port
    open: true,
  },
  build: {
    outDir: 'dist',
    sourcemap: false,  // Disable in production
  },
})
```

## Tailwind v4 CSS Configuration (Reference)

Replace `tailwind.config.js` with CSS-first config:

```css
/* src/styles/tailwind.css */
@import "tailwindcss";

@theme {
  /* Custom colors from current tailwind.config.js */
  --color-gray-750: #2b3544;
  --color-gray-850: #1a2231;
  --color-gray-950: #0f1521;
  --color-blue-350: #7dabf8;
  --color-blue-450: #4f85e6;
  --color-blue-550: #3a72d6;
  --color-blue-650: #2a5db8;
  --color-blue-750: #214a98;
  --color-blue-850: #183878;

  /* Custom fonts */
  --font-sans: "Inter", ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, "Noto Sans", sans-serif;
  --font-mono: "JetBrains Mono", ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, "Liberation Mono", "Courier New", monospace;

  /* Custom spacing */
  --spacing-72: 18rem;
  --spacing-84: 21rem;
  --spacing-96: 24rem;
  --spacing-128: 32rem;
}
```

## Key Breaking Changes to Address

### Tailwind v3 -> v4 (verified from official upgrade guide)

| Change | Impact on This Project | Action |
|--------|----------------------|--------|
| `@tailwind base/components/utilities` removed | `tailwind.css` uses all three | Replace with `@import "tailwindcss"` |
| `bg-opacity-*` removed | Check all files for opacity utilities | Use `bg-black/50` slash syntax instead |
| `shadow-sm` renamed to `shadow-xs` | Check all shadow usage | Update class names |
| `rounded-sm` renamed to `rounded-xs` | Check all border-radius usage | Update class names |
| `blur-sm` renamed to `blur-xs` | Check all blur usage | Update class names |
| `ring` (3px) renamed to `ring-3` | Check all ring usage | Update class names |
| `outline-none` renamed to `outline-hidden` | Check all outline usage | Update class names |
| `!` important prefix moves to suffix | Check all `!` usage | `!flex` becomes `flex!` |
| `tailwind.config.js` not auto-detected | Large config exists | Migrate to `@theme` in CSS or use `@config` directive |
| `@layer utilities` becomes `@utility` | Custom utilities in tailwind.css | Refactor `@layer utilities` blocks |
| `@layer components` becomes `@layer` | Custom components in tailwind.css | Keep `@layer components` (still works) or use `@utility` |
| Default border color changed | Could affect all borders | Add base layer override if needed |
| Variant stacking order reversed | Check compound variants like `first:*:pt-0` | Reverse order: `*:first:pt-0` |

### React 18 -> 19 (verified from official upgrade guide)

| Change | Impact on This Project | Action |
|--------|----------------------|--------|
| `propTypes` removed for functions | Not used (no propTypes in codebase) | None |
| `defaultProps` removed for functions | Not used | None |
| String refs removed | Not used (no class components) | None |
| `forwardRef` deprecated | Not used in current codebase | None |
| `<Context.Provider>` deprecated | Used in `AppContext.js` | Change to `<AppContext value={...}>` |
| `act` import moved | Testing files | Import from `react` not `react-dom/test-utils` |
| `useRef` requires argument | Not using TypeScript, low risk | Add `null` argument if issues arise |
| `process.env.NODE_ENV` checks | Used in api.js | Already compatible, but move to `import.meta.env.MODE` with Vite |

**Automated migration:** Run `npx codemod@latest react/19/migration-recipe` to handle most changes automatically.

**Automated Tailwind migration:** Run `npx @tailwindcss/upgrade` to handle class renames, config migration, and dependency updates. Requires Node.js 20+.

## Sources

- React 19 release blog: https://react.dev/blog/2024/12/05/react-19 (verified via WebFetch, **HIGH** confidence)
- React 19.2 release: https://react.dev/blog (verified via WebFetch, confirms React 19.2 released Oct 2025, **HIGH** confidence)
- React 19 upgrade guide: https://react.dev/blog/2024/04/25/react-19-upgrade-guide (verified via WebFetch, **HIGH** confidence)
- Tailwind CSS v4.2 installation: https://tailwindcss.com/docs/installation (verified via WebFetch, **HIGH** confidence)
- Tailwind CSS v3->v4 upgrade guide: https://tailwindcss.com/docs/upgrade-guide (verified via WebFetch, **HIGH** confidence)
- Vite version: Training data indicates Vite 6.x current as of early 2026 (**MEDIUM** confidence -- verify exact version with `npm view vite version`)
- shadcn/ui: Training data indicates React 19 + Tailwind v4 support (**MEDIUM** confidence -- verify at https://ui.shadcn.com/docs/installation/vite)
- Vitest version: Training data indicates v3.x current (**MEDIUM** confidence -- verify with `npm view vitest version`)
