# Architecture Patterns

**Domain:** Frontend modernization (CRA + React 18 + Tailwind v3 -> Vite + React 19 + Tailwind v4 + shadcn/ui)
**Researched:** 2026-03-12

## Current Architecture (What Exists Today)

```
frontend/
  src/
    index.js              <- CRA entry point (ReactDOM.createRoot)
    App.js                <- React Router v6 (createBrowserRouter), Layout wrapper
    context/
      AppContext.js        <- Global state: theme, recentMolecules, apiConfig, loading, error
    hooks/
      useApiCall.js        <- Generic API call wrapper (loading/error states)
      useMolecule.js       <- Domain hook (SMILES validation, descriptors, depiction)
    services/
      api.js               <- Axios instance (REACT_APP_API_URL env var)
      chemService.js       <- /chem endpoint calls
      convertService.js    <- /convert endpoint calls
      depictService.js     <- /depict endpoint calls
      toolsService.js      <- /tools endpoint calls
      ocsrService.js       <- /ocsr endpoint calls
    components/
      common/              <- Header, Footer, MoleculeCard, SMILESInput, etc. (9 files)
      chem/                <- 14 view components (DescriptorsView, ClassyfireView, etc.)
      convert/             <- 3 view components (FormatConversionView, Mol2D/3DView)
      depict/              <- 6 view components (Depict2DMultiView, StructureVisualizerView, etc.)
      tools/               <- 4 view components (InChIView, SugarRemovalView, etc.)
      ocsr/                <- 1 view component (OCRView)
    pages/
      HomePage.js, ChemPage.js, ConvertPage.js, DepictPage.js,
      ToolsPage.js, OCSRPage.js, AboutPage.js, TermsOfService.js, PrivacyPolicy.js
    styles/
      tailwind.css         <- @tailwind directives + CSS variables + custom components/utilities
      animations.css       <- Custom keyframe animations (orbit, fade, slide, etc.)
    utils/
      inchiUtils.js        <- InChI parsing utilities
      rinchiUtils.js       <- RInChI parsing utilities
  public/
    index.html             <- CRA HTML template (with %PUBLIC_URL%)
  tailwind.config.js       <- Tailwind v3 JS config (darkMode: 'class', custom colors, fonts)
  postcss.config.js        <- PostCSS config (tailwindcss + autoprefixer plugins)
  package.json             <- CRA scripts (react-scripts), dependencies
```

### Current Data Flow

```
User Input
  |
  v
Page Component (e.g., ChemPage)
  |
  v
View Component (e.g., DescriptorsView)
  |-- uses SMILESInput (common) for molecule input
  |-- calls service function (e.g., chemService.calculateDescriptors)
  |
  v
Service Layer (services/chemService.js)
  |-- uses Axios instance (services/api.js)
  |-- reads process.env.REACT_APP_API_URL
  |
  v
FastAPI Backend (/latest/chem/descriptors)
  |
  v
Response rendered by View Component
```

### Current Component Hierarchy

```
AppProvider (context: isDarkMode, recentMolecules, apiConfig, isLoading, globalError)
  |
  RouterProvider
    |
    Layout (flex col, bg-gray-900)
      |
      +-- Header (navigation, theme toggle)
      +-- <Outlet /> (page content)
      |     |
      |     Page (e.g., ChemPage)
      |       |
      |       View Components (e.g., DescriptorsView, ClassyfireView)
      |         |
      |         Common Components (SMILESInput, MoleculeCard, LoadingScreen)
      |
      +-- Footer
```

### Current Styling Architecture

Three files control styling today:

1. **`tailwind.config.js`** -- JS-based theme (darkMode: 'class', custom colors gray-750/850/950, custom blue shades, Inter + JetBrains Mono fonts, custom animations, @tailwindcss/forms + @tailwindcss/typography plugins)
2. **`postcss.config.js`** -- PostCSS pipeline (tailwindcss + autoprefixer)
3. **`styles/tailwind.css`** -- CSS variables for light/dark themes in `:root`/`.dark` blocks, custom component classes (.btn, .card, .form-group, .glass, .text-gradient), background patterns (dots, grid, mesh, noise)

The existing CSS variable approach in `tailwind.css` with `:root` and `.dark` selectors is already structurally similar to what shadcn/ui uses. The migration formalizes this pattern.

---

## Recommended Architecture (Target State)

### Directory Structure After Migration

```
frontend/
  index.html                  <- Vite entry point (moved from public/)
  vite.config.ts              <- Vite config + @tailwindcss/vite plugin
  components.json             <- shadcn/ui CLI configuration
  tsconfig.json               <- TypeScript config (path aliases: @/ -> ./src)
  src/
    main.tsx                  <- Vite entry point (replaces index.js)
    App.tsx                   <- React Router, Layout wrapper
    app.css                   <- SINGLE CSS file: Tailwind v4 config + theme + shadcn vars
    context/
      AppContext.tsx           <- Global state (theme, recentMolecules, apiConfig)
    hooks/
      useApiCall.ts            <- Generic API call wrapper
      useMolecule.ts           <- Domain hook (SMILES validation, descriptors)
    services/
      api.ts                   <- Axios instance (import.meta.env.VITE_API_URL)
      chemService.ts           <- /chem endpoint calls
      convertService.ts        <- /convert endpoint calls
      depictService.ts         <- /depict endpoint calls
      toolsService.ts          <- /tools endpoint calls
      ocsrService.ts           <- /ocsr endpoint calls
    components/
      ui/                      <- NEW: shadcn/ui primitives (owned source, not node_modules)
        button.tsx
        card.tsx
        input.tsx
        textarea.tsx
        select.tsx
        tabs.tsx
        dialog.tsx
        dropdown-menu.tsx
        badge.tsx
        separator.tsx
        skeleton.tsx
        tooltip.tsx
        sheet.tsx
        sonner.tsx
        toggle.tsx
        scroll-area.tsx
      common/                  <- App-specific shared components (composed from ui/)
        Header.tsx
        Footer.tsx
        MoleculeCard.tsx
        HighlightedMoleculeCard.tsx
        SMILESInput.tsx
        SMILESDisplay.tsx
        MolFileUpload.tsx
        Navigation.tsx
        LoadingScreen.tsx
        ThemeToggle.tsx        <- NEW: Dark/light mode toggle
      chem/                    <- 14 view components (unchanged structure)
      convert/                 <- 3 view components
      depict/                  <- 6 view components
      tools/                   <- 4 view components
      ocsr/                    <- 1 view component
    pages/                     <- 9 page-level route components (unchanged structure)
    utils/
      inchiUtils.ts
      rinchiUtils.ts
    lib/
      utils.ts                 <- NEW: cn() = clsx + tailwind-merge
```

**What was removed:**
- `tailwind.config.js` -- replaced by `@theme` directive in `app.css`
- `postcss.config.js` -- replaced by `@tailwindcss/vite` plugin
- `styles/tailwind.css` -- merged into `app.css`
- `styles/animations.css` -- merged into `app.css`
- `styles/` directory -- eliminated entirely

### Component Boundaries

| Component Layer | Responsibility | Communicates With | Ownership |
|----------------|---------------|-------------------|-----------|
| **ui/** (shadcn primitives) | Accessible UI primitives (Radix UI underneath). Button, Card, Input, Dialog, etc. Stateless, generic. | Nothing -- receives props only | shadcn CLI generates; project owns source |
| **common/** | App-specific shared components. Header, Footer, SMILESInput, MoleculeCard. Composed from ui/ primitives. | ui/ components, AppContext, services (indirectly) | Project team |
| **{feature}/** (chem, convert, depict, tools, ocsr) | Feature-specific view components. Each handles one API interaction. | common/ components, services/, hooks/ | Project team |
| **pages/** | Route-level orchestrators. Compose feature views, handle URL params. | {feature}/ components, React Router | Project team |
| **services/** | API communication layer. Axios calls to FastAPI backend. No UI knowledge. | api.ts (Axios instance), backend API | Project team |
| **hooks/** | Reusable stateful logic. API call state, molecule operations. | services/, context/ | Project team |
| **context/** | Global application state. Theme, recent molecules, API config. | localStorage (persistence) | Project team |
| **lib/utils.ts** | Utility function `cn()` for merging Tailwind classes. | clsx, tailwind-merge | shadcn convention |

### Data Flow (Post-Migration)

The data flow **does not change**. The change is in _how_ the view layer is rendered: custom Tailwind utility soup is replaced by composable shadcn/ui primitives with consistent styling.

```
User Input
  |
  v
Page Component (e.g., ChemPage.tsx)
  |
  v
View Component (e.g., DescriptorsView.tsx)
  |-- uses SMILESInput (common/) which wraps Input (ui/)
  |-- uses Card, Button, Badge, Skeleton (ui/) for layout
  |-- calls service function (chemService.calculateDescriptors)
  |
  v
Service Layer (services/chemService.ts)
  |-- uses Axios instance (services/api.ts)
  |-- reads import.meta.env.VITE_API_URL
  |
  v
FastAPI Backend (/latest/chem/descriptors)
  |
  v
Response rendered by View Component using ui/ primitives
```

---

## Key Architectural Changes: Before vs After

### 1. Build System: CRA -> Vite

| Aspect | CRA (Before) | Vite (After) |
|--------|-------------|--------------|
| Entry HTML | `public/index.html` (webpack processed) | `index.html` at project root |
| JS entry | `src/index.js` (CRA convention) | `src/main.tsx` (referenced via `<script type="module">` in index.html) |
| Dev server | Webpack Dev Server (slow cold start) | Vite native ESM dev server (instant HMR) |
| Env vars | `REACT_APP_*`, `process.env.REACT_APP_*` | `VITE_*`, `import.meta.env.VITE_*` |
| Config | Hidden in react-scripts | Explicit `vite.config.ts` |
| Build output | `build/` directory | `dist/` directory |
| Plugin system | None (must eject) | Vite plugin array in config |

### 2. CSS Architecture: Tailwind v3 JS Config -> Tailwind v4 CSS-First

This is the most architecturally significant change. Three separate config files collapse into one CSS file.

| Aspect | Tailwind v3 (Before) | Tailwind v4 (After) |
|--------|---------------------|---------------------|
| Config file | `tailwind.config.js` (JavaScript) | `src/app.css` (`@theme` directive) |
| Import syntax | `@tailwind base; @tailwind components; @tailwind utilities;` | `@import "tailwindcss";` (single import) |
| Dark mode config | `darkMode: 'class'` in JS config | `@custom-variant dark (&:where(.dark, .dark *));` in CSS |
| Custom colors | `theme.extend.colors` in JS object | `@theme { --color-gray-750: #2b3544; }` in CSS |
| Custom fonts | `theme.extend.fontFamily` in JS | `@theme { --font-sans: "Inter", ...; }` in CSS |
| Custom animations | `theme.extend.animation` in JS | `@theme { --animate-pulse-slow: pulse 3s ...; }` in CSS |
| Custom utilities | `@layer utilities { .glass { ... } }` | `@utility glass { ... }` |
| PostCSS pipeline | `tailwindcss` + `autoprefixer` plugins | `@tailwindcss/vite` plugin (no PostCSS) |
| TW plugins | `require('@tailwindcss/forms')` | Not needed (shadcn/ui provides form components) |
| Content scanning | `content: ["./src/**/*.{js,jsx}"]` | Automatic (Vite plugin detects sources) |
| Important prefix | `!flex !bg-red-500` (start) | `flex! bg-red-500!` (end) |
| Default border color | `gray-200` | `currentColor` |

**Confidence:** HIGH -- all changes verified against official Tailwind CSS v4 upgrade guide and @theme documentation (fetched from tailwindcss.com).

### 3. Component System: Custom CSS -> shadcn/ui

| Aspect | Before | After |
|--------|--------|-------|
| Button styling | `.btn-primary` class or inline `px-4 py-2 rounded-md bg-blue-600 hover:bg-blue-700...` | `<Button variant="default">` from `ui/button.tsx` |
| Card styling | `.card` class in tailwind.css | `<Card><CardHeader><CardContent>` from `ui/card.tsx` |
| Form inputs | `.form-group`, `.form-label` classes | `<Input>`, `<Textarea>`, `<Select>` from ui/ |
| Loading states | Custom `LoadingScreen.jsx` | `<Skeleton>` from ui/ for content shimmer effects |
| Dropdown/select | `@headlessui/react` | shadcn/ui Select/DropdownMenu (Radix UI underneath) |
| Toast/notifications | None | `<Sonner>` from ui/ |
| Icons | Mixed: `@fortawesome/*` + `react-icons` | `lucide-react` (shadcn default icon set) |
| Class merging | None (class conflicts possible) | `cn()` from `lib/utils.ts` (clsx + tailwind-merge) |

**Key principle:** shadcn/ui is NOT an npm package dependency. The CLI copies component source files into `src/components/ui/`. You own the code and can modify it. The `components.json` file configures paths and styling conventions.

### 4. Dark Mode Architecture

**Before (current):**
```
tailwind.config.js: darkMode: 'class'
AppContext.js: manages isDarkMode state + localStorage persistence
document.documentElement.classList.add("dark") for activation
tailwind.css: CSS variables in :root {} and .dark {} blocks
Components: use dark: variant classes (dark:bg-slate-800, etc.)
```

**After (target):**
```
app.css: @custom-variant dark (&:where(.dark, .dark *));
app.css: shadcn/ui CSS variables in :root {} and .dark {} blocks
AppContext.tsx: SAME pattern (isDarkMode + localStorage) -- it works, keep it
Components: use semantic classes (bg-background, text-foreground, bg-card)
Components: still can use dark: variant when needed
```

The dark mode _mechanism_ stays the same (class on `<html>`). What changes is the _theming system_ -- shadcn/ui's semantic CSS variables (`--background`, `--foreground`, `--primary`, `--card`, `--muted`, etc.) replace the ad-hoc custom variables in `tailwind.css`.

This means most components STOP needing `dark:` variants because `bg-background` already resolves differently in light vs dark mode.

**Confidence for @custom-variant syntax:** HIGH -- verified from official Tailwind CSS v4 dark mode docs.

### 5. The `app.css` File: New Configuration Hub

Three config files (`tailwind.config.js`, `postcss.config.js`, `styles/tailwind.css`) collapse into one. Here is the target structure with actual syntax from Tailwind v4 docs:

```css
/* src/app.css */

/* 1. Import Tailwind (replaces three @tailwind directives) */
@import "tailwindcss";

/* 2. Dark mode (replaces darkMode: 'class' in JS config) */
@custom-variant dark (&:where(.dark, .dark *));

/* 3. Theme tokens (replaces tailwind.config.js theme.extend) */
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

  /* Fonts from current config */
  --font-sans: "Inter", ui-sans-serif, system-ui, -apple-system,
    BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial,
    "Noto Sans", sans-serif;
  --font-mono: "JetBrains Mono", ui-monospace, SFMono-Regular, Menlo,
    Monaco, Consolas, "Liberation Mono", "Courier New", monospace;

  /* Custom animations */
  --animate-pulse-slow: pulse 3s cubic-bezier(0.4, 0, 0.6, 1) infinite;
  --animate-pulse-fast: pulse 1s cubic-bezier(0.4, 0, 0.6, 1) infinite;
  --animate-bounce-slow: bounce 2s infinite;

  /* Custom spacing */
  --spacing-128: 32rem;

  /* Custom box shadows */
  --shadow-inner-lg: inset 0 2px 4px 0 rgba(0, 0, 0, 0.2);
  --shadow-blue-glow: 0 0 15px rgba(66, 153, 225, 0.5);

  /* Custom border radius */
  --radius-3xl: 3rem;
}

/* 4. shadcn/ui CSS variables (generated by npx shadcn@latest init) */
@layer base {
  :root {
    --background: 210 40% 96%;
    --foreground: 222 47% 11%;
    --card: 0 0% 100%;
    --card-foreground: 222 47% 11%;
    --primary: 199 89% 48%;        /* Sky-500 matching current theme */
    --primary-foreground: 0 0% 100%;
    --secondary: 210 40% 96%;
    --muted: 210 40% 96%;
    --muted-foreground: 215 16% 47%;
    --accent: 210 40% 96%;
    --border: 214 32% 91%;
    --input: 214 32% 91%;
    --ring: 199 89% 48%;
    --radius: 0.5rem;
  }
  .dark {
    --background: 222 47% 6%;
    --foreground: 210 40% 98%;
    --card: 217 33% 17%;
    --card-foreground: 210 40% 98%;
    --primary: 199 89% 48%;
    --primary-foreground: 0 0% 100%;
    --secondary: 217 33% 17%;
    --muted: 217 33% 17%;
    --muted-foreground: 215 20% 65%;
    --accent: 217 33% 17%;
    --border: 217 33% 25%;
    --input: 217 33% 25%;
    --ring: 199 89% 48%;
  }

  body {
    @apply min-h-screen antialiased bg-background text-foreground;
    font-family: var(--font-sans);
  }
}

/* 5. Custom utilities (replaces @layer utilities from old tailwind.css) */
@utility glass {
  @apply bg-white/70 backdrop-blur-lg border border-slate-300/50;
  @apply dark:bg-slate-800/70 dark:backdrop-blur-xl dark:border-slate-700/50;
}

@utility text-gradient {
  @apply text-transparent bg-clip-text bg-gradient-to-r from-sky-600 to-indigo-600;
  @apply dark:from-sky-400 dark:to-cyan-400;
}

/* 6. Animations (migrated from animations.css) */
@utility animate-orbit-1 {
  animation: orbit-1 4s linear infinite;
}
/* ... remaining orbit/fade/slide utilities ... */

@keyframes orbit-1 {
  0% { transform: translate(-50%, 0) rotate(0deg) translateY(-28px) rotate(0deg); }
  100% { transform: translate(-50%, 0) rotate(360deg) translateY(-28px) rotate(-360deg); }
}
/* ... remaining keyframes ... */
```

**Confidence:** HIGH for `@import "tailwindcss"`, `@custom-variant`, `@theme` syntax -- all verified from official docs. MEDIUM for exact shadcn/ui CSS variable values -- these will be generated by `npx shadcn@latest init` and fine-tuned visually.

---

## Vite Config Structure (Target)

```typescript
// vite.config.ts
import { defineConfig } from "vite"
import react from "@vitejs/plugin-react"
import tailwindcss from "@tailwindcss/vite"
import path from "path"

export default defineConfig({
  plugins: [
    react(),
    tailwindcss(),
  ],
  resolve: {
    alias: {
      "@": path.resolve(__dirname, "./src"),
    },
  },
  server: {
    port: 3000,        // Match current CRA port for Docker/dev compatibility
    host: "0.0.0.0",   // Required for Docker
  },
  build: {
    outDir: "dist",
    sourcemap: true,
  },
})
```

---

## Migration Order (Build Dependencies)

Each phase depends on the previous one. This order is non-negotiable.

### Phase 1: CRA -> Vite (Build System)
**Must come first.** Everything else requires a working Vite build.

| Step | Action | Why |
|------|--------|-----|
| 1a | Install `vite`, `@vitejs/plugin-react`; create `vite.config.ts` | Build tool replacement |
| 1b | Move `public/index.html` to project root | Vite convention |
| 1c | Update `index.html`: remove `%PUBLIC_URL%`, add `<script type="module" src="/src/main.jsx">` | Vite entry point |
| 1d | Rename `src/index.js` -> `src/main.jsx` | Vite convention |
| 1e | Replace `process.env.REACT_APP_*` with `import.meta.env.VITE_*` in `api.js`, `AppContext.js` | Env var format change |
| 1f | Update `.env.development`, `.env.production`: rename env var keys | Match Vite prefix |
| 1g | Update `package.json` scripts: `vite` / `vite build` / `vite preview` | Build commands |
| 1h | Remove `react-scripts` from dependencies | CRA removal |
| 1i | Update `Dockerfile.frontend`: output `dist/` not `build/` | Docker build fix |
| 1j | **VERIFY:** app builds, dev server works, all pages load, API calls succeed | Gate |

### Phase 2: Tailwind v3 -> Tailwind v4 (CSS Architecture)
**Must come before shadcn/ui.** The shadcn/ui init CLI reads Tailwind version.

| Step | Action | Why |
|------|--------|-----|
| 2a | Install `tailwindcss@4`, `@tailwindcss/vite` | New Tailwind |
| 2b | Remove `tailwindcss@3`, `autoprefixer`, `postcss`, `@tailwindcss/forms`, `@tailwindcss/typography` | Old deps |
| 2c | Add `tailwindcss()` plugin to `vite.config.ts` | Vite integration |
| 2d | Delete `postcss.config.js` | Replaced by Vite plugin |
| 2e | Create `src/app.css` with `@import "tailwindcss"` + `@custom-variant dark` + `@theme` | CSS-first config |
| 2f | Migrate `tailwind.config.js` theme values to `@theme` directive | Config format change |
| 2g | Migrate `styles/tailwind.css` custom CSS into `app.css` | Consolidate |
| 2h | Migrate `styles/animations.css` keyframes into `app.css` | Consolidate |
| 2i | Delete `tailwind.config.js`, `postcss.config.js`, `styles/` directory | Cleanup |
| 2j | Run `npx @tailwindcss/upgrade` to auto-fix renamed utilities in components | Automated fixes |
| 2k | Fix remaining breaking changes: `!` prefix -> suffix, border-color defaults | Manual fixes |
| 2l | Update `main.jsx` to import `./app.css` | Entry point |
| 2m | **VERIFY:** all pages render, dark mode works, custom colors work, animations work | Gate |

### Phase 3: shadcn/ui Integration (Component Primitives)
**Depends on Vite + Tailwind v4.**

| Step | Action | Why |
|------|--------|-----|
| 3a | Set up path aliases: `@/` -> `./src` in `tsconfig.json` + `vite.config.ts` | shadcn/ui requires @ imports |
| 3b | Install: `tailwind-merge`, `clsx`, `class-variance-authority`, `lucide-react` | shadcn dependencies |
| 3c | Run `npx shadcn@latest init` | Generates `components.json`, `lib/utils.ts`, updates CSS |
| 3d | Add core: `npx shadcn@latest add button card input textarea select label` | Essential primitives |
| 3e | Add layout: `npx shadcn@latest add tabs dialog dropdown-menu sheet separator` | Navigation/layout |
| 3f | Add feedback: `npx shadcn@latest add badge skeleton tooltip sonner scroll-area` | Status/loading |
| 3g | **VERIFY:** shadcn/ui components render in both light/dark modes | Gate |

### Phase 4: Common Component Refactor
**Depends on shadcn/ui primitives.**

| Step | Action | Why |
|------|--------|-----|
| 4a | Create `ThemeToggle.tsx` using Button + lucide Sun/Moon | New component |
| 4b | Refactor `Header.tsx` to use shadcn Button, NavigationMenu | Replace hand-rolled nav |
| 4c | Refactor `Footer.tsx` | Consistent styling |
| 4d | Refactor `SMILESInput.tsx` to use Input + Button | Core interaction component |
| 4e | Refactor `MoleculeCard.tsx` to use Card, Badge, Button | Most-used display component |
| 4f | Refactor `MolFileUpload.tsx` | File upload UX |
| 4g | Replace `@headlessui/react` with shadcn equivalents | Remove old dep |
| 4h | Replace `@fortawesome/*` + `react-icons` with `lucide-react` | Unify icons |
| 4i | **VERIFY:** all common components work in both modes | Gate |

### Phase 5: Feature View Migration
**Depends on common components being stable.**

| Step | Action | Why |
|------|--------|-----|
| 5a | Migrate `chem/` (14 views) -- use ui/ primitives | Largest section |
| 5b | Migrate `depict/` (6 views) | Visual-heavy section |
| 5c | Migrate `convert/` (3 views) | Medium section |
| 5d | Migrate `tools/` (4 views) | Medium section |
| 5e | Migrate `ocsr/` (1 view) | Smallest section |
| 5f | Migrate pages/ (9 pages) | Layout consistency |
| 5g | **VERIFY:** all feature views function correctly | Gate |

### Phase 6: React 19 + TypeScript + Optimization
**Can be done after core migration is stable.**

| Step | Action | Why |
|------|--------|-----|
| 6a | Upgrade React 18 -> 19 | Performance, new hooks |
| 6b | Rename .js/.jsx -> .ts/.tsx | Type safety |
| 6c | Add types for service responses, context, hooks | Gradual typing |
| 6d | Add React.lazy() + Suspense for page code splitting | Performance |
| 6e | Remove unused deps: `@headlessui/react`, `@fortawesome/*`, `react-icons`, `web-vitals` | Cleanup |
| 6f | Remove `package.json` overrides (no longer needed) | Cleanup |
| 6g | Final `npm audit` -- target zero vulnerabilities | Security |

---

## Patterns to Follow

### Pattern 1: shadcn/ui Composition
**What:** Build app-specific components by composing shadcn/ui primitives. Never put app logic in `ui/` files.
**When:** Any shared or feature component.

```tsx
// GOOD: common/SMILESInput.tsx composes ui/ primitives
import { Input } from "@/components/ui/input"
import { Button } from "@/components/ui/button"
import { Label } from "@/components/ui/label"

export function SMILESInput({ value, onChange, onSubmit, isLoading }) {
  return (
    <div className="flex gap-2">
      <div className="flex-1">
        <Label htmlFor="smiles">SMILES</Label>
        <Input id="smiles" placeholder="Enter SMILES..." value={value}
          onChange={(e) => onChange(e.target.value)} className="font-mono" />
      </div>
      <Button onClick={onSubmit} disabled={isLoading}>
        {isLoading ? "Loading..." : "Submit"}
      </Button>
    </div>
  )
}
```

### Pattern 2: cn() for Class Merging
**What:** Always use `cn()` from `lib/utils.ts` for conditional/merged Tailwind classes.
**When:** Any component accepting `className` props or having conditional styles.

```tsx
import { cn } from "@/lib/utils"

export function MoleculeCard({ className, isHighlighted, children }) {
  return (
    <Card className={cn(
      "transition-shadow",
      isHighlighted && "ring-2 ring-primary",
      className
    )}>
      {children}
    </Card>
  )
}
```

### Pattern 3: Semantic Color Variables
**What:** Use `bg-background`, `text-foreground`, `bg-card`, `text-muted-foreground` instead of hardcoded colors.
**When:** Always. This makes dark mode automatic.

```tsx
// GOOD: Adapts to dark mode via CSS variables
<div className="bg-background text-foreground border-border">

// AVOID: Requires dark: variant for every color
<div className="bg-white text-gray-900 dark:bg-gray-900 dark:text-white">
```

### Pattern 4: Service Layer Unchanged
**What:** Only env var access changes in the service layer.
**When:** Phase 1 (Vite migration).

```typescript
// ONLY this line changes in api.ts:
const API_URL = import.meta.env.VITE_API_URL || "https://dev.api.naturalproducts.net/latest";
```

Axios interceptors, service function signatures, error handling -- all unchanged.

---

## Anti-Patterns to Avoid

### Anti-Pattern 1: Modifying ui/ Files for App Logic
**What:** Adding domain logic (API calls, state management) directly into `components/ui/button.tsx`.
**Why bad:** Blocks future shadcn/ui updates. Breaks primitive/composition boundary.
**Instead:** Compose in `common/` or feature directories.

### Anti-Pattern 2: Mixing Old and New Styling Systems
**What:** Some components using `.card` CSS class, others using `<Card>`, others using raw Tailwind for the same pattern.
**Why bad:** Inconsistent theming. Maintenance nightmare.
**Instead:** Once a shadcn/ui component exists, migrate ALL uses of that pattern. Remove the old CSS class.

### Anti-Pattern 3: Keeping PostCSS Config with Vite
**What:** Using `@tailwindcss/postcss` instead of `@tailwindcss/vite`.
**Why bad:** Slower, more complex. Vite plugin is purpose-built.
**Instead:** Use `@tailwindcss/vite` plugin. Delete `postcss.config.js`.

### Anti-Pattern 4: npm install shadcn-ui
**What:** Trying to install shadcn/ui as a package dependency.
**Why bad:** It is not a package. It is a code generator CLI.
**Instead:** `npx shadcn@latest add <component>` copies source files into your project.

### Anti-Pattern 5: Multiple Icon Libraries
**What:** Keeping `@fortawesome/react-fontawesome`, `react-icons`, AND `lucide-react`.
**Why bad:** Inconsistent sizing/weight. Triple bundle cost.
**Instead:** Consolidate to `lucide-react`. It covers all common icons.

### Anti-Pattern 6: Big-Bang TypeScript Migration
**What:** Renaming all files to .tsx and enabling strict mode simultaneously.
**Why bad:** Hundreds of type errors at once. Cannot verify other changes work.
**Instead:** Set `"strict": false` initially. Rename files. Gradually enable strict checks.

---

## Critical File Mapping

| Current File | Target File | Change |
|-------------|------------|--------|
| `public/index.html` | `index.html` (root) | Moved, simplified |
| `src/index.js` | `src/main.tsx` | Renamed, typed |
| `src/App.js` | `src/App.tsx` | Renamed, typed |
| `tailwind.config.js` | `src/app.css` (`@theme`) | JS -> CSS |
| `postcss.config.js` | DELETED | Replaced by Vite plugin |
| `src/styles/tailwind.css` | `src/app.css` | Merged |
| `src/styles/animations.css` | `src/app.css` | Merged |
| N/A | `vite.config.ts` | New |
| N/A | `components.json` | New (shadcn config) |
| N/A | `src/lib/utils.ts` | New (cn function) |
| N/A | `src/components/ui/*.tsx` | New (shadcn primitives) |
| `.env.development` | `.env.development` | `REACT_APP_*` -> `VITE_*` |
| `.env.production` | `.env.production` | `REACT_APP_*` -> `VITE_*` |

## Docker Build Impact

Minimal but required change:

```dockerfile
# Output directory changes from build/ to dist/
COPY --from=build /app/dist /usr/share/nginx/html
# (was: COPY --from=build /app/build /usr/share/nginx/html)
```

The nginx config, multi-stage build pattern, and port mapping remain identical.

---

## Sources

- Tailwind CSS v4 Upgrade Guide: https://tailwindcss.com/docs/upgrade-guide (HIGH confidence, fetched from official docs)
- Tailwind CSS v4 Dark Mode: https://tailwindcss.com/docs/dark-mode (HIGH confidence, fetched from official docs)
- Tailwind CSS v4 @theme Directive: https://tailwindcss.com/docs/theme (HIGH confidence, fetched from official docs)
- Tailwind CSS v4 Compatibility: https://tailwindcss.com/docs/compatibility (HIGH confidence, fetched from official docs)
- shadcn/ui Architecture: https://ui.shadcn.com/docs (MEDIUM confidence -- training data; WebFetch blocked; verify during `npx shadcn@latest init`)
- Vite Project Structure: https://vite.dev/guide/ (MEDIUM confidence -- training data; WebFetch blocked; verify during setup)
- React 19: https://react.dev/blog/2024/12/05/react-19 (MEDIUM confidence, training data)
- Current codebase: Direct file inspection of all source files (HIGH confidence)
