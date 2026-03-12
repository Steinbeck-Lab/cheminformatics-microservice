# Phase 4: Component System + Dark Mode - Research

**Researched:** 2026-03-12
**Domain:** shadcn/ui component system, CSS variable theming, dark mode, icon migration
**Confidence:** HIGH

## Summary

Phase 4 replaces all hand-rolled UI primitives (buttons, cards, form inputs, dialogs) with shadcn/ui components, migrates the CSS variable system from custom names to shadcn/ui's standard OKLCH-based variables, swaps all icons from react-icons (Heroicons + FontAwesome) to lucide-react, and adds system preference detection with FOUC prevention for dark mode.

The codebase has **46 files** importing from react-icons (43 from `react-icons/hi`, 1 from `react-icons/hi2`, 3 from `react-icons/fa`), plus **1 file** using `@fortawesome/react-fontawesome`. There are **290 native HTML element instances** (`<button>`, `<input>`, `<select>`, `<textarea>`) across 38 files that need evaluation for shadcn/ui component replacement. The existing CSS variable system uses 15 custom variables (`--bg-primary`, `--text-primary`, etc.) referenced in 7 files (110 occurrences) that must be remapped to shadcn/ui's standard naming. Four tabbed pages (ChemPage, DepictPage, ConvertPage, ToolsPage) use custom `LayoutGroup`-based animated tab indicators.

**Primary recommendation:** Initialize shadcn/ui with `npx shadcn@latest init`, use OKLCH color format (the current shadcn/ui default for Tailwind v4), keep custom tab implementations with LayoutGroup animations rather than replacing with shadcn Tabs (the existing animation system is deeply integrated and higher quality), and use a blocking `<script>` in index.html for FOUC prevention.

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions
- Full replacement: remove all custom --bg-primary, --text-primary, --border-primary, etc. CSS variables
- Adopt shadcn/ui's standard variable naming (--background, --foreground, --card, --primary, --muted, --border, --input, --ring, etc.)
- Map current sky/slate color palette into the shadcn/ui variable system -- the app should look the same, just with standardized variable names underneath
- Install clsx + tailwind-merge, create lib/utils.ts with standard cn() helper function
- Remove all btn-*, card-*, card-header, card-body, card-footer, form-group, form-label @utility classes from tailwind.css -- shadcn/ui components replace them
- Update custom utilities (glass, text-gradient) to reference new shadcn/ui CSS variables
- Remove @tailwindcss/forms plugin (shadcn/ui Input/Select/Textarea handle their own styling)
- Keep @tailwindcss/typography plugin (needed for prose content on About, Privacy, Terms pages)
- Install shadcn/ui components: Button, Card, Dialog, Sheet, Input, Select, Textarea, Tabs
- Migration approach: component-by-component across ALL 9 pages (e.g., all buttons first, then all cards, then all inputs)
- Wrap shadcn/ui components inside motion.div wrappers to preserve existing page transitions, stagger animations, and spring physics
- 4 pages use tabbed interfaces: ChemPage, DepictPage, ConvertPage, ToolsPage -- currently use custom tab implementations with motion LayoutGroup for animated indicators
- Full swap from react-icons to lucide-react in one pass across all 45 files
- Replace ALL icon families (HiOutline Heroicons + FontAwesome) with lucide equivalents
- Uninstall react-icons entirely after migration -- one icon library for the whole app
- Standardize icon sizes: 16px (inline/small), 20px (buttons/nav), 24px (headers/features)
- Keep the existing animated pill toggle in Header (spring physics, LayoutGroup animation)
- Swap HiOutlineSun/HiOutlineMoon icons for lucide Sun/Moon icons
- Keep the mobile theme toggle in the mobile menu
- First-visit only: auto-detect system preference via prefers-color-scheme when no localStorage value exists
- Once the user manually toggles, their choice overrides system preference permanently
- Change default from hardcoded dark mode to system preference for new visitors
- Add inline blocking script in index.html that reads localStorage and sets .dark class before React renders -- eliminates flash of wrong theme (FOUC)

### Claude's Discretion
- Color format for CSS variables (HSL vs OKLCH) -- pick best for shadcn/ui compatibility and Tailwind v4
- Tab component approach: shadcn Tabs with motion indicator, shadcn Tabs as-is, or keep current custom tabs -- evaluate per page
- lucide icon stroke weight (default 2px vs lighter 1.5px) -- match current design aesthetic
- Any additional shadcn/ui components beyond the required set if clearly needed during migration

### Deferred Ideas (OUT OF SCOPE)
- Bento grid layouts for pages -- Phase 5 (UX) or Phase 6 (from Phase 3 backlog)
- Liquid glass visual effect -- Phase 6 (from Phase 3 backlog)
- Claymorphism/glassmorphism combination styles -- Phase 6 (from Phase 3 backlog)
- Toast notifications (Sonner) -- Phase 5 (LOAD-02)
- Skeleton loading states -- Phase 5 (LOAD-01)
- Command palette (Cmd+K) -- Phase 6 (POWER-01)
</user_constraints>

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|-----------------|
| COMP-01 | shadcn/ui installed and configured for Vite + Tailwind v4 | Standard Stack section: full init process, components.json, path aliases, cn() utility |
| COMP-02 | Core buttons replaced with shadcn/ui Button component | Architecture Patterns: Button migration strategy, 290 native elements across 38 files |
| COMP-03 | Cards and content containers replaced with shadcn/ui Card | Architecture Patterns: Card migration, only 2 files use card utility class but many use card-like patterns |
| COMP-04 | Modals and dialogs replaced with shadcn/ui Dialog/Sheet | Architecture Patterns: Dialog/Sheet for mobile menu and any modal patterns |
| COMP-05 | Form inputs replaced with shadcn/ui Input/Select/Textarea | Architecture Patterns: Form input migration, removal of @tailwindcss/forms |
| COMP-06 | Icons migrated from react-icons to lucide-react | Standard Stack: lucide-react icon mapping, 46 files with react-icons imports |
| COMP-07 | Consistent component usage across all 9 pages | Architecture Patterns: component-by-component migration ensures consistency |
| THEME-01 | Dark mode toggle available in UI (header/nav) | Architecture Patterns: existing pill toggle preserved, icon swap only |
| THEME-02 | System preference auto-detection (prefers-color-scheme) | Code Examples: FOUC script and AppContext system preference detection |
| THEME-03 | Theme preference persisted in localStorage | Code Examples: existing localStorage persistence enhanced with system preference |
| THEME-04 | CSS variable-based theming using shadcn/ui approach | Standard Stack: OKLCH variables, @theme inline directive, full variable mapping |
| THEME-05 | All pages and components render correctly in both light and dark modes | Architecture Patterns: sky/slate palette mapped to shadcn variables, visual verification |
</phase_requirements>

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| shadcn/ui (CLI) | latest (@latest) | Component scaffolding and code generation | Official CLI generates Tailwind v4 + React 19 compatible components |
| radix-ui | unified package | Headless UI primitives (Dialog, Select, Tabs, etc.) | shadcn/ui Feb 2026 uses unified radix-ui package instead of individual @radix-ui/react-* |
| lucide-react | ^0.577.0 | Icon library (1000+ tree-shakable SVG icons) | Default icon library for shadcn/ui, replaces react-icons |
| class-variance-authority | latest | Component variant management | Required by shadcn/ui Button and other variant-based components |
| clsx | latest | Conditional className utility | Part of standard cn() helper |
| tailwind-merge | latest | Tailwind class conflict resolution | Part of standard cn() helper, prevents duplicate/conflicting utility classes |
| tw-animate-css | latest | CSS-first animation utilities for Tailwind v4 | Replaces tailwindcss-animate plugin, required by shadcn/ui components |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| @tailwindcss/typography | ^0.5.19 | Prose styling for text-heavy pages | Already installed; keep for About, Privacy, Terms pages |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| OKLCH colors | HSL colors | OKLCH is the shadcn/ui default for Tailwind v4, perceptually uniform, better for dark mode. HSL was the v3 approach. Use OKLCH. |
| shadcn Tabs | Keep custom LayoutGroup tabs | Custom tabs have animated indicators with spring physics deeply integrated into 4 pages. Shadcn Tabs would require reimplementing all animation logic. Keep custom tabs. |
| next-themes | Custom ThemeProvider | next-themes is Next.js focused. shadcn/ui provides a Vite-specific ThemeProvider pattern. Use the existing AppContext with enhancements. |
| Default stroke weight (2px) | 1.5px stroke weight | Heroicons use 1.5px stroke; lucide default is 2px. Use default 2px -- it matches the shadcn/ui design language and provides better visibility. |

### Packages to Remove
| Package | Reason |
|---------|--------|
| react-icons | Replaced by lucide-react |
| @fortawesome/react-fontawesome | Replaced by lucide-react |
| @fortawesome/free-solid-svg-icons | Replaced by lucide-react |
| @tailwindcss/forms | shadcn/ui Input/Select/Textarea handle their own styling |

**Installation:**
```bash
cd frontend
npm install lucide-react class-variance-authority clsx tailwind-merge tw-animate-css
npm install shadcn
npx shadcn@latest init
npx shadcn@latest add button card dialog sheet input select textarea tabs
npm uninstall react-icons @fortawesome/react-fontawesome @fortawesome/free-solid-svg-icons @tailwindcss/forms
```

## Architecture Patterns

### Recommended Project Structure
```
frontend/src/
  components/
    ui/               # shadcn/ui generated components (Button, Card, Dialog, etc.)
      button.tsx
      card.tsx
      dialog.tsx
      sheet.tsx
      input.tsx
      select.tsx
      textarea.tsx
      tabs.tsx
    common/           # Existing shared components (Header, Footer, Navigation, etc.)
    chem/             # Chemistry tool views
    convert/          # Conversion tool views
    depict/           # Depiction tool views
    ocsr/             # OCSR views
    tools/            # Tool views
  lib/
    utils.ts          # cn() helper function
  context/
    AppContext.tsx     # Enhanced with system preference detection
  styles/
    tailwind.css      # Rewritten with shadcn/ui CSS variables
```

### Pattern 1: shadcn/ui Initialization for Vite + Tailwind v4

**What:** Configure shadcn/ui for an existing Vite + Tailwind v4 + TypeScript project
**When to use:** First step of Phase 4

**Setup requirements:**

1. Add path alias to `tsconfig.app.json`:
```json
{
  "compilerOptions": {
    "baseUrl": ".",
    "paths": {
      "@/*": ["./src/*"]
    }
  }
}
```

2. Add resolve alias to `vite.config.mts`:
```typescript
import path from "path";
// ... existing imports
export default defineConfig({
  // ... existing config
  resolve: {
    alias: {
      "@": path.resolve(__dirname, "./src"),
    },
  },
});
```

3. Create `components.json` in frontend/:
```json
{
  "$schema": "https://ui.shadcn.com/schema.json",
  "style": "new-york",
  "rsc": false,
  "tsx": true,
  "tailwind": {
    "config": "",
    "css": "src/styles/tailwind.css",
    "baseColor": "slate",
    "cssVariables": true,
    "prefix": ""
  },
  "aliases": {
    "components": "@/components",
    "utils": "@/lib/utils",
    "ui": "@/components/ui",
    "lib": "@/lib",
    "hooks": "@/hooks"
  },
  "iconLibrary": "lucide"
}
```

4. Create `src/lib/utils.ts`:
```typescript
import { clsx, type ClassValue } from "clsx";
import { twMerge } from "tailwind-merge";

export function cn(...inputs: ClassValue[]) {
  return twMerge(clsx(inputs));
}
```

### Pattern 2: CSS Variable Migration (Custom to shadcn/ui Standard)

**What:** Replace 15 custom CSS variables with shadcn/ui's OKLCH-based standard variables
**When to use:** During tailwind.css rewrite

**Current variables to remove:**
- `--bg-primary`, `--bg-secondary` -> `--background`, `--card`
- `--text-primary`, `--text-secondary` -> `--foreground`, `--muted-foreground`
- `--text-accent`, `--text-accent-hover` -> `--primary`
- `--border-primary`, `--border-secondary` -> `--border`
- `--link-primary`, `--link-hover` -> `--primary`
- `--input-bg`, `--input-border`, `--input-text`, `--input-focus-ring` -> `--input`, `--ring`
- `--placeholder-color` -> handled by shadcn Input component
- `--shadow-color-rgb` -> remove (not needed with shadcn approach)

**New tailwind.css structure:**
```css
@import "tailwindcss";
@import "tw-animate-css";

@plugin "@tailwindcss/typography";

@custom-variant dark (&:where(.dark, .dark *));

:root {
  --radius: 0.625rem;
  --background: oklch(0.97 0.002 250);        /* slate-100 equivalent */
  --foreground: oklch(0.15 0.015 260);         /* slate-900 equivalent */
  --card: oklch(1 0 0);                         /* white */
  --card-foreground: oklch(0.15 0.015 260);    /* slate-900 */
  --popover: oklch(1 0 0);
  --popover-foreground: oklch(0.15 0.015 260);
  --primary: oklch(0.55 0.15 235);             /* sky-600 equivalent */
  --primary-foreground: oklch(1 0 0);           /* white */
  --secondary: oklch(0.93 0.005 260);          /* slate-200 equivalent */
  --secondary-foreground: oklch(0.15 0.015 260);
  --muted: oklch(0.93 0.005 260);
  --muted-foreground: oklch(0.45 0.015 260);   /* slate-600 equivalent */
  --accent: oklch(0.93 0.005 260);
  --accent-foreground: oklch(0.15 0.015 260);
  --destructive: oklch(0.577 0.245 27.325);
  --border: oklch(0.82 0.01 260);              /* slate-300 equivalent */
  --input: oklch(0.82 0.01 260);
  --ring: oklch(0.55 0.15 235);                /* sky-500 equivalent */
}

.dark {
  --background: oklch(0.16 0.02 250);          /* custom gray-950 #0f1521 */
  --foreground: oklch(0.93 0.005 260);         /* slate-100 */
  --card: oklch(0.28 0.02 250);                /* slate-800 equivalent */
  --card-foreground: oklch(0.93 0.005 260);
  --popover: oklch(0.28 0.02 250);
  --popover-foreground: oklch(0.93 0.005 260);
  --primary: oklch(0.65 0.17 230);             /* sky-400 equivalent */
  --primary-foreground: oklch(0.15 0.015 260);
  --secondary: oklch(0.35 0.015 255);          /* slate-700 equivalent */
  --secondary-foreground: oklch(1 0 0);
  --muted: oklch(0.35 0.015 255);
  --muted-foreground: oklch(0.7 0.01 255);     /* slate-300 */
  --accent: oklch(0.35 0.015 255);
  --accent-foreground: oklch(1 0 0);
  --destructive: oklch(0.577 0.245 27.325);
  --border: oklch(0.35 0.015 255);             /* slate-700 */
  --input: oklch(0.35 0.015 255);
  --ring: oklch(0.65 0.17 230);                /* sky-400 */
}

@theme inline {
  --color-background: var(--background);
  --color-foreground: var(--foreground);
  --color-card: var(--card);
  --color-card-foreground: var(--card-foreground);
  --color-popover: var(--popover);
  --color-popover-foreground: var(--popover-foreground);
  --color-primary: var(--primary);
  --color-primary-foreground: var(--primary-foreground);
  --color-secondary: var(--secondary);
  --color-secondary-foreground: var(--secondary-foreground);
  --color-muted: var(--muted);
  --color-muted-foreground: var(--muted-foreground);
  --color-accent: var(--accent);
  --color-accent-foreground: var(--accent-foreground);
  --color-destructive: var(--destructive);
  --color-border: var(--border);
  --color-input: var(--input);
  --color-ring: var(--ring);
  --radius: 0.625rem;
}
```

**Critical note:** OKLCH values above are approximate mappings of the existing sky/slate palette. The exact values must be computed from the actual hex/RGB values used by Tailwind's slate and sky color scales. Use a converter tool or compute from Tailwind v4's built-in OKLCH values for slate-100, slate-900, sky-600, etc.

### Pattern 3: Component Wrapping with motion/react

**What:** Wrap shadcn/ui components with motion wrappers to preserve existing animations
**When to use:** When replacing native elements with shadcn/ui components in animated contexts

**Example -- Button in animated context:**
```typescript
import { Button } from "@/components/ui/button";
import { motion } from "motion/react";

// Option A: motion wrapper around Button
<motion.div whileHover={{ scale: 1.02 }} whileTap={{ scale: 0.98 }}>
  <Button variant="default" size="default">
    Submit
  </Button>
</motion.div>

// Option B: motion component with asChild (for Radix slot composition)
// Note: Only works if Button supports asChild prop
```

**Example -- Card in staggered list:**
```typescript
import { Card, CardHeader, CardContent, CardFooter } from "@/components/ui/card";

<motion.div variants={itemVariants}>
  <Card>
    <CardHeader>Title</CardHeader>
    <CardContent>Content</CardContent>
    <CardFooter>Footer</CardFooter>
  </Card>
</motion.div>
```

### Pattern 4: Tab Pages -- Keep Custom Implementation

**What:** Preserve existing LayoutGroup-based tab animations on 4 pages
**When to use:** ChemPage, DepictPage, ConvertPage, ToolsPage

**Recommendation:** Do NOT replace with shadcn/ui Tabs component. Reasons:
1. The existing tab implementation uses `LayoutGroup` with spring-physics animated indicators (stiffness 500, damping 30)
2. Tab pages have complex custom rendering: mobile dropdown, animated tab switching with AnimatePresence, per-tab icons and descriptions
3. ChemPage has a category-based sidebar (not standard tabs) with 14 tools organized into categories
4. The shadcn Tabs component provides basic tab functionality but would require re-implementing all animation logic from scratch
5. The visual quality of the current implementation exceeds what shadcn Tabs offers out of the box

**Action:** Update the tab styling to use shadcn/ui CSS variables for colors but keep the structural implementation. Replace `isDarkMode ? "dark-class" : "light-class"` patterns with the semantic CSS variable approach where possible.

### Pattern 5: FOUC Prevention Script

**What:** Inline blocking script in index.html to apply theme before React renders
**When to use:** Added to index.html `<head>` section

```html
<script>
  (function() {
    var theme = localStorage.getItem('darkMode');
    if (theme === 'true') {
      document.documentElement.classList.add('dark');
    } else if (theme === null) {
      // No saved preference: detect system preference
      if (window.matchMedia('(prefers-color-scheme: dark)').matches) {
        document.documentElement.classList.add('dark');
      }
    }
    // If theme === 'false', no class added (light mode default)
  })();
</script>
```

**Critical notes:**
- Must be in `<head>`, before any CSS or JS loads
- Must be render-blocking (no `async` or `defer`)
- Must use ES5 syntax only (no arrow functions, no const/let) for maximum browser compatibility
- Keep it minimal -- this blocks rendering

### Anti-Patterns to Avoid
- **Mixing old and new component styles:** Never have a page with both hand-rolled `<button className="btn btn-primary">` and `<Button variant="default">`. Complete all buttons across all pages before moving to cards.
- **Adding shadcn ThemeProvider alongside existing AppContext:** Do NOT create a separate ThemeProvider. Enhance the existing AppContext with system preference detection. The app already has a working theme system.
- **Using `forwardRef` in new components:** React 19 + latest shadcn/ui have removed forwardRef. Components accept ref as a regular prop.
- **Importing from @radix-ui/react-*:** Use the unified `radix-ui` package. shadcn CLI (Feb 2026+) generates imports from `radix-ui`.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| CSS class merging | Custom conditional class logic | `cn()` from clsx + tailwind-merge | Handles Tailwind class conflicts (e.g., bg-red-500 vs bg-blue-500) correctly |
| Component variants | Custom prop-to-class mappings | class-variance-authority (cva) | Type-safe variant system used by all shadcn/ui components |
| Accessible dialogs | Custom modal with portal + focus trap | shadcn/ui Dialog (Radix primitive) | Focus trapping, escape handling, aria attributes, scroll locking |
| Accessible select | Custom dropdown with state management | shadcn/ui Select (Radix primitive) | Keyboard navigation, typeahead, ARIA, portal positioning |
| Mobile slide-out menu | Custom AnimatePresence menu | shadcn/ui Sheet | Side overlay with proper focus management, backdrop, escape handling |
| Icon set consistency | Mixing multiple icon libraries | lucide-react exclusively | Tree-shakable, consistent stroke weight, type-safe, shadcn default |

**Key insight:** shadcn/ui components are not a dependency -- they are source code copied into your project. You own and can modify them. The value is in the accessibility primitives (Radix), not the styling (which you control via CSS variables).

## Common Pitfalls

### Pitfall 1: OKLCH Color Value Accuracy
**What goes wrong:** Approximate OKLCH values don't match the intended sky/slate palette, causing subtle but noticeable visual differences
**Why it happens:** Manual OKLCH approximation without conversion tools
**How to avoid:** Use Tailwind v4's built-in color values. Tailwind v4 internally uses OKLCH. Extract the exact OKLCH values from Tailwind's source for slate-100, slate-900, sky-600, sky-400, etc.
**Warning signs:** Colors look "off" compared to pre-migration -- slightly more or less saturated

### Pitfall 2: CSS Variable Scope with @theme inline
**What goes wrong:** Variables defined outside `@theme inline` are not available as Tailwind utility classes (e.g., `bg-background` doesn't work)
**Why it happens:** Tailwind v4 requires explicit registration of CSS variables in `@theme inline` for utility class generation
**How to avoid:** Every CSS variable used in utility classes MUST be registered in `@theme inline` with the `--color-` prefix
**Warning signs:** Build succeeds but utility classes don't apply any styles

### Pitfall 3: Path Alias Not Working in Tests
**What goes wrong:** Vitest can't resolve `@/components/ui/button` imports
**Why it happens:** Vitest uses its own module resolution, separate from Vite's resolve.alias
**How to avoid:** Add the resolve alias to the `test` section of vite.config.mts or ensure Vitest inherits Vite config properly (it should with the current setup since test config is inside vite.config.mts)
**Warning signs:** Tests fail with "Cannot find module @/..."

### Pitfall 4: Removing @tailwindcss/forms Too Early
**What goes wrong:** Form inputs lose all styling before shadcn/ui Input components are in place
**Why it happens:** @tailwindcss/forms provides default styling for native form elements
**How to avoid:** Remove @tailwindcss/forms ONLY AFTER all native form inputs are replaced with shadcn/ui components. This should be one of the last steps.
**Warning signs:** Inputs appear as unstyled browser defaults

### Pitfall 5: LayoutGroup ID Conflicts
**What goes wrong:** Tab indicator animations break or animate to wrong positions after component migration
**Why it happens:** Wrapping existing LayoutGroup-animated elements in new component wrappers changes the layout tree
**How to avoid:** Test tab animations after each change. Keep LayoutGroup usage exactly as-is -- do not wrap tab indicators in additional layers
**Warning signs:** Tab indicator "jumps" instead of smoothly animating

### Pitfall 6: Icon Import Name Mismatches
**What goes wrong:** Build errors because lucide-react uses different naming conventions than react-icons
**Why it happens:** react-icons uses `HiOutlineBeaker` (PascalCase with family prefix), lucide-react uses `Beaker` or `FlaskConical` (PascalCase without prefix)
**How to avoid:** Create a complete icon mapping table before starting migration. Verify each icon exists in lucide-react's catalog. Some icons have different names (e.g., HiOutlineBeaker -> FlaskConical, HiOutlineCollection -> LayoutList)
**Warning signs:** TypeScript errors on import, or visually wrong icons after swap

### Pitfall 7: react-icons/hi vs react-icons/hi2 Confusion
**What goes wrong:** Icon mapping assumes all icons come from the same set
**Why it happens:** HomePage imports from `react-icons/hi2` (Heroicons v2) while all other files import from `react-icons/hi` (Heroicons v1). The icon names differ between versions.
**How to avoid:** Check each file's import source. hi2 icons have different names (e.g., `HiOutlineArrowsRightLeft` in hi2 vs `HiOutlineSwitchHorizontal` in hi)
**Warning signs:** Icon looks different than expected after migration

## Code Examples

### cn() Helper Function
```typescript
// Source: shadcn/ui official docs - Manual Installation
// File: src/lib/utils.ts
import { clsx, type ClassValue } from "clsx";
import { twMerge } from "tailwind-merge";

export function cn(...inputs: ClassValue[]) {
  return twMerge(clsx(inputs));
}
```

### shadcn/ui Button Variant Usage
```typescript
// Source: shadcn/ui docs - Button component
import { Button } from "@/components/ui/button";

// Variants: default, destructive, outline, secondary, ghost, link
// Sizes: default, sm, lg, icon
<Button variant="default">Primary Action</Button>
<Button variant="secondary">Secondary Action</Button>
<Button variant="destructive">Delete</Button>
<Button variant="outline" size="sm">Small Outline</Button>
<Button variant="ghost" size="icon"><Sun className="h-5 w-5" /></Button>
```

### System Preference Detection in AppContext
```typescript
// Enhanced AppContext initialization
useEffect(() => {
  const savedDarkMode = localStorage.getItem("darkMode");
  if (savedDarkMode !== null) {
    // User has a saved preference -- use it
    setIsDarkMode(savedDarkMode === "true");
  } else {
    // No saved preference -- detect system preference
    const prefersDark = window.matchMedia("(prefers-color-scheme: dark)").matches;
    setIsDarkMode(prefersDark);
  }
  // ... rest of initialization
}, []);
```

### Icon Migration Example
```typescript
// BEFORE (react-icons/hi):
import { HiOutlineBeaker, HiOutlineDownload, HiOutlineCheck } from "react-icons/hi";
<HiOutlineBeaker className="h-5 w-5" />

// AFTER (lucide-react):
import { FlaskConical, Download, Check } from "lucide-react";
<FlaskConical className="h-5 w-5" />
// Note: lucide-react default strokeWidth is 2, size defaults to 24px
// Use className="h-5 w-5" (20px) for button/nav context
```

### Icon Mapping Reference (Most Used)

| react-icons | lucide-react | Notes |
|-------------|-------------|-------|
| HiOutlineSun | Sun | Theme toggle |
| HiOutlineMoon | Moon | Theme toggle |
| HiOutlineMenu | Menu | Mobile menu |
| HiOutlineX | X | Close button |
| HiOutlineHome | House | Navigation |
| HiOutlineBeaker | FlaskConical | Chemistry tool |
| HiOutlineRefresh | RefreshCw | Refresh/standardize |
| HiOutlineDownload | Download | Download button |
| HiOutlineCheck | Check | Success/check |
| HiOutlineClipboard | Clipboard | Copy to clipboard |
| HiOutlineCode | Code | Code/HOSE codes |
| HiOutlineSearch | Search | Search/lookup |
| HiOutlineExternalLink | ExternalLink | External links |
| HiOutlineEye | Eye | View/preview |
| HiOutlineCube | Box | 3D/cube icon |
| HiOutlineViewGrid | LayoutGrid | Grid view |
| HiOutlinePencil | Pencil | Edit/draw |
| HiOutlineDocumentDuplicate | Copy | Duplicate/copy |
| HiOutlineDatabase | Database | Database/data |
| HiOutlineFingerPrint | Fingerprint | Identity/unique |
| HiOutlineTag | Tag | Label/classify |
| HiOutlineFilter | Filter | Filter/all filters |
| HiOutlineCollection | LayoutList | Collection/list |
| HiOutlineTemplate | Layout | Template/layout |
| HiOutlineChartSquareBar | BarChart3 | Chart/stats |
| HiOutlineGlobeAlt | Globe | Global/lookup |
| HiOutlineLightningBolt | Zap | Quick/lightning |
| HiOutlineExclamationCircle | AlertCircle | Warning/error |
| HiOutlineInformationCircle | Info | Information |
| HiOutlineCheckCircle | CheckCircle | Success circle |
| HiOutlineDocumentSearch | FileSearch | Document search |
| HiOutlineDocumentText | FileText | Document text |
| HiOutlinePuzzle | Puzzle | Puzzle/structure gen |
| HiOutlineSwitchHorizontal | ArrowLeftRight | Switch/convert |
| HiOutlineAdjustments | SlidersHorizontal | Settings/adjust |
| HiOutlineCamera | Camera | Camera/OCSR |
| HiOutlinePresentationChartLine | LineChart | Presentation/chart |
| HiUserGroup | Users | User group |
| HiChevronDown | ChevronDown | Dropdown indicator |
| HiOutlineArrowRight | ArrowRight | Navigation arrow |
| HiOutlineArrowsExpand | Maximize2 | Expand/fullscreen |
| HiOutlineCubeTransparent | BoxSelect | Transparent cube |
| HiOutlineArrowsRightLeft (hi2) | ArrowLeftRight | Bidirectional arrows |
| HiOutlineWrenchScrewdriver (hi2) | Wrench | Tools/wrench |
| FaGithub | Github | GitHub link |
| FaBook | BookOpen | Documentation |
| FaCode | Code | Source code |
| FaFlask | FlaskConical | Lab/flask |
| FaConnectdevelop | Network | Connected/network |
| FaAtom | Atom | Atom icon |
| FaUniversity | Landmark | University/institution |
| FaCoffee | Coffee | Coffee/support |
| FaShieldAlt | Shield | Security/shield |
| FaFileContract | FileCheck | Terms/contract |
| faAtom (FontAwesome) | Atom | Atom (used in StructureVisualizerView) |
| faCircle (FontAwesome) | Circle | Dot indicator |

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Individual @radix-ui/react-* packages | Unified `radix-ui` package | Feb 2026 | Single dependency, cleaner imports |
| tailwindcss-animate plugin | tw-animate-css (CSS-only) | Feb 2025 (Tailwind v4) | No JS plugin, CSS-first approach |
| HSL color format in CSS vars | OKLCH color format | Feb 2025 (shadcn v4 update) | Perceptually uniform, better dark mode |
| forwardRef in components | Direct ref prop (React 19) | 2025 | Simpler component code |
| `hsl(var(--primary))` in utility | Direct `var(--primary)` with @theme inline | Tailwind v4 | Simpler variable access |
| tailwind.config.js | CSS-first @theme directive | Tailwind v4 | No JS config file needed |

**Deprecated/outdated:**
- `tailwindcss-animate`: Replaced by tw-animate-css for Tailwind v4
- Individual `@radix-ui/react-*` packages: Use unified `radix-ui` package
- `forwardRef` in React 19: Use ref as regular prop
- `style: "default"` in shadcn components.json: Deprecated, use `"new-york"`

## Open Questions

1. **Exact OKLCH values for sky/slate palette**
   - What we know: Tailwind v4 uses OKLCH internally for all its color values. The mapping needs to preserve the exact sky-600/slate-900 look.
   - What's unclear: The precise OKLCH values that correspond to this project's specific palette choices (including custom gray-950 #0f1521).
   - Recommendation: During implementation, use a hex-to-oklch converter for #0f1521 and cross-reference Tailwind v4's source for standard slate/sky OKLCH values. The `npx shadcn@latest init` with baseColor "slate" will generate reasonable defaults.

2. **shadcn CLI behavior with existing tailwind.css**
   - What we know: `npx shadcn@latest init` modifies the CSS file specified in components.json
   - What's unclear: Whether it will merge with or overwrite the existing custom utilities (glass, text-gradient, bg-gradient-radial, transform-style-3d)
   - Recommendation: Back up tailwind.css before running init. The custom @utility blocks should be added back after init. Alternatively, run init with a temporary CSS path and manually merge.

3. **Content Security Policy for inline scripts**
   - What we know: index.html has a CSP with `script-src 'self' 'unsafe-inline'` which allows inline scripts
   - What's unclear: Whether the FOUC prevention script will work correctly with the existing CSP
   - Recommendation: The CSP already allows 'unsafe-inline', so the blocking script will work. No CSP changes needed.

## Validation Architecture

### Test Framework
| Property | Value |
|----------|-------|
| Framework | Vitest 4.1.0 + Testing Library |
| Config file | vite.config.mts (test section) |
| Quick run command | `cd frontend && npm test` |
| Full suite command | `cd frontend && npm run test:coverage` |

### Phase Requirements -> Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| COMP-01 | shadcn/ui installed, cn() utility works | unit | `cd frontend && npx vitest run src/__tests__/utils.test.ts -x` | No -- Wave 0 |
| COMP-02 | Button component renders with variants | unit | `cd frontend && npx vitest run src/__tests__/components/button.test.tsx -x` | No -- Wave 0 |
| COMP-03 | Card component renders correctly | unit | `cd frontend && npx vitest run src/__tests__/components/card.test.tsx -x` | No -- Wave 0 |
| COMP-04 | Dialog/Sheet opens and closes | unit | `cd frontend && npx vitest run src/__tests__/components/dialog.test.tsx -x` | No -- Wave 0 |
| COMP-05 | Input/Select/Textarea render and accept input | unit | `cd frontend && npx vitest run src/__tests__/components/input.test.tsx -x` | No -- Wave 0 |
| COMP-06 | lucide-react icons render (no react-icons imports remain) | smoke | `cd frontend && npx vitest run src/__tests__/icon-migration.test.ts -x` | No -- Wave 0 |
| COMP-07 | All 9 pages render without errors | smoke | `cd frontend && npx vitest run src/__tests__/App.test.tsx -x` | Yes (existing) |
| THEME-01 | Dark mode toggle changes isDarkMode state | unit | `cd frontend && npx vitest run src/__tests__/context/theme.test.tsx -x` | No -- Wave 0 |
| THEME-02 | System preference detection works | unit | `cd frontend && npx vitest run src/__tests__/context/theme.test.tsx -x` | No -- Wave 0 |
| THEME-03 | localStorage persistence works | unit | `cd frontend && npx vitest run src/__tests__/context/theme.test.tsx -x` | No -- Wave 0 |
| THEME-04 | CSS variables change between light/dark | manual-only | Visual inspection -- CSS variables are applied by browser | N/A |
| THEME-05 | Pages render correctly in both modes | manual-only | Visual inspection of all 9 pages in both themes | N/A |

### Sampling Rate
- **Per task commit:** `cd frontend && npm test`
- **Per wave merge:** `cd frontend && npm run test:coverage`
- **Phase gate:** Full suite green + `npm run build` succeeds + `npm run typecheck` succeeds

### Wave 0 Gaps
- [ ] `src/__tests__/utils.test.ts` -- covers COMP-01 (cn() utility)
- [ ] `src/__tests__/components/button.test.tsx` -- covers COMP-02
- [ ] `src/__tests__/components/card.test.tsx` -- covers COMP-03
- [ ] `src/__tests__/components/dialog.test.tsx` -- covers COMP-04
- [ ] `src/__tests__/components/input.test.tsx` -- covers COMP-05
- [ ] `src/__tests__/icon-migration.test.ts` -- covers COMP-06 (grep for react-icons imports)
- [ ] `src/__tests__/context/theme.test.tsx` -- covers THEME-01, THEME-02, THEME-03

## Sources

### Primary (HIGH confidence)
- [shadcn/ui Manual Installation](https://ui.shadcn.com/docs/installation/manual) - Full dependency list, cn() utility, CSS variables, components.json schema
- [shadcn/ui Vite Installation](https://ui.shadcn.com/docs/installation/vite) - Vite-specific setup steps
- [shadcn/ui Tailwind v4](https://ui.shadcn.com/docs/tailwind-v4) - OKLCH migration, @theme inline, tw-animate-css
- [shadcn/ui Theming](https://ui.shadcn.com/docs/theming) - Complete CSS variable list, dark mode convention
- [shadcn/ui Dark Mode Vite](https://ui.shadcn.com/docs/dark-mode/vite) - ThemeProvider, useTheme hook, system preference detection
- [shadcn/ui components.json](https://ui.shadcn.com/docs/components-json) - Full schema, all properties, Tailwind v4 settings
- [shadcn/ui Unified Radix Package (Feb 2026)](https://ui.shadcn.com/docs/changelog/2026-02-radix-ui) - Single radix-ui dependency
- [lucide-react npm](https://www.npmjs.com/package/lucide-react) - v0.577.0, 1000+ icons

### Secondary (MEDIUM confidence)
- [Lucide Icons catalog](https://lucide.dev/icons/) - Icon name verification, used for mapping table
- [tw-animate-css npm](https://www.npmjs.com/package/tw-animate-css) - CSS-first animation replacement
- [FOUC prevention approaches](https://dev.to/gaisdav/how-to-prevent-theme-flash-in-a-react-instant-dark-mode-switching-o20) - Blocking script pattern

### Tertiary (LOW confidence)
- Icon mapping table: Derived from lucide-react catalog browsing and naming pattern matching. Individual icon name mappings should be verified against lucide.dev/icons/ during implementation.

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - All libraries verified via official docs and npm
- Architecture: HIGH - shadcn/ui official docs provide complete Vite + Tailwind v4 setup
- CSS variables: MEDIUM - OKLCH values are approximate; exact mapping needs tool-based conversion
- Icon mapping: MEDIUM - Most mappings verified via lucide.dev; some edge cases (HiOutlineCollection, HiOutlineCubeTransparent) need runtime verification
- FOUC prevention: HIGH - Standard pattern well-documented across multiple sources
- Tab approach: HIGH - Strong recommendation to keep custom tabs based on codebase analysis

**Research date:** 2026-03-12
**Valid until:** 2026-04-12 (shadcn/ui and lucide-react are fast-moving but core APIs stable)
