# Phase 6: Animations + Power Features - Research

**Researched:** 2026-03-13
**Domain:** React animations (motion/react), command palette (cmdk), inline SMILES preview, breadcrumb navigation
**Confidence:** HIGH

## Summary

Phase 6 adds polish animations (page transitions, button micro-interactions, focus rings, tab content transitions, staggered lists, layout animations) and power user features (Cmd+K command palette, inline SMILES preview, breadcrumb navigation). The project already has motion/react v12.36.0 with AnimatePresence, LayoutGroup, and spring physics used extensively across 24+ files. The primary challenge is page transitions with React Router v6's `createBrowserRouter` + `Outlet` pattern, which requires a custom `AnimatedOutlet` component using `useOutlet` + `cloneElement`. Button micro-interactions are CSS-only (per user decision) and partially exist in the Button component already (`hover:scale-[1.03]` and `active:scale-[0.97]`). The command palette uses shadcn/ui's Command component (built on cmdk) which integrates cleanly with the existing shadcn/ui setup.

**Primary recommendation:** Use the established motion/react patterns for all JS-driven animations, CSS transitions for button micro-interactions (augmenting existing clay-btn styles), shadcn/ui Command component for the palette, and the existing `/latest/depict/2D` endpoint with 500ms debounce for SMILES preview.

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions
- **Page Transitions (ANIM-01):** Cross-fade with subtle vertical rise: opacity 0->1 + translateY(8px->0) on enter, reverse on exit. Overlap timing (AnimatePresence mode="sync"). Duration ~250ms using existing spring physics (stiffness 100, damping 20). Wraps around the Suspense/Outlet in Layout component.
- **Tab Content Transitions (ANIM-04):** Subtle cross-fade (~150ms) when switching tools within tab-based pages (ChemPage, DepictPage, ConvertPage, ToolsPage). Tab indicator already animates via layoutId -- content area adds matching fade transition. Quick and lightweight, no directional sliding.
- **Button Micro-interactions (ANIM-02):** Subtle hover: translateY(-2px) + box-shadow with soft glow matching page gradient palette (opacity ~20%). Subtle press: scale(0.97) + translateY(0). Implemented via CSS transitions (200ms ease hover, 100ms ease press) -- NOT motion.div wrappers, for performance across many buttons. Applied globally to interactive buttons.
- **Staggered List Animations (ANIM-05):** Results lists and data displays animate in with staggered timing. Use motion/react for consistency.
- **Layout Animations (ANIM-06):** Content state changes (expand/collapse, show/hide) animate layout shifts. Use motion/react layout animations.
- **Command Palette (POWER-01):** Cmd+K (Mac) / Ctrl+K (Windows). Four searchable scopes: pages (7), individual tools (~21), example molecules (6), recent molecules (from AppContext). Glass-styled modal with grouped/categorized results. Keyboard navigation (arrow keys + Enter). Fuzzy search across all categories.
- **Inline SMILES Preview (POWER-02):** 2D structure preview card below SMILES input, ~120x120px, rendered via `/latest/depict/2D` with 500ms debounce. Glass shimmer skeleton while loading. Invalid SMILES: hide silently. Integrated into SMILESInput component.
- **Breadcrumb Navigation (POWER-03):** Below glass header, above page content. Two-level: Section > Tool. First level clickable. Only on tool pages (Chem, Convert, Depict, Tools). Desktop only (hidden <768px). Tool name cross-fades (~150ms) when switching tabs.

### Claude's Discretion
- Focus ring animation specifics (ANIM-03)
- Staggered list animation timing and style (ANIM-05)
- Layout animation details for expand/collapse (ANIM-06)
- Command palette library choice (cmdk, custom, or shadcn/ui command component)
- Command palette keyboard navigation details
- Exact glass styling values for command palette and SMILES preview
- Breadcrumb separator character and typography
- Performance optimization (will-change, GPU layers, animation cleanup)

### Deferred Ideas (OUT OF SCOPE)
None -- discussion stayed within phase scope.
</user_constraints>

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|-----------------|
| ANIM-01 | Smooth page transitions when navigating between routes | AnimatedOutlet pattern with useOutlet + cloneElement + AnimatePresence mode="sync"; motion.div with key=location.pathname |
| ANIM-02 | Micro-interactions on buttons (press scale, hover glow/lift) | CSS transition utilities in tailwind.css; augment existing clay-btn with translateY + box-shadow glow |
| ANIM-03 | Focus ring animations on interactive elements | CSS transition on outline/ring properties; extend existing focus-visible patterns in button.tsx |
| ANIM-04 | View transitions when switching tabs within tool pages | AnimatePresence mode="wait" already exists in all 4 tab pages; reduce to simple cross-fade (~150ms) |
| ANIM-05 | Staggered list animations for results and data displays | motion/react staggerChildren variant pattern; already used in ChemPage sidebar |
| ANIM-06 | Layout animations on content state changes | motion/react layout prop + AnimatePresence for expand/collapse; LayoutGroup for coordinated animations |
| POWER-01 | Cmd+K / Ctrl+K command palette for quick tool navigation | shadcn/ui Command component (wraps cmdk); install via CLI; centralized navigation data registry |
| POWER-02 | Inline molecule structure preview when entering SMILES | depictService.get2DDepictionUrl with small dimensions + debounce hook; GlassSkeleton for loading |
| POWER-03 | Breadcrumb navigation showing current location in tool hierarchy | Custom Breadcrumbs component using useLocation + useParams; route-to-section mapping object |
</phase_requirements>

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| motion | ^12.36.0 | Page transitions, staggered lists, layout animations | Already installed and used in 24+ files; AnimatePresence, LayoutGroup, layoutId patterns established |
| cmdk | ^1.1.1 | Command palette search/keyboard navigation engine | Powers shadcn/ui Command; used by Linear, Raycast, Vercel; composable, accessible, fast |
| radix-ui | ^1.4.3 | Dialog primitive for command palette overlay | Already installed; shadcn/ui Dialog already uses it |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| shadcn/ui Command | (generated) | Pre-styled command palette component | Install via `npx shadcn@latest add command`; wraps cmdk with project styling |
| axios (via api.ts) | existing | SMILES preview API calls | Already configured with base URL, interceptors, error handling |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| shadcn/ui Command (cmdk) | Custom command palette | cmdk handles fuzzy search, keyboard nav, a11y out of the box; custom would take 3-5x longer |
| CSS transitions for buttons | motion.div wrappers | CSS transitions are more performant for dozens of buttons; motion.div adds JS overhead per instance |
| Backend /depict/2D for preview | Client-side RDKit.js | Backend already exists and works; RDKit.js would add ~15MB to bundle |

**Installation:**
```bash
cd frontend
npx shadcn@latest add command
npm install cmdk
```

Note: The shadcn CLI will install cmdk as a dependency and generate the Command component files. If the CLI creates files in an incorrect directory (known issue from Phase 4), manually move to `src/components/ui/`.

## Architecture Patterns

### Recommended New Component Structure
```
src/
  components/
    common/
      Breadcrumbs.tsx          # NEW: breadcrumb navigation
      CommandPalette.tsx        # NEW: Cmd+K command palette
      SMILESPreview.tsx         # NEW: inline 2D structure preview
    ui/
      command.tsx               # NEW: shadcn/ui generated command component
  hooks/
    useDebounce.ts              # NEW: debounce hook for SMILES preview
  data/
    navigationRegistry.ts       # NEW: centralized tool/page data for palette + breadcrumbs
```

### Pattern 1: AnimatedOutlet for Page Transitions (ANIM-01)
**What:** Custom component that captures the outlet element via `useOutlet()` and clones it with `key={location.pathname}` so AnimatePresence can track route changes.
**When to use:** In the Layout component, replacing the bare `<Outlet />`.
**Why needed:** React Router v6's `<Outlet />` renders an OutletContext.Provider as its direct child, which AnimatePresence cannot track. The useOutlet + cloneElement pattern gives AnimatePresence a keyed direct child to animate.

```typescript
// AnimatedOutlet component
import { useOutlet, useLocation } from "react-router-dom";
import { AnimatePresence, motion } from "motion/react";
import React from "react";

const pageTransition = {
  type: "spring" as const,
  stiffness: 100,
  damping: 20,
};

export function AnimatedOutlet() {
  const location = useLocation();
  const element = useOutlet();

  return (
    <AnimatePresence mode="sync">
      {element && (
        <motion.div
          key={location.pathname}
          initial={{ opacity: 0, y: 8 }}
          animate={{ opacity: 1, y: 0 }}
          exit={{ opacity: 0, y: -8 }}
          transition={pageTransition}
        >
          {React.cloneElement(element, { key: location.pathname })}
        </motion.div>
      )}
    </AnimatePresence>
  );
}
```

**Critical integration point:** The Layout component in App.tsx currently has:
```tsx
<Suspense fallback={<RouteLoadingFallback />}>
  <Outlet />
</Suspense>
```
This must become:
```tsx
<Suspense fallback={<RouteLoadingFallback />}>
  <AnimatedOutlet />
</Suspense>
```

### Pattern 2: CSS Button Micro-interactions (ANIM-02)
**What:** CSS utility classes for hover lift + glow and press scale, applied globally via button component.
**When to use:** For all interactive buttons (except ghost/link variants which are intentionally flat).
**Why CSS not motion.div:** User explicitly decided CSS for performance -- dozens of buttons rendered simultaneously.

```css
/* In tailwind.css -- augment existing clay-btn */
@utility btn-hover-lift {
  &:hover {
    transform: translateY(-2px);
    box-shadow:
      0 4px 12px -2px var(--btn-glow-color, rgb(56 189 248 / 0.2)),
      /* existing clay-btn shadows preserved */
      inset 0 -3px 6px 0 rgb(0 0 0 / 0.08),
      inset 0 2px 3px 0 rgb(255 255 255 / 0.4),
      0 6px 16px -4px rgb(0 0 0 / 0.12);
    transition: transform 200ms ease, box-shadow 200ms ease;
  }
  &:active {
    transform: translateY(0) scale(0.97);
    transition: transform 100ms ease, box-shadow 100ms ease;
  }
}
```

**Existing state:** Button component already has `hover:scale-[1.03]` and `active:scale-[0.97]` in the base cva class. The task is to REPLACE these with the specified translateY(-2px) + glow hover and scale(0.97) press. The ghost and link variants already override with `hover:scale-100 active:scale-100`.

### Pattern 3: Navigation Data Registry (for POWER-01 + POWER-03)
**What:** Centralized data source mapping all pages and tools with their routes, names, icons, and parent sections.
**When to use:** Command palette entries + breadcrumb section/tool mapping.

```typescript
// src/data/navigationRegistry.ts
export interface NavEntry {
  id: string;
  name: string;
  path: string;
  section?: string;       // parent section name for breadcrumbs
  sectionPath?: string;   // parent section route
  icon?: React.ComponentType;
  keywords?: string[];    // extra search terms for fuzzy match
}

export const pages: NavEntry[] = [
  { id: "home", name: "Home", path: "/home", icon: House },
  { id: "chem", name: "Chemical Analysis", path: "/chem", icon: FlaskConical },
  // ... 7 total
];

export const tools: NavEntry[] = [
  { id: "stereoisomers", name: "Stereoisomers", path: "/chem/stereoisomers",
    section: "Chemical Analysis", sectionPath: "/chem", icon: Copy },
  { id: "hosecodes", name: "HOSE Codes", path: "/chem/hosecodes",
    section: "Chemical Analysis", sectionPath: "/chem", icon: Code },
  // ... ~25 total across all 4 tab pages
];

export const exampleMolecules = [
  { name: "Caffeine", smiles: "CN1C=NC2=C1C(=O)N(C(=O)N2C)C" },
  // ... 6 total (from EXAMPLE_MOLECULES in SMILESInput.tsx)
];
```

### Pattern 4: Debounced SMILES Preview (POWER-02)
**What:** Custom hook + preview component that calls the existing depict API with debounce.
**When to use:** Inside SMILESInput component, below the input field.

```typescript
// src/hooks/useDebounce.ts
export function useDebounce<T>(value: T, delay: number): T {
  const [debouncedValue, setDebouncedValue] = useState(value);
  useEffect(() => {
    const timer = setTimeout(() => setDebouncedValue(value), delay);
    return () => clearTimeout(timer);
  }, [value, delay]);
  return debouncedValue;
}

// SMILESPreview usage in SMILESInput.tsx
const debouncedSmiles = useDebounce(value, 500);
// Use get2DDepictionUrl(debouncedSmiles, { width: 240, height: 240 })
// to construct an <img> src URL -- no axios call needed, just URL construction
```

**Key insight:** The `get2DDepictionUrl()` function in depictService.ts constructs a URL that can be used directly as an `<img src>`. This avoids needing an axios call + state management. The browser handles the image fetch. On error (invalid SMILES), the img onError handler hides the preview silently.

### Pattern 5: Command Palette with shadcn/ui Command (POWER-01)
**What:** Modal command palette triggered by Cmd+K, using shadcn/ui Command component.
**When to use:** Global component rendered in Layout, always available.

```typescript
// CommandPalette.tsx -- key structure
<CommandDialog open={open} onOpenChange={setOpen}>
  <CommandInput placeholder="Search pages, tools, molecules..." />
  <CommandList>
    <CommandEmpty>No results found.</CommandEmpty>
    <CommandGroup heading="Pages">
      {filteredPages.map(page => (
        <CommandItem key={page.id} onSelect={() => navigate(page.path)}>
          <page.icon className="mr-2 h-4 w-4" />
          {page.name}
        </CommandItem>
      ))}
    </CommandGroup>
    <CommandSeparator />
    <CommandGroup heading="Tools">
      {/* ~25 tool entries */}
    </CommandGroup>
    <CommandSeparator />
    <CommandGroup heading="Molecules">
      {/* 6 examples + recent from AppContext */}
    </CommandGroup>
  </CommandList>
</CommandDialog>
```

**Glass styling:** Override CommandDialog's content with glass-bold + backdrop-blur classes to match the design system. The shadcn Dialog already has animation support (fade-in/zoom-in).

### Anti-Patterns to Avoid
- **Wrapping every button in motion.div:** The user explicitly decided CSS transitions for buttons. Do NOT use motion.div wrappers for hover/press effects.
- **Using Outlet directly with AnimatePresence:** AnimatePresence cannot track Outlet's internal Provider. MUST use useOutlet + cloneElement pattern.
- **Heavy page transition animations:** User specified subtle 250ms cross-fade with 8px rise, not dramatic slides or zooms. Keep it restrained.
- **Fetching SMILES preview via axios:** Use get2DDepictionUrl() to construct a URL and use `<img src>` -- simpler, cacheable by browser.
- **Separate route data in each component:** Centralize in navigationRegistry.ts to avoid duplication between command palette, breadcrumbs, and page components.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Command palette fuzzy search | Custom fuzzy search algorithm | cmdk (via shadcn/ui Command) | cmdk handles search scoring, keyboard navigation, a11y, grouping |
| Keyboard shortcut detection | Manual keydown listeners | useEffect + single global listener | Single listener with meta/ctrl detection is sufficient; no library needed |
| Debounce | Custom debounce function | Simple useDebounce hook (5 lines) | Too small for a library; too common to inline everywhere |
| Route-aware breadcrumbs | react-router breadcrumb library | Custom component with useLocation + registry | Registry approach is simpler than configuring a library for 4 routes |

**Key insight:** cmdk/shadcn Command handles the hard parts (fuzzy search, keyboard navigation, screen reader announcements, focus management). The rest are simple enough to hand-write.

## Common Pitfalls

### Pitfall 1: AnimatePresence Not Detecting Route Changes
**What goes wrong:** Wrapping `<Outlet />` in AnimatePresence and motion.div with `key={location.pathname}` -- the animation never fires because Outlet re-renders its Provider child rather than unmounting/remounting it.
**Why it happens:** React Router v6's Outlet renders an OutletContext.Provider, not the page component directly. AnimatePresence only tracks direct children keys.
**How to avoid:** Use `useOutlet()` hook to get the element, then `React.cloneElement(element, { key: location.pathname })` as a direct child of AnimatePresence.
**Warning signs:** Pages swap instantly without any fade, or exit animation never plays.

### Pitfall 2: Suspense Boundary Interfering with Page Transitions
**What goes wrong:** The Suspense boundary around the AnimatedOutlet catches the lazy-loaded page's loading state, causing a flash of the fallback during transitions.
**Why it happens:** When navigating to a new lazy-loaded page, React suspends and shows the fallback before the page loads, breaking the smooth exit animation.
**How to avoid:** Keep the Suspense boundary OUTSIDE the AnimatePresence. The RouteLoadingFallback shows for the initial load; subsequent navigations benefit from prefetching. Alternatively, ensure the Suspense fallback is also wrapped in a motion.div with the same animation.
**Warning signs:** Brief flash of loading skeleton between page transitions.

### Pitfall 3: Button Glow Color Hardcoded Instead of Per-Page
**What goes wrong:** All buttons glow the same color regardless of which page they're on.
**Why it happens:** Hover glow should match per-page gradient palette from Phase 4.1.
**How to avoid:** Use a CSS custom property `--btn-glow-color` that each page's gradient context sets. Or use a simpler approach: a single neutral glow color (sky/primary) that works everywhere.
**Warning signs:** Glow color clashes with page gradient palette.

### Pitfall 4: SMILES Preview Flickering on Every Keystroke
**What goes wrong:** Preview image flashes/flickers as user types because each keystroke triggers a new API URL.
**Why it happens:** Without debounce, every character change creates a new image URL, causing the browser to start fetching a new image.
**How to avoid:** Use 500ms debounce (per user decision). Only update the img src after debounce settles. Show GlassSkeleton while loading, hide on error.
**Warning signs:** Image rapidly appears/disappears/reloads while typing.

### Pitfall 5: Command Palette Steals Focus from Inputs
**What goes wrong:** Pressing Cmd+K while typing in a SMILES input opens the palette but the input loses its text or cursor position.
**Why it happens:** Focus moves to the command palette search input.
**How to avoid:** This is expected behavior -- same as VS Code/Linear. The palette closes on Escape, returning focus to the previous element. Ensure proper focus restoration.
**Warning signs:** Users complain about losing work when accidentally hitting Cmd+K.

### Pitfall 6: shadcn CLI Creates Files in Wrong Directory
**What goes wrong:** The shadcn CLI generates command.tsx in a literal `@/` directory instead of `src/components/ui/`.
**Why it happens:** Known issue from Phase 4 -- the CLI doesn't always resolve the alias correctly.
**How to avoid:** After running `npx shadcn@latest add command`, verify the file landed in `src/components/ui/command.tsx`. If not, move it manually.
**Warning signs:** Import errors after running the CLI.

### Pitfall 7: AnimatePresence mode="sync" vs mode="wait"
**What goes wrong:** Using mode="wait" for page transitions causes a visible gap (old page fully exits, delay, new page enters). Using mode="sync" overlaps them but can cause layout issues with two pages rendered simultaneously.
**Why it happens:** User specified mode="sync" for overlap timing (cross-fade). With mode="sync", both pages exist in the DOM simultaneously during the transition.
**How to avoid:** Use `position: absolute` or similar technique on the exiting page to prevent it from affecting layout during overlap. Or use `mode="popLayout"` which removes the exiting element from layout flow.
**Warning signs:** Two pages stacked vertically during transition instead of overlapping.

## Code Examples

### Existing Animation Patterns (for consistency)

**Spring physics already established (Phase 2):**
```typescript
// Tab indicators: fast, snappy
{ type: "spring", stiffness: 500, damping: 30 }

// Page entrances: smooth, gentle
{ type: "spring", stiffness: 100, damping: 20 }

// Nav pill: medium
{ type: "spring", stiffness: 350, damping: 30 }
```

**Staggered container (already in ChemPage):**
```typescript
const sidebarStaggerContainer = {
  hidden: { opacity: 0 },
  visible: {
    opacity: 1,
    transition: { staggerChildren: 0.05, delayChildren: 0.1 },
  },
};
```

**AnimatePresence with mode="wait" (already in 4 tab pages):**
```typescript
// ChemPage, ConvertPage, DepictPage, ToolsPage all have:
<AnimatePresence mode="wait">
  <motion.div
    key={activeTabId}
    initial="hidden"
    animate="visible"
    exit="hidden"  // or "exit"
    variants={contentVariants}
  >
    <ActiveComponent />
  </motion.div>
</AnimatePresence>
```

### SMILES Preview URL Construction (no fetch needed)
```typescript
// Source: frontend/src/services/depictService.ts
import { get2DDepictionUrl } from "@/services/depictService";

// Construct URL for small preview image
const previewUrl = get2DDepictionUrl(debouncedSmiles, {
  width: 240,
  height: 240,
  toolkit: "rdkit",
});
// Use as: <img src={previewUrl} onError={hidePreview} />
```

### Global Keyboard Shortcut for Command Palette
```typescript
useEffect(() => {
  const handleKeyDown = (e: KeyboardEvent) => {
    if ((e.metaKey || e.ctrlKey) && e.key === "k") {
      e.preventDefault();
      setOpen(prev => !prev);
    }
  };
  document.addEventListener("keydown", handleKeyDown);
  return () => document.removeEventListener("keydown", handleKeyDown);
}, []);
```

### Breadcrumb Route Mapping
```typescript
// Derive breadcrumb from current location
const location = useLocation();
const { toolId } = useParams();
const segments = location.pathname.split("/").filter(Boolean);

// Map route prefix to section name
const sectionMap: Record<string, { name: string; path: string }> = {
  chem: { name: "Chemical Analysis", path: "/chem" },
  convert: { name: "Format Conversion", path: "/convert" },
  depict: { name: "Depiction", path: "/depict" },
  tools: { name: "Tools", path: "/tools" },
};
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| framer-motion package | motion/react package | 2024 (v11+) | Already migrated in Phase 2 |
| AnimatePresence + TransitionGroup | AnimatePresence alone | Stable | Simpler API, built-in exit animations |
| Custom fuzzy search | cmdk library | 2022+ | Battle-tested, accessible, composable |
| Heavy animation libs (GSAP) | motion/react + CSS transitions | Ongoing | motion/react handles 95% of needs; CSS for micro-interactions |

**Deprecated/outdated:**
- `framer-motion` import path: Project already uses `motion/react` (migrated in Phase 2)
- `react-cmdk` library: Last updated 3 years ago; use `cmdk` by Paco Coursey instead
- React Router v5 `<TransitionGroup>` patterns: v6 requires useOutlet approach

## Open Questions

1. **AnimatePresence mode="sync" Layout Overlap**
   - What we know: User wants mode="sync" for overlap cross-fade. Two pages will exist simultaneously.
   - What's unclear: Whether absolute positioning during exit will cause layout shifts on mobile.
   - Recommendation: Test with `mode="popLayout"` first (removes exiting element from flow). Fall back to `mode="sync"` with absolute positioning if needed.

2. **Button Glow Color Strategy**
   - What we know: User wants glow matching page gradient palette at ~20% opacity.
   - What's unclear: Whether per-page CSS custom property (`--btn-glow-color`) is needed or if a single primary/sky glow works universally.
   - Recommendation: Start with primary color glow (`var(--primary)` at 20% opacity). Simpler, consistent, and matches all page palettes. Per-page override can be added later if needed.

3. **Command Palette: Molecule Selection Navigation**
   - What we know: Selecting a molecule navigates to /chem with SMILES pre-filled.
   - What's unclear: How to pre-fill SMILES in the target page's input (no global SMILES state currently).
   - Recommendation: Use URL search params (`/chem?smiles=...`) or a short-lived AppContext state. URL params are more robust and shareable.

## Validation Architecture

### Test Framework
| Property | Value |
|----------|-------|
| Framework | Vitest 4.1.0 + @testing-library/react 16.3.2 |
| Config file | `frontend/vite.config.mts` (test section) |
| Quick run command | `cd frontend && npx vitest run --reporter=verbose` |
| Full suite command | `cd frontend && npx vitest run` |

### Phase Requirements -> Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| ANIM-01 | Page transitions animate on route change | integration | `cd frontend && npx vitest run src/__tests__/components/page-transitions.test.tsx -x` | No -- Wave 0 |
| ANIM-02 | Buttons have hover lift + press scale CSS | unit | `cd frontend && npx vitest run src/__tests__/components/button-animations.test.tsx -x` | No -- Wave 0 |
| ANIM-03 | Focus rings animate on interactive elements | unit | `cd frontend && npx vitest run src/__tests__/components/focus-ring.test.tsx -x` | No -- Wave 0 |
| ANIM-04 | Tab content cross-fades on switch | manual-only | N/A (animation timing not testable in jsdom) | N/A |
| ANIM-05 | Results lists stagger in | manual-only | N/A (stagger timing not testable in jsdom) | N/A |
| ANIM-06 | Layout changes animate | manual-only | N/A (layout animation not testable in jsdom) | N/A |
| POWER-01 | Cmd+K opens command palette; search + navigate | integration | `cd frontend && npx vitest run src/__tests__/components/command-palette.test.tsx -x` | No -- Wave 0 |
| POWER-02 | SMILES input shows inline 2D preview | integration | `cd frontend && npx vitest run src/__tests__/components/smiles-preview.test.tsx -x` | No -- Wave 0 |
| POWER-03 | Breadcrumbs show section > tool on tool pages | unit | `cd frontend && npx vitest run src/__tests__/components/breadcrumbs.test.tsx -x` | No -- Wave 0 |

### Sampling Rate
- **Per task commit:** `cd frontend && npx vitest run --reporter=verbose`
- **Per wave merge:** `cd frontend && npx vitest run`
- **Phase gate:** Full suite green before `/gsd:verify-work`

### Wave 0 Gaps
- [ ] `src/__tests__/components/command-palette.test.tsx` -- covers POWER-01 (open/close, search, navigation)
- [ ] `src/__tests__/components/smiles-preview.test.tsx` -- covers POWER-02 (debounce, show/hide, error handling)
- [ ] `src/__tests__/components/breadcrumbs.test.tsx` -- covers POWER-03 (section/tool display, visibility rules)
- [ ] `src/__tests__/components/button-animations.test.tsx` -- covers ANIM-02 (CSS class presence)
- [ ] `src/__tests__/components/page-transitions.test.tsx` -- covers ANIM-01 (AnimatedOutlet renders with key)

Note: ANIM-04, ANIM-05, ANIM-06 are animation-timing tests that cannot be meaningfully verified in jsdom. These should be verified visually during manual testing.

## Existing Code Inventory

### Files That Need Modification
| File | Change | Requirement |
|------|--------|-------------|
| `src/App.tsx` Layout component | Replace `<Outlet />` with `<AnimatedOutlet />`; add `<CommandPalette />` | ANIM-01, POWER-01 |
| `src/styles/tailwind.css` | Add button hover/press utility; focus ring animation | ANIM-02, ANIM-03 |
| `src/components/ui/button.tsx` | Replace scale hover/press with translateY + glow pattern | ANIM-02 |
| `src/components/common/SMILESInput.tsx` | Add SMILESPreview below input field | POWER-02 |
| `src/pages/ChemPage.tsx` | Simplify tab content animation to cross-fade; add stagger to results | ANIM-04, ANIM-05 |
| `src/pages/ConvertPage.tsx` | Simplify tab content animation to cross-fade | ANIM-04 |
| `src/pages/DepictPage.tsx` | Simplify tab content animation to cross-fade | ANIM-04 |
| `src/pages/ToolsPage.tsx` | Simplify tab content animation to cross-fade | ANIM-04 |

### New Files
| File | Purpose | Requirement |
|------|---------|-------------|
| `src/components/common/AnimatedOutlet.tsx` | Page transition wrapper | ANIM-01 |
| `src/components/common/CommandPalette.tsx` | Command palette modal | POWER-01 |
| `src/components/common/SMILESPreview.tsx` | Inline 2D structure preview | POWER-02 |
| `src/components/common/Breadcrumbs.tsx` | Breadcrumb navigation | POWER-03 |
| `src/components/ui/command.tsx` | shadcn/ui Command component | POWER-01 |
| `src/hooks/useDebounce.ts` | Debounce hook | POWER-02 |
| `src/data/navigationRegistry.ts` | Centralized page/tool/molecule data | POWER-01, POWER-03 |

### Complete Tool Inventory (for Command Palette)
**Pages (7):** Home, Chemical Analysis, Format Conversion, Depiction, Tools, OCSR, About

**Chem tools (14):** Stereoisomers, HOSE Codes, Standardize, Functional Groups, Tautomers, Fix Radicals, Descriptors, NP-likeness, Similarity, Structure Finder, Check Structure, All Filters, COCONUT Preprocessing, ClassyFire

**Convert tools (3):** Format Conversion, 2D Coordinates, 3D Coordinates

**Depict tools (4):** Structure Explorer, 2D Depiction, 3D Depiction, Draw a Structure

**Tools tools (4):** Sugar Detection, Structure Generation, InChI Converter, RInChI Converter

**Example molecules (6):** Caffeine, Aspirin, Sucrose, Cholesterol, Paracetamol, Ibuprofen

**Recent molecules:** From AppContext.recentMolecules (up to 10, filtered to last 24h)

**Total command palette entries:** ~38 static + up to 10 dynamic recent

## Sources

### Primary (HIGH confidence)
- Project codebase analysis: App.tsx, Navigation.tsx, ChemPage.tsx, ConvertPage.tsx, DepictPage.tsx, ToolsPage.tsx, SMILESInput.tsx, AppContext.tsx, button.tsx, tailwind.css, depictService.ts, api.ts
- [shadcn/ui Command docs](https://ui.shadcn.com/docs/components/radix/command) -- component API, installation
- [cmdk npm](https://www.npmjs.com/package/cmdk) -- v1.1.1, used by Linear/Raycast/Vercel

### Secondary (MEDIUM confidence)
- [React Router + AnimatePresence pattern](https://github.com/remix-run/react-router/discussions/8604) -- useOutlet + cloneElement approach
- [AnimatePresence page transitions guide](https://dev.to/joserfelix/page-transitions-in-react-1c8g) -- mode="wait" vs layout patterns
- [Framer Motion AnimatePresence + Outlet](https://medium.com/@antonio.falcescu/animating-react-pages-with-react-router-dom-outlet-and-framer-motion-animatepresence-bd5438b3433b) -- detailed implementation pattern

### Tertiary (LOW confidence)
- None -- all findings verified against project code or official sources.

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- motion/react and shadcn/ui already in project; cmdk is the established standard
- Architecture: HIGH -- patterns verified against existing codebase structure and React Router v6 docs
- Pitfalls: HIGH -- AnimatePresence/Outlet pitfall is well-documented; other pitfalls from Phase 4 experience
- Code examples: HIGH -- based on actual project code and verified library APIs

**Research date:** 2026-03-13
**Valid until:** 2026-04-13 (stable libraries, no upcoming breaking changes)
