# Phase 5: Loading States + UX - Research

**Researched:** 2026-03-13
**Domain:** React UX patterns (loading states, toast notifications, error handling, responsive layout, code splitting)
**Confidence:** HIGH

## Summary

Phase 5 addresses the application's feedback layer and performance. Currently, all 21 tool view components use a full-page `LoadingScreen` overlay during API calls, 3 files use `alert()` for error feedback, and all 9 page components are eagerly imported in `App.tsx`, creating a 1.27MB main bundle. The application already has a glassmorphism design system established in Phase 4.1, so the new loading/error/toast components must integrate cohesively with glass styling.

The codebase is well-structured for this work: the Axios instance (`services/api.ts`) already has interceptor stubs ready for global error handling, the glass utility classes provide the foundation for skeleton shimmer styling, and the existing per-component loading patterns (each tool view manages its own `loading`/`error` state) make the refactoring mechanical rather than architectural.

**Primary recommendation:** Add Sonner via shadcn/ui CLI, create reusable `GlassSkeleton` and `GlassErrorCard` components, refactor all 21 tool views from `LoadingScreen` to per-component skeletons, wire Axios response interceptor for network/5xx toasts, implement React.lazy code splitting for all pages, and fix responsive breakpoints.

<user_constraints>

## User Constraints (from CONTEXT.md)

### Locked Decisions
- Glass shimmer skeletons: frosted glass appearance with a translucent highlight sweep animation (left-to-right shimmer)
- Content-shaped skeletons: mirror the actual result layout (molecule image placeholder, text line blocks, property rows)
- Replace the full-page LoadingScreen overlay everywhere -- all loading becomes per-component glass skeletons in the output area
- Submit buttons show inline spinner icon + text change (e.g., "Analyzing..." with spinning loader), disabled during loading
- Input area stays usable/visible during loading -- only the output area shows skeletons
- Sonner library for toast system
- Position: bottom-right corner
- Glass-styled toasts: backdrop-blur + translucent background matching glassmorphism design system
- Auto-dismiss: 3 seconds for success toasts, 5 seconds for error toasts
- Toast scope: toasts for side-effects only (downloads, copy-to-clipboard, network errors)
- Inline feedback for primary content: API results appear in output panel, validation errors shown near input, error + retry button in output area
- Replace all 3 existing alert() calls (MoleculeCard, HighlightedMoleculeCard, NPlikenessView)
- Glass error cards: frosted glass background with red/amber accent left border, error icon, user-friendly message, and Retry button
- User-friendly error messages only -- no status codes or raw error text shown to user. Technical details logged to browser console
- Retry button scrolls to / focuses the input area so user can modify and resubmit
- Two-layer error handling: global Axios interceptor catches network/5xx errors (shows toast), per-component handling for 4xx/business errors (shows inline glass error card)
- Error messages are contextual: "Could not analyze this molecule", "Service temporarily unavailable" -- not generic "Something went wrong"
- Stacked mobile layout: input above output vertically on mobile (<1024px), resizable panels disabled
- Tab strip becomes horizontally scrollable on mobile
- Keep existing Sheet hamburger menu -- polish tap targets to 44px minimum, comfortable spacing
- Remove 3D CaffeineMolecule on mobile (<768px) for GPU performance -- show glass hero card with title + CTA only
- Feature cards stack vertically in single column on mobile
- Glass effects reduced on mobile: backdrop-blur 12px instead of 24px (already established in Phase 4.1)
- Standard Tailwind breakpoints (sm:640, md:768, lg:1024, xl:1280) -- no custom breakpoints

### Claude's Discretion
- Code splitting strategy: React.lazy + Suspense for route-level splitting, Suspense fallback component design
- Navigation wayfinding improvements (UX-03): active nav highlighting, visual indicators for current page
- Lighthouse performance optimization targets and techniques
- Empty state design for tools with no results yet
- Exact shimmer animation CSS implementation (keyframes, gradient technique)
- Which components get content-shaped skeletons vs simpler glass block skeletons

### Deferred Ideas (OUT OF SCOPE)
None -- discussion stayed within phase scope

</user_constraints>

<phase_requirements>

## Phase Requirements

| ID | Description | Research Support |
|----|-------------|-----------------|
| LOAD-01 | Skeleton loading states replace full-page spinners during API calls | GlassSkeleton component with shimmer animation; refactor 21 LoadingScreen usages to per-component skeletons |
| LOAD-02 | Toast notification system (Sonner) replaces alert() calls | Sonner v2.0.7 via shadcn/ui CLI; glass-styled Toaster; replace 3 alert() calls |
| LOAD-03 | Per-component error states with clear messages and retry actions | GlassErrorCard component with contextual messages; retry focuses input area |
| LOAD-04 | Loading indicators for individual API operations (not full-page) | Submit button inline spinner + text change; output-area-only skeleton |
| UX-01 | Responsive layout works well on mobile and tablet viewports | Stacked layout <1024px; scrollable tab strip; 44px tap targets; hide 3D on mobile |
| UX-02 | Code splitting with React.lazy for route-level lazy loading | React.lazy all 9 pages in App.tsx; manualChunks for vendor splitting |
| UX-03 | Improved navigation structure and visual wayfinding | Active nav highlighting already exists via LayoutGroup; enhance breadcrumb or current-page indicator |
| UX-04 | Faster initial page load (measurable improvement via Lighthouse) | Route splitting + vendor chunking should reduce main bundle from 1.27MB to ~400KB |

</phase_requirements>

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| sonner | 2.0.7 | Toast notifications | Locked decision; 12K GitHub stars; shadcn/ui integration; no hook setup needed |
| react (existing) | 19.2.4 | UI framework | Already installed |
| motion (existing) | 12.36.0 | Animations for skeleton shimmer, error card entrance | Already installed |
| axios (existing) | 1.13.5 | HTTP client with interceptor support | Already installed, interceptor stubs ready |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| react-error-boundary | latest | Error boundary for lazy-loaded routes | Wrap Suspense boundaries for error recovery |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| sonner | react-hot-toast | Sonner is locked decision; has better shadcn integration |
| Custom skeleton | react-loading-skeleton | Custom matches glassmorphism design; no extra dependency |
| react-error-boundary | Manual class component | react-error-boundary provides hooks + retry; minimal footprint |

**Installation:**
```bash
cd frontend && npm install sonner react-error-boundary
```

Note: `sonner` can also be added via `npx shadcn@latest add sonner` which creates the `Toaster` wrapper component at `src/components/ui/sonner.tsx`. This is the recommended approach since the project already uses shadcn/ui.

## Architecture Patterns

### New Components Structure
```
src/
  components/
    feedback/
      GlassSkeleton.tsx        # Base shimmer skeleton with glass styling
      GlassErrorCard.tsx        # Error card with retry button
      ToolSkeleton.tsx          # Content-shaped skeleton for tool output areas
      EmptyState.tsx            # Empty state for tools with no results yet
    ui/
      sonner.tsx                # shadcn/ui Toaster wrapper (glass-styled)
      skeleton.tsx              # shadcn/ui base Skeleton (if needed)
```

### Pattern 1: Glass Shimmer Skeleton
**What:** A translucent frosted-glass skeleton with a left-to-right shimmer sweep
**When to use:** Replace LoadingScreen in all 21 tool view components
**Example:**
```typescript
// GlassSkeleton.tsx
import { cn } from "@/lib/utils";

interface GlassSkeletonProps {
  className?: string;
  variant?: "text" | "image" | "card" | "row";
}

export const GlassSkeleton = ({ className, variant = "text" }: GlassSkeletonProps) => (
  <div
    className={cn(
      "relative overflow-hidden rounded-lg",
      "bg-white/10 dark:bg-slate-700/20",
      "backdrop-blur-sm",
      "animate-glass-shimmer",
      {
        "h-4 w-full": variant === "text",
        "h-48 w-full": variant === "image",
        "h-32 w-full": variant === "card",
        "h-6 w-3/4": variant === "row",
      },
      className
    )}
  >
    {/* Shimmer overlay gradient */}
    <div className="absolute inset-0 -translate-x-full animate-shimmer-sweep
      bg-gradient-to-r from-transparent via-white/20 dark:via-white/5 to-transparent" />
  </div>
);
```

### Pattern 2: Two-Layer Error Handling
**What:** Global Axios interceptor for network/5xx errors (toast), per-component for 4xx (inline card)
**When to use:** All API calls through the service layer
**Example:**
```typescript
// services/api.ts - Add to existing response interceptor
import { toast } from "sonner";

api.interceptors.response.use(
  (response) => response,
  (error) => {
    const status = error?.response?.status;

    // Network errors (no response) or 5xx server errors -> toast
    if (!error.response || (status && status >= 500)) {
      toast.error("Service temporarily unavailable", {
        description: "Please try again in a moment.",
        duration: 5000,
      });
    }
    // Rate limiting -> toast
    if (status === 429) {
      toast.error("Too many requests", {
        description: "Please wait before trying again.",
        duration: 5000,
      });
    }

    // 4xx errors are NOT toasted — they propagate to per-component error handling
    // Log technical details to console only
    if (import.meta.env.DEV) {
      console.error("API Error:", status, error?.response?.data);
    }

    return Promise.reject(error);
  }
);
```

### Pattern 3: Route-Level Code Splitting
**What:** React.lazy + Suspense for all page imports in App.tsx
**When to use:** All 9 page-level route components
**Example:**
```typescript
// App.tsx
import React, { lazy, Suspense } from "react";

const HomePage = lazy(() => import("./pages/HomePage"));
const ChemPage = lazy(() => import("./pages/ChemPage"));
const ConvertPage = lazy(() => import("./pages/ConvertPage"));
const DepictPage = lazy(() => import("./pages/DepictPage"));
const ToolsPage = lazy(() => import("./pages/ToolsPage"));
const OCSRPage = lazy(() => import("./pages/OCSRPage"));
const AboutPage = lazy(() => import("./pages/AboutPage"));
const TermsOfService = lazy(() => import("./pages/TermsOfService"));
const PrivacyPolicy = lazy(() => import("./pages/PrivacyPolicy"));

// In Layout component
const Layout = () => (
  <div className="flex flex-col min-h-screen bg-background text-foreground">
    <Header />
    <main className="grow">
      <Suspense fallback={<RouteLoadingFallback />}>
        <Outlet />
      </Suspense>
    </main>
    <Footer />
    <Toaster />
    <ComparisonTray />
    <ComparisonView />
  </div>
);
```

### Pattern 4: Tool View Refactoring
**What:** Replace LoadingScreen + error div with GlassSkeleton + GlassErrorCard
**When to use:** All 21 tool view components follow this same pattern
**Before:**
```tsx
{loading && <LoadingScreen text="Calculating..." />}
{error && !loading && (
  <div className="p-4 rounded-md bg-red-50 ...">
    <AlertCircle ... />
    <span>{error}</span>
  </div>
)}
```
**After:**
```tsx
{/* Submit button with inline spinner */}
<Button disabled={loading}>
  {loading ? (
    <>
      <Loader2 className="mr-2 h-4 w-4 animate-spin" />
      Analyzing...
    </>
  ) : (
    <>
      <Calculator className="mr-2 h-4 w-4" />
      Analyze Molecule
    </>
  )}
</Button>

{/* Output area: skeleton while loading, error card on failure */}
{loading && <ToolSkeleton variant="descriptors" />}
{error && !loading && (
  <GlassErrorCard
    message="Could not analyze this molecule"
    onRetry={() => inputRef.current?.focus()}
  />
)}
```

### Anti-Patterns to Avoid
- **Full-page loading overlay:** Never use `LoadingScreen` or `fixed inset-0` loader for API calls. Always use per-component skeletons in the output area only.
- **Raw error messages to users:** Never show `error.message` or status codes. Map to user-friendly messages; log details to console.
- **Global loading state:** Do NOT use `AppContext.isLoading` for per-tool operations. Each tool view manages its own loading state (they already do this correctly).
- **Toast for primary content errors:** Toasts are for side-effects only. API result errors use inline GlassErrorCard.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Toast notifications | Custom toast with portal + animations | Sonner | Accessibility, stacking, auto-dismiss, mobile swipe — all handled |
| Error boundary | Manual try/catch in render | react-error-boundary | Provides reset/retry, works with Suspense, declarative |
| Shimmer animation | JS-based animation loop | CSS @keyframes + Tailwind @theme | GPU-accelerated, no JS overhead, consistent timing |
| Code splitting | Manual webpack/rollup config | React.lazy + dynamic import() | Vite handles chunking automatically with dynamic imports |
| Vendor chunk splitting | Complex rollup config | Vite manualChunks | Separate react, motion, lucide into cacheable chunks |

**Key insight:** The glass shimmer skeleton is pure CSS -- no animation library needed. The shimmer is a translating gradient pseudo-element with `transform: translateX()` animation, fully GPU-composited.

## Common Pitfalls

### Pitfall 1: Sonner Styling Specificity
**What goes wrong:** Sonner's default styles override custom classes. The project uses glassmorphism styling that needs to override Sonner defaults.
**Why it happens:** Sonner applies inline styles with higher specificity than Tailwind utility classes.
**How to avoid:** Use `unstyled: true` on the Toaster OR use `!important` via Tailwind's `!` prefix in `classNames` prop. Alternatively, use `toastOptions.className` with the glass utility.
**Warning signs:** Toasts don't have backdrop-blur or appear fully opaque when they should be translucent.

### Pitfall 2: Suspense + createBrowserRouter Compatibility
**What goes wrong:** React Router's `createBrowserRouter` with `createRoutesFromElements` expects components directly, not lazy elements.
**Why it happens:** The `element` prop needs the actual React element, and lazy components must be wrapped in Suspense.
**How to avoid:** Wrap the `<Outlet />` in Layout with a single `<Suspense>` boundary, OR use the route-level `lazy` property with React Router's built-in lazy loading.
**Warning signs:** "A component suspended while rendering, but no fallback UI was specified."

### Pitfall 3: Layout Shift from Skeletons
**What goes wrong:** Skeleton placeholders are different sizes than actual content, causing CLS (Cumulative Layout Shift).
**Why it happens:** Skeleton dimensions don't match final content dimensions.
**How to avoid:** Use fixed heights matching the output area, or use `min-h-[X]` on the skeleton container matching the typical result height.
**Warning signs:** Content "jumps" when loading completes.

### Pitfall 4: Axios Interceptor Toast Storms
**What goes wrong:** Multiple failed API calls trigger many toasts simultaneously, flooding the screen.
**Why it happens:** No deduplication of error toasts for concurrent failures.
**How to avoid:** Use Sonner's `id` parameter for network error toasts (e.g., `toast.error("...", { id: "network-error" })`). Same ID replaces instead of stacking.
**Warning signs:** 5+ toasts appearing simultaneously during a network outage.

### Pitfall 5: CaffeineMolecule3D on Mobile
**What goes wrong:** Three.js/WebGL causes GPU perf issues, battery drain, and potentially crashes on mobile devices.
**Why it happens:** The 3D scene with bloom post-processing is heavy for mobile GPUs.
**How to avoid:** Already lazy-loaded, but context says to hide completely on mobile (<768px). Use `window.matchMedia` or a responsive hook, and conditionally render `null` instead of the component.
**Warning signs:** Poor Lighthouse performance score on mobile, janky scroll.

### Pitfall 6: Tailwind v4 animate-pulse May Need Explicit Keyframes
**What goes wrong:** The `animate-pulse` utility may not work out of the box in Tailwind v4 if the keyframes aren't defined.
**Why it happens:** Tailwind v4 does not auto-include all keyframes; some need explicit `@keyframes` definition.
**How to avoid:** Define custom shimmer keyframes in `tailwind.css` under `@theme` for animation, or under `@keyframes` block. Test that `animate-pulse` works. The project already defines custom animations in `@theme` so follow that pattern.
**Warning signs:** Skeleton appears as a static gray block with no animation.

## Code Examples

### CSS: Glass Shimmer Animation
```css
/* Add to tailwind.css */
@theme {
  /* Custom shimmer animation for skeleton loading */
  --animate-shimmer: shimmer 2s ease-in-out infinite;
}

@keyframes shimmer {
  0% {
    transform: translateX(-100%);
  }
  100% {
    transform: translateX(100%);
  }
}
```

### Sonner Toaster with Glass Styling
```typescript
// src/components/ui/sonner.tsx
import { Toaster as SonnerToaster } from "sonner";
import { useAppContext } from "@/context/AppContext";

export function Toaster() {
  const { isDarkMode } = useAppContext();

  return (
    <SonnerToaster
      position="bottom-right"
      theme={isDarkMode ? "dark" : "light"}
      toastOptions={{
        duration: 3000,
        classNames: {
          toast:
            "!backdrop-blur-xl !bg-white/70 dark:!bg-slate-800/70 !border !border-white/20 dark:!border-slate-600/30 !shadow-lg",
          error: "!border-l-4 !border-l-red-500",
          success: "!border-l-4 !border-l-green-500",
          description: "!text-muted-foreground",
        },
      }}
    />
  );
}
```

### GlassErrorCard Component
```typescript
// src/components/feedback/GlassErrorCard.tsx
import { AlertTriangle, RotateCcw } from "lucide-react";
import { Button } from "@/components/ui/button";
import { cn } from "@/lib/utils";

interface GlassErrorCardProps {
  message: string;
  onRetry?: () => void;
  className?: string;
}

export const GlassErrorCard = ({ message, onRetry, className }: GlassErrorCardProps) => (
  <div
    className={cn(
      "glass-bold rounded-xl p-6",
      "border-l-4 border-l-red-500 dark:border-l-red-400",
      className
    )}
    role="alert"
  >
    <div className="flex items-start gap-3">
      <AlertTriangle className="h-5 w-5 text-red-500 dark:text-red-400 shrink-0 mt-0.5" />
      <div className="flex-1">
        <p className="text-sm text-foreground font-medium">{message}</p>
        {onRetry && (
          <Button
            variant="ghost"
            size="sm"
            onClick={onRetry}
            className="mt-3 text-primary hover:text-primary/80"
          >
            <RotateCcw className="mr-2 h-4 w-4" />
            Try again
          </Button>
        )}
      </div>
    </div>
  </div>
);
```

### Vite Manual Chunks Configuration
```typescript
// vite.config.mts build.rollupOptions addition
build: {
  sourcemap: "hidden",
  outDir: "dist",
  rollupOptions: {
    output: {
      manualChunks: {
        "vendor-react": ["react", "react-dom", "react-router-dom"],
        "vendor-motion": ["motion"],
        "vendor-ui": ["lucide-react", "radix-ui", "class-variance-authority", "clsx", "tailwind-merge"],
      },
    },
  },
},
```

### Error Message Mapping
```typescript
// src/lib/error-messages.ts
const ERROR_MESSAGES: Record<string, Record<string, string>> = {
  chem: {
    default: "Could not analyze this molecule",
    parse: "This molecule could not be parsed. Please check the SMILES notation",
    timeout: "Analysis is taking too long. Please try a simpler molecule",
  },
  convert: {
    default: "Could not convert this structure",
    parse: "The input format could not be recognized",
  },
  depict: {
    default: "Could not generate the depiction",
  },
  tools: {
    default: "This tool encountered an error",
  },
  ocsr: {
    default: "Could not recognize the structure from this image",
  },
  network: {
    default: "Service temporarily unavailable",
    timeout: "The request timed out. Please try again",
    rate_limit: "Too many requests. Please wait a moment",
  },
};

export function getErrorMessage(domain: string, error: unknown): string {
  const messages = ERROR_MESSAGES[domain] || ERROR_MESSAGES.network;
  if (error instanceof Error) {
    if (error.message.includes("timeout")) return messages.timeout || messages.default;
    if (error.message.includes("parse") || error.message.includes("invalid"))
      return messages.parse || messages.default;
  }
  return messages.default;
}
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Full-page loading spinners | Per-component skeleton placeholders | Industry standard since 2020 | Better perceived performance, no content flash |
| Browser alert() | Toast notification libraries | Widely adopted 2018+ | Non-blocking, styled, auto-dismiss |
| Eager route imports | React.lazy + Suspense | React 16.6+ (2018), stable concurrent in React 18+ | Reduced initial bundle, faster TTI |
| Single monolithic bundle | Vendor chunk splitting | Vite/Rollup default behavior with dynamic imports | Better caching, parallel loading |

**Deprecated/outdated:**
- `LoadingScreen.tsx`: The full-page overlay component will be replaced entirely. Can be deleted after all 21 usages are migrated.
- `AppContext.isLoading/globalError`: The global loading/error state is unused by tool views (they each manage local state). Consider removing in Phase 5 cleanup.

## Current Bundle Analysis

| Chunk | Size | Gzipped | Notes |
|-------|------|---------|-------|
| index (main) | 1,268 KB | 345 KB | All pages + all deps. Too large. |
| CaffeineMolecule3D | 955 KB | 256 KB | Already lazy-loaded (Three.js) |
| jszip.min | 98 KB | 30 KB | Already separate (dynamic import) |
| CSS | 214 KB | 27 KB | Single CSS bundle |

**After code splitting target:**
| Chunk | Est. Size | Improvement |
|-------|-----------|-------------|
| vendor-react | ~150 KB | Cacheable, rarely changes |
| vendor-motion | ~100 KB | Cacheable |
| vendor-ui | ~80 KB | Cacheable |
| app-core (layout, routing) | ~50 KB | Entry point, fast load |
| Each page chunk | 15-50 KB | Loaded on demand |

## Scope of Changes

### Files to Create (new)
- `src/components/feedback/GlassSkeleton.tsx` -- Base glass shimmer skeleton
- `src/components/feedback/GlassErrorCard.tsx` -- Error card with retry
- `src/components/feedback/ToolSkeleton.tsx` -- Content-shaped skeletons per tool type
- `src/components/feedback/EmptyState.tsx` -- Empty state for tools
- `src/components/feedback/RouteLoadingFallback.tsx` -- Suspense fallback for route loading
- `src/components/ui/sonner.tsx` -- Glass-styled Toaster wrapper
- `src/lib/error-messages.ts` -- Contextual error message mapping

### Files to Modify (21 tool views)
All replace `LoadingScreen` import + usage with skeleton + error card:
- `src/components/chem/AllFiltersView.tsx`
- `src/components/chem/ClassyfireView.tsx`
- `src/components/chem/CoconutPreProcessingView.tsx`
- `src/components/chem/DescriptorsView.tsx`
- `src/components/chem/ErtlFunctionalGroupView.tsx`
- `src/components/chem/FixRadicalsView.tsx`
- `src/components/chem/HOSECodeView.tsx`
- `src/components/chem/NPlikenessView.tsx`
- `src/components/chem/PubChemLookupView.tsx`
- `src/components/chem/StandardizeView.tsx`
- `src/components/chem/StandardizedTautomerView.tsx`
- `src/components/chem/StereoisomersView.tsx`
- `src/components/chem/StructureErrorView.tsx`
- `src/components/chem/TanimotoView.tsx`
- `src/components/convert/FormatConversionView.tsx`
- `src/components/convert/Mol2DView.tsx`
- `src/components/convert/Mol3DView.tsx`
- `src/components/depict/Depict2DMultiView.tsx`
- `src/components/depict/StructureVisualizerView.tsx`
- `src/components/tools/StructureGenView.tsx`
- `src/components/tools/SugarRemovalView.tsx`

### Files to Modify (alert() replacement)
- `src/components/common/MoleculeCard.tsx` -- line 184: `alert("Failed to download...")` -> `toast.error(...)`
- `src/components/common/HighlightedMoleculeCard.tsx` -- line 315: same pattern
- `src/components/chem/NPlikenessView.tsx` -- line 156: `alert(info text)` -> convert to tooltip or info dialog

### Files to Modify (infrastructure)
- `src/App.tsx` -- React.lazy imports, Suspense wrapper, Toaster component
- `src/services/api.ts` -- Axios response interceptor for toast notifications
- `src/styles/tailwind.css` -- Shimmer keyframes, glass-skeleton utility
- `vite.config.mts` -- manualChunks for vendor splitting
- `src/pages/HomePage.tsx` -- Hide CaffeineMolecule3D on mobile (<768px)

### Files to Potentially Delete
- `src/components/common/LoadingScreen.tsx` -- After all 21 usages are migrated

## Responsive Strategy

### Breakpoint Usage (Standard Tailwind)
| Breakpoint | Width | Layout Behavior |
|------------|-------|----------------|
| Default (mobile) | <640px | Single column, stacked I/O, scrollable tabs, no 3D, glass-mobile blur |
| sm | 640px+ | Minor spacing adjustments |
| md | 768px+ | CaffeineMolecule3D appears, sidebar shows on tool pages |
| lg | 1024px+ | Side-by-side input/output panels, full desktop nav |
| xl | 1280px+ | Full labels in nav, wider panels |

### Mobile-Specific Changes
1. **HomePage:** Hide `<CaffeineMolecule3D />` below md (768px). Conditional render via `useMediaQuery` or CSS `hidden md:block`.
2. **Tool Pages (ChemPage, ToolsPage, etc.):** Already have `isMobile` state. Ensure sidebar is fully hidden, content fills width.
3. **Tab Strips:** Add `overflow-x-auto` + `flex-nowrap` + `scrollbar-hide` for horizontal scrolling.
4. **Touch Targets:** Audit all buttons/links for minimum 44x44px. Particularly sidebar nav items and action buttons.

## Open Questions

1. **react-error-boundary version**
   - What we know: It works with React 19 (peer dep allows >=16.13)
   - What's unclear: Whether the latest version has any React 19 specific optimizations
   - Recommendation: Install latest, test at implementation time

2. **Lighthouse baseline measurement**
   - What we know: The current main bundle is 1.27MB (345KB gzipped). After splitting, entry should be ~250-400KB gzipped.
   - What's unclear: Exact Lighthouse scores before/after (need to run against deployed or local build)
   - Recommendation: Run Lighthouse on current build as first task, record baseline, compare after code splitting

3. **NPlikenessView info alert**
   - What we know: Uses `alert()` to show info about NP-likeness score ranges (not an error)
   - What's unclear: Whether this should become a tooltip, popover, or dialog
   - Recommendation: Convert to a glass-styled popover (shadcn/ui Popover) or simply expand the existing "About" info box -- this is NOT a toast scenario

## Validation Architecture

### Test Framework
| Property | Value |
|----------|-------|
| Framework | Vitest 4.1.0 |
| Config file | `frontend/vite.config.mts` (test section) |
| Quick run command | `cd frontend && npx vitest run` |
| Full suite command | `cd frontend && npx vitest run --coverage` |

### Phase Requirements -> Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| LOAD-01 | GlassSkeleton renders with shimmer animation class | unit | `cd frontend && npx vitest run src/__tests__/components/glass-skeleton.test.tsx -t "shimmer"` | Wave 0 |
| LOAD-02 | Toaster component renders, toast.error fires | unit | `cd frontend && npx vitest run src/__tests__/components/sonner.test.tsx` | Wave 0 |
| LOAD-03 | GlassErrorCard renders message and retry button | unit | `cd frontend && npx vitest run src/__tests__/components/glass-error-card.test.tsx` | Wave 0 |
| LOAD-04 | Submit button shows spinner when loading=true | unit | `cd frontend && npx vitest run src/__tests__/components/tool-loading.test.tsx` | Wave 0 |
| UX-01 | Responsive layout renders without horizontal overflow | manual-only | Manual viewport testing at 375px/768px | N/A |
| UX-02 | React.lazy pages load on navigation | unit | `cd frontend && npx vitest run src/__tests__/App.test.tsx -t "lazy"` | Wave 0 |
| UX-03 | Active nav link highlighted | unit | `cd frontend && npx vitest run src/__tests__/components/navigation.test.tsx` | Wave 0 |
| UX-04 | Build produces multiple chunks (not single monolith) | smoke | `cd frontend && npm run build 2>&1 \| grep -c "dist/assets/"` | N/A |

### Sampling Rate
- **Per task commit:** `cd frontend && npx vitest run`
- **Per wave merge:** `cd frontend && npx vitest run --coverage`
- **Phase gate:** Full suite green before `/gsd:verify-work`

### Wave 0 Gaps
- [ ] `src/__tests__/components/glass-skeleton.test.tsx` -- covers LOAD-01
- [ ] `src/__tests__/components/sonner.test.tsx` -- covers LOAD-02
- [ ] `src/__tests__/components/glass-error-card.test.tsx` -- covers LOAD-03
- [ ] `src/__tests__/components/tool-loading.test.tsx` -- covers LOAD-04
- [ ] `src/__tests__/components/navigation.test.tsx` -- covers UX-03

## Sources

### Primary (HIGH confidence)
- Sonner official documentation (sonner.emilkowal.ski) -- API reference, Toaster props, toast types, styling
- shadcn/ui Sonner integration (ui.shadcn.com/docs/components/radix/sonner) -- shadcn CLI setup
- shadcn/ui Skeleton (ui.shadcn.com/docs/components/radix/skeleton) -- Base skeleton pattern
- React official docs (react.dev/reference/react/Suspense) -- Suspense API
- Vite build output analysis -- Direct measurement of current 1.27MB bundle

### Secondary (MEDIUM confidence)
- Axios interceptor patterns -- Verified against axios official docs
- Tailwind v4 custom animation keyframes -- Verified against project's existing `@theme` pattern
- React Router + React.lazy compatibility -- Verified against current `createBrowserRouter` setup

### Tertiary (LOW confidence)
- react-error-boundary React 19 compatibility -- Not explicitly verified, inferred from peer deps
- Exact post-split bundle size estimates -- Estimates based on library sizes, actual results will vary

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - Sonner is locked decision with verified API; all other tools already in project
- Architecture: HIGH - Patterns verified against existing codebase structure; mechanical refactoring
- Pitfalls: HIGH - Based on direct code inspection (found 21 LoadingScreen usages, 3 alert() calls, verified bundle sizes)
- Code splitting: MEDIUM - React.lazy + Vite works, but exact chunk sizes are estimates

**Research date:** 2026-03-13
**Valid until:** 2026-04-13 (stable libraries, low churn)
