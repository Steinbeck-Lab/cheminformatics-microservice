# Feature Landscape

**Domain:** Cheminformatics web tool frontend modernization
**Researched:** 2026-03-12
**Overall confidence:** MEDIUM-HIGH (verified React 19, Tailwind v4, upgrade paths via official docs; component patterns based on training data + codebase analysis)

## Table Stakes

Features users expect from a modern, polished science web tool in 2025-2026. Missing any of these makes the product feel dated or unfinished.

| Feature | Why Expected | Complexity | Notes |
|---------|--------------|------------|-------|
| **All 9 pages work identically** | Users don't care about internal tooling changes | Low | Zero functional regressions during migration. Verify every page after each migration phase. |
| **Dark/light mode with persistence** | Every modern dev/science tool has it. Already partially implemented in AppContext with localStorage. | Low | Enhance: add `prefers-color-scheme` system detection, CSS variable-based theming via shadcn/ui. Existing `isDarkMode` + `toggleDarkMode` + `document.documentElement.classList` logic is the right pattern. |
| **Fast dev server startup** | Developers expect sub-second HMR after moving off CRA | Low | Vite provides this out of the box. CRA's webpack takes 10-30s cold start; Vite takes <1s. |
| **Consistent component system** | Current UI has hand-rolled Tailwind components with varying styles. Modern tools use a design system. | Medium | Replace ad-hoc components with shadcn/ui: Button, Card, Dialog, Input, Select, Tabs, Badge, Tooltip. |
| **Accessible focus management & keyboard navigation** | WCAG 2.1 AA expected at universities. Current components lack focus trapping, arrow-key nav. | Medium | shadcn/ui (Radix UI) handles this by default. Focus rings, keyboard nav, ARIA attributes included. |
| **Loading states with skeleton screens** | Current LoadingScreen is a full-page overlay. Modern pattern is inline skeleton placeholders. | Low | shadcn/ui Skeleton component for cards, tables, form results. Keep molecule animation for initial load only. |
| **Toast notifications for async feedback** | Copy-to-clipboard, download complete, API errors need non-blocking feedback. | Low | shadcn/ui Sonner (toast). Replace `alert()` calls and inline "Copied!" states. Auto-dismiss, stackable. |
| **Responsive layout that works on tablets** | Scientists use iPads in labs. Current responsive exists but needs verification after Tailwind v4 class renames. | Medium | Test at 768px (iPad) and 1024px (laptop). Touch-friendly tap targets (min 44px). |
| **Smooth page transitions** | framer-motion already installed. Page transitions between routes are currently abrupt. | Low | `AnimatePresence` with `mode="wait"` around route outlet. Fast transitions (200-300ms). |
| **Proper error states per component** | Each API call result area should have clear error states with retry. | Low | shadcn/ui Alert component for inline errors. "Something went wrong" + retry button. |
| **Monospace display for chemical identifiers** | SMILES, InChI, InChIKey, HOSE codes are code-like strings needing monospace + horizontal scroll. | Low | Already partially done. Standardize across all components. |

## Differentiators

Features that set this product apart. Not expected, but make it feel premium.

| Feature | Value Proposition | Complexity | Notes |
|---------|-------------------|------------|-------|
| **Command palette (Cmd+K)** | Power users jump between tools instantly. No other cheminformatics web tool has this. | Medium | shadcn/ui `Command` component (cmdk). Index all tools/pages. Fuzzy search. Highest-impact single feature. |
| **Animated molecule structure cards** | Current MoleculeCard is static. Hover effects, loading shimmer, expand/collapse. | Low | framer-motion `whileHover`, `layoutId` for shared-element transitions. |
| **Animated tab/view transitions** | Crossfade/slide when switching tools within a page (e.g., Descriptors to Tanimoto). | Low | `AnimatePresence` with `mode="wait"` around active view. Fade + slide-up (200ms). |
| **Micro-interactions on forms** | Button press scale-down, submit spinner->checkmark, focus glow. | Low | framer-motion `whileTap={{ scale: 0.98 }}`. Loading state on submit buttons. |
| **Inline molecule preview on SMILES input** | Live 2D preview when typing SMILES. Debounced API call (500ms). | Medium | Call `/depict/2D` endpoint. Small preview (64x64) next to input. Loading skeleton while fetching. |
| **Breadcrumb trail for nested tool views** | When in ChemPage > Descriptors, show orientation breadcrumb. | Low | shadcn/ui Breadcrumb component. Generate from route params. |
| **Resizable input/output panels** | Side-by-side input form and results with draggable divider on wide screens. | Medium | shadcn/ui ResizablePanelGroup (wraps react-resizable-panels). |

## Anti-Features

Features to explicitly NOT build during this migration.

| Anti-Feature | Why Avoid | What to Do Instead |
|--------------|-----------|-------------------|
| **TypeScript migration** | Scope creep. All 55+ files are .js/.jsx. Converting adds weeks with zero user-visible benefit. | Keep .js/.jsx. TypeScript is NOT required for shadcn/ui or any target stack component. |
| **Server-side rendering (SSR)** | Client-side SPA. No SEO requirements. Backend is separate FastAPI. | Stay with client-side Vite SPA. |
| **State management library** | AppContext has 5 state values. No complex state. Redux/Zustand is overkill. | Keep AppContext as-is. |
| **User accounts / authentication** | Backend is stateless, no auth. Would require backend changes (out of scope). | Use localStorage for personalization. |
| **Molecule drawing editor** | JSME/Ketcher integration is a full project on its own. | Keep current approach. Future milestone if needed. |
| **Internationalization (i18n)** | Chemistry notation is universal. Minimal UI text. English-standard domain. | Keep English only. |
| **React Router v7 migration** | v6 is stable and current. v7 is a different paradigm (framework-mode). Unnecessary risk. | Stay on react-router-dom v6. |
| **Custom theme builder** | Dark + light sufficient. shadcn/ui provides complete theming system. | Use shadcn/ui CSS variable theming. |
| **PWA / offline support** | Requires active API for all features. Offline mode non-functional. | Skip entirely. |
| **PDF report generation** | Large scope for minimal value. Browser print works. | Add print-friendly CSS if needed. |

## Feature Dependencies

```
Vite migration ---------> React 19 upgrade --------> Tailwind v4 migration --> shadcn/ui integration
                                                                                      |
                                    +-------------------------------------------------+
                                    |
                                    v
                          Component replacement (Button, Card, Input, Select, Tabs, Toast, Dialog)
                                    |
                                    v
                          Dark mode polish (system preference, CSS variables)
                                    |
                                    v
                          Page-by-page UI polish + transitions + micro-interactions
                                    |
                                    v
                          Differentiators (Command palette, inline preview, resizable panels)
```

## MVP Recommendation

Prioritize in this order for maximum impact with minimum effort:

### Phase 1: Foundation (must complete first)
1. **Vite migration** -- unblocks everything, fixes CRA deprecation and security vulns
2. **React 19 upgrade** -- automated codemod handles most changes
3. **Tailwind v4 migration** -- automated upgrade tool + manual class fixes
4. **shadcn/ui integration** -- component system foundation + dark mode theming

### Phase 2: Core Polish
1. **Dark mode with system preference** -- enhance existing implementation
2. **Toast notifications (Sonner)** -- replace alert() and inline states
3. **Skeleton loading states** -- replace full-page loader
4. **Page transitions** -- AnimatePresence on routes
5. **Component-by-component replacement** -- swap hand-rolled for shadcn/ui

### Phase 3: Differentiators
1. **Command palette (Cmd+K)** -- highest-impact differentiator
2. **Tab/view transitions** -- low effort, big polish
3. **Micro-interactions** -- low effort, premium feel
4. **Inline molecule preview** -- science-specific delight
5. **Breadcrumb navigation** -- orientation for nested views

### Defer
- **Result comparison mode**: High complexity, Phase 4+ at earliest
- **Persistent workspace**: Medium complexity, build incrementally
- **Resizable panels**: Nice-to-have, not critical

## Sources

- Existing codebase analysis (55+ source files, package.json, tailwind.config.js, CSS files)
- React 19 features: https://react.dev/blog/2024/12/05/react-19 (verified via WebFetch)
- Tailwind v4 features: https://tailwindcss.com/docs/installation (verified via WebFetch, v4.2)
- shadcn/ui patterns: Training data (MEDIUM confidence -- verify at https://ui.shadcn.com/docs)
