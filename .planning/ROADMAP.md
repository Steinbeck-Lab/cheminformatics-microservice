# Roadmap: Cheminformatics Microservice Frontend Overhaul

## Overview

This roadmap modernizes the Cheminformatics Microservice frontend from a deprecated CRA + React 18 + Tailwind v3 stack to a production-grade Vite + React 19 + Tailwind v4 + shadcn/ui stack with dark mode and polished UX. The migration follows a strict dependency chain: Vite must come first (Tailwind v4 requires Vite plugin, shadcn CLI assumes Vite), then React 19 (independent upgrade, but required before shadcn/ui), then Tailwind v4 (CSS variable foundation for theming), then shadcn/ui + dark mode (component system and theming together), and finally UX and animation polish layered on top of the stable component system.

## Phases

**Phase Numbering:**
- Integer phases (1, 2, 3): Planned milestone work
- Decimal phases (2.1, 2.2): Urgent insertions (marked with INSERTED)

Decimal phases appear between their surrounding integers in numeric order.

- [x] **Phase 1: Vite Migration** - Replace deprecated CRA with Vite, resolve all security vulnerabilities, update dependencies (completed 2026-03-12)
- [x] **Phase 2: React 19 + TypeScript 5** - Upgrade React 18 to 19, TypeScript 4.9 to 5.x, verify framer-motion compatibility (completed 2026-03-12)
- [x] **Phase 3: Tailwind v4 Migration** - Upgrade Tailwind CSS v3 to v4 with CSS-first config and class rename resolution (completed 2026-03-12)
- [ ] **Phase 4: Component System + Dark Mode** - Integrate shadcn/ui components across all pages with CSS variable theming and dark mode (gap closure in progress)
- [ ] **Phase 5: Loading States + UX** - Add skeleton loaders, toast notifications, error states, responsive improvements, and performance optimization (gap closure in progress)
- [ ] **Phase 6: Animations + Power Features** - Polish with smooth transitions, micro-interactions, command palette, and inline SMILES preview

## Phase Details

### Phase 1: Vite Migration
**Goal**: The frontend builds and serves on Vite with zero security vulnerabilities and all dependencies at latest stable versions
**Depends on**: Nothing (first phase)
**Requirements**: BUILD-01, BUILD-02, BUILD-03, BUILD-04, BUILD-05, BUILD-06
**Success Criteria** (what must be TRUE):
  1. Running `npm run dev` starts the app via Vite and all 9 pages load correctly in the browser
  2. Running `npm run build` produces a production bundle with no errors or warnings
  3. `npm audit` reports zero vulnerabilities and package.json contains no `overrides` section
  4. Docker frontend build (Dockerfile.frontend) succeeds and serves the app correctly
  5. All environment variables use the VITE_* prefix and the app communicates with the backend API
**Plans**: 3 plans

Plans:
- [x] 01-01-PLAN.md -- Core Vite setup: uninstall CRA, install Vite stack, migrate env vars, update deps
- [x] 01-02-PLAN.md -- Developer tooling: ESLint v9 flat config, Vitest, husky, CI workflow
- [x] 01-03-PLAN.md -- Docker update and full migration verification (9-page visual check)

### Phase 2: React 19 + TypeScript 5
**Goal**: The application runs on React 19 and TypeScript 5 with all existing functionality intact, and animations use the modern motion package
**Depends on**: Phase 1
**Requirements**: FRAME-01, FRAME-03, FRAME-05
**Success Criteria** (what must be TRUE):
  1. React 19.x is installed and the app renders all 9 pages without console errors or warnings
  2. TypeScript 5.x compiles the project with no type errors
  3. All animations work correctly using the motion package (page transitions, component animations)
**Plans**: 2 plans

Plans:
- [x] 02-01-PLAN.md -- React 19 upgrade, TypeScript conversion (61 files .js/.jsx to .ts/.tsx), tsconfig, type definitions, ESLint TypeScript support
- [x] 02-02-PLAN.md -- Motion package migration (framer-motion to motion/react), animation enhancements, visual verification

### Phase 3: Tailwind v4 Migration
**Goal**: All styling runs on Tailwind CSS v4 with CSS-first configuration and no visual regressions
**Depends on**: Phase 2
**Requirements**: FRAME-02, FRAME-04
**Success Criteria** (what must be TRUE):
  1. Tailwind v4 is installed with the @tailwindcss/vite plugin and the old postcss/tailwind JS config files are removed
  2. All 9 pages render with correct styling -- no missing colors, broken layouts, or invisible elements from silent class renames
  3. CSS variables in tailwind.css use Tailwind v4 @theme block syntax instead of the v3 theme() function
**Plans**: 2 plans

Plans:
- [x] 03-01-PLAN.md -- Upgrade tool + package swap, CSS-first config migration (@theme, @plugin, @custom-variant), theme() replacement, class renames
- [x] 03-02-PLAN.md -- Edge case audit (transform-none, ring defaults, bg-opacity merging), final sweep, 9-page visual verification

### Phase 4: Component System + Dark Mode
**Goal**: All UI elements use shadcn/ui components with a unified CSS variable theming system and functional dark mode toggle
**Depends on**: Phase 3
**Requirements**: COMP-01, COMP-02, COMP-03, COMP-04, COMP-05, COMP-06, COMP-07, THEME-01, THEME-02, THEME-03, THEME-04, THEME-05
**Success Criteria** (what must be TRUE):
  1. shadcn/ui is installed and configured, with a components/ui/ directory containing Button, Card, Dialog, Input, Select, Textarea, and Sheet components
  2. All 9 pages use shadcn/ui components consistently -- no mix of hand-rolled and shadcn/ui buttons, cards, or form inputs
  3. A dark mode toggle is visible in the header/navigation and switching it changes the entire UI between light and dark themes
  4. Theme preference persists across page reloads (localStorage) and auto-detects system preference on first visit
  5. Icons throughout the app use lucide-react instead of react-icons
**Plans**: 5 plans

Plans:
- [x] 04-01-PLAN.md -- shadcn/ui foundation: install deps, path aliases, cn() utility, CSS variable rewrite to OKLCH, FOUC prevention, system preference detection
- [x] 04-02-PLAN.md -- Icon migration: replace all react-icons and @fortawesome with lucide-react across 43+ files, uninstall old packages
- [x] 04-03-PLAN.md -- Component migration: replace all buttons with shadcn Button, cards with Card, form inputs with Input/Select/Textarea
- [x] 04-04-PLAN.md -- Final pass: mobile menu to Sheet, isDarkMode ternary cleanup, consistency sweep, visual verification (both modes)
- [ ] 04-05-PLAN.md -- Gap closure: migrate remaining hand-rolled buttons and cards in OCSRPage, HomePage, AboutPage, and Footer to shadcn/ui

### Phase 04.1: Visual Design System — Glassmorphism, Liquid Glass & Bento Grids (INSERTED)

**Goal**: All pages use a cohesive glassmorphism + liquid glass fusion aesthetic with bento grid layouts, resizable panels, animated molecule cards, and result comparison mode
**Depends on**: Phase 4
**Requirements**: VISUAL-01, VISUAL-02, VISUAL-03, VISUAL-04, VISUAL-05, VISUAL-06
**Success Criteria** (what must be TRUE):
  1. Cards, containers, and panels across all 9 pages use glassmorphism styling (frosted glass backgrounds, subtle borders, backdrop blur)
  2. Tool pages use bento grid layouts with varied, editorial-style sizing instead of uniform grids
  3. Molecule cards animate on hover with shimmer/glow effects and support expand/collapse
  4. Wide-screen tool pages have resizable input/output panels with a draggable divider
  5. Users can compare two molecule results side-by-side in a comparison mode
  6. Claymorphism accents (soft shadows, raised surfaces) blend with glass effects for depth
**Plans**: 7 plans

Plans:
- [x] 04.1-00-PLAN.md -- Wave 0: Test stub files for all 6 VISUAL requirements (vitest scaffolds)
- [x] 04.1-01-PLAN.md -- Foundation: CSS glass/clay/gradient utilities, per-page gradient mesh config, bento layout config, useMediaQuery hook, GradientMesh and LiquidGlassFilter components
- [x] 04.1-02-PLAN.md -- BentoGrid component + ChemPage migration from sidebar tabs to bento grid (14 tools)
- [x] 04.1-03-PLAN.md -- Remaining bento pages: ConvertPage, DepictPage, ToolsPage migration + HomePage glass hero and bento feature cards
- [x] 04.1-04-PLAN.md -- Resizable panels (shadcn Resizable + ResizableToolPanel), MoleculeCard glow/lift/expand enhancements, BentoCard-to-ResizableToolPanel wiring
- [x] 04.1-05-PLAN.md -- Comparison mode: ComparisonContext, ComparisonTray, ComparisonView, AddToCompareButton, App.tsx wiring, ChemPage compare integration
- [x] 04.1-06-PLAN.md -- Glass sweep: Header/Footer glass treatment, OCSR/About/Privacy/Terms gradient mesh + glass containers, refraction removed and tool pages reverted to tabs per user feedback

### Phase 5: Loading States + UX
**Goal**: The application provides clear feedback during all async operations and works well across all viewport sizes with measurably faster load times
**Depends on**: Phase 4
**Requirements**: LOAD-01, LOAD-02, LOAD-03, LOAD-04, UX-01, UX-02, UX-03, UX-04
**Success Criteria** (what must be TRUE):
  1. API calls show skeleton loading placeholders instead of full-page spinners, and individual operations show per-component loading indicators
  2. Success and error feedback uses toast notifications (Sonner) instead of browser alert() calls
  3. Failed API calls show inline error messages with a visible retry button on the affected component
  4. The layout is usable on mobile (375px) and tablet (768px) viewports with no horizontal overflow or overlapping elements
  5. Lighthouse performance score improves measurably over the Phase 4 baseline due to route-level code splitting and lazy loading
**Plans**: 4 plans

Plans:
- [x] 05-01-PLAN.md -- Foundation: feedback components (GlassSkeleton, GlassErrorCard, ToolSkeleton, EmptyState), Sonner toast system, error message mapping, Axios interceptor
- [x] 05-02-PLAN.md -- Tool view refactoring: replace LoadingScreen in all 21 tool views with glass skeletons + error cards, replace alert() calls, delete LoadingScreen
- [x] 05-03-PLAN.md -- Code splitting (React.lazy + vendor chunks), responsive fixes (mobile 3D, scrollable tabs, tap targets), navigation wayfinding
- [ ] 05-04-PLAN.md -- Gap closure: replace NPlikenessView alert() with info expansion box, add Loader2 spinner to all tool view submit buttons

### Phase 6: Animations + Power Features
**Goal**: The UI feels premium and fluid with smooth transitions throughout, and power users can navigate instantly via command palette
**Depends on**: Phase 5
**Requirements**: ANIM-01, ANIM-02, ANIM-03, ANIM-04, ANIM-05, ANIM-06, POWER-01, POWER-02, POWER-03
**Success Criteria** (what must be TRUE):
  1. Navigating between pages shows smooth animated transitions (not abrupt content swaps)
  2. Buttons respond to hover with a lift/glow effect and to press with a scale-down animation
  3. Results lists and data displays animate in with staggered timing, and tab/section changes animate content layout shifts
  4. Pressing Cmd+K (Mac) or Ctrl+K (Windows) opens a command palette that allows quick navigation to any tool page
  5. Typing a SMILES string in any input field shows an inline 2D structure preview of the molecule
**Plans**: TBD

Plans:
- [ ] 06-01: TBD
- [ ] 06-02: TBD
- [ ] 06-03: TBD

## Progress

**Execution Order:**
Phases execute in numeric order: 1 -> 2 -> 3 -> 4 -> 4.1 -> 5 -> 6

| Phase | Plans Complete | Status | Completed |
|-------|----------------|--------|-----------|
| 1. Vite Migration | 3/3 | Complete   | 2026-03-12 |
| 2. React 19 + TypeScript 5 | 2/2 | Complete | 2026-03-12 |
| 3. Tailwind v4 Migration | 2/2 | Complete   | 2026-03-12 |
| 4. Component System + Dark Mode | 4/5 | Gap closure | 2026-03-13 |
| 4.1 Visual Design System | 7/7 | Complete | 2026-03-13 |
| 5. Loading States + UX | 3/4 | Gap closure | - |
| 6. Animations + Power Features | 0/3 | Not started | - |
