# Phase 6: Animations + Power Features - Context

**Gathered:** 2026-03-13
**Status:** Ready for planning

<domain>
## Phase Boundary

Polish the UI with smooth page transitions, button micro-interactions, staggered list animations, and tab content transitions. Add power user features: Cmd+K command palette for quick navigation to pages/tools/molecules, inline 2D SMILES structure preview below inputs, and breadcrumb navigation showing current section and tool. All animations must feel premium and fluid while maintaining the "elegance over complexity" principle established in Phase 4.1.

</domain>

<decisions>
## Implementation Decisions

### Page Transitions (ANIM-01)
- Cross-fade with subtle vertical rise: opacity 0→1 + translateY(8px→0) on enter, reverse on exit
- Overlap timing (AnimatePresence mode="sync"): old page fades out while new page fades in simultaneously
- Duration ~250ms using existing spring physics (stiffness 100, damping 20)
- Wraps around the Suspense/Outlet in Layout component

### Tab Content Transitions (ANIM-04)
- Subtle cross-fade (~150ms) when switching tools within tab-based pages (ChemPage, DepictPage, ConvertPage, ToolsPage)
- Tab indicator already animates via layoutId — content area adds matching fade transition
- Quick and lightweight, no directional sliding

### Button Micro-interactions (ANIM-02)
- Subtle hover: translateY(-2px) + box-shadow with soft glow matching page gradient palette (opacity ~20%)
- Subtle press: scale(0.97) + translateY(0)
- Implemented via CSS transitions (200ms ease hover, 100ms ease press) — NOT motion.div wrappers, for performance across many buttons
- Applied globally to interactive buttons

### Focus Ring Animations (ANIM-03)
- Claude's discretion on exact implementation
- Should feel consistent with the subtle button micro-interactions

### Staggered List Animations (ANIM-05)
- Results lists and data displays animate in with staggered timing
- Claude's discretion on exact stagger delay and animation style
- Should use motion/react for consistency with existing animation patterns

### Layout Animations (ANIM-06)
- Content state changes (expand/collapse, show/hide) animate layout shifts
- Claude's discretion on exact implementation
- Should use motion/react layout animations

### Command Palette (POWER-01)
- Triggered by Cmd+K (Mac) / Ctrl+K (Windows)
- Four searchable scopes: pages (7), individual tools (~21), example molecules (6), recent molecules (from AppContext)
- Glass-styled modal: backdrop-blur + translucent background matching glassmorphism design system
- Categorized results grouped under headers: Pages, Tools, Molecules, Recent
- Fuzzy search across all categories
- Selecting a molecule navigates to /chem with that SMILES pre-filled in the input
- Selecting a tool deep-links to the correct page + tab (e.g., /chem/descriptors)
- Keyboard navigation (arrow keys + Enter) for power users

### Inline SMILES Preview (POWER-02)
- 2D structure preview card appears below the SMILES input field when valid SMILES is entered
- Glass-styled card, ~120x120px structure image
- Rendered via backend API call to /latest/depict/2D with 500ms debounce after typing stops
- Glass shimmer skeleton shown while loading
- Structure image only — no molecule name or formula (keep it simple)
- Invalid SMILES or API errors: preview card hidden silently (no error messages — user is still typing)
- Integrated into the SMILESInput component used across all 23 tool views

### Breadcrumb Navigation (POWER-03)
- Appears below the glass header, above page content
- Two-level hierarchy: Section > Tool (e.g., "Chemical Analysis > Descriptors")
- First level (section name) is clickable — navigates to page with default tool
- Only shows on tool pages (Chem, Convert, Depict, Tools) — hidden on Home, About, OCSR, Privacy, Terms
- Desktop only (hidden on mobile <768px) — tab bar already provides context on mobile
- Tool name cross-fades (~150ms) when switching tabs; section name stays static

### Claude's Discretion
- Focus ring animation specifics (ANIM-03)
- Staggered list animation timing and style (ANIM-05)
- Layout animation details for expand/collapse (ANIM-06)
- Command palette library choice (cmdk, custom, or shadcn/ui command component)
- Command palette keyboard navigation details
- Exact glass styling values for command palette and SMILES preview
- Breadcrumb separator character and typography
- Performance optimization (will-change, GPU layers, animation cleanup)

</decisions>

<specifics>
## Specific Ideas

- "Elegance over complexity" — carried from Phase 4.1, all animations should be subtle and refined
- Button micro-interactions via CSS transitions, not motion.div, for performance across dozens of buttons
- Command palette grouped by category like VS Code / Linear — familiar pattern for power users
- SMILES preview uses the existing /latest/depict/2D backend endpoint — no new dependencies needed
- Page transitions use existing spring physics from Phase 2 (stiffness 100, damping 20)
- Breadcrumbs are a wayfinding enhancement, not a primary navigation mechanism

</specifics>

<code_context>
## Existing Code Insights

### Reusable Assets
- motion/react: Already used in 24 files with AnimatePresence, LayoutGroup, layoutId, spring physics
- Navigation.tsx: Has LayoutGroup + layoutId for animated active pill — pattern for page identity
- SMILESInput.tsx: Shared across 23 tool views — integration point for inline preview
- AppContext.tsx: Has recentMolecules array — feeds into command palette recent section
- RouteLoadingFallback.tsx: Existing Suspense fallback — page transitions wrap around this
- GlassSkeleton: Existing shimmer skeleton component — reuse for SMILES preview loading state
- Glass @utility in tailwind.css: Existing glassmorphism utilities — apply to command palette + preview card

### Established Patterns
- Spring physics: stiffness 500/damping 30 (tabs), stiffness 100/damping 20 (pages) — Phase 2
- motion.div wrapper pattern for component animations
- AnimatePresence for enter/exit animations (used in Footer, ComparisonTray)
- Clay-interactive pattern on buttons (soft shadows, raised feel) — Phase 4.1
- Dark mode via .dark class on documentElement
- cn() utility for conditional class merging
- Per-page gradient palettes from Phase 4.1

### Integration Points
- App.tsx Layout component: Wrap Outlet with AnimatePresence for page transitions
- All 4 tab-based pages: Add AnimatePresence to tab content area for cross-fade
- tailwind.css: Add button hover/press animation utilities
- SMILESInput.tsx: Add debounced preview below input field
- New components: CommandPalette, Breadcrumbs, SMILESPreview
- App.tsx: Add CommandPalette to Layout, add global Cmd+K listener
- Navigation data (navLinks in Navigation.tsx): Reuse for command palette page entries

</code_context>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 06-animations-power-features*
*Context gathered: 2026-03-13*
