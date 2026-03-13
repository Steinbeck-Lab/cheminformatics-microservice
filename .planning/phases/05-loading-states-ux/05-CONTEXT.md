# Phase 5: Loading States + UX - Context

**Gathered:** 2026-03-13
**Status:** Ready for planning

<domain>
## Phase Boundary

Add skeleton loading states, toast notifications (Sonner), inline error states with retry, and per-component loading indicators across all tool pages. Make the layout responsive on mobile (375px) and tablet (768px) with no horizontal overflow. Improve load performance via route-level code splitting and lazy loading. Replace all `alert()` calls and full-page LoadingScreen with modern alternatives.

</domain>

<decisions>
## Implementation Decisions

### Skeleton & Loading Style
- Glass shimmer skeletons: frosted glass appearance with a translucent highlight sweep animation (left-to-right shimmer)
- Content-shaped skeletons: mirror the actual result layout (molecule image placeholder, text line blocks, property rows)
- Replace the full-page LoadingScreen overlay everywhere — all loading becomes per-component glass skeletons in the output area
- Submit buttons show inline spinner icon + text change (e.g., "Analyzing..." with spinning loader), disabled during loading
- Input area stays usable/visible during loading — only the output area shows skeletons

### Toast Notifications
- Sonner library for toast system
- Position: bottom-right corner
- Glass-styled toasts: backdrop-blur + translucent background matching glassmorphism design system
- Auto-dismiss: 3 seconds for success toasts, 5 seconds for error toasts
- Toast scope: toasts for side-effects only (downloads, copy-to-clipboard, network errors)
- Inline feedback for primary content: API results appear in output panel, validation errors shown near input, error + retry button in output area
- Replace all 3 existing alert() calls (MoleculeCard, HighlightedMoleculeCard, NPlikenessView)

### Error & Retry UX
- Glass error cards: frosted glass background with red/amber accent left border, error icon, user-friendly message, and Retry button
- User-friendly error messages only — no status codes or raw error text shown to user. Technical details logged to browser console
- Retry button scrolls to / focuses the input area so user can modify and resubmit
- Two-layer error handling: global Axios interceptor catches network/5xx errors (shows toast), per-component handling for 4xx/business errors (shows inline glass error card)
- Error messages are contextual: "Could not analyze this molecule", "Service temporarily unavailable" — not generic "Something went wrong"

### Responsive Adaptation
- Stacked mobile layout: input above output vertically on mobile (<1024px), resizable panels disabled
- Tab strip becomes horizontally scrollable on mobile
- Keep existing Sheet hamburger menu — polish tap targets to 44px minimum, comfortable spacing
- Remove 3D CaffeineMolecule on mobile (<768px) for GPU performance — show glass hero card with title + CTA only
- Feature cards stack vertically in single column on mobile
- Glass effects reduced on mobile: backdrop-blur 12px instead of 24px (already established in Phase 4.1)
- Standard Tailwind breakpoints (sm:640, md:768, lg:1024, xl:1280) — no custom breakpoints

### Claude's Discretion
- Code splitting strategy: React.lazy + Suspense for route-level splitting, Suspense fallback component design
- Navigation wayfinding improvements (UX-03): active nav highlighting, visual indicators for current page
- Lighthouse performance optimization targets and techniques
- Empty state design for tools with no results yet
- Exact shimmer animation CSS implementation (keyframes, gradient technique)
- Which components get content-shaped skeletons vs simpler glass block skeletons

</decisions>

<specifics>
## Specific Ideas

- "Accuracy first, no band-aids" — production quality carried from Phase 1
- Glass shimmer skeletons should feel like frosted glass "loading" — cohesive with the glassmorphism design language
- Error cards use left accent border (red/amber) similar to VS Code warning panels — clear but not alarming
- Toast notifications should feel like glass floating notifications — backdrop-blur matches the header/footer glass style
- Mobile experience should feel intentionally designed, not just "desktop squished down"

</specifics>

<code_context>
## Existing Code Insights

### Reusable Assets
- `LoadingScreen.tsx`: Full-page molecule animation overlay — to be replaced everywhere with per-component skeletons, then potentially deleted
- `glass` @utility in tailwind.css: Existing glassmorphism utility (blur + translucent bg) — skeleton shimmer builds on this
- shadcn/ui components (Button, Card, Input, Select, Tabs): Already established from Phase 4
- motion/react: Animation library for skeleton shimmer, error card entrance, toast transitions
- AppContext.tsx: Has global isLoading state — may need refactoring for per-component loading
- Axios service layer (services/): Integration point for global error interceptor

### Established Patterns
- Glass aesthetic: backdrop-blur 24px (desktop), 12px (mobile), per-page gradient palettes (Phase 4.1)
- Dark mode via .dark class on document.documentElement
- cn() utility from lib/utils.ts for conditional class merging
- motion.div wrapper pattern for shadcn/ui component animations
- ResizablePanelGroup with autoSaveId for localStorage persistence (Phase 4.1)

### Integration Points
- App.tsx: Route-level code splitting with React.lazy + Suspense for all 9 pages (currently eagerly imported)
- services/*.ts: Add global Axios response interceptor for network/5xx error toasts
- All tool view components: Replace loading/error patterns with skeleton + glass error card
- tailwind.css: Add skeleton shimmer animation utilities (glass-shimmer keyframes)
- package.json: Add sonner dependency
- 3 files with alert(): MoleculeCard.tsx, HighlightedMoleculeCard.tsx, NPlikenessView.tsx — replace with toast.error()

</code_context>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 05-loading-states-ux*
*Context gathered: 2026-03-13*
