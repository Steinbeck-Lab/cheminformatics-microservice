# Phase 4: Component System + Dark Mode - Context

**Gathered:** 2026-03-12
**Status:** Ready for planning

<domain>
## Phase Boundary

Replace all hand-rolled UI elements (buttons, cards, dialogs, form inputs, tabs) with shadcn/ui components across all 9 pages. Migrate the theming system from custom CSS variables (--bg-primary, etc.) to shadcn/ui's standard HSL-based variable system. Migrate all icons from react-icons to lucide-react. Add system preference auto-detection for dark mode and prevent theme flash on page load. Maintain all existing motion/react animations by wrapping shadcn/ui components.

</domain>

<decisions>
## Implementation Decisions

### CSS Variable Theming
- Full replacement: remove all custom --bg-primary, --text-primary, --border-primary, etc. CSS variables
- Adopt shadcn/ui's standard variable naming (--background, --foreground, --card, --primary, --muted, --border, --input, --ring, etc.)
- Map current sky/slate color palette into the shadcn/ui variable system — the app should look the same, just with standardized variable names underneath
- Install clsx + tailwind-merge, create lib/utils.ts with standard cn() helper function
- Remove all btn-*, card-*, card-header, card-body, card-footer, form-group, form-label @utility classes from tailwind.css — shadcn/ui components replace them
- Update custom utilities (glass, text-gradient) to reference new shadcn/ui CSS variables
- Remove @tailwindcss/forms plugin (shadcn/ui Input/Select/Textarea handle their own styling)
- Keep @tailwindcss/typography plugin (needed for prose content on About, Privacy, Terms pages)

### Component Set
- Install shadcn/ui components: Button, Card, Dialog, Sheet, Input, Select, Textarea, Tabs
- Migration approach: component-by-component across ALL 9 pages (e.g., all buttons first, then all cards, then all inputs)
- Ensures consistency — no page has mixed old/new component styles at any point
- Wrap shadcn/ui components inside motion.div wrappers to preserve existing page transitions, stagger animations, and spring physics
- The shadcn component handles styling and accessibility; motion handles animation

### Tab Implementation
- 4 pages use tabbed interfaces: ChemPage, DepictPage, ConvertPage, ToolsPage
- Currently use custom tab implementations with motion LayoutGroup for animated indicators

### Icon Migration
- Full swap from react-icons to lucide-react in one pass across all 45 files
- Replace ALL icon families (HiOutline Heroicons + FontAwesome) with lucide equivalents
- Uninstall react-icons entirely after migration — one icon library for the whole app
- Standardize icon sizes: 16px (inline/small), 20px (buttons/nav), 24px (headers/features)

### Dark Mode Toggle
- Keep the existing animated pill toggle in Header (spring physics, LayoutGroup animation)
- Swap HiOutlineSun/HiOutlineMoon icons for lucide Sun/Moon icons
- Keep the mobile theme toggle in the mobile menu

### System Preference Detection
- First-visit only: auto-detect system preference via prefers-color-scheme when no localStorage value exists
- Once the user manually toggles, their choice overrides system preference permanently
- Change default from hardcoded dark mode to system preference for new visitors
- Add inline blocking script in index.html that reads localStorage and sets .dark class before React renders — eliminates flash of wrong theme (FOUC)

### Claude's Discretion
- Color format for CSS variables (HSL vs OKLCH) — pick best for shadcn/ui compatibility and Tailwind v4
- Tab component approach: shadcn Tabs with motion indicator, shadcn Tabs as-is, or keep current custom tabs — evaluate per page
- lucide icon stroke weight (default 2px vs lighter 1.5px) — match current design aesthetic
- Any additional shadcn/ui components beyond the required set if clearly needed during migration

</decisions>

<specifics>
## Specific Ideas

- "Avoid band-aids, shortcuts — this is a production system, accuracy first" (carried from Phase 1)
- "Make everything future proof" (carried from Phase 2)
- Component-by-component migration ensures zero inconsistency across pages at any point
- Pill toggle is a signature UI element — preserve its animation quality while swapping to lucide icons
- FOUC prevention is critical for professional feel — blocking script before React renders

</specifics>

<code_context>
## Existing Code Insights

### Reusable Assets
- AppContext (context/AppContext.tsx): Already has isDarkMode, toggleDarkMode, localStorage persistence — needs system preference detection added
- Header (components/common/Header.tsx): Animated pill toggle with spring physics — swap icons, keep animation
- Navigation (components/common/Navigation.tsx): Desktop + mobile nav — component styling will be updated
- tailwind.css: @custom-variant dark already configured, base styles and @utility classes to be replaced

### Established Patterns
- Dark mode via .dark class on document.documentElement (class-based, already works)
- motion/react used in 15+ files: motion.div, AnimatePresence, LayoutGroup, useScroll, useTransform
- Tab-based pages use custom LayoutGroup indicators for animated tab switching
- CSS variables in :root/.dark blocks for light/dark theming
- TypeScript strict mode with React 19 — shadcn/ui components will integrate cleanly

### Integration Points
- tailwind.css: Major rewrite — remove old variables and utility classes, add shadcn/ui variables
- index.html: Add inline theme detection script before React bundle
- AppContext.tsx: Add prefers-color-scheme detection, change default from dark to system preference
- package.json: Add lucide-react, clsx, tailwind-merge, @radix-ui/* peer deps; remove react-icons
- 45 files: Icon import replacements (react-icons → lucide-react)
- 9 pages: Component replacements (hand-rolled → shadcn/ui)
- components/ui/: New directory for shadcn/ui component files

</code_context>

<deferred>
## Deferred Ideas

- Bento grid layouts for pages — Phase 5 (UX) or Phase 6 (from Phase 3 backlog)
- Liquid glass visual effect — Phase 6 (from Phase 3 backlog)
- Claymorphism/glassmorphism combination styles — Phase 6 (from Phase 3 backlog)
- Toast notifications (Sonner) — Phase 5 (LOAD-02)
- Skeleton loading states — Phase 5 (LOAD-01)
- Command palette (Cmd+K) — Phase 6 (POWER-01)

</deferred>

---

*Phase: 04-component-system-dark-mode*
*Context gathered: 2026-03-12*
