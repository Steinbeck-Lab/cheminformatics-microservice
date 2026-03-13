# Requirements: Cheminformatics Microservice Frontend Overhaul

**Defined:** 2026-03-12
**Core Value:** The frontend must feel noticeably modern, fast, and polished — while never breaking existing functionality. Zero security alerts, production-grade quality.

## v1 Requirements

Requirements for this milestone. Each maps to roadmap phases.

### Build & Infrastructure

- [x] **BUILD-01**: Frontend builds and runs on Vite instead of Create React App (react-scripts)
- [x] **BUILD-02**: All dependabot security vulnerabilities are resolved (zero alerts)
- [x] **BUILD-03**: All package.json overrides are eliminated (proper dependency upgrades)
- [x] **BUILD-04**: All dependencies updated to latest stable versions
- [x] **BUILD-05**: Docker frontend build (Dockerfile.frontend) works with Vite
- [x] **BUILD-06**: Environment variables migrated from REACT_APP_* to VITE_* prefix

### Framework Upgrades

- [x] **FRAME-01**: React upgraded from 18.2 to 19.x (latest stable)
- [x] **FRAME-02**: Tailwind CSS upgraded from v3.3 to v4.x with CSS-first config
- [x] **FRAME-03**: TypeScript upgraded from 4.9.5 to 5.x (latest stable)
- [x] **FRAME-04**: All Tailwind v3 class renames and breaking changes resolved
- [x] **FRAME-05**: Framer Motion verified compatible with React 19

### Component System

- [x] **COMP-01**: shadcn/ui installed and configured for Vite + Tailwind v4
- [x] **COMP-02**: Core buttons replaced with shadcn/ui Button component
- [x] **COMP-03**: Cards and content containers replaced with shadcn/ui Card
- [x] **COMP-04**: Modals and dialogs replaced with shadcn/ui Dialog/Sheet
- [x] **COMP-05**: Form inputs replaced with shadcn/ui Input/Select/Textarea
- [x] **COMP-06**: Icons migrated from react-icons to lucide-react
- [x] **COMP-07**: Consistent component usage across all 9 pages

### Theming & Dark Mode

- [x] **THEME-01**: Dark mode toggle available in UI (header/nav)
- [x] **THEME-02**: System preference auto-detection (prefers-color-scheme)
- [x] **THEME-03**: Theme preference persisted in localStorage
- [x] **THEME-04**: CSS variable-based theming using shadcn/ui approach
- [x] **THEME-05**: All pages and components render correctly in both light and dark modes

### Visual Design System

- [x] **VISUAL-01**: Glassmorphism + liquid glass fusion aesthetic applied to cards, containers, and panels across all pages
- [x] **VISUAL-02**: Glass-styled tabbed navigation for tool pages (tab-based UX with glassmorphism aesthetic, replacing the originally planned bento grid layouts per user preference)
- [x] **VISUAL-03**: Animated molecule structure cards with hover shimmer, glow effects, and expand/collapse
- [x] **VISUAL-04**: Resizable input/output panels with draggable divider on wide screens
- [x] **VISUAL-05**: Result comparison mode for side-by-side molecule analysis
- [x] **VISUAL-06**: Claymorphism accents blended with glass effects for depth and tactile feel

### Loading & Feedback

- [x] **LOAD-01**: Skeleton loading states replace full-page spinners during API calls
- [x] **LOAD-02**: Toast notification system (Sonner) replaces alert() calls
- [x] **LOAD-03**: Per-component error states with clear messages and retry actions
- [x] **LOAD-04**: Loading indicators for individual API operations (not full-page)

### Animations & Polish

- [x] **ANIM-01**: Smooth page transitions when navigating between routes
- [x] **ANIM-02**: Micro-interactions on buttons (press scale, hover glow/lift)
- [x] **ANIM-03**: Focus ring animations on interactive elements
- [x] **ANIM-04**: View transitions when switching tabs within tool pages
- [x] **ANIM-05**: Staggered list animations for results and data displays
- [x] **ANIM-06**: Layout animations on content state changes (expand/collapse, show/hide)

### Power User Features

- [ ] **POWER-01**: Cmd+K / Ctrl+K command palette for quick tool navigation
- [ ] **POWER-02**: Inline molecule structure preview when entering SMILES strings
- [ ] **POWER-03**: Breadcrumb navigation showing current location in tool hierarchy

### UX Improvements

- [x] **UX-01**: Responsive layout works well on mobile and tablet viewports
- [x] **UX-02**: Code splitting with React.lazy for route-level lazy loading
- [x] **UX-03**: Improved navigation structure and visual wayfinding
- [x] **UX-04**: Faster initial page load (measurable improvement via Lighthouse)

## v2 Requirements

Deferred to future release. Tracked but not in current roadmap.

### Advanced Features

- **ADV-01**: Interactive 3D molecule viewer (3Dmol.js or similar)
- **ADV-02**: Result comparison mode (side-by-side molecule comparison)
- **ADV-03**: Persistent workspace (save and resume sessions)
- **ADV-04**: Export results to PDF/CSV
- **ADV-05**: Keyboard shortcuts for common operations beyond Cmd+K

## Out of Scope

| Feature | Reason |
|---------|--------|
| Backend API changes | This is frontend-only; backend is stable |
| New API endpoints | Not adding new backend functionality |
| Mobile native app | Web-first, responsive is sufficient for now |
| Complete visual rebrand | Evolving current design, not replacing it |
| Real-time collaboration | Complexity too high for this milestone |
| User accounts/auth | Not part of current architecture |
| Internationalization (i18n) | Not requested, can add later |

## Traceability

Which phases cover which requirements. Updated during roadmap creation.

| Requirement | Phase | Status |
|-------------|-------|--------|
| BUILD-01 | Phase 1 | Complete |
| BUILD-02 | Phase 1 | Complete |
| BUILD-03 | Phase 1 | Complete |
| BUILD-04 | Phase 1 | Complete |
| BUILD-05 | Phase 1 | Complete |
| BUILD-06 | Phase 1 | Complete |
| FRAME-01 | Phase 2 | Complete |
| FRAME-02 | Phase 3 | Complete |
| FRAME-03 | Phase 2 | Complete |
| FRAME-04 | Phase 3 | Complete |
| FRAME-05 | Phase 2 | Complete |
| COMP-01 | Phase 4 | Complete |
| COMP-02 | Phase 4 | Complete |
| COMP-03 | Phase 4 | Complete |
| COMP-04 | Phase 4 | Complete |
| COMP-05 | Phase 4 | Complete |
| COMP-06 | Phase 4 | Complete |
| COMP-07 | Phase 4 | Complete |
| THEME-01 | Phase 4 | Complete |
| THEME-02 | Phase 4 | Complete |
| THEME-03 | Phase 4 | Complete |
| THEME-04 | Phase 4 | Complete |
| THEME-05 | Phase 4 | Complete |
| VISUAL-01 | Phase 4.1 | Complete |
| VISUAL-02 | Phase 4.1 | Complete |
| VISUAL-03 | Phase 4.1 | Complete |
| VISUAL-04 | Phase 4.1 | Complete |
| VISUAL-05 | Phase 4.1 | Complete |
| VISUAL-06 | Phase 4.1 | Complete |
| LOAD-01 | Phase 5 | Complete |
| LOAD-02 | Phase 5 | Complete |
| LOAD-03 | Phase 5 | Complete |
| LOAD-04 | Phase 5 | Complete |
| ANIM-01 | Phase 6 | Complete |
| ANIM-02 | Phase 6 | Complete |
| ANIM-03 | Phase 6 | Complete |
| ANIM-04 | Phase 6 | Complete |
| ANIM-05 | Phase 6 | Complete |
| ANIM-06 | Phase 6 | Complete |
| POWER-01 | Phase 6 | Pending |
| POWER-02 | Phase 6 | Pending |
| POWER-03 | Phase 6 | Pending |
| UX-01 | Phase 5 | Complete |
| UX-02 | Phase 5 | Complete |
| UX-03 | Phase 5 | Complete |
| UX-04 | Phase 5 | Complete |

**Coverage:**
- v1 requirements: 46 total
- Mapped to phases: 46
- Unmapped: 0

---
*Requirements defined: 2026-03-12*
*Last updated: 2026-03-12 after roadmap creation*
