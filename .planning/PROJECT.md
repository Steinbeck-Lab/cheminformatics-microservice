# Cheminformatics Microservice — Frontend Overhaul

## What This Is

A comprehensive frontend modernization of the Cheminformatics Microservice UI. The existing React + Tailwind frontend serves as a web interface for cheminformatics services (structure analysis, format conversion, 2D/3D depiction, descriptor generation, OCSR). This project upgrades the entire frontend stack, resolves all security vulnerabilities, and polishes the UI/UX while preserving current functionality.

## Core Value

The frontend must feel noticeably modern, fast, and polished — while never breaking existing functionality. Zero security alerts, production-grade quality.

## Requirements

### Validated

- ✓ React SPA with client-side routing — existing
- ✓ Pages: Home, Chem, Convert, Depict, OCSR, Tools, About, Privacy, Terms — existing
- ✓ Component architecture per feature domain (chem, convert, depict, ocsr, tools, common) — existing
- ✓ API service layer communicating with FastAPI backend — existing
- ✓ Tailwind CSS styling — existing
- ✓ Framer Motion animations — existing
- ✓ File upload (react-dropzone) — existing
- ✓ Axios HTTP client — existing
- ✓ Prettier formatting — existing

### Active

- [ ] Migrate from Create React App (react-scripts) to Vite
- [ ] Upgrade React 18 → React 19
- [ ] Upgrade Tailwind CSS v3 → v4
- [ ] Integrate shadcn/ui component library
- [ ] Add dark mode with user toggle (light/dark)
- [ ] Fix all dependabot security vulnerabilities
- [ ] Update all dependencies to latest stable versions
- [ ] Modernize layout and visual design (evolve current look, don't reinvent)
- [ ] Smooth animations and micro-interactions throughout
- [ ] Improve navigation and UX workflows
- [ ] Responsive design improvements
- [ ] Faster load and interaction times (code splitting, lazy loading)
- [ ] Production-grade build optimization

### Out of Scope

- Backend API changes — this is frontend-only
- Adding new API endpoints or features
- Changing the backend Python code
- Mobile native app
- Complete visual rebrand — evolving the current design, not replacing it

## Context

- **Current stack:** React 18.2, react-scripts 5.0.1 (CRA, deprecated), Tailwind 3.3.5, TypeScript 4.9.5, framer-motion 12.6.2
- **Current structure:** 9 pages, 6 component directories (chem, common, convert, depict, ocsr, tools), service layer for API calls
- **CRA is deprecated:** react-scripts is no longer maintained, Vite is the recommended replacement
- **Existing overrides in package.json:** nth-check, postcss, svgo, underscore, serialize-javascript — these are security patches that should be resolved by proper upgrades
- **Backend is stable:** FastAPI backend with RDKit/CDK/OpenBabel — not being touched
- **The frontend communicates via `REACT_APP_API_URL` env var** — will need to change to `VITE_API_URL` after Vite migration
- **Docker builds exist** for frontend — Dockerfile.frontend will need updating for Vite

## Constraints

- **No breaking changes:** All existing pages and functionality must continue working after each phase
- **Backend untouched:** Only frontend code changes
- **Env var migration:** CRA's `REACT_APP_*` prefix must migrate to Vite's `VITE_*` prefix
- **Docker compatibility:** Frontend Docker build must continue working
- **CI compatibility:** GitHub Actions test workflow must continue passing
- **Design evolution:** Keep the current design language, make it better — not a complete redesign

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| CRA → Vite | CRA is deprecated, Vite is faster and actively maintained | — Pending |
| React 18 → 19 | Latest stable with performance improvements and new hooks | — Pending |
| Tailwind v3 → v4 | CSS-first config, better performance, modern features | — Pending |
| Add shadcn/ui | Tailwind-native, customizable, consistent component system | — Pending |
| Dark mode with toggle | User requested, modern expectation | — Pending |
| Evolve design, don't replace | User wants improvement, not reinvention — don't break what works | — Pending |

---
*Last updated: 2026-03-12 after initialization*
