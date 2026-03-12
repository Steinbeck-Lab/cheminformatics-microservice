# Phase 1: Vite Migration - Context

**Gathered:** 2026-03-12
**Status:** Ready for planning

<domain>
## Phase Boundary

Replace deprecated Create React App (react-scripts) with Vite. Resolve all security vulnerabilities (zero npm audit alerts). Update all dependencies to latest stable versions. Eliminate all package.json overrides via proper upgrades. Update Docker build for Vite. Migrate environment variables from REACT_APP_* to VITE_* prefix. Set up frontend CI pipeline.

</domain>

<decisions>
## Implementation Decisions

### Testing Framework
- Use Vitest (Vite-native, ESM-first, shares vite.config.js)
- Keep @testing-library/react and @testing-library/jest-dom as DOM testing utilities (swap jest-dom matchers for vitest-compatible version)
- Use jsdom as the test environment (complete browser emulation, matches CRA behavior)
- Configure @vitest/coverage-v8 with coverage thresholds for local/CI enforcement
- Do NOT integrate with Codecov for frontend — no external coverage reporting

### Dependency Cleanup
- Keep @headlessui/react until Phase 4 (shadcn/ui replaces it there)
- If a dependency still pulls vulnerable transitive deps after upgrading, replace the parent dependency — no overrides, no band-aids
- Defer icon library consolidation to Phase 4 (COMP-06) — keep both FontAwesome and react-icons working for now
- Add npm audit to CI as a vulnerability gate
- Upgrade TypeScript from 4.9.5 to 5.x in this phase (needed for Vite compatibility, avoids half-upgraded state)
- Audit dompurify and jszip for actual usage — remove if unused
- Stay with npm as package manager (no switch to pnpm)

### Claude's Discretion — Dependencies
- Icon cleanup timing: remove web-vitals (CRA-specific) now; defer icon library consolidation to Phase 4
- Lockfile integrity: evaluate whether lockfile-lint adds meaningful security value for this project's threat model and decide accordingly
- Path aliases: evaluate codebase import depth and decide whether @ alias for src/ adds enough value to justify the config

### Docker & Node Version
- Upgrade Docker build stage from Node 18-alpine to Node 22 LTS (22-alpine)
- Use Vite default output directory `dist/` (update Dockerfile COPY path from /app/build to /app/dist)
- Keep multi-stage Docker build pattern (Node build stage + nginx:alpine serve stage)
- Keep build-time environment variable injection via --build-arg (same pattern as current, just VITE_API_URL instead of REACT_APP_API_URL)

### ESLint Configuration
- Use ESLint v9+ flat config format (eslint.config.js)
- Preset: eslint-plugin-react + eslint-plugin-react-hooks + eslint-plugin-jsx-a11y (standard React + accessibility rules)
- Add eslint-config-prettier to disable formatting rules that conflict with Prettier
- Add ESLint to CI as a separate check step
- Set up husky + lint-staged for pre-commit hooks (ESLint + Prettier on staged files)

### Vite Plugin Selection
- Use @vitejs/plugin-react-swc (SWC-based, 20x faster transforms — no custom Babel plugins needed)
- Add vite-plugin-checker for in-dev TypeScript and ESLint error overlay in browser
- Skip vite-plugin-svgr (no SVG component imports in codebase — not needed)

### Import Resolution
- Migrate all bare imports to explicit relative paths (standard Node/ESM resolution, no CRA-isms)
- No src/ root alias to maintain CRA-like bare imports

### Dev Server Configuration
- Keep port 3000 (consistent with docker-compose, docs, and team muscle memory)

### Claude's Discretion — Dev Server
- API proxy: evaluate whether Vite server.proxy for /v1/* and /latest/* to localhost:8000 improves the dev experience given existing CORS config, and decide accordingly

### index.html & Public Assets
- Move index.html from public/ to frontend/ root (Vite requirement)
- Replace %PUBLIC_URL%/ with / in all paths
- Add Vite entry script: `<script type="module" src="/src/index.js">`
- Keep all meta tags, Matomo analytics, 3Dmol.js CDN, CSP tag intact

### Claude's Discretion — index.html
- CSP meta tag: decide the best approach for making connect-src environment-aware (dev needs localhost:8000, prod needs the strict list)
- Public assets: audit and clean up if warranted, or keep as-is if all are used

### Environment Files
- Rename all REACT_APP_* vars to VITE_* (REACT_APP_API_URL -> VITE_API_URL)
- Replace all process.env.REACT_APP_* references with import.meta.env.VITE_* in source files (5 files identified)
- GENERATE_SOURCEMAP handled by vite.config.js build.sourcemap setting instead of env var
- Create committed .env with defaults (VITE_API_URL=http://localhost:8000/latest) for local dev
- Keep .env.production with VITE_API_URL=https://api.naturalproducts.net/latest
- Update .env.local.example to use new VITE_ prefix

### Build Optimization
- Use Vite/Rollup default chunk splitting (no manual vendor chunks)
- Enable hidden source maps (build.sourcemap: 'hidden') for error tracking without exposing source
- Add rollup-plugin-visualizer as devDep with `npm run build:analyze` script for bundle size tracking

### CI Workflow
- Create new frontend-test.yml workflow file (separate from backend test.yml)
- Steps: ESLint check, Prettier format check, Vitest tests, npm audit, production build verification
- Trigger on push/PR only when files in frontend/ change (paths filter)

</decisions>

<specifics>
## Specific Ideas

- "Avoid band-aids, shortcuts, or specific quick fixes. This is a production system, so accuracy is the top priority."
- Clean break from CRA patterns — follow Vite conventions wherever possible (dist/ output, import.meta.env, flat config)
- All existing 9 pages must continue to load and function correctly after migration

</specifics>

<code_context>
## Existing Code Insights

### Reusable Assets
- No SVG component imports (ReactComponent pattern) — simplifies migration, no SVGR needed
- CSS imports are straightforward (tailwind.css + animations.css in index.js)
- All images referenced via URLs or public/ directory (no JS image imports)

### Established Patterns
- React Router v6 with createBrowserRouter + RouterProvider (modern pattern, Vite-compatible)
- AppContext provider wrapping RouterProvider (context architecture)
- Service layer in services/api.js centralizes API calls via Axios
- Framer Motion animations used throughout (motion.div, AnimatePresence)
- Tailwind CSS v3 with extensive custom theme config (custom colors, spacing, shadows, fonts)

### Integration Points
- 5 files reference REACT_APP_API_URL: AppContext.js, api.js, MoleculeCard.jsx, HighlightedMoleculeCard.jsx, StandardizeView.jsx
- Dockerfile: ARG REACT_APP_API_URL, COPY /app/build to nginx
- docker-compose: maps frontend to port 3000
- package.json scripts: start/build/test all use react-scripts
- eslintConfig in package.json: extends react-app (CRA-specific)
- tailwind.config.js: module.exports format with require() for plugins
- postcss.config.js: module.exports format

</code_context>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 01-vite-migration*
*Context gathered: 2026-03-12*
