# Phase 2: React 19 + TypeScript 5 - Context

**Gathered:** 2026-03-12
**Status:** Ready for planning

<domain>
## Phase Boundary

Upgrade React 18.2 to React 19.2.4, convert the entire JavaScript codebase to TypeScript (.ts/.tsx), migrate from framer-motion to the motion package with enhanced animations, and verify all 9 pages render correctly with no console errors or type errors.

</domain>

<decisions>
## Implementation Decisions

### TypeScript Adoption
- Convert ALL source files from .js/.jsx to .ts/.tsx (atomic conversion in one commit)
- Add tsconfig.json with proper React 19 + Vite configuration
- Remove the treat-js-as-jsx Vite plugin (no longer needed with .tsx files)
- Convert vite.config.mjs to vite.config.ts (Vite natively supports TS config)
- Update lint-staged globs from .js/.jsx to .ts/.tsx
- Update ESLint config for TypeScript support (future-proof setup)
- Update Prettier config to target .ts/.tsx files
- Update Vitest setup files from .js to .ts

### React 19 Upgrade
- Upgrade to React 19.2.4 (latest stable, includes CVE-2025-55182 patch)
- Upgrade react-dom to 19.2.4 to match
- Remove @headlessui/react (dead dependency — not imported anywhere in source code)
- Upgrade @testing-library/react and all test dependencies for React 19 compatibility
- No forwardRef, defaultProps, or propTypes in codebase — clean migration path
- ReactDOM.createRoot already in use — no changes needed for root rendering

### Motion Package Migration
- Migrate from `framer-motion` (v12.6.2) to `motion` package (latest v12.35.2)
- Update all 15 files: change imports from `framer-motion` to `motion/react`
- This is the industry standard going forward — framer-motion is maintained for backward compatibility only
- Enhance ALL animation areas:
  - Page transitions: smoother, more polished transitions between 9 pages
  - Scroll animations: enhance useScroll/useTransform effects (AboutPage, Footer)
  - Layout animations: improve LayoutGroup transitions on tab-based pages (ChemPage, DepictPage, ConvertPage, ToolsPage)
  - Spring physics: upgrade spring configs for more natural, physics-based motion

### File Extensions & Config
- Atomic conversion: rename all ~50 .js/.jsx files to .ts/.tsx in single operation
- Update all import paths to use new extensions
- Full dev config alignment: lint-staged, ESLint, Prettier, Vitest, build scripts all target .ts/.tsx
- Remove esbuild .js loader config from vite.config (no longer needed)

### Claude's Discretion
- TypeScript strictness level: evaluate codebase complexity and Phase 4 shadcn/ui requirements, pick appropriate strictness
- Typing approach per file: use explicit interfaces for API service layer and complex props, inferred types where sufficient
- React 19 new API adoption: selectively adopt new APIs (use hook, Actions, useOptimistic) where they genuinely improve code quality, don't refactor for the sake of it
- ESLint TypeScript plugin selection: decide between @typescript-eslint rules and simpler approaches based on what's future-proof
- Animation enhancement scope: improve existing animations using Motion 12 features where quality improves without risk

</decisions>

<specifics>
## Specific Ideas

- "Use all up to date LTS software for all" — latest stable versions are non-negotiable
- "Make smooth cool animation and do not break anything — should work super smooth in all browsers"
- "Make everything future proof" — decisions should align with where the ecosystem is heading
- Motion package is the recommended choice per official docs and industry adoption (30.7k GitHub stars, 3.6M weekly npm downloads)
- React 19.2.4 chosen over earlier versions to include the critical security patch (even though RCE only affects Server Components, not SPAs)
- Phase 1 decision carried forward: "Avoid band-aids, shortcuts, or specific quick fixes. This is a production system, so accuracy is the top priority."

</specifics>

<code_context>
## Existing Code Insights

### Reusable Assets
- `useMolecule` hook: custom hook with useState/useCallback/useEffect — will need TypeScript types for molecule data shape
- `useApiCall` hook: generic API call wrapper — good candidate for generic TypeScript typing
- `AppContext`: createContext/useContext pattern — needs typed context value
- `services/api.js`: Axios-based service layer — needs typed API response interfaces

### Established Patterns
- React Router v6 with createBrowserRouter + RouterProvider (modern pattern, works with React 19)
- Framer Motion used in 15 files: motion.div, AnimatePresence, LayoutGroup, useScroll, useTransform, useAnimation, useInView
- Tailwind CSS v3 with custom theme (stays as-is until Phase 3)
- Component structure: pages/, components/{chem,common,convert,depict,ocsr,tools}/, hooks/, services/, context/, utils/
- No forwardRef, defaultProps, propTypes, createRef, or React.Children usage anywhere

### Integration Points
- vite.config.mjs: needs conversion to .ts, removal of treat-js-as-jsx plugin and esbuild .js loader
- eslint.config.mjs: needs TypeScript parser/plugin additions and .ts/.tsx glob updates
- package.json lint-staged: globs need updating from .js/.jsx to .ts/.tsx
- Vitest setup: ./src/__tests__/setup.js needs .ts rename
- index.jsx (entry point): needs .tsx rename and corresponding index.html script src update

</code_context>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 02-react-19-typescript-5*
*Context gathered: 2026-03-12*
