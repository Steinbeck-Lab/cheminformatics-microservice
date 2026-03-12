# Research Summary: Cheminformatics Microservice Frontend Modernization

**Domain:** Frontend stack migration (CRA -> Vite, React 18 -> 19, Tailwind v3 -> v4, + shadcn/ui)
**Researched:** 2026-03-12
**Overall confidence:** MEDIUM-HIGH

## Executive Summary

The Cheminformatics Microservice frontend is a well-structured React 18 SPA with 9 pages, 27+ view components, dark mode toggle, framer-motion animations, and an Axios-based service layer calling a FastAPI backend. The codebase is functional but built on a deprecated foundation: Create React App (react-scripts 5.0.1, unmaintained), Tailwind CSS v3 with JS-based configuration, and React 18 (two major versions behind).

The target stack -- Vite + React 19 + Tailwind v4 + shadcn/ui -- is well-established as of early 2026. React 19.2 (released Oct 2025) is mature and stable. Tailwind CSS v4.2 is current with a CSS-first configuration approach and a dedicated Vite plugin (`@tailwindcss/vite`). The ecosystem has converged: Vite is the default build tool for React projects, Tailwind v4 is the standard utility CSS framework, and shadcn/ui is the dominant Tailwind-native component system.

The migration is straightforward for this codebase because it avoids most React 19 breaking changes: no class components, no string refs, no propTypes, no legacy context patterns, already using `createRoot`. The Tailwind v4 migration is more involved due to 12+ breaking changes (class renames, config format change, syntax changes), but an official automated upgrade tool handles most of it. The riskiest area is the existing CSS variable-based theming in `tailwind.css`, which uses Tailwind v3's `theme()` function and must be rewritten for v4's `@theme` block syntax.

This project does NOT require TypeScript migration (all 55+ files are .js/.jsx), state management libraries (AppContext is sufficient), or React Router upgrades (v6 is current). The scope is tightly focused: modernize the build tool, upgrade frameworks, add a component system, polish dark mode.

## Key Findings

**Stack:** Vite 6.x + React 19.2 + Tailwind v4.2 + shadcn/ui (Radix UI + lucide-react). Replaces CRA + React 18 + Tailwind v3 + hand-rolled components. All security vulnerability overrides in package.json are eliminated by this upgrade.

**Architecture:** Existing page->component->service architecture is sound. Only addition is `components/ui/` directory for shadcn/ui primitives and `lib/utils.js` for the `cn()` class merging utility. No structural changes needed.

**Critical pitfall:** Tailwind v4 class renames are SILENT failures -- old classes produce no CSS with no error. Must run `npx @tailwindcss/upgrade` and visually verify every page. The existing `theme()` function calls in CSS variables will break and need manual rewriting.

## Implications for Roadmap

Based on research, suggested phase structure:

1. **CRA to Vite Migration** - Build tool replacement
   - Addresses: CRA deprecation, security vulnerabilities, slow builds (30s -> 5s)
   - Avoids: Env var prefix errors by systematically replacing all `REACT_APP_*` with `VITE_*`
   - Includes: Dockerfile update (node:20, VITE_*, /app/dist), entry point rename

2. **React 18 to 19 Upgrade** - Framework upgrade
   - Addresses: Access to React 19 APIs (useActionState, ref-as-prop, Context as provider)
   - Avoids: Breaking changes by running automated codemod first
   - Includes: @testing-library/react v16 upgrade

3. **Tailwind v3 to v4 Migration** - Styling framework upgrade
   - Addresses: CSS-first config, Vite plugin integration, modern utility names
   - Avoids: Silent class rename failures by running automated upgrade tool + visual verification
   - Includes: Rewriting tailwind.css theme variables, removing postcss.config.js

4. **shadcn/ui Integration + Dark Mode** - Component system + theming
   - Addresses: Consistent components, accessibility (Radix UI), CSS variable theming
   - Avoids: Theming conflicts by establishing single CSS variable strategy before component replacement
   - Includes: Core components (Button, Card, Input, Select, Tabs, Dialog, Toast, Skeleton)

5. **UI Polish + Differentiators** - Animations, transitions, power-user features
   - Addresses: Premium feel, smooth navigation, Command palette (Cmd+K)
   - Avoids: Scope creep by deferring complex features (result comparison, persistent workspace)

**Phase ordering rationale:**
- Vite must come first: Tailwind v4's `@tailwindcss/vite` plugin requires Vite; shadcn CLI assumes Vite project
- React 19 before Tailwind v4: Independent upgrade, but must be done before shadcn/ui (which assumes React 19)
- Tailwind v4 before shadcn/ui: shadcn/ui's latest component templates use Tailwind v4 CSS variable syntax
- Dark mode with shadcn/ui: shadcn/ui's theming system IS the dark mode solution
- Polish last: Layered on top of working, themed components

**Research flags for phases:**
- Phase 1: Standard CRA-to-Vite migration. Well-documented. LOW research risk.
- Phase 2: Standard React upgrade. Automated codemod available. LOW research risk.
- Phase 3: MEDIUM research risk -- `theme()` function replacement needs manual work; verify upgrade tool handles CSS variable definitions
- Phase 4: MEDIUM research risk -- shadcn/ui CLI behavior with JavaScript (not TypeScript) projects should be tested; verify Tailwind v4 support
- Phase 5: LOW risk -- additive features using established component patterns

## Confidence Assessment

| Area | Confidence | Notes |
|------|------------|-------|
| Stack | HIGH | React 19.2, Tailwind v4.2 verified via official docs. Vite 6.x from training data (verify exact version). |
| Features | MEDIUM-HIGH | Table stakes verified via codebase analysis. Differentiator value based on proven patterns. |
| Architecture | HIGH | Existing structure is sound. Modifications are minimal and well-documented (shadcn/ui conventions). |
| Pitfalls | HIGH | Tailwind v4 and React 19 pitfalls verified from official upgrade guides. CRA-to-Vite pitfalls well-known. |

## Gaps to Address

- Exact Vite version should be verified at implementation time (`npm view vite version`)
- shadcn/ui CLI support for JavaScript (vs TypeScript) projects needs verification
- framer-motion v12 + React 19 compatibility should be verified (likely fine, but untested)
- Icon mapping from react-icons/hi to lucide-react needs a detailed table before migration
- The `theme()` function usage in existing CSS needs a concrete replacement plan for each variable
