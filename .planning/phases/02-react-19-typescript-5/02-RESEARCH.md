# Phase 2: React 19 + TypeScript 5 - Research

**Researched:** 2026-03-12
**Domain:** React 19 upgrade, JavaScript-to-TypeScript conversion, Motion package migration
**Confidence:** HIGH

## Summary

This phase involves three interconnected upgrades: (1) React 18.3.1 to React 19.2.4, (2) converting 61 JavaScript source files to TypeScript (.ts/.tsx), and (3) migrating from `framer-motion` to the `motion` package. The codebase is in excellent shape for all three upgrades -- no forwardRef, defaultProps, propTypes, or legacy context usage exists, ReactDOM.createRoot is already in use, and the framer-motion API is identical in the `motion` package. The `@headlessui/react` v1.7.19 is listed in package.json but is NOT imported anywhere in source code -- it is dead weight that should be removed rather than upgraded.

The TypeScript conversion is the largest task by volume (61 files, ~28,000 lines), but the codebase follows consistent patterns (service layer with JSDoc, hooks with consistent return shapes, context with a single provider) that make typing straightforward. The motion package migration is purely an import path change from `"framer-motion"` to `"motion/react"` across 15 files.

**Primary recommendation:** Execute in three waves: (1) TypeScript infrastructure + file renames, (2) React 19 + dependency upgrades, (3) motion package migration + animation enhancements. Each wave can be verified independently.

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions
- Convert ALL source files from .js/.jsx to .ts/.tsx (atomic conversion in one commit)
- Add tsconfig.json with proper React 19 + Vite configuration
- Remove the treat-js-as-jsx Vite plugin (no longer needed with .tsx files)
- Convert vite.config.mjs to vite.config.ts (Vite natively supports TS config)
- Update lint-staged globs from .js/.jsx to .ts/.tsx
- Update ESLint config for TypeScript support (future-proof setup)
- Update Prettier config to target .ts/.tsx files
- Update Vitest setup files from .js to .ts
- Upgrade to React 19.2.4 (latest stable, includes CVE-2025-55182 patch)
- Upgrade react-dom to 19.2.4 to match
- Upgrade @headlessui/react from v1.7.17 to v2.x for React 19 compatibility
- Upgrade @testing-library/react and all test dependencies for React 19 compatibility
- No forwardRef, defaultProps, or propTypes in codebase -- clean migration path
- ReactDOM.createRoot already in use -- no changes needed for root rendering
- Migrate from `framer-motion` (v12.6.2) to `motion` package (latest v12.35.2)
- Update all 15 files: change imports from `framer-motion` to `motion/react`
- Enhance ALL animation areas (page transitions, scroll, layout, spring physics)
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

### Deferred Ideas (OUT OF SCOPE)
None -- discussion stayed within phase scope
</user_constraints>

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|-----------------|
| FRAME-01 | React upgraded from 18.2 to 19.x (latest stable) | React 19.2.4 identified; no breaking patterns in codebase; createRoot already used; no forwardRef/defaultProps/propTypes; react-router-dom v6 compatible |
| FRAME-03 | TypeScript upgraded from 4.9.5 to 5.x (latest stable) | TypeScript 5.9.3 already in devDependencies; needs tsconfig.json, file renames, type definitions; @types/react@19 and @types/react-dom@19 required |
| FRAME-05 | Framer Motion verified compatible with React 19 | Motion v12.35.2 supports React 18.2+; migration is import path change only (`"framer-motion"` to `"motion/react"`); API is identical |
</phase_requirements>

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| react | 19.2.4 | UI framework | Latest stable with CVE-2025-55182 patch |
| react-dom | 19.2.4 | DOM rendering | Must match react version |
| typescript | ^5.9.3 | Type system | Already in devDependencies, latest stable |
| motion | ^12.35.2 | Animation library | Official successor to framer-motion, same API |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| @types/react | ^19.2.14 | React type definitions | Required for TypeScript + React 19 |
| @types/react-dom | ^19.2.14 | ReactDOM type definitions | Required for TypeScript + React 19 |
| @typescript-eslint/eslint-plugin | ^8.x | ESLint TypeScript rules | For .ts/.tsx linting in ESLint flat config |
| @typescript-eslint/parser | ^8.x | TypeScript parser for ESLint | Required by typescript-eslint plugin |
| @testing-library/react | ^16.3.2 | React test utilities | Already installed, supports React 19 |

### Remove (Dead Dependencies)
| Library | Reason |
|---------|--------|
| @headlessui/react | Not imported in any source file -- dead dependency |
| framer-motion | Replaced by motion package |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| @typescript-eslint full config | typescript-eslint (unified package) | The `typescript-eslint` npm package is the modern unified approach; it bundles parser + plugin. Use this instead of separate installs. |
| Strict TypeScript | Relaxed TypeScript | Strict is better for Phase 4 (shadcn/ui expects strict). Start strict, suppress individual cases with `// @ts-expect-error` where needed. |

**Installation:**
```bash
cd frontend
npm install react@^19.2.4 react-dom@^19.2.4 motion@^12.35.2
npm install -D @types/react@^19 @types/react-dom@^19 typescript-eslint@^8
npm uninstall framer-motion @headlessui/react
```

## Architecture Patterns

### Recommended TypeScript Project Structure
```
frontend/
├── tsconfig.json              # Root config (references app + node)
├── tsconfig.app.json          # App-specific config (src files)
├── tsconfig.node.json         # Node config (vite.config.ts)
├── vite.config.ts             # Converted from .mjs
├── eslint.config.mjs          # Updated with TypeScript parser/plugin
├── src/
│   ├── index.tsx              # Entry point (renamed from .jsx)
│   ├── App.tsx                # Root component (renamed from .js)
│   ├── types/                 # Shared type definitions
│   │   ├── api.ts             # API response types
│   │   ├── molecule.ts        # Molecule data types
│   │   └── global.d.ts        # Global declarations (3Dmol, etc.)
│   ├── context/
│   │   └── AppContext.tsx      # Typed context with explicit interface
│   ├── hooks/
│   │   ├── useApiCall.ts      # Generic typed hook
│   │   └── useMolecule.ts     # Typed molecule hook
│   ├── services/
│   │   ├── api.ts             # Typed Axios instance
│   │   ├── chemService.ts     # Typed service functions
│   │   ├── convertService.ts
│   │   ├── depictService.ts
│   │   ├── ocsrService.ts
│   │   └── toolsService.ts
│   ├── pages/                 # All .tsx
│   ├── components/            # All .tsx
│   └── __tests__/
│       ├── setup.ts           # Renamed from .js
│       └── App.test.tsx       # Renamed from .jsx
```

### Pattern 1: TypeScript Config for Vite + React 19

**What:** Three-file tsconfig structure matching Vite's recommended setup.
**When to use:** All Vite + React + TypeScript projects.

**tsconfig.json (root):**
```json
{
  "references": [
    { "path": "./tsconfig.app.json" },
    { "path": "./tsconfig.node.json" }
  ],
  "files": []
}
```

**tsconfig.app.json:**
```json
{
  "compilerOptions": {
    "target": "ES2020",
    "useDefineForClassFields": true,
    "lib": ["ES2020", "DOM", "DOM.Iterable"],
    "module": "ESNext",
    "skipLibCheck": true,
    "moduleResolution": "bundler",
    "allowImportingTsExtensions": true,
    "isolatedModules": true,
    "moduleDetection": "force",
    "noEmit": true,
    "jsx": "react-jsx",
    "strict": true,
    "noUnusedLocals": true,
    "noUnusedParameters": true,
    "noFallthroughCasesInSwitch": true,
    "noUncheckedIndexedAccess": true,
    "forceConsistentCasingInFileNames": true,
    "resolveJsonModule": true,
    "esModuleInterop": true
  },
  "include": ["src"]
}
```

**tsconfig.node.json:**
```json
{
  "compilerOptions": {
    "target": "ES2022",
    "lib": ["ES2023"],
    "module": "ESNext",
    "skipLibCheck": true,
    "moduleResolution": "bundler",
    "allowImportingTsExtensions": true,
    "isolatedModules": true,
    "moduleDetection": "force",
    "noEmit": true,
    "strict": true,
    "noUnusedLocals": true,
    "noUnusedParameters": true,
    "noFallthroughCasesInSwitch": true
  },
  "include": ["vite.config.ts"]
}
```

### Pattern 2: Typed Context Pattern

**What:** Strongly typed React context with explicit interface.
**When to use:** AppContext.tsx conversion.

```typescript
// Source: Derived from existing AppContext.js patterns
interface ApiConfig {
  baseUrl: string;
  timeout: number;
}

interface RecentMolecule {
  smiles: string;
  name?: string;
  timestamp: string;
  _savedAt: number;
}

interface AppContextValue {
  isDarkMode: boolean;
  toggleDarkMode: () => void;
  recentMolecules: RecentMolecule[];
  addRecentMolecule: (molecule: Omit<RecentMolecule, "_savedAt">) => void;
  clearRecentMolecules: () => void;
  apiConfig: ApiConfig;
  updateApiConfig: (newConfig: Partial<ApiConfig>) => void;
  isLoading: boolean;
  setIsLoading: React.Dispatch<React.SetStateAction<boolean>>;
  globalError: string | null;
  setGlobalError: React.Dispatch<React.SetStateAction<string | null>>;
}

const AppContext = createContext<AppContextValue | undefined>(undefined);

export const useAppContext = (): AppContextValue => {
  const context = useContext(AppContext);
  if (!context) {
    throw new Error("useAppContext must be used within an AppProvider");
  }
  return context;
};
```

### Pattern 3: Typed Service Layer

**What:** TypeScript interfaces for API services.
**When to use:** All service file conversions.

```typescript
// src/services/chemService.ts
import api from "./api";

const CHEM_URL = "/chem";

export interface StructureErrorResult {
  messages: string;
  standardized?: {
    smi: string;
  };
}

export const checkStructureErrors = async (
  smiles: string,
  fix = false
): Promise<StructureErrorResult> => {
  const response = await api.get<StructureErrorResult>(`${CHEM_URL}/errors`, {
    params: { smiles, fix },
  });
  return response.data;
};
```

### Pattern 4: Generic Typed Hook

**What:** TypeScript generics for the useApiCall hook.
**When to use:** useApiCall.ts conversion.

```typescript
// src/hooks/useApiCall.ts
interface UseApiCallOptions<T> {
  immediate?: boolean;
  initialArgs?: unknown[];
  onSuccess?: (data: T) => void;
  onError?: (error: Error) => void;
}

interface UseApiCallReturn<T> {
  call: (...args: unknown[]) => Promise<T>;
  data: T | null;
  isLoading: boolean;
  error: string | null;
  lastCallTime: Date | null;
  reset: () => void;
}

function useApiCall<T>(
  apiFunction: (...args: unknown[]) => Promise<T>,
  options: UseApiCallOptions<T> = {}
): UseApiCallReturn<T> {
  // ... implementation
}
```

### Pattern 5: Motion Import Migration

**What:** Import path change from framer-motion to motion.
**When to use:** All 15 files using framer-motion.

```typescript
// BEFORE (framer-motion)
import { motion, AnimatePresence, LayoutGroup } from "framer-motion";

// AFTER (motion package)
import { motion, AnimatePresence, LayoutGroup } from "motion/react";
```

The API is identical -- no prop or component changes needed.

### Pattern 6: Component Props Typing

**What:** Props interface pattern for React components.
**When to use:** All component .tsx files.

```typescript
// For components with simple props, use inline or interface
interface NavigationProps {
  isMobile?: boolean;
  closeMenu?: () => void;
  isAnimated?: boolean;
  menuItemVariants?: Record<string, unknown>;
}

const Navigation: React.FC<NavigationProps> = ({
  isMobile = false,
  closeMenu = () => {},
}) => {
  // ...
};
```

### Anti-Patterns to Avoid
- **`any` as escape hatch:** Use `unknown` instead of `any` where types are unclear. If truly needed, use `// @ts-expect-error` with explanation.
- **Over-typing:** Don't type what TypeScript can infer. Function return types, simple variable assignments, and map/filter operations rarely need explicit types.
- **Barrel exports with re-export all:** Avoid `export * from './types'` for type files -- use named exports for tree-shaking.
- **Typing third-party globals as `any`:** The 3Dmol library loaded via CDN script tag needs a proper global declaration in `global.d.ts`.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| TypeScript ESLint config | Custom parser/plugin setup | `typescript-eslint` unified package | Single package bundles parser + plugin + configs; use `tseslint.configs.recommended` |
| tsconfig from scratch | Custom compiler options | Vite's three-file reference pattern | Matches Vite's build expectations; separates app vs. node environments |
| React 19 type codemods | Manual type fixes | `npx types-react-codemod@latest preset-19 ./src` | Automates JSX namespace changes and deprecated type updates |
| API response types | Manual interface writing | Derive types from existing JSDoc annotations | All service functions have JSDoc; use those as type spec |
| Motion migration | File-by-file manual edits | Find-and-replace `from "framer-motion"` with `from "motion/react"` | Import path is the only change; API is identical |

**Key insight:** This phase is primarily mechanical transformation (rename, re-type, re-import) rather than architectural change. Automation and batch operations are critical to avoid introducing bugs through inconsistency.

## Common Pitfalls

### Pitfall 1: index.html Script Source Not Updated
**What goes wrong:** After renaming `index.jsx` to `index.tsx`, the app fails to load because `index.html` still references `/src/index.jsx`.
**Why it happens:** The HTML file is outside the rename batch and easy to forget.
**How to avoid:** Update `<script type="module" src="/src/index.tsx"></script>` in `index.html` as the very first step after file renames.
**Warning signs:** Blank white page in browser, 404 in network tab for index.jsx.

### Pitfall 2: TypeScript Strict Mode Overwhelm
**What goes wrong:** Enabling `strict: true` on a 28,000-line JS codebase produces hundreds of errors, blocking progress.
**Why it happens:** Strict mode enables 7+ individual checks simultaneously (strictNullChecks, noImplicitAny, etc.).
**How to avoid:** Start with strict mode ON but address errors file by file. The codebase already uses consistent patterns (JSDoc, default parameters) that minimize surprises. Focus on the most impactful strict checks: `strictNullChecks` and `noImplicitAny`.
**Warning signs:** More than 50 type errors on first compile -- if this happens, consider temporarily relaxing `noUnusedLocals` and `noUnusedParameters` (they can be warnings in ESLint instead).

### Pitfall 3: React Router Future Flags Conflict
**What goes wrong:** After upgrading to React 19, the `v7_startTransition` and `v7_relativeSplatPath` future flags on createBrowserRouter may cause warnings or behavioral changes.
**Why it happens:** React 19 has built-in support for transitions via `startTransition`; react-router-dom v6 future flags may conflict.
**How to avoid:** Keep react-router-dom at v6.30.3 (already installed). The `v7_startTransition` flag is designed for forward compatibility and should work with React 19. If warnings appear, the flags can be safely removed since React 19 natively supports transitions.
**Warning signs:** Console warnings about `startTransition` being called outside a transition.

### Pitfall 4: @types/react Version Mismatch
**What goes wrong:** Installing `@types/react@^18` with React 19 causes type errors, particularly around `ref` handling and JSX namespace changes.
**Why it happens:** React 19 changed the JSX type namespace and ref behavior; @types/react@18 doesn't know about these changes.
**How to avoid:** Install `@types/react@^19` and `@types/react-dom@^19`. Run `npx types-react-codemod@latest preset-19 ./src` if any existing type patterns are affected.
**Warning signs:** Type errors mentioning `JSX.Element`, `React.FC` children prop, or `ref` type issues.

### Pitfall 5: Vitest Setup File Path Not Updated
**What goes wrong:** Tests fail to run because vite.config.ts still references `./src/__tests__/setup.js` but the file has been renamed to `.ts`.
**Why it happens:** The `setupFiles` path in vite.config test configuration is a string, not validated at compile time.
**How to avoid:** Update `setupFiles: "./src/__tests__/setup.ts"` in vite.config.ts during the config conversion step.
**Warning signs:** "Cannot find module" error when running `npm test`.

### Pitfall 6: DOMPurify and Third-Party Types
**What goes wrong:** Libraries like `dompurify`, `jszip`, and `react-dropzone` may not have type definitions or have outdated ones.
**Why it happens:** Some npm packages include their own types; others need `@types/` packages.
**How to avoid:** Check each dependency: `dompurify` has `@types/dompurify`; `jszip` ships its own types; `react-dropzone` ships its own types; `axios` ships its own types. Only install @types packages where needed.
**Warning signs:** "Could not find a declaration file for module" TypeScript errors.

### Pitfall 7: 3Dmol Global Variable Not Declared
**What goes wrong:** TypeScript errors on `$3Dmol` usage because it's loaded via CDN `<script>` tag, not imported.
**Why it happens:** Globals from CDN scripts have no type information.
**How to avoid:** Create `src/types/global.d.ts` with `declare const $3Dmol: any;` (or a more specific type if available). Also declare `import.meta.env` types for Vite environment variables.
**Warning signs:** "Cannot find name '$3Dmol'" TypeScript error.

### Pitfall 8: CSS Module Type Declarations Missing
**What goes wrong:** TypeScript errors when importing `.css` files.
**Why it happens:** TypeScript doesn't know how to handle CSS imports by default.
**How to avoid:** Vite's `vite/client` types (referenced via `/// <reference types="vite/client" />` in a `.d.ts` file or included in tsconfig) handle `.css`, `.svg`, and other asset imports. Ensure `vite/client` is in the types.
**Warning signs:** "Cannot find module './styles/tailwind.css'" error.

## Code Examples

### vite.config.ts (Converted)
```typescript
/// <reference types="vitest/config" />
import { defineConfig } from "vite";
import react from "@vitejs/plugin-react-swc";

export default defineConfig({
  plugins: [
    react({
      jsxRuntime: "automatic",
    }),
    // treat-js-as-jsx plugin REMOVED -- no longer needed with .tsx files
  ],
  server: {
    port: 3000,
    proxy: {
      "/v1": {
        target: "http://localhost:8000",
        changeOrigin: true,
      },
      "/latest": {
        target: "http://localhost:8000",
        changeOrigin: true,
      },
    },
  },
  // optimizeDeps.esbuildOptions.loader REMOVED -- no longer needed
  build: {
    sourcemap: "hidden",
    outDir: "dist",
  },
  test: {
    globals: true,
    environment: "jsdom",
    setupFiles: "./src/__tests__/setup.ts",
    css: true,
    coverage: {
      provider: "v8",
      reporter: ["text", "lcov"],
    },
  },
});
```

### ESLint Config with TypeScript (eslint.config.mjs)
```javascript
import js from "@eslint/js";
import globals from "globals";
import reactPlugin from "eslint-plugin-react";
import reactHooks from "eslint-plugin-react-hooks";
import jsxA11y from "eslint-plugin-jsx-a11y";
import prettierConfig from "eslint-config-prettier";
import tseslint from "typescript-eslint";

export default [
  js.configs.recommended,
  ...tseslint.configs.recommended,
  {
    files: ["**/*.{ts,tsx}"],
    plugins: {
      react: reactPlugin,
      "react-hooks": reactHooks,
      "jsx-a11y": jsxA11y,
    },
    languageOptions: {
      ecmaVersion: "latest",
      sourceType: "module",
      globals: {
        ...globals.browser,
        ...globals.es2021,
      },
      parserOptions: {
        ecmaFeatures: {
          jsx: true,
        },
      },
    },
    settings: {
      react: {
        version: "detect",
      },
    },
    rules: {
      ...reactPlugin.configs.recommended.rules,
      ...reactPlugin.configs["jsx-runtime"].rules,
      "react-hooks/rules-of-hooks": "error",
      "react-hooks/exhaustive-deps": "warn",
      ...jsxA11y.flatConfigs.recommended.rules,
      "react/prop-types": "off",
      "react/no-unescaped-entities": "warn",
      "react/no-unknown-property": ["error", { ignore: ["jsx", "global"] }],
      "@typescript-eslint/no-unused-vars": [
        "warn",
        { argsIgnorePattern: "^_", varsIgnorePattern: "^_" },
      ],
      "no-unused-vars": "off", // Replaced by @typescript-eslint version
      "no-case-declarations": "warn",
      "jsx-a11y/click-events-have-key-events": "warn",
      "jsx-a11y/no-static-element-interactions": "warn",
      "jsx-a11y/no-noninteractive-element-interactions": "warn",
      "jsx-a11y/interactive-supports-focus": "warn",
      "jsx-a11y/label-has-associated-control": "warn",
      "jsx-a11y/no-autofocus": "warn",
    },
  },
  prettierConfig,
  {
    ignores: ["dist/", "build/", "node_modules/", "public/standalone/"],
  },
];
```

### Global Type Declarations (src/types/global.d.ts)
```typescript
/// <reference types="vite/client" />

// 3Dmol.js loaded via CDN script tag in index.html
declare const $3Dmol: {
  createViewer: (
    element: HTMLElement,
    config?: Record<string, unknown>
  ) => unknown;
  [key: string]: unknown;
};

// Vite environment variables
interface ImportMetaEnv {
  readonly VITE_API_URL: string;
}

interface ImportMeta {
  readonly env: ImportMetaEnv;
}
```

### Typed Axios Instance (src/services/api.ts)
```typescript
import axios, { AxiosInstance, InternalAxiosRequestConfig, AxiosResponse } from "axios";

const API_URL: string =
  import.meta.env.VITE_API_URL || "https://dev.api.naturalproducts.net/latest";

const api: AxiosInstance = axios.create({
  baseURL: API_URL,
  headers: {
    "Content-Type": "application/json",
  },
  timeout: 30000,
});

api.interceptors.request.use(
  (config: InternalAxiosRequestConfig) => config,
  (error: unknown) => Promise.reject(error)
);

api.interceptors.response.use(
  (response: AxiosResponse) => response,
  (error: unknown) => {
    if (import.meta.env.DEV && error instanceof Error) {
      console.error("API Error:", error.message);
    }
    return Promise.reject(error);
  }
);

export default api;
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| `framer-motion` package | `motion` package | Late 2024 | Import path change only; `framer-motion` still maintained but `motion` is the future |
| `@types/react@18` | `@types/react@19` | December 2024 | JSX namespace moved into module; ref as prop; children no longer implicit in FC |
| `React.FC<Props>` | `function Component(props: Props)` | React 19 convention | FC no longer includes implicit children; plain functions preferred |
| forwardRef wrapper | ref as regular prop | React 19 | forwardRef still works but is deprecated; ref is now a regular prop |
| defaultProps on functions | ES6 default parameters | React 19 | defaultProps on function components deprecated |
| Separate eslint parser + plugin | `typescript-eslint` unified package | 2024 | Single import: `import tseslint from "typescript-eslint"` |

**Deprecated/outdated:**
- `React.FC` with implicit children -- children must be explicitly declared in React 19 types
- `framer-motion` import path -- use `motion/react` for new projects
- `.eslintrc` format -- ESLint 9 uses flat config exclusively (already migrated in Phase 1)

## TypeScript Strictness Recommendation (Claude's Discretion)

**Recommendation: Enable `strict: true`.**

Rationale:
1. The codebase is well-structured with consistent patterns -- strict mode should produce manageable error counts
2. Phase 4 (shadcn/ui) expects strict TypeScript -- setting it up now avoids a second migration
3. The service layer already has JSDoc types that map cleanly to TypeScript interfaces
4. React 19 + @types/react@19 work best with strict mode enabled
5. If error count is overwhelming (>100 on first compile), `noUnusedLocals` and `noUnusedParameters` can be temporarily set to false in tsconfig and enforced via ESLint warnings instead

## Motion Animation Enhancement Recommendations (Claude's Discretion)

Existing animations are already sophisticated. Enhancements should be targeted, not wholesale rewrites:

1. **Spring physics upgrade:** Replace hardcoded `duration` transitions with spring-based transitions where physics-based motion improves feel (e.g., sidebar pill indicators already use springs well -- extend this pattern)
2. **useScroll/useTransform on AboutPage:** These are already in use; consider adding `useMotionValueEvent` for scroll-triggered state changes
3. **LayoutGroup animations:** Already well-implemented on ChemPage, DepictPage, ConvertPage, ToolsPage for tab indicators; keep as-is
4. **Page transitions:** Consider adding subtle page-level `AnimatePresence` wrapper in the Layout component for route transitions (currently absent -- pages animate independently)
5. **DO NOT change:** Footer canvas particles, Header mobile menu animations, or LoadingScreen -- these work well and are fragile to modify

## React 19 New API Adoption (Claude's Discretion)

**Recommendation: Minimal adoption in this phase.**

- `use()` hook: NOT recommended for this phase. The codebase uses `useApiCall` custom hook which handles loading/error states well. Introducing `use()` would require rethinking the data fetching pattern -- defer to a future phase.
- `useOptimistic()`: NOT applicable -- no optimistic UI patterns exist in this codebase.
- `useFormStatus()` / `useActionState()`: NOT applicable -- no server actions or form actions patterns exist (this is a client-side SPA).
- `ref` as prop: No action needed since no forwardRef exists in the codebase.

## Dependency Compatibility Matrix

| Package | Current | Target | React 19 Compatible | Notes |
|---------|---------|--------|---------------------|-------|
| react | 18.3.1 | 19.2.4 | Yes (is React 19) | |
| react-dom | 18.3.1 | 19.2.4 | Yes (is React 19) | |
| react-router-dom | 6.30.3 | 6.30.3 (keep) | Yes | Designed to bridge React 18/19 |
| framer-motion | 12.36.0 | REMOVE | N/A | Replace with motion package |
| motion | N/A | ^12.35.2 | Yes (supports 18.2+) | |
| @headlessui/react | 1.7.19 | REMOVE | N/A | Not imported anywhere in source |
| axios | 1.13.5 | 1.13.5 (keep) | N/A (no React dep) | Ships own types |
| dompurify | 3.3.1 | 3.3.1 (keep) | N/A (no React dep) | Needs @types/dompurify |
| jszip | 3.10.1 | 3.10.1 (keep) | N/A (no React dep) | Ships own types |
| react-dropzone | 14.2.3 | 14.2.3 (keep) | Yes | Ships own types |
| react-icons | 4.12.0 | 4.12.0 (keep) | Yes | Stateless SVG components |
| @fortawesome/* | 6.7.2/0.2.2 | Keep | Yes | |
| @testing-library/react | 16.3.2 | 16.3.2 (keep) | Yes | Already supports React 19 |
| @testing-library/jest-dom | 6.9.1 | 6.9.1 (keep) | N/A | |
| @testing-library/user-event | 14.6.1 | 14.6.1 (keep) | Yes | |
| @vitejs/plugin-react-swc | 4.3.0 | 4.3.0 (keep) | N/A | Vite plugin, not React-dependent |
| vite | 6.4.1 | 6.4.1 (keep) | N/A | |
| vitest | 4.1.0 | 4.1.0 (keep) | N/A | |
| eslint | 9.39.4 | 9.39.4 (keep) | N/A | |
| typescript | 5.9.3 | 5.9.3 (keep) | N/A | Already installed |

## File Rename Inventory

Total files to rename: 61

**By category:**
- Pages (9): All .jsx -> .tsx
- Components/common (7): All .jsx -> .tsx
- Components/chem (14): All .jsx -> .tsx
- Components/convert (3): All .jsx -> .tsx
- Components/depict (6): All .jsx -> .tsx
- Components/ocsr (1): .jsx -> .tsx
- Components/tools (4): All .jsx -> .tsx
- Services (6): All .js -> .ts
- Hooks (2): All .js -> .ts
- Context (1): .js -> .tsx
- Utils (2): .js -> .ts
- App (1): .js -> .tsx
- Test files (2): .js -> .ts and .jsx -> .tsx

**Critical rename notes:**
- `src/index.jsx` -> `src/index.tsx` -- also update `index.html` script src
- `src/App.js` -> `src/App.tsx` -- also update import in index.tsx
- `src/__tests__/setup.js` -> `src/__tests__/setup.ts` -- also update vite.config.ts setupFiles

## Open Questions

1. **tailwind.config.js content scanning**
   - What we know: `tailwind.config.js` content array already includes `"./src/**/*.{js,jsx,ts,tsx}"` -- it will scan .ts/.tsx files correctly
   - What's unclear: Nothing -- this is already handled
   - Recommendation: No changes needed to tailwind.config.js

2. **Potential test behavior changes with React 19**
   - What we know: @testing-library/react v16.3.2 supports React 19, but React 19 changed Suspense behavior
   - What's unclear: The existing smoke test is minimal (one test), so risk is low; but if more tests are added, `act()` wrapping requirements may be stricter
   - Recommendation: Run existing test after upgrade; if it passes, no changes needed. If Suspense-related issues appear, wrap renders in `await act(async () => { ... })`

3. **postcss.config.js format**
   - What we know: Currently `module.exports = { ... }` (CJS format), which works because there is no `"type": "module"` in package.json
   - What's unclear: Whether converting to .ts is worth it given Phase 3 (Tailwind v4) will likely change the PostCSS setup
   - Recommendation: Leave postcss.config.js as-is -- it's CJS by design and will be addressed in Phase 3

## Validation Architecture

### Test Framework
| Property | Value |
|----------|-------|
| Framework | Vitest 4.1.0 |
| Config file | `vite.config.ts` (test section) |
| Quick run command | `cd frontend && npm test` |
| Full suite command | `cd frontend && npm run test:coverage` |

### Phase Requirements -> Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| FRAME-01 | React 19.x installed, app renders all 9 pages | smoke | `cd frontend && npm test` | Yes (App.test.tsx) |
| FRAME-01 | No console errors/warnings | manual | Browser DevTools console check | N/A |
| FRAME-03 | TypeScript compiles with no errors | build | `cd frontend && npx tsc --noEmit` | N/A (tsconfig creates this) |
| FRAME-03 | ESLint passes on .ts/.tsx | lint | `cd frontend && npm run lint` | N/A |
| FRAME-05 | Motion animations work correctly | manual | Visual verification of page transitions, layout animations | N/A |
| FRAME-05 | Build succeeds with motion package | build | `cd frontend && npm run build` | N/A |

### Sampling Rate
- **Per task commit:** `cd frontend && npx tsc --noEmit && npm test && npm run build`
- **Per wave merge:** `cd frontend && npx tsc --noEmit && npm test && npm run build && npm run lint`
- **Phase gate:** Full suite green + visual verification of all 9 pages + DevTools console clean

### Wave 0 Gaps
- [ ] `tsconfig.json` + `tsconfig.app.json` + `tsconfig.node.json` -- must be created before any TS compilation
- [ ] `src/types/global.d.ts` -- must exist for 3Dmol, Vite env vars, CSS module imports
- [ ] `@types/react@^19` + `@types/react-dom@^19` -- must be installed before TS compilation
- [ ] Update `package.json` scripts: lint/format globs from `.js/.jsx` to `.ts/.tsx`

## Sources

### Primary (HIGH confidence)
- React 19 Upgrade Guide: https://react.dev/blog/2024/04/25/react-19-upgrade-guide
- React v19 Release: https://react.dev/blog/2024/12/05/react-19
- Motion Upgrade Guide: https://motion.dev/docs/react-upgrade-guide
- Motion for React docs: https://motion.dev/docs/react
- typescript-eslint Getting Started: https://typescript-eslint.io/getting-started/
- Vite TypeScript support: https://vite.dev/guide/features.html#typescript

### Secondary (MEDIUM confidence)
- @headlessui/react v2 discussion: https://github.com/tailwindlabs/headlessui/discussions/3354 -- confirmed React 19 compatibility in v2.2.1+, but moot since package is unused
- @testing-library/react releases: https://github.com/testing-library/react-testing-library/releases -- v16 supports React 19
- motion npm: https://www.npmjs.com/package/motion -- v12.35.2 latest, supports React 18.2+
- @types/react npm: https://www.npmjs.com/package/@types/react -- v19.2.14 latest for React 19

### Tertiary (LOW confidence)
- None -- all critical findings verified with primary sources

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - Versions verified via npm registry and official docs; compatibility matrix cross-referenced
- Architecture: HIGH - TypeScript config patterns from Vite official docs; code examples derived from actual codebase analysis
- Pitfalls: HIGH - Each pitfall identified from direct codebase inspection (e.g., index.html script src, setup.js path, 3Dmol global)
- Motion migration: HIGH - Official docs confirm identical API with import path change only
- React 19 compatibility: HIGH - No breaking patterns (forwardRef, defaultProps, propTypes) found in codebase

**Research date:** 2026-03-12
**Valid until:** 2026-04-12 (stable ecosystem, no imminent breaking changes expected)
