# Phase 1: Vite Migration - Research

**Researched:** 2026-03-12
**Domain:** Build toolchain migration (CRA to Vite), dependency cleanup, ESLint/testing setup, Docker
**Confidence:** HIGH

## Summary

This phase replaces Create React App (react-scripts 5.0.1) with Vite as the build tool for the frontend. CRA is officially sunset and the source of all 16 npm audit vulnerabilities (9 low, 1 moderate, 6 high). Removing react-scripts eliminates the entire vulnerability chain. The project has zero existing frontend tests, a straightforward CSS/JS structure with no SVG component imports, and 5 files referencing `REACT_APP_*` environment variables plus 3 files using `process.env.PUBLIC_URL`.

The migration is well-understood territory with no exotic patterns in this codebase. The main complexity lies in: (1) correctly wiring environment variables (8 files need changes), (2) moving index.html and updating asset paths, (3) setting up ESLint v9 flat config from scratch (replacing CRA's built-in config), (4) establishing Vitest as the test framework, and (5) updating the Docker build.

**Primary recommendation:** Use Vite 6.4.1 (latest 6.x) rather than Vite 7 or 8. Vite 8 just released with a Rolldown migration and plugin ecosystem compatibility is not yet confirmed. Vite 6 is battle-tested and all plugins explicitly declare support.

<user_constraints>
## User Constraints (from CONTEXT.md)

### Locked Decisions
- Testing: Vitest with @testing-library/react, @testing-library/jest-dom (vitest-compatible), jsdom environment, @vitest/coverage-v8. No frontend Codecov integration.
- Dependency cleanup: Keep @headlessui/react until Phase 4. Replace parent deps if transitive vulns persist (no overrides). Defer icon consolidation to Phase 4. Add npm audit to CI. Upgrade TypeScript to 5.x. Audit dompurify and jszip for actual usage. Stay with npm (no pnpm).
- Docker: Upgrade to Node 22-alpine. Output dir dist/. Keep multi-stage build. Build-time env injection via --build-arg with VITE_API_URL.
- ESLint: v9+ flat config (eslint.config.js). Plugins: eslint-plugin-react + eslint-plugin-react-hooks + eslint-plugin-jsx-a11y + eslint-config-prettier. ESLint in CI. Husky + lint-staged for pre-commit.
- Vite plugins: @vitejs/plugin-react-swc (SWC-based). vite-plugin-checker for dev overlay. Skip vite-plugin-svgr.
- Import resolution: Explicit relative paths (no src/ alias, no CRA bare imports).
- Dev server: Keep port 3000.
- index.html: Move from public/ to frontend/ root. Replace %PUBLIC_URL%/ with /. Add `<script type="module" src="/src/index.js">`. Keep meta tags, Matomo, 3Dmol.js, CSP.
- Env files: REACT_APP_* -> VITE_*. process.env.REACT_APP_* -> import.meta.env.VITE_*. Committed .env with defaults. Keep .env.production. Update .env.local.example.
- Build optimization: Default chunk splitting. Hidden source maps. rollup-plugin-visualizer as devDep with build:analyze script.
- CI: New frontend-test.yml workflow. Steps: ESLint, Prettier, Vitest, npm audit, build verification. Trigger on frontend/ path changes.

### Claude's Discretion
- Remove web-vitals now (CRA-specific, confirmed unused in codebase).
- Evaluate lockfile-lint for security value and decide.
- Evaluate whether @ path alias for src/ adds value.
- Evaluate Vite server.proxy for API proxying given existing CORS config.
- Decide CSP meta tag approach for environment-aware connect-src.
- Audit and clean up public assets if warranted.

### Deferred Ideas (OUT OF SCOPE)
None -- discussion stayed within phase scope.
</user_constraints>

<phase_requirements>
## Phase Requirements

| ID | Description | Research Support |
|----|-------------|-----------------|
| BUILD-01 | Frontend builds and runs on Vite instead of CRA | Vite 6.4.1 + @vitejs/plugin-react-swc setup, vite.config.js pattern, index.html migration, env variable mapping |
| BUILD-02 | All dependabot security vulnerabilities resolved (zero alerts) | Removing react-scripts eliminates all 16 current vulnerabilities. Fresh dependency tree with Vite has zero known vulns |
| BUILD-03 | All package.json overrides eliminated | Current overrides (nth-check, postcss, svgo, underscore, serialize-javascript) all trace to react-scripts transitive deps. Removing react-scripts removes the need for all overrides |
| BUILD-04 | All dependencies updated to latest stable versions | Version matrix provided below for all deps. TypeScript 4.9.5 -> 5.x upgrade included |
| BUILD-05 | Docker frontend build works with Vite | Node 22-alpine, dist/ output path, VITE_API_URL build-arg pattern |
| BUILD-06 | Environment variables migrated from REACT_APP_* to VITE_* | 8 files identified needing changes (5 REACT_APP_API_URL + 3 PUBLIC_URL + env files + Dockerfile + docker-compose) |
</phase_requirements>

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| vite | 6.4.1 | Build tool, dev server, HMR | Latest stable 6.x; all plugins explicitly support ^6. Vite 7/8 have breaking changes and incomplete plugin ecosystem support |
| @vitejs/plugin-react-swc | 4.2.3 | React JSX transform via SWC | 20x faster than Babel-based plugin. Peer dep: vite ^4-^7. No custom Babel needed for this codebase |
| vite-plugin-checker | 0.11.0 | TypeScript + ESLint dev overlay | Shows TS/ESLint errors in browser overlay during dev. v0.11.0 supports eslint >=7 and vite >=5.4.20 (v0.12.0 requires eslint >=9.39.1 -- too restrictive) |
| typescript | 5.9.3 | Type checking | Latest stable. Required for Vite compatibility (current 4.9.5 is too old) |

### Testing
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| vitest | 4.1.0 | Test runner | Vite-native, shares vite.config.js, ESM-first. Peer dep: vite ^6-^8 |
| @vitest/coverage-v8 | 4.1.0 | Coverage provider | V8-based coverage for Vitest |
| @testing-library/react | 16.3.2 | React DOM testing | Standard React testing utilities. Supports React 18 |
| @testing-library/jest-dom | 6.9.1 | DOM assertion matchers | Vitest-compatible via setup file. Provides toBeInTheDocument() etc. |
| @testing-library/user-event | 14.6.1 | User interaction simulation | Standard companion to testing-library/react |
| jsdom | (bundled) | Browser environment | Used by Vitest as test environment. Peer dep of vitest, not installed separately |

### Linting & Formatting
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| eslint | 9.39.4 | Linter | Latest stable 9.x. Supports flat config natively |
| eslint-plugin-react | 7.37.5 | React-specific lint rules | Exports flat configs: flat.recommended, flat['jsx-runtime'] |
| eslint-plugin-react-hooks | 7.0.1 | Hooks rules | v7 has proper flat config support |
| eslint-plugin-jsx-a11y | 6.10.2 | Accessibility rules | Standard a11y linting for JSX |
| eslint-config-prettier | 10.1.8 | Disable ESLint formatting rules | Prevents ESLint/Prettier conflicts |
| globals | 17.4.0 | Browser/node globals for ESLint | Required for flat config: `globals.browser` |
| prettier | 3.6.2 | Code formatter | Already in project, keep current |
| husky | 9.1.7 | Git hooks manager | Runs lint-staged on pre-commit |
| lint-staged | 16.3.3 | Run linters on staged files | ESLint + Prettier on staged .js/.jsx files |

### Build Utilities
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| rollup-plugin-visualizer | 7.0.1 | Bundle analysis | devDep only. Generates interactive treemap for bundle size tracking |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Vite 6.4.1 | Vite 8.0.0 | Vite 8 uses Rolldown (new bundler). @vitejs/plugin-react-swc peer deps don't include ^8 yet. Too fresh for production migration |
| Vite 6.4.1 | Vite 7.3.1 | Viable but 6.x is more battle-tested and plugin-react-swc explicitly supports ^6. Marginal benefit |
| @vitejs/plugin-react-swc | @vitejs/plugin-react | Babel-based, slower. No benefit since no custom Babel plugins needed |
| vite-plugin-checker 0.11.0 | vite-plugin-checker 0.12.0 | 0.12.0 requires eslint >= 9.39.1 specifically. Works but unnecessarily strict peer dep |
| ESLint 9.x | ESLint 10.x | ESLint 10.0.3 exists but is brand new. 9.x has proven ecosystem compatibility |
| jsdom | happy-dom | happy-dom is faster but less complete browser emulation. jsdom matches CRA behavior (user decision) |

**Installation (full command):**
```bash
cd frontend

# Remove CRA
npm uninstall react-scripts web-vitals

# Install Vite core
npm install -D vite@^6.4.1 @vitejs/plugin-react-swc@^4.2.3 vite-plugin-checker@^0.11.0

# Install TypeScript (upgrade from devDep 4.9.5)
npm install -D typescript@^5.9.3

# Install testing
npm install -D vitest@^4.1.0 @vitest/coverage-v8@^4.1.0 jsdom
npm install -D @testing-library/react@^16.3.2 @testing-library/jest-dom@^6.9.1 @testing-library/user-event@^14.6.1

# Install ESLint + plugins
npm install -D eslint@^9.39.4 eslint-plugin-react@^7.37.5 eslint-plugin-react-hooks@^7.0.1 eslint-plugin-jsx-a11y@^6.10.2 eslint-config-prettier@^10.1.8 globals@^17.4.0

# Install git hooks
npm install -D husky@^9.1.7 lint-staged@^16.3.3

# Install build utilities
npm install -D rollup-plugin-visualizer@^7.0.1

# Move testing-library deps from dependencies to devDependencies
npm uninstall @testing-library/jest-dom @testing-library/react @testing-library/user-event
# (they get re-added as devDeps by the install -D commands above)
```

## Architecture Patterns

### Recommended Project Structure After Migration
```
frontend/
  index.html                   # Moved from public/ (Vite requirement)
  vite.config.js               # Vite configuration
  eslint.config.mjs            # ESLint v9 flat config (.mjs to force ESM)
  .prettierrc                  # Prettier config (may already exist or use package.json)
  .env                         # Committed defaults (VITE_API_URL=http://localhost:8000/latest)
  .env.production              # Production overrides
  .env.local.example           # Template for local overrides
  .husky/
    pre-commit                 # Runs lint-staged
  public/                      # Static assets (images, manifest, standalone/)
    favicon.ico
    img/
    standalone/                # Embedded standalone apps (InChI viewer etc.)
    manifest.json
  src/
    index.js                   # Entry point (unchanged)
    App.js                     # Router setup (unchanged)
    components/                # React components
    context/                   # AppContext
    hooks/                     # Custom hooks
    pages/                     # Page components
    services/                  # API service layer
    styles/                    # CSS files
    utils/                     # Utility functions
    __tests__/                 # Test files (new)
      setup.js                 # Vitest setup (jest-dom matchers)
```

### Pattern 1: Vite Configuration
**What:** Central vite.config.js with React SWC plugin, checker plugin, dev server, and build options
**When to use:** Single config file for all Vite behavior
**Example:**
```javascript
// frontend/vite.config.js
import { defineConfig } from 'vite';
import react from '@vitejs/plugin-react-swc';
import checker from 'vite-plugin-checker';

export default defineConfig({
  plugins: [
    react(),
    checker({
      typescript: true,
      eslint: {
        lintCommand: 'eslint "./src/**/*.{js,jsx}"',
        useFlatConfig: true,
      },
    }),
  ],
  server: {
    port: 3000,
    // Optional: API proxy (Claude's discretion)
    // proxy: {
    //   '/v1': 'http://localhost:8000',
    //   '/latest': 'http://localhost:8000',
    // },
  },
  build: {
    sourcemap: 'hidden',
    outDir: 'dist',
  },
  test: {
    globals: true,
    environment: 'jsdom',
    setupFiles: './src/__tests__/setup.js',
    css: true,
  },
});
```

### Pattern 2: Environment Variable Migration
**What:** Replace all process.env references with Vite equivalents
**When to use:** Every file that reads environment variables

| CRA Pattern | Vite Pattern | Notes |
|-------------|-------------|-------|
| `process.env.REACT_APP_API_URL` | `import.meta.env.VITE_API_URL` | 5 files: api.js, AppContext.js, MoleculeCard.jsx, HighlightedMoleculeCard.jsx, StandardizeView.jsx |
| `process.env.PUBLIC_URL` | `import.meta.env.BASE_URL` | 3 files: StructureDrawView.jsx, RInChIView.jsx, InChIView.jsx. Note: BASE_URL defaults to '/' |
| `process.env.NODE_ENV` | `import.meta.env.MODE` | 1 file: api.js (error interceptor check). Or use `import.meta.env.PROD` / `import.meta.env.DEV` booleans |
| `GENERATE_SOURCEMAP=false` | `build.sourcemap` in vite.config.js | No longer an env var |

### Pattern 3: ESLint v9 Flat Config
**What:** Modern ESLint configuration using the flat config format
**When to use:** eslint.config.mjs at frontend root (use .mjs to force ESM without requiring "type": "module" in package.json)
**Example:**
```javascript
// frontend/eslint.config.mjs
import js from '@eslint/js';
import globals from 'globals';
import reactPlugin from 'eslint-plugin-react';
import reactHooks from 'eslint-plugin-react-hooks';
import jsxA11y from 'eslint-plugin-jsx-a11y';
import prettierConfig from 'eslint-config-prettier';

export default [
  js.configs.recommended,
  {
    files: ['**/*.{js,jsx}'],
    plugins: {
      react: reactPlugin,
      'react-hooks': reactHooks,
      'jsx-a11y': jsxA11y,
    },
    languageOptions: {
      ecmaVersion: 'latest',
      sourceType: 'module',
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
        version: 'detect',
      },
    },
    rules: {
      ...reactPlugin.configs.recommended.rules,
      ...reactPlugin.configs['jsx-runtime'].rules,
      'react-hooks/rules-of-hooks': 'error',
      'react-hooks/exhaustive-deps': 'warn',
      ...jsxA11y.configs.recommended.rules,
    },
  },
  prettierConfig,
  {
    ignores: ['dist/', 'build/', 'node_modules/', 'public/standalone/'],
  },
];
```

### Pattern 4: Vitest Setup File
**What:** Configure jest-dom matchers for Vitest
**Example:**
```javascript
// frontend/src/__tests__/setup.js
import '@testing-library/jest-dom/vitest';
import { cleanup } from '@testing-library/react';
import { afterEach } from 'vitest';

afterEach(() => {
  cleanup();
});
```

### Pattern 5: index.html Migration
**What:** Move and transform index.html from CRA format to Vite format
**Key changes:**
```html
<!-- BEFORE (public/index.html - CRA) -->
<link rel="icon" href="%PUBLIC_URL%/favicon.ico" />
<link rel="apple-touch-icon" href="%PUBLIC_URL%/img/logo_small.png" />
<link rel="manifest" href="%PUBLIC_URL%/manifest.json" />
<!-- No script tag - CRA injects it -->

<!-- AFTER (frontend/index.html - Vite) -->
<link rel="icon" href="/favicon.ico" />
<link rel="apple-touch-icon" href="/img/logo_small.png" />
<link rel="manifest" href="/manifest.json" />
<!-- Explicit module script entry point (before closing </body>) -->
<script type="module" src="/src/index.js"></script>
```

### Anti-Patterns to Avoid
- **Using `process.env` anywhere in frontend code:** Vite does NOT polyfill `process`. Any remaining `process.env` reference will cause a runtime error ("process is not defined"). Every occurrence must be migrated.
- **Importing CSS modules with `.module.css` when they don't exist:** This codebase uses plain CSS + Tailwind, no CSS modules. Don't add unnecessary complexity.
- **Adding `"type": "module"` to package.json prematurely:** Not needed for Vite. Using `.mjs` extension for ESLint config avoids this entirely. Vite processes its own config files as ESM regardless.
- **Keeping browserslist config:** CRA used browserslist for build targets. Vite uses `build.target` in vite.config.js (defaults to modern browsers). The browserslist section in package.json should be removed.
- **Keeping eslintConfig in package.json:** The CRA-specific `"eslintConfig": { "extends": ["react-app", "react-app/jest"] }` must be removed entirely. Replaced by eslint.config.mjs.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| React JSX transform | Custom Babel config | @vitejs/plugin-react-swc | SWC handles JSX, Fast Refresh, and HMR automatically |
| TypeScript checking in dev | Manual tsc --watch | vite-plugin-checker | Runs in worker thread, shows errors in browser overlay |
| ESLint checking in dev | Manual eslint --watch | vite-plugin-checker (eslint option) | Same worker thread, same overlay |
| Test setup for React | Custom test runner wiring | Vitest with `environment: 'jsdom'` | Shares Vite config, understands Vite transforms |
| Pre-commit linting | Custom git hook scripts | husky + lint-staged | Standard, well-maintained, runs only on staged files |
| Bundle analysis | Custom size tracking | rollup-plugin-visualizer | Interactive treemap, gzip/brotli size display |

**Key insight:** Vite's plugin ecosystem handles all the build complexity that CRA previously hid. The migration is primarily about replacing CRA's hidden config with explicit, visible config files.

## Common Pitfalls

### Pitfall 1: Missing process.env References
**What goes wrong:** App crashes at runtime with "process is not defined" error
**Why it happens:** Vite does NOT inject a process polyfill. CRA did this automatically.
**How to avoid:** Search for ALL `process.env` references across the entire codebase, not just the obvious REACT_APP_ ones. This project has 9 occurrences across 8 files:
- `process.env.REACT_APP_API_URL` (5 files: api.js, AppContext.js, MoleculeCard.jsx, HighlightedMoleculeCard.jsx, StandardizeView.jsx)
- `process.env.PUBLIC_URL` (3 files: StructureDrawView.jsx, RInChIView.jsx, InChIView.jsx)
- `process.env.NODE_ENV` (1 file: api.js error interceptor)
**Warning signs:** Grep for `process\.env` in src/ -- any match is a migration miss.

### Pitfall 2: CSP Meta Tag Blocking Local Dev
**What goes wrong:** Browser blocks API requests during local development because Content-Security-Policy connect-src doesn't include localhost:8000
**Why it happens:** The current CSP meta tag hardcodes `connect-src 'self' https://dev.api.naturalproducts.net https://api.naturalproducts.net https://matomo.nfdi4chem.de` -- no localhost.
**How to avoid:** Use Vite server.proxy so API calls go to 'self' and CSP is satisfied. Alternatively, add localhost:8000 to connect-src for dev builds.
**Warning signs:** API calls fail in dev mode with CSP violation in browser console.

### Pitfall 3: Standalone Embedded Apps Using PUBLIC_URL
**What goes wrong:** Iframe src for standalone InChI/RInChI/StructureDraw editors breaks
**Why it happens:** Three components use `process.env.PUBLIC_URL` to construct iframe src paths to `/standalone/index.html`. The standalone directory is in public/ and contains a full pre-built app (with its own index.html, static/, asset-manifest.json).
**How to avoid:** Replace with `/standalone/index.html` directly (since BASE_URL defaults to '/') or use `import.meta.env.BASE_URL`.
**Warning signs:** Embedded editors don't load; blank iframes.

### Pitfall 4: Docker Build Cache Invalidation
**What goes wrong:** Docker builds are slow because npm ci runs on every code change
**Why it happens:** If COPY . . comes before npm ci, every source change invalidates the npm cache layer.
**How to avoid:** Keep the existing multi-stage pattern: COPY package*.json first, then npm ci, then COPY source. Vite doesn't change this pattern.
**Warning signs:** Docker builds taking 2+ minutes when they should take 30 seconds for code-only changes.

### Pitfall 5: Forgetting to Update docker-compose.yml
**What goes wrong:** docker-compose still passes REACT_APP_API_URL, which Vite ignores
**Why it happens:** docker-compose.yml line 37: `REACT_APP_API_URL=http://localhost:8000/latest` is a CRA env var name
**How to avoid:** Update to `VITE_API_URL` in docker-compose.yml AND in the Dockerfile ARG/ENV declarations.
**Warning signs:** App loads but API calls go to the hardcoded fallback URL instead of the Docker-configured one.

### Pitfall 6: Vitest Not Finding Tests
**What goes wrong:** `npm run test` passes with "0 tests found" -- gives false confidence
**Why it happens:** Vitest uses include patterns (default: `**/*.{test,spec}.{js,mjs,cjs,ts,mts,cts,jsx,tsx}`) and the project has zero tests currently.
**How to avoid:** Create at least one smoke test (e.g., App renders without crashing) as part of the migration to verify the test pipeline works end-to-end.
**Warning signs:** CI passes with 0 tests -- always suspicious.

### Pitfall 7: ESLint Flat Config Import Syntax
**What goes wrong:** `require is not defined` error when running ESLint
**Why it happens:** eslint.config.js uses ESM import syntax by default. If the project doesn't have `"type": "module"` in package.json, Node.js treats .js files as CommonJS.
**How to avoid:** Use `eslint.config.mjs` (force ESM) instead of `.js`. This avoids needing to modify package.json.
**Warning signs:** ESLint fails to start with syntax errors about import/export.

### Pitfall 8: Testing Library Version Mismatch
**What goes wrong:** `toBeInTheDocument` is not a function
**Why it happens:** @testing-library/jest-dom v6+ exports vitest-specific matchers at `@testing-library/jest-dom/vitest`. Importing from the base path may not auto-extend vitest's expect.
**How to avoid:** In setup file, import `@testing-library/jest-dom/vitest` (not just `@testing-library/jest-dom`).
**Warning signs:** Test assertions on DOM elements fail with "not a function" errors.

## Code Examples

### Complete vite.config.js
```javascript
// Source: Vite official docs + verified plugin APIs
import { defineConfig } from 'vite';
import react from '@vitejs/plugin-react-swc';
import checker from 'vite-plugin-checker';

export default defineConfig({
  plugins: [
    react(),
    checker({
      typescript: true,
      eslint: {
        lintCommand: 'eslint "./src/**/*.{js,jsx}"',
        useFlatConfig: true,
      },
    }),
  ],
  server: {
    port: 3000,
    open: false,
  },
  build: {
    sourcemap: 'hidden',
    outDir: 'dist',
  },
});
```

### Updated package.json scripts
```json
{
  "scripts": {
    "dev": "vite",
    "build": "vite build",
    "preview": "vite preview",
    "test": "vitest run",
    "test:watch": "vitest",
    "test:coverage": "vitest run --coverage",
    "lint": "eslint \"src/**/*.{js,jsx}\"",
    "lint:fix": "eslint \"src/**/*.{js,jsx}\" --fix",
    "format": "prettier --write \"src/**/*.{js,jsx,json,css,md}\"",
    "format:check": "prettier --check \"src/**/*.{js,jsx,json,css,md}\"",
    "build:analyze": "vite build",
    "prepare": "husky"
  }
}
```

Note on `build:analyze`: rollup-plugin-visualizer should be added to vite.config.js plugins conditionally or as a separate config:
```javascript
// In vite.config.js -- add to plugins array for analysis builds
import { visualizer } from 'rollup-plugin-visualizer';

// In plugins:
visualizer({ open: true, filename: 'stats.html' }),
```

A practical approach: use an environment variable to enable it:
```javascript
...(process.env.ANALYZE && [visualizer({ open: true, filename: 'stats.html' })]),
```
Then: `ANALYZE=true npm run build`

### Updated Dockerfile
```dockerfile
# Stage 1: Build
FROM node:22-alpine AS build
WORKDIR /app
COPY package*.json ./
RUN npm ci --silent
COPY . .
ARG VITE_API_URL=http://localhost:8000/latest
ENV VITE_API_URL=$VITE_API_URL
RUN npm run build

# Stage 2: Serve
FROM nginx:alpine
COPY --from=build /app/dist /usr/share/nginx/html
# nginx config for React Router SPA support
RUN echo 'server { \
    listen 80; \
    server_name localhost; \
    root /usr/share/nginx/html; \
    index index.html; \
    location / { \
        try_files $uri $uri/ /index.html; \
    } \
    location ~* \.(js|css|png|jpg|jpeg|gif|ico|svg)$ { \
        expires 1y; \
        add_header Cache-Control "public, max-age=31536000, immutable"; \
    } \
}' > /etc/nginx/conf.d/default.conf
EXPOSE 80
CMD ["nginx", "-g", "daemon off;"]
```

### Environment Variable Migration Map
```javascript
// BEFORE (CRA)
const API_URL = process.env.REACT_APP_API_URL || "https://dev.api.naturalproducts.net/latest";

// AFTER (Vite)
const API_URL = import.meta.env.VITE_API_URL || "https://dev.api.naturalproducts.net/latest";
```

```javascript
// BEFORE (CRA) -- checking environment
if (process.env.NODE_ENV !== "production") { ... }

// AFTER (Vite) -- use boolean helpers
if (import.meta.env.DEV) { ... }
// OR: if (import.meta.env.MODE !== "production") { ... }
```

```javascript
// BEFORE (CRA) -- PUBLIC_URL for iframe
src={`${process.env.PUBLIC_URL}/standalone/index.html`}

// AFTER (Vite) -- just use absolute path (base is '/')
src="/standalone/index.html"
// OR with BASE_URL: src={`${import.meta.env.BASE_URL}standalone/index.html`}
```

### Vitest Configuration (inline in vite.config.js)
```javascript
// frontend/vite.config.js
/// <reference types="vitest/config" />
import { defineConfig } from 'vite';
import react from '@vitejs/plugin-react-swc';

export default defineConfig({
  plugins: [react()],
  test: {
    globals: true,
    environment: 'jsdom',
    setupFiles: './src/__tests__/setup.js',
    css: true,
    coverage: {
      provider: 'v8',
      reporter: ['text', 'lcov'],
    },
  },
});
```

### Smoke Test Example
```javascript
// frontend/src/__tests__/App.test.jsx
import { render } from '@testing-library/react';
import { describe, it, expect } from 'vitest';
import App from '../App';

describe('App', () => {
  it('renders without crashing', () => {
    const { container } = render(<App />);
    expect(container).toBeTruthy();
  });
});
```

### lint-staged Configuration (in package.json)
```json
{
  "lint-staged": {
    "src/**/*.{js,jsx}": [
      "eslint --fix",
      "prettier --write"
    ],
    "src/**/*.{json,css,md}": [
      "prettier --write"
    ]
  }
}
```

### .husky/pre-commit
```bash
npx lint-staged
```

### CI Workflow Example
```yaml
# .github/workflows/frontend-test.yml
name: frontend-test

on:
  pull_request:
    branches: [main, development]
    paths:
      - 'frontend/**'
  push:
    branches: [main, development]
    paths:
      - 'frontend/**'

jobs:
  test:
    runs-on: ubuntu-latest
    defaults:
      run:
        working-directory: frontend
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-node@v4
        with:
          node-version: '22'
          cache: 'npm'
          cache-dependency-path: frontend/package-lock.json
      - run: npm ci
      - name: ESLint
        run: npm run lint
      - name: Prettier
        run: npm run format:check
      - name: Tests
        run: npm run test
      - name: Security audit
        run: npm audit --audit-level=low
      - name: Production build
        run: npm run build
```

## Dependency Audit Results

### dompurify: KEEP (actively used)
- Used in 2 components for sanitizing HTML/SVG content before rendering
- Security-critical: sanitizes server-returned content using DOMPurify.sanitize()
- **Verdict:** Essential. Keep.

### jszip: KEEP (actively used)
- Used in 2 files: `depictService.js` and `Depict2DMultiView.jsx`
- Dynamically imported (`await import('jszip')`) for multi-molecule ZIP download feature
- **Verdict:** Actively used for batch download. Keep.

### web-vitals: REMOVE
- Zero references in src/ (confirmed via grep)
- CRA-specific package for performance metrics reporting
- **Verdict:** Remove. Not used.

### Overrides Analysis
All 5 overrides trace back to react-scripts transitive dependencies:
- `nth-check`: via react-scripts -> css-select
- `postcss`: via react-scripts -> css-minimizer-webpack-plugin (already in devDeps independently for Tailwind)
- `svgo`: via react-scripts -> css-minimizer-webpack-plugin -> svgo
- `underscore`: via react-scripts -> workbox-build
- `serialize-javascript`: via react-scripts -> css-minimizer-webpack-plugin, workbox-build

Removing react-scripts eliminates all need for overrides. The independent postcss devDep (^8.4.31) remains for Tailwind and has no vulnerabilities.

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| CRA (react-scripts) | Vite | CRA deprecated Feb 2025 | Faster builds, ESM-native, no hidden config |
| ESLint .eslintrc (legacy config) | eslint.config.js (flat config) | ESLint 9 (Apr 2024) | Simpler, composable, no extends chains |
| Jest | Vitest | Vitest 1.0 (Dec 2023) | Shares Vite config, faster, ESM-native |
| process.env.REACT_APP_* | import.meta.env.VITE_* | Vite convention | ESM standard, no Node.js polyfill needed |
| Webpack + Babel | Vite + SWC | Industry shift 2023-2025 | 10-20x faster dev server startup, instant HMR |
| eslint-config-react-app | eslint-plugin-react + hooks + a11y | CRA deprecated | Explicit plugin control, no hidden rules |

**Deprecated/outdated:**
- `react-scripts`: Officially sunset. No security patches.
- `web-vitals`: CRA boilerplate. Not needed without CRA's reportWebVitals().
- `eslintConfig.extends: ["react-app"]`: CRA-specific config. Replaced by explicit plugins.
- `browserslist` in package.json: CRA used this for Babel targets. Vite uses `build.target` config option.

## Open Questions

1. **Vite server.proxy vs CORS (Claude's Discretion)**
   - What we know: Backend already has `ALLOWED_ORIGINS=*` (default) CORS config. Frontend works without proxy today.
   - What's unclear: Whether a Vite proxy would improve DX by avoiding CSP connect-src issues in dev
   - Recommendation: Add proxy for `/v1` and `/latest` to localhost:8000. This solves the CSP issue in dev and is standard Vite practice. Zero impact on production (proxy only applies to dev server).

2. **CSP Meta Tag Environment Awareness (Claude's Discretion)**
   - What we know: Current CSP hardcodes production URLs in connect-src. Doesn't include localhost.
   - What's unclear: Best approach -- build-time injection? Separate dev/prod? Proxy workaround?
   - Recommendation: If using Vite proxy (question 1), CSP is satisfied in dev because API calls go to 'self'. For production, keep current CSP. No changes needed to the CSP tag itself.

3. **Path Alias @ for src/ (Claude's Discretion)**
   - What we know: Codebase uses relative paths (`../../components/common/Header`). Imports are 2-3 levels deep max.
   - What's unclear: Whether @ alias improves readability enough to justify config overhead
   - Recommendation: Skip. Import depth is manageable (max 3 levels). Adding alias requires vite.config.js resolve.alias + jsconfig/tsconfig paths + ESLint resolver config. Not worth it for this codebase.

4. **lockfile-lint (Claude's Discretion)**
   - What we know: Protects against registry confusion attacks and modified lockfiles
   - What's unclear: Whether this project's threat model warrants it (public repo, npm-only deps)
   - Recommendation: Skip. The project uses npm exclusively (no mixed registries), is open source, and CI runs `npm ci` which validates lockfile integrity. lockfile-lint adds marginal value here.

5. **eslint.config.js vs eslint.config.mjs**
   - What we know: Without `"type": "module"` in package.json, .js files are CJS and ESM imports fail
   - What's unclear: Whether adding "type": "module" to package.json has side effects on other tools
   - Recommendation: Use `eslint.config.mjs` to force ESM without touching package.json. Safest approach. vite-plugin-checker supports .mjs config files.

## Validation Architecture

### Test Framework
| Property | Value |
|----------|-------|
| Framework | Vitest 4.1.0 + @testing-library/react 16.3.2 |
| Config file | vite.config.js (test section) -- Wave 0 creates this |
| Quick run command | `cd frontend && npx vitest run` |
| Full suite command | `cd frontend && npx vitest run --coverage` |

### Phase Requirements to Test Map
| Req ID | Behavior | Test Type | Automated Command | File Exists? |
|--------|----------|-----------|-------------------|-------------|
| BUILD-01 | Vite dev server starts and serves app | smoke | `cd frontend && npx vite build` (build succeeds = Vite works) | N/A -- build verification |
| BUILD-02 | Zero npm audit vulnerabilities | security | `cd frontend && npm audit --audit-level=low` | N/A -- audit command |
| BUILD-03 | No overrides in package.json | manual-only | Inspect package.json for `"overrides"` key | N/A -- file inspection |
| BUILD-04 | Dependencies at latest stable | manual-only | `cd frontend && npm outdated` | N/A -- npm command |
| BUILD-05 | Docker build succeeds | integration | `docker build -t test-frontend frontend/` | N/A -- Docker command |
| BUILD-06 | Env vars use VITE_* prefix | unit | `cd frontend && npx vitest run src/__tests__/env.test.js -x` | Wave 0 |

### Sampling Rate
- **Per task commit:** `cd frontend && npx vite build && npm audit --audit-level=low`
- **Per wave merge:** `cd frontend && npx vitest run --coverage && npx vite build`
- **Phase gate:** Full suite green + `npm audit` clean + Docker build succeeds

### Wave 0 Gaps
- [ ] `frontend/vite.config.js` -- Vite + Vitest configuration
- [ ] `frontend/src/__tests__/setup.js` -- jest-dom matcher setup for Vitest
- [ ] `frontend/src/__tests__/App.test.jsx` -- smoke test: App renders
- [ ] `frontend/src/__tests__/env.test.js` -- verify env vars use import.meta.env pattern
- [ ] `frontend/eslint.config.mjs` -- ESLint v9 flat config
- [ ] Framework install: `npm install -D vitest @vitest/coverage-v8 jsdom` -- if not yet installed

## Sources

### Primary (HIGH confidence)
- npm registry (`npm view <pkg> version/peerDependencies`) -- verified all version numbers and peer dependency compatibility on 2026-03-12
- [Vite official docs: Env Variables and Modes](https://vite.dev/guide/env-and-mode) -- environment variable patterns
- [Vite official releases](https://vite.dev/releases) -- version history, Vite 6/7/8 status
- [ESLint official docs: Configuration Files](https://eslint.org/docs/latest/use/configure/configuration-files) -- flat config format
- [Vitest Configuration](https://vitest.dev/config/) -- test configuration options

### Secondary (MEDIUM confidence)
- [eslint-plugin-react GitHub](https://github.com/jsx-eslint/eslint-plugin-react) -- flat config exports (flat.recommended, flat['jsx-runtime'])
- [eslint-plugin-react-hooks flat config PR](https://github.com/facebook/react/pull/30774) -- confirmed flat config support in v5.2+/v7
- [vite-plugin-checker docs](https://vite-plugin-checker.netlify.app/) -- ESLint checker configuration, useFlatConfig option
- [Vite 8 Beta announcement](https://vite.dev/blog/announcing-vite8-beta) -- Rolldown migration details, plugin compatibility
- [CRA to Vite migration guide (DEV Community)](https://dev.to/solitrix02/goodbye-cra-hello-vite-a-developers-2026-survival-guide-for-migration-2a9f) -- general migration patterns
- [Vitest + React Testing Library guide (DEV Community)](https://dev.to/cristiansifuentes/mastering-vitest-react-testing-library-fixing-beforeeach-tobeinthedocument-and-jsdom-2379) -- setup patterns and gotchas
- [Husky getting started](https://typicode.github.io/husky/get-started.html) -- husky 9 setup
- [plugin-react-swc Vite 8 support issue](https://github.com/vitejs/vite-plugin-react/issues/1012) -- confirms peer deps not yet updated for Vite 8

### Tertiary (LOW confidence)
- None -- all findings verified against npm registry or official docs

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH -- all versions verified via npm registry, peer deps cross-checked
- Architecture: HIGH -- patterns derived from official Vite docs and confirmed by codebase analysis
- Pitfalls: HIGH -- each pitfall identified by examining actual codebase files and verified against Vite behavior
- Dependency audit: HIGH -- grep-verified usage of dompurify (3 usages), jszip (2 usages), web-vitals (0 usages)
- ESLint setup: MEDIUM -- flat config is well-documented but plugin interop can have edge cases (mitigated by using .mjs extension)

**Research date:** 2026-03-12
**Valid until:** 2026-04-12 (stable domain, 30-day window reasonable)
