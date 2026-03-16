# Codebase Structure

**Analysis Date:** 2026-03-12

## Directory Layout

```
cheminformatics-microservice/
├── app/                           # FastAPI application (Python backend)
│   ├── __init__.py                # Package metadata
│   ├── main.py                    # FastAPI app entry point, versioning, middleware setup
│   ├── exception_handlers.py       # Custom exception handler for InvalidInputException
│   ├── limiter.py                 # slowapi rate limiter configuration
│   ├── routers/                   # API endpoint handlers (one per feature area)
│   │   ├── __init__.py            # (empty)
│   │   ├── chem.py                # /chem/* — chemistry analysis endpoints
│   │   ├── converters.py          # /convert/* — molecular format conversion endpoints
│   │   ├── depict.py              # /depict/* — 2D/3D structure visualization endpoints
│   │   ├── tools.py               # /tools/* — sugar removal, structure generation endpoints
│   │   └── ocsr.py                # /ocsr/* — optical structure recognition (conditional, INCLUDE_OCSR env var)
│   ├── schemas/                   # Pydantic request/response models per feature area
│   │   ├── __init__.py            # Exports common schemas (HealthCheck)
│   │   ├── chem_schema.py         # Response models for /chem endpoints
│   │   ├── converters_schema.py   # Response models for /convert endpoints
│   │   ├── depict_schema.py       # Response models for /depict endpoints
│   │   ├── tools_schema.py        # Response models for /tools endpoints
│   │   ├── ocsr_schema.py         # Response models for /ocsr endpoints
│   │   ├── error.py               # Error response models (ErrorResponse, BadRequestModel)
│   │   ├── healthcheck.py         # HealthCheck model (used by all /health endpoints)
│   │   ├── classyfire.py          # ClassyFire classification API models
│   │   ├── pubchem_schema.py      # PubChem retrieval response model
│   │   ├── chemblstandardizer.py  # ChEMBL standardization response models
│   │   ├── coconut.py             # COCONUT descriptor/preprocessing models
│   │   └── msg.py                 # General message schema
│   ├── modules/                   # Chemistry logic and toolkit wrappers
│   │   ├── toolkits/              # Toolkit abstraction layer (3 parallel engines)
│   │   │   ├── helpers.py         # parse_input(), parse_SMILES() — input validation & dispatch
│   │   │   ├── rdkit_wrapper.py   # RDKit functions (descriptors, depiction, filters, 2D/3D conformers)
│   │   │   ├── cdk_wrapper.py     # CDK via JPype (canonical SMILES, InChI, HOSE codes, depiction)
│   │   │   └── openbabel_wrapper.py # OpenBabel functions (canonical SMILES, InChI)
│   │   ├── cdk_depict/            # Modular CDK depiction feature plugins
│   │   │   ├── __init__.py        # Feature registration/initialization
│   │   │   ├── abbreviations.py   # Functional group abbreviations
│   │   │   ├── annotations.py     # Molecule annotations/labels
│   │   │   ├── aromatic_display.py # Aromatic ring rendering
│   │   │   ├── dative_bonds.py    # Dative/coordinate bond rendering
│   │   │   ├── hydrogen_display.py # Hydrogen atom visibility control
│   │   │   ├── mdl_hilite.py      # MDL highlight color parsing
│   │   │   ├── multicenter_bonds.py # Multicenter bond handling
│   │   │   ├── radical_perception.py # Radical/unpaired electron rendering
│   │   │   ├── reaction_depiction.py # Reaction arrow/layout
│   │   │   ├── advanced_controls.py # Advanced rendering options
│   │   │   ├── cxsmiles_parser.py # CXSMILES extension parsing
│   │   │   ├── style_presets.py   # Rendering style templates
│   │   │   └── data/              # SVG assets, style templates, data files
│   │   ├── all_descriptors.py     # Molecular descriptor calculation (RDKit-based)
│   │   ├── depiction.py           # RDKit-based 2D depiction (generates PNG/SVG)
│   │   ├── depiction_enhanced.py  # CDK-based 2D depiction with plugin system
│   │   ├── converters/            # (implied) Format conversion utilities
│   │   ├── coconut/               # COCONUT natural product database integration
│   │   │   ├── descriptors.py     # COCONUT descriptor calculation
│   │   │   └── preprocess.py      # COCONUT data preprocessing
│   │   ├── tools/                 # Domain-specific chemistry tools
│   │   │   ├── sugar_removal.py   # SugarRemovalUtility (linear/circular sugar removal)
│   │   │   └── surge.py           # SURGE structure generation
│   │   ├── classyfire.py          # ClassyFire REST API client
│   │   ├── decimer.py             # DECIMER optical structure recognition (TensorFlow model)
│   │   ├── fix_radicals.py        # Radical atom fixing utility
│   │   ├── npscorer.py            # Natural Product likeliness scoring
│   │   └── pubchem_retrieve.py    # PubChem REST API client
│   └── templates/                 # Jinja2 HTML templates (if any)
│       ├── mol.html               # Molecule viewer template
│       └── style.css              # Template styling
├── frontend/                      # React 18 SPA (TypeScript-compatible JavaScript)
│   ├── src/
│   │   ├── App.js                 # Root component, router setup (React Router v6)
│   │   ├── index.js               # React DOM mount point
│   │   ├── components/            # Reusable UI components
│   │   │   ├── chem/              # Chemistry analysis page components
│   │   │   ├── common/            # Shared components (Header, Footer, common UI)
│   │   │   ├── convert/           # Format conversion page components
│   │   │   ├── depict/            # Structure visualization page components
│   │   │   ├── ocsr/              # OCSR/image input page components
│   │   │   └── tools/             # Chemistry tools page components
│   │   ├── pages/                 # Page-level components (full page routes)
│   │   │   ├── HomePage.js
│   │   │   ├── ChemPage.js
│   │   │   ├── ConvertPage.js
│   │   │   ├── DepictPage.js
│   │   │   ├── ToolsPage.js
│   │   │   ├── OCSRPage.js
│   │   │   ├── AboutPage.js
│   │   │   ├── TermsOfService.js
│   │   │   └── PrivacyPolicy.js
│   │   ├── services/              # API communication layer
│   │   │   ├── api.js             # Axios instance with base URL from REACT_APP_API_URL
│   │   │   ├── chemService.js     # API calls for /chem endpoints
│   │   │   ├── convertService.js  # API calls for /convert endpoints
│   │   │   ├── depictService.js   # API calls for /depict endpoints
│   │   │   ├── toolsService.js    # API calls for /tools endpoints
│   │   │   └── ocsrService.js     # API calls for /ocsr endpoints
│   │   ├── context/               # React context for global state (AppContext)
│   │   ├── hooks/                 # Custom React hooks
│   │   ├── utils/                 # Utility functions (parsing, formatting, etc.)
│   │   └── styles/                # Global CSS/Tailwind configuration
│   ├── public/                    # Static assets (favicon, manifest, images)
│   ├── package.json               # npm dependencies, scripts, build config
│   ├── tailwind.config.js         # Tailwind CSS configuration
│   ├── postcss.config.js          # PostCSS configuration
│   └── Dockerfile                 # Docker image for frontend (nginx)
├── tests/                         # Pytest test suite (88% coverage)
│   ├── test_chem.py               # Tests for /chem endpoints
│   ├── test_converters.py         # Tests for /convert endpoints
│   ├── test_depict.py             # Tests for /depict endpoints
│   ├── test_depiction_enhanced.py # Tests for CDK depiction module
│   ├── test_cdk_*.py              # Granular tests for each CDK depict plugin
│   ├── test_tools.py              # Tests for /tools endpoints
│   ├── test_cdk_*.py              # Tests for CDK toolkit wrapper functions
│   ├── test_ocsr.py               # Tests for /ocsr endpoints (conditional)
│   ├── test_classyfire.py         # Tests for ClassyFire integration
│   ├── test_helper.py             # Tests for input parsing helpers
│   ├── test_limiter.py            # Tests for rate limiting
│   ├── test_main.py               # Tests for main app setup
│   └── [test_data]                # Test fixtures (ketcher.cdx, *.png images)
├── docs/                          # Sphinx/Markdown documentation
│   ├── index.md / index.rst       # Documentation home
│   ├── installation.md            # Setup instructions
│   ├── architecture.md            # System architecture overview
│   ├── chem.md                    # /chem endpoint documentation
│   ├── conversions.md             # /convert endpoint documentation
│   ├── depict.md                  # /depict endpoint documentation
│   ├── tools.md                   # /tools endpoint documentation
│   ├── ocsr.md                    # /ocsr endpoint documentation
│   ├── api-examples.md            # Example API calls
│   ├── docker.md                  # Docker deployment guide
│   ├── cluster-deployment.md      # Kubernetes deployment guide
│   ├── scaling.md                 # Scaling strategies
│   ├── sustainability.md          # Project sustainability
│   ├── introduction.md            # Project overview
│   ├── contributors.md            # Contributors list
│   ├── issues.md                  # Known issues
│   ├── python_modules.rst         # Python API reference
│   └── public/                    # Documentation assets (logos, diagrams)
├── .github/                       # GitHub Actions CI/CD
│   └── workflows/
│       └── test.yml               # Test and lint workflow
├── ops/                           # Operations scripts and configs
│   ├── docker-compose.yml         # Multi-container orchestration
│   ├── docker-compose-dev.yml     # Development variant
│   ├── docker-compose-prod.yml    # Production variant
│   ├── docker-compose.lite.yml    # Lite variant (no OCSR)
│   └── [deployment scripts]       # Zero-downtime deploy scripts
├── Dockerfile                     # Main backend container image
├── Dockerfile.lite                # Lite backend (no TensorFlow/DECIMER)
├── docker-compose.yml             # Docker Compose (api:8000, frontend:3000, prometheus:9090, grafana:3001)
├── cm-dashboard.json              # Grafana dashboard configuration
├── requirements.txt               # Python dependencies (full)
├── requirements_lite.txt          # Python dependencies (lite, no OCSR)
├── package.json                   # Frontend npm dependencies
├── package-lock.json              # Frontend dependency lock
├── .pre-commit-config.yaml        # Pre-commit hook setup (black, flake8, reorder-python-imports, pydocstringformatter)
├── .readthedocs.yml               # ReadTheDocs build configuration
├── .gitignore                     # Git ignore rules
├── .env.example                   # Example env var template (NEVER read actual .env files)
├── .env.development               # Development environment config (tracked, no secrets)
├── .env.production                # Production environment config (tracked, no secrets)
├── CHANGELOG.md                   # Versioned release notes (release-please managed)
├── CITATION.cff                   # Citation metadata
├── CLAUDE.md                      # Claude Code instructions (this file's project context)
├── README.md                      # Project overview, quickstart
├── LICENSE                        # License (CC BY 4.0)
├── Makefile / make.bat            # Documentation build helpers
└── .planning/                     # GSD planning documents (generated)
    └── codebase/                  # Architecture/structure documentation
        ├── ARCHITECTURE.md        # (this file's counterpart)
        └── STRUCTURE.md           # (this file)
```

## Directory Purposes

**`app/`:**
- Purpose: Main Python FastAPI backend application
- Contains: HTTP routers, Pydantic schemas, chemistry modules, exception handling, middleware
- Key files: `main.py` (entry point), `exception_handlers.py`, `limiter.py`

**`app/routers/`:**
- Purpose: API endpoint implementations grouped by feature area
- Contains: APIRouter instances with request handlers
- Key files: `chem.py` (70+ endpoints), `converters.py`, `depict.py`, `tools.py`, `ocsr.py`
- Pattern: Each file defines a router with `/health` endpoint and domain-specific endpoints

**`app/schemas/`:**
- Purpose: Pydantic request/response models for automatic validation and OpenAPI documentation
- Contains: Response models (no request models in body — most use query params)
- Key files: `chem_schema.py` (20+ response models), `converters_schema.py`, `depict_schema.py`, `error.py` (error response models)
- Pattern: Organized by router; each endpoint response gets a dedicated model

**`app/modules/`:**
- Purpose: Chemistry logic layer — functions that perform actual computations
- Contains: Descriptor calculation, depiction, format conversion, external API clients
- Key files: `all_descriptors.py`, `depiction.py`, `depiction_enhanced.py`
- Pattern: Pure functions that take molecule objects and return data; no HTTP concerns

**`app/modules/toolkits/`:**
- Purpose: Abstraction layer for 3 parallel chemistry engines
- Contains: Wrapper functions for RDKit, CDK (via JPype), OpenBabel
- Key files: `rdkit_wrapper.py` (37KB, ~800 lines), `cdk_wrapper.py` (24KB), `helpers.py` (parse logic)
- Pattern: One wrapper per toolkit; functions return toolkit-specific molecule objects or computed results

**`app/modules/cdk_depict/`:**
- Purpose: Modular CDK depiction feature plugins
- Contains: 13 specialized depiction feature files (abbreviations, annotations, aromatic display, etc.)
- Key files: `abbreviations.py`, `annotations.py`, `hydrogen_display.py`, `style_presets.py`
- Pattern: Each feature is independent; imported and composed by `depiction_enhanced.py`

**`frontend/src/`:**
- Purpose: React 18 single-page application (SPA) for user interaction
- Contains: Components, pages, services (API clients), global context, utilities
- Key files: `App.js` (router setup), `services/` (API layer), `components/` (reusable UI)
- Pattern: Component structure mirrors API structure (chem/, convert/, depict/, tools/, ocsr/)

**`frontend/src/services/`:**
- Purpose: API communication layer (Axios-based HTTP clients)
- Contains: Service functions that call backend endpoints
- Key files: `api.js` (Axios config with base URL), `chemService.js`, `depictService.js`, etc.
- Pattern: One service file per API router; functions map to endpoints

**`tests/`:**
- Purpose: Pytest test suite (88% coverage)
- Contains: Unit tests for modules, integration tests for endpoints
- Key files: `test_chem.py` (54KB, 700+ lines), `test_converters.py`, `test_depict.py`
- Pattern: Test files mirror structure (test_<router>.py); CDK depict features have granular test_cdk_<feature>.py files

**`docs/`:**
- Purpose: Sphinx/Markdown documentation for end users and developers
- Contains: Installation guide, API reference, deployment guide, examples
- Key files: `index.md`, `architecture.md`, `chem.md`, `docker.md`
- Pattern: Markdown for user-facing docs; RST for Python API autodoc

**`.github/workflows/`:**
- Purpose: GitHub Actions CI/CD automation
- Contains: Test, lint, build, publish workflows
- Key files: `test.yml` (runs flake8, pytest, uploads coverage to Codecov)

**`ops/`:**
- Purpose: Deployment and operations orchestration
- Contains: Docker Compose variants, zero-downtime deploy scripts
- Key files: `docker-compose.yml` (all services), `docker-compose.lite.yml` (OCSR disabled)

## Key File Locations

**Entry Points:**
- `app/main.py`: FastAPI app initialization, router registration, middleware setup
- `frontend/src/index.js`: React DOM mount point
- `frontend/src/App.js`: React Router configuration, layout setup
- `Dockerfile`: Backend container build
- `frontend/Dockerfile`: Frontend container build
- `docker-compose.yml`: Full stack orchestration (API, frontend, Prometheus, Grafana)

**Configuration:**
- `app/main.py`: Environment variable parsing (ALLOWED_ORIGINS, INCLUDE_OCSR, RELEASE_VERSION, HOMEPAGE_URL, CMS_INTERNAL_AUTH_TOKEN)
- `.env.example`: Template for environment variables (database URI, API keys, etc. — NOT tracked if secrets)
- `.env.development`, `.env.production`: Non-secret config (tracked in git)
- `requirements.txt`, `requirements_lite.txt`: Python dependencies
- `package.json`: Frontend dependencies
- `tailwind.config.js`, `postcss.config.js`: Frontend build config
- `.pre-commit-config.yaml`: Git hook configuration (black, flake8, reorder-python-imports)

**Core Logic:**
- `app/modules/toolkits/helpers.py`: Input parsing and toolkit dispatch (`parse_input()`)
- `app/modules/toolkits/rdkit_wrapper.py`: RDKit functions (descriptors, depiction, conformers, filters)
- `app/modules/toolkits/cdk_wrapper.py`: CDK functions (canonical SMILES, InChI, depiction)
- `app/modules/all_descriptors.py`: Descriptor calculation
- `app/modules/depiction.py`: RDKit depiction
- `app/modules/depiction_enhanced.py`: CDK depiction with plugins
- `app/modules/tools/sugar_removal.py`: Natural product preprocessing (SugarRemovalUtility)

**Testing:**
- `tests/test_chem.py`: Tests for chemistry analysis endpoints
- `tests/test_converters.py`: Tests for format conversion
- `tests/test_depict.py`: Tests for depiction endpoints
- `tests/test_helper.py`: Tests for input parsing
- `conftest.py`: (if exists) Shared pytest fixtures and configuration

**Error Handling & Security:**
- `app/exception_handlers.py`: Custom exception handler for InvalidInputException
- `app/limiter.py`: Rate limiting configuration and SSRF validation
- `app/routers/ocsr.py:_validate_url()`: SSRF protection for external image fetching

## Naming Conventions

**Files:**
- Routers: `<domain>.py` (chem.py, converters.py, depict.py) — lowercase, underscore-separated
- Schemas: `<domain>_schema.py` (chem_schema.py, converters_schema.py) — lowercase, underscore-separated
- Modules: `<feature>.py` (all_descriptors.py, depiction.py) — lowercase, underscore-separated
- Tests: `test_<domain>.py` (test_chem.py, test_converters.py) — pytest convention

**Directories:**
- Features: `<domain>/` (routers/, schemas/, modules/) — lowercase, underscore-separated
- Packages: `app/`, `frontend/src/` — project organization
- Toolkit wrappers: `toolkits/` — indicates toolkit abstraction layer

**Functions:**
- Toolkit wrappers: `get_<result>_<toolkit>()` (e.g., `get_tanimoto_similarity_rdkit()`, `get_CDK_HOSE_codes()`)
- Schema response models: `Generate<Result>Response` or `<Result>Response` (e.g., `GenerateDescriptorsResponse`, `HealthCheck`)
- API endpoints: `@router.get("/path")` with handler function named `<operation>()` (e.g., `def properties()`)

**Variables:**
- Molecule objects: `mol`, `molecule`, `rdkit_mol`, `cdk_mol`
- Strings: `smiles`, `inchi`, `inchi_key`, `input`
- Collections: `molecules` (plural), `descriptors` (dict), `results` (list)
- Booleans: `is_valid`, `include_<feature>`, `standardize`

## Where to Add New Code

**New Chemistry Analysis Endpoint:**
1. Add endpoint handler in `app/routers/chem.py` with `@router.get()` or `@router.post()`
2. Create response schema in `app/schemas/chem_schema.py` as `Generate<Feature>Response`
3. Implement logic function in `app/modules/<feature>.py` (or extend existing module)
4. If logic needs toolkit dispatch, call `parse_input()` from `app/modules/toolkits/helpers.py`, then call toolkit-specific function in wrappers
5. Add tests in `tests/test_chem.py` following existing test pattern

**New Format Converter:**
1. Add endpoint in `app/routers/converters.py`
2. Add schema in `app/schemas/converters_schema.py`
3. Implement conversion logic in `app/modules/toolkits/<toolkit>_wrapper.py` (e.g., `get_InChI_cdk()`)
4. Use toolkit dispatch in `helpers.py` if needed
5. Add tests in `tests/test_converters.py`

**New Depiction Feature:**
1. For RDKit: add function to `app/modules/depiction.py`
2. For CDK: create `app/modules/cdk_depict/<feature>.py` plugin file following existing patterns (see abbreviations.py, annotations.py)
3. Register feature in `app/modules/cdk_depict/__init__.py`
4. Call feature in `app/modules/depiction_enhanced.py` composition
5. Add endpoint in `app/routers/depict.py` if new endpoint needed
6. Add schema in `app/schemas/depict_schema.py`
7. Add granular test file `tests/test_cdk_<feature>.py`

**New Frontend Page:**
1. Create page component in `frontend/src/pages/<FeatureName>Page.js`
2. Create feature components in `frontend/src/components/<feature>/`
3. Add service functions in `frontend/src/services/<feature>Service.js` that use `frontend/src/services/api.js` Axios instance
4. Add route in `frontend/src/App.js` router configuration
5. Add navigation link in `frontend/src/components/common/Header.js` if needed

**New Toolkit Wrapper:**
1. Create `app/modules/toolkits/<toolkit>_wrapper.py`
2. Add import in `app/modules/toolkits/helpers.py:parse_SMILES()`
3. Add case to `framework` parameter dispatch logic
4. Implement wrapper functions following naming convention `get_<result>_<toolkit>()`
5. Add tests in `tests/test_helper.py` and toolkit-specific test files

**New Tool (Sugar Removal, Structure Generation, etc.):**
1. Create module in `app/modules/tools/<tool>.py` (e.g., `sugar_removal.py`)
2. Add endpoint in `app/routers/tools.py`
3. Add schema in `app/schemas/tools_schema.py`
4. Add tests in `tests/test_tools.py`

## Special Directories

**`app/templates/`:**
- Purpose: Jinja2 HTML templates for server-side rendering (if any)
- Generated: No
- Committed: Yes
- Contains: mol.html, style.css

**`frontend/public/`:**
- Purpose: Static assets served by React (favicon, manifest, logo)
- Generated: No (build output goes to `frontend/build/`)
- Committed: Yes (source files)
- Contents: favicon.ico, manifest.json, index.html, images

**`frontend/build/`:**
- Purpose: Production build output from `npm run build`
- Generated: Yes (by npm build script)
- Committed: No (in .gitignore)
- Contents: Minified JS/CSS, HTML, static assets

**`docs/public/`:**
- Purpose: Documentation assets (logos, diagrams, screenshots)
- Generated: No
- Committed: Yes
- Contents: PNG, SVG images used in markdown

**`htmlcov/`:**
- Purpose: HTML coverage report output from `pytest --cov --cov-report=html`
- Generated: Yes
- Committed: No (in .gitignore)

**`prometheus_data/`:**
- Purpose: Prometheus time-series database volume (from docker-compose)
- Generated: Yes (by Prometheus container)
- Committed: No (in .gitignore)

**`.ipynb_checkpoints/`:**
- Purpose: Jupyter notebook checkpoints
- Generated: Yes
- Committed: No (in .gitignore)

**`.venv/`:**
- Purpose: Python virtual environment (local development)
- Generated: Yes (by `python -m venv .venv`)
- Committed: No (in .gitignore)

