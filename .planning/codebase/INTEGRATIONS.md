# External Integrations

**Analysis Date:** 2026-03-12

## APIs & External Services

**Chemical Structure Classification:**
- ClassyFire API - Classifies chemical compounds by structure
  - SDK/Client: `httpx` (async HTTP client in `app/modules/classyfire.py`)
  - Endpoint: `http://classyfire.wishartlab.com/queries/?format=json`
  - Timeout: 300 seconds (5 min) for classification + 60s connection timeout
  - Async interface: `classify(smiles)` returns classification response
  - Used by: OCSR router for compound classification

**Chemical Data Retrieval:**
- PubChem API - Retrieves chemical compound data and properties
  - SDK/Client: `requests` library with retry strategy in `app/modules/pubchem_retrieve.py`
  - Base URL: `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound`
  - Auth: None (public API)
  - Retry strategy: 3 retries with 0.5s backoff for status codes [429, 500, 502, 503, 504]
  - Timeout: 10 seconds (configurable in PubChemClient)
  - Cache: 128-entry LRU cache for queries
  - Used by: `chem` router for compound information endpoints

**Structure Generation:**
- SURGE Tool - Structure generation and assembly
  - Type: External binary tool (downloaded from GitHub releases at runtime)
  - Location: `/usr/bin/surge` (in Docker)
  - Version: v2.0 for Linux x86_64
  - Downloaded by: Dockerfile build process
  - Used by: `app/modules/tools/surge.py` for structure generation tasks

## Data Storage

**Databases:**
- Not detected - No persistent database integration. Stateless API design.

**File Storage:**
- Local filesystem only - Temporary files for image processing/conversion
  - Image uploads handled via FastAPI `UploadFile`
  - Processed images in temp directories (GIF→PNG conversion, PDFs)
  - No external S3 or cloud storage detected

**Caching:**
- LRU cache only - In-memory caching via Python `@lru_cache` decorator
  - PubChem client: 128-entry cache (configurable)
  - No Redis or Memcached integration

## Authentication & Identity

**Auth Provider:**
- Custom rate-limit based authentication (not full identity auth)
  - Implementation: Internal auth token via `CMS_INTERNAL_AUTH_TOKEN` env var
  - Location: `app/limiter.py`
  - Purpose: Bypass rate limiting for internal services
  - Mechanism: Requests with `x-internal-auth` header matching token share single rate-limit bucket
  - No user login/authentication system

**Rate Limiting:**
- slowapi - Request rate limiting library
  - Config: `app/limiter.py` with `slowapi.middleware.SlowAPIMiddleware`
  - Strategy: Token bucket with internal auth bypass
  - Per-router health endpoints can have independent rate limits

## Monitoring & Observability

**Error Tracking:**
- None detected - Error logging only via Python logger to stdout

**Logs:**
- Python logging to stdout
  - Format: `"%(asctime)s - %(name)s - %(levelname)s - %(message)s"`
  - Loggers: Named per module (e.g., `pubchem_client`, `ocsr`)
  - SSRF validation logging in OCSR router (`app/routers/ocsr.py`)
  - Async logging support for error cases

**Metrics & Monitoring:**
- Prometheus via prometheus-fastapi-instrumentator
  - Metrics endpoint: `/metrics`
  - Instrumentation: Automatic request/response metrics
  - Config: `Instrumentator().instrument(app).expose(app)` in `app/main.py`
- Grafana dashboards
  - Dashboard config: `cm-dashboard.json`
  - Admin password: `GRAFANA_ADMIN_PASSWORD` env var (required)
  - Data source: Prometheus (localhost:9090)
  - Port: 3001 (mapped from container 3000)
- Prometheus scrape config
  - Location: `prometheus_data/prometheus.yml`
  - Scrapes metrics from API and other services

## CI/CD & Deployment

**Hosting:**
- Docker Hub - `nfdi4chem/cheminformatics-microservice`
  - Full variant: `nfdi4chem/cheminformatics-microservice:latest`
  - Lite variant: `nfdi4chem/cheminformatics-microservice:latest-lite`
  - Frontend: `nfdi4chem/cheminformatics-microservice-frontend:latest`

**CI Pipeline:**
- GitHub Actions
  - Test workflow (`.github/workflows/test.yml`):
    - Trigger: PR/push to main or development branches
    - Python 3.10 on ubuntu-latest
    - Steps: Lint (flake8), test (pytest), coverage upload (Codecov)
    - Coverage token: `CODECOV_TOKEN` secret
    - Skips version files (CHANGELOG.md, package.json)
  - Production build workflow (`.github/workflows/prod-build.yml`):
    - Trigger: GitHub release published
    - Builds 3 Docker images: full, lite, frontend
    - Pushes to Docker Hub with version tags and `:latest`
    - Uses Docker credentials: `DOCKER_USERNAME`, `DOCKER_PASSWORD` secrets
  - Development build workflow (`.github/workflows/dev-build.yml`):
    - Pushes development builds to Docker Hub
  - Documentation deploy (`.github/workflows/deploy-doc.yml`):
    - Builds and deploys docs via VitePress
  - Release automation (`.github/workflows/release-please.yml`):
    - Auto-generates changelog and creates release PRs

**Codecov Integration:**
- Service: Codecov
  - Token: `CODECOV_TOKEN` secret in GitHub Actions
  - Reports: XML coverage format from pytest-cov
  - Flags: "service" for backend tests

## Environment Configuration

**Required env vars:**
- `GRAFANA_ADMIN_PASSWORD` - Grafana admin password (REQUIRED, no default, docker-compose will fail without it)
- `CMS_INTERNAL_AUTH_TOKEN` - Rate limit bypass token (32+ chars recommended)

**Optional env vars:**
- `INCLUDE_OCSR` (default "true") - Enable/disable OCSR module and TensorFlow loading
- `ALLOWED_ORIGINS` (default "*") - CORS allowed origins
- `RELEASE_VERSION` (default "1.0") - API version in OpenAPI schema
- `HOMEPAGE_URL` (default "/latest/docs") - Root path redirect
- `WORKERS` (default 2) - Uvicorn worker count
- `NODE_ENV` (default "development") - Frontend environment (development|production)
- `REACT_APP_API_URL` (default "http://localhost:8000/latest") - Frontend API base URL

**Secrets location:**
- `.env` file (git-ignored, not committed)
- `.env.development` - Development environment template
- `.env.production` - Production environment template
- `.env.example` - Template with all available variables
- Docker secrets: GitHub Actions secrets for DOCKER_USERNAME, DOCKER_PASSWORD, CODECOV_TOKEN

## Webhooks & Callbacks

**Incoming:**
- None detected - API is request-response only, no webhook listeners

**Outgoing:**
- GitHub Actions webhooks only (release-please for automated PRs)
- No third-party webhook notifications

## Health Checks

**Router-Level Health Endpoints:**
- `/v1/chem/health` - Chemical analysis module health
- `/v1/converters/health` - Format converters health
- `/v1/depict/health` - Visualization module health
- `/v1/tools/health` - Utilities health
- `/v1/ocsr/health` - OCSR module health (if INCLUDE_OCSR=true)
- `/latest/health` - Main API health (redirects to latest version)

**Docker Health Check:**
- Command: `curl -f http://localhost:8000/latest/chem/health || exit 1`
- Interval: 90 seconds
- Timeout: 10 seconds
- Retries: 20
- Start period: 60 seconds

## Security

**URL Validation (SSRF Protection):**
- OCSR router implements `_validate_url()` in `app/routers/ocsr.py`
  - Blocks private/internal IP addresses (127.0.0.1, 10.x.x.x, 172.16-31.x.x, 192.168.x.x)
  - Blocks localhost and .local domains
  - Allows only http/https schemes
  - Used for image URL processing endpoints

**Security Headers:**
- Custom middleware in `app/main.py` adds:
  - `X-Content-Type-Options: nosniff`
  - `X-Frame-Options: DENY`
  - `X-XSS-Protection: 1; mode=block`
  - `Referrer-Policy: strict-origin-when-cross-origin`

**Error Handling:**
- Generic error messages returned to clients
  - Details logged server-side only via Python logger
  - Invalid input exceptions return 422 (Unprocessable Entity)
  - Custom exception handler: `app/exception_handlers.py`

---

*Integration audit: 2026-03-12*
