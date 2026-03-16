# Architecture

**Analysis Date:** 2026-03-12

## Pattern Overview

**Overall:** Multi-router FastAPI microservice with parallel chemistry toolkit abstraction layer.

**Key Characteristics:**
- FastAPI-based REST API with versioned routing via `VersionedFastAPI` (routes served at `/v1/...` and `/latest/...`)
- Input validation and parsing abstraction in `helpers.parse_input()` dispatches to framework-specific toolkit wrappers
- Three chemistry toolkits (RDKit, CDK via JPype, OpenBabel) run in parallel with framework parameter-based selection
- Pydantic schema-based request/response validation per endpoint
- Custom exception handler for invalid input (`InvalidInputException` → 422 response)
- Rate limiting via `slowapi` with internal auth token bypass support
- Middleware chain for security headers, CORS, and monitoring instrumentation

## Layers

**HTTP/API Layer:**
- Purpose: Handle incoming HTTP requests, parse query parameters/body, validate with Pydantic schemas, return JSON responses
- Location: `app/routers/` (chem.py, converters.py, depict.py, tools.py, ocsr.py)
- Contains: APIRouter instances with endpoint handlers decorated with `@router.get()` / `@router.post()`
- Depends on: Pydantic schemas (`app/schemas/`), module functions, limiter, exception handlers
- Used by: FastAPI app (`app/main.py`), clients (frontend SPA, REST clients)

**Validation & Parsing Layer:**
- Purpose: Accept user input (SMILES, molecular formula, file uploads), validate format, normalize, dispatch to appropriate toolkit
- Location: `app/modules/toolkits/helpers.py:parse_input()` and `parse_SMILES()`
- Contains: `parse_input()` dispatches on `framework` parameter; falls back to CDK if RDKit fails
- Depends on: Toolkit wrappers (rdkit_wrapper, cdk_wrapper, openbabel_wrapper)
- Used by: All routers for molecule input processing

**Toolkit Abstraction Layer:**
- Purpose: Provide a unified interface to three chemistry engines; handle framework selection and cross-toolkit conversions
- Location: `app/modules/toolkits/` (rdkit_wrapper.py, cdk_wrapper.py, openbabel_wrapper.py)
- Contains: Functions like `get_CDK_IAtomContainer()`, `Chem.MolFromSmiles()` (RDKit), `get_ob_mol()` (OpenBabel); JPype JVM initialization for CDK
- Depends on: RDKit (Python), CDK JARs (loaded via JPype + pystow), OpenBabel (compiled C++ library)
- Used by: Validation layer, descriptor/depiction/conversion modules

**Chemistry Logic Layer:**
- Purpose: Perform chemical computations (descriptors, depictions, conversions, analysis)
- Location: `app/modules/` (all_descriptors.py, depiction.py, depiction_enhanced.py, depiction_cdk_depict/, classyfire.py, coconut/, decimer.py, fix_radicals.py, npscorer.py, pubchem_retrieve.py, tools/)
- Contains: Domain-specific functions (e.g., `get_tanimoto_similarity()`, `get_rdkit_depiction()`, `get_cdk_depiction()`, descriptor calculation, OCSR via DECIMER)
- Depends on: Toolkit wrappers, external services (PubChem API, ClassyFire API), TensorFlow (DECIMER, conditionally loaded via `INCLUDE_OCSR` env var)
- Used by: Routers

**Cross-Cutting Concerns:**
- Purpose: Rate limiting, error handling, monitoring, security
- Location: `app/limiter.py` (rate limiting), `app/exception_handlers.py` (custom exceptions), `app/main.py` (middleware chain, security headers)
- Contains: SlowAPI limiter with custom key function (`custom_rate_limit_key_func()`), `InvalidInputException`, `SecurityHeadersMiddleware`, CORS middleware, Prometheus instrumentation
- Used by: All layers

## Data Flow

**Standard Endpoint Data Flow:**

1. HTTP request arrives at FastAPI router endpoint (e.g., `/chem/properties`)
2. Pydantic auto-validates query/body parameters against schema in `app/schemas/*.py`
3. If validation fails, FastAPI returns 422; if invalid input is raised, caught by `InvalidInputException` handler → 422 with sanitized message
4. Router calls `parse_input(smiles, framework)` from `app/modules/toolkits/helpers.py`
5. `parse_input()` dispatches to `parse_SMILES()` → calls toolkit wrapper based on `framework` param:
   - `framework="rdkit"`: calls RDKit's `Chem.MolFromSmiles()`; if fails and molecule has R groups, tries CDK fallback
   - `framework="cdk"`: calls `get_CDK_IAtomContainer()`
   - `framework="openbabel"`: calls `get_ob_mol()`
6. Parsed molecule object returned to router
7. Router passes molecule to chemistry logic function (e.g., `get_tanimoto_similarity()`, `get_rdkit_depiction()`)
8. Logic function performs computation using toolkit-specific APIs
9. Result serialized to Pydantic response schema
10. FastAPI returns 200 + JSON response

**Example: `/depict/2D` endpoint:**
1. Request: `GET /v1/depict/2D?smiles=CC&depictionType=fixed`
2. Parse SMILES via `parse_input("CC", framework="rdkit")`
3. Call `get_cdk_depiction()` or `get_rdkit_depiction()` based on `depictionType`
4. Depiction module imports CDK-specific depiction modules from `app/modules/cdk_depict/`
5. Returns SVG string or image bytes in `Depict2DResponse` schema
6. Frontend receives and renders image

**Error Flow:**
- Pydantic validation errors → FastAPI 422 response
- Invalid SMILES/input → `parse_input()` raises `InvalidInputException(name, value)` → caught by `input_exception_handler()` → 422 with message `{"detail": "Invalid input for 'smiles'"}`
- Toolkit/external service failures → caught, logged, returned as HTTP error (400, 500, etc.) with generic message to client

**State Management:**
- Stateless HTTP: each request is independent
- Rate limiter maintains per-IP request counters (Redis if distributed)
- OCSR model (DECIMER) loaded once at startup if `INCLUDE_OCSR=true`
- CDK JARs downloaded once via `pystow` on first startup

## Key Abstractions

**Molecule Object:**
- Purpose: Represents a chemical structure in toolkit-agnostic way before dispatch to specific toolkit
- Examples: `rdkit.Chem.Mol`, `jpype._jproxy.Jpype...` (CDK IAtomContainer), `openbabel.openbabel.OBMol`
- Pattern: Each toolkit provides its own molecule object type; parsed in `helpers.py` and passed to toolkit-specific functions in modules

**Toolkit Wrapper Functions:**
- Purpose: Encapsulate toolkit-specific APIs into reusable functions
- Examples: `rdkit_wrapper.get_tanimoto_similarity_rdkit()`, `cdk_wrapper.get_CDK_HOSE_codes()`, `openbabel_wrapper.get_ob_InChI()`
- Pattern: Function name includes toolkit name; takes toolkit-specific molecule object; returns data structures or raises exceptions

**Depiction Modules:**
- Purpose: Generate 2D/3D visual representations (SVG, PNG)
- Examples: `depiction.py:get_rdkit_depiction()` (RDKit-based 2D), `depiction_enhanced.py:get_cdk_depiction()` (CDK with feature plugins from `cdk_depict/`)
- Pattern: Modular feature plugins in `cdk_depict/` (abbreviations, annotations, aromatic_display, etc.) composed by `depiction_enhanced.py`

**Response Schemas:**
- Purpose: Define and validate API response structure via Pydantic
- Examples: `chem_schema.GenerateDescriptorsResponse`, `depict_schema.Depict2DResponse`, `converters_schema.GenerateInChIResponse`
- Pattern: One model per endpoint response, organized by router in `app/schemas/`

## Entry Points

**HTTP Entry Point:**
- Location: `app/main.py`
- Triggers: FastAPI startup (e.g., `uvicorn app.main:app --reload`)
- Responsibilities: Create base FastAPI app, include all routers, wrap with VersionedFastAPI, register middleware (security headers, CORS, rate limiting), register exception handlers, expose `/health` endpoint

**Routers (Secondary Entry Points):**
- Location: `app/routers/chem.py`, `converters.py`, `depict.py`, `tools.py`, `ocsr.py` (conditional on `INCLUDE_OCSR`)
- Triggers: HTTP requests to `/v1/<prefix>/...` or `/latest/<prefix>/...`
- Responsibilities: Define endpoint handlers, call validation logic, dispatch to modules, return serialized responses

**Python REPL/Script Entry Point:**
- Direct import: `from app.modules.toolkits.rdkit_wrapper import get_properties`
- Used in tests and notebooks for unit testing modules in isolation

## Error Handling

**Strategy:** Exceptions are caught and logged server-side; generic messages returned to clients. Input validation errors map to 422; business logic errors to 400/500.

**Patterns:**

1. **Pydantic Validation Error**: Query/body param fails schema validation → FastAPI auto-returns 422 with validation detail
   - Example: Missing required query param `smiles` → 422 with "field required"

2. **InvalidInputException**: Custom exception raised in `parse_input()` for malformed SMILES/input
   - Raised in: `app/modules/toolkits/helpers.py:parse_SMILES()`
   - Caught by: `input_exception_handler()` in `app/exception_handlers.py`
   - Response: 422 with message `{"detail": "Invalid input for '{field}'"}`
   - Logged: Server-side as warning with actual value

3. **HTTPException**: Manual error responses in router endpoints
   - Example: `raise HTTPException(status_code=400, detail="Invalid framework")`
   - Used for business logic errors (unsupported parameter values, external service failures)

4. **Unhandled Exceptions**: Caught by FastAPI default handler
   - Logged: Server-side traceback
   - Response: 500 with generic message (details hidden from client)

5. **Rate Limit Exceeded**: `slowapi.errors.RateLimitExceeded`
   - Caught by: `rate_limit_exceeded_handler()` in `app/limiter.py`
   - Response: 429 "Too Many Requests"
   - Bypass: Requests with header `x-internal-auth: {CMS_INTERNAL_AUTH_TOKEN}` exempt

## Cross-Cutting Concerns

**Logging:**
- Framework: Python `logging` module
- Pattern: Each module defines `logger = logging.getLogger(__name__)`
- Usage: Log warnings in exception handlers, errors in toolkit wrappers, info in endpoints
- Example: `logger.warning("Invalid input for %s: %s", exc.name, exc.value)` in `exception_handlers.py`

**Validation:**
- Mechanism: Pydantic BaseModel schemas with Field descriptors
- Input validation: Query/body parameters validated automatically by FastAPI before handler called
- Business validation: Custom logic in endpoints (e.g., `_validate_url()` in ocsr.py for SSRF protection)
- Invalid input exception: Raised in parsing layer, caught by central exception handler

**Authentication & Authorization:**
- Mechanism: `x-internal-auth` header matching `CMS_INTERNAL_AUTH_TOKEN` env var
- Purpose: Exempt internal service-to-service calls from rate limiting
- Pattern: `_is_authenticated()` in `app/limiter.py` checks header; `custom_rate_limit_key_func()` exempts authenticated requests by giving them unique key per request
- No role-based access control; all endpoints public or protected by token

**Rate Limiting:**
- Framework: `slowapi` (wrapper around `requests-limiter`)
- Pattern: Decorated endpoints with `@limiter.limit("10/minute")`
- Key function: `custom_rate_limit_key_func()` — returns IP address for unauthenticated, unique key for authenticated
- Registered: Applied via `SlowAPIMiddleware` in main.py

**Security Headers:**
- Middleware: `SecurityHeadersMiddleware` in `app/main.py`
- Headers added:
  - `X-Content-Type-Options: nosniff` — prevent MIME type sniffing
  - `X-Frame-Options: DENY` — prevent clickjacking
  - `X-XSS-Protection: 1; mode=block` — enable browser XSS filter
  - `Referrer-Policy: strict-origin-when-cross-origin` — control referrer info

**CORS:**
- Framework: `CORSMiddleware` from `fastapi.middleware.cors`
- Configuration: Origins from `ALLOWED_ORIGINS` env var (default `"*"`); if `*`, credentials disabled; otherwise allowed
- Usage: Frontend SPA can call API from different origin

**Monitoring & Observability:**
- Prometheus metrics: `prometheus-fastapi-instrumentator` auto-instruments all endpoints
- Metrics exposed at: `/metrics` endpoint (FastAPI default)
- Grafana dashboard: `cm-dashboard.json` in root (requires `GRAFANA_ADMIN_PASSWORD` env var)
- Structured logging: Logs sent to stdout; can be collected by container orchestration

**SSRF Protection:**
- Function: `_validate_url()` in `app/routers/ocsr.py`
- Blocks: Private/loopback IP addresses (127.0.0.1, 192.168.x.x, etc.), DNS names (localhost, *.local), non-HTTP(S) schemes
- Used in: OCSR endpoint that fetches images from user-provided URLs

