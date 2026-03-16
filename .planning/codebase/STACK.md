# Technology Stack

**Analysis Date:** 2026-03-12

## Languages

**Primary:**
- Python 3.11 - Backend API and cheminformatics logic (`app/`)
- JavaScript/TypeScript 4.9.5 - Frontend React application (`frontend/src/`)

**Secondary:**
- Java - Chemistry Toolkit via CDK (Java-based, bridged via JPype)

## Runtime

**Environment:**
- Python 3.11 (via Conda miniconda3:24.1.2-0 in Docker)
- Node.js 18.20.8 (frontend)
- Java 11 OpenJDK (for CDK toolkit)

**Package Manager:**
- pip (Python) - primary dependency manager
- npm 10.8.2 (frontend)
- Conda (development environment)

## Frameworks

**Core Backend:**
- FastAPI >=0.110.3 - REST API framework with automatic OpenAPI/Swagger docs
- fastapi-versioning >=0.10.0 - API versioning at `/v{major}/` and `/latest/` paths
- Uvicorn with Gunicorn - ASGI server with configurable workers (default 2)

**Frontend:**
- React 18.2.0 - UI framework
- React Router v6.30.3 - Client-side routing
- React Scripts 5.0.1 - Development tooling and build pipeline

**Testing:**
- pytest - Python test runner
- pytest-asyncio - Async test support
- pytest-cov - Coverage reporting
- @testing-library/react, @testing-library/jest-dom - Frontend testing utilities

**Build/Dev:**
- black v22.6.0 - Python code formatter
- flake8 - Python linter with custom ignores (E402, E501, W503, E203; F401 ignored in `__init__.py`)
- reorder-python-imports - Import ordering tool
- pydocstringformatter - Docstring formatting
- prettier 3.6.2 - Frontend code formatter
- pre-commit - Git hooks for code quality
- Tailwind CSS 3.3.5 - Utility-first CSS framework (frontend)

## Key Dependencies

**Chemistry Toolkits:**
- rdkit - RDKit chemistry toolkit (primary, Python-native)
- jpype1 ==1.4.1 - Java/Python bridge for CDK toolkit via JPype JVM
- openbabel-wheel - OpenBabel chemistry toolkit (wheel distribution)
- pystow >=0.4.9 - Automatic download/caching of Java JAR files (CDK, SugarRemovalUtility, OPSIN)

**Deep Learning (OCSR):**
- tensorflow ==2.15.1 - TensorFlow for DECIMER deep learning model
- Keras-Preprocessing ==1.1.2 - TensorFlow preprocessing utilities
- decimer ==2.7.1 - DECIMER model package (optional via INCLUDE_OCSR env var)
- decimer-image-segmentation - Image segmentation for OCSR (installed from git)

**Image Processing:**
- opencv-python ==4.8.1.78 - OpenCV for image manipulation
- pillow ==12.1.1 - PIL for image handling
- pillow-heif ==0.10.0 - HEIF format support (iPhone photos)
- pdf2image ==1.16.2 - PDF to image conversion
- scikit-image - Scientific image processing
- matplotlib >=3.4.3 - Plotting and visualization

**HTTP & API Clients:**
- httpx >=0.24.1 - Async HTTP client for ClassyFire and OCSR endpoints
- requests - Synchronous HTTP client for PubChem API
- websockets >=14.0 - WebSocket support

**Data Processing:**
- pandas - Data manipulation
- chembl_structure_pipeline - ChEMBL chemistry data processing
- mapchiral - Chiral compound utilities
- selfies >=2.1.1 - SELFIES molecular representation library
- unicodedata2 ==15.0.0 - Unicode data utilities

**Feature Extraction:**
- efficientnet - EfficientNet neural network models
- imantics ==0.1.12 - Instance segmentation utilities

**Monitoring & Rate Limiting:**
- prometheus-fastapi-instrumentator - Prometheus metrics instrumentation for FastAPI
- slowapi - Rate limiting library (uses `slowapi.middleware.SlowAPIMiddleware`)

**Utilities:**
- fastapi-pagination ==0.10.0 - Pagination support for list endpoints
- python-multipart - Form data and file upload parsing
- pydantic - Data validation and serialization
- jinja2 - Template engine
- PyMuPDF (pymupdf) - PDF manipulation

**Structure Generation:**
- HOSE_code_generator (from git) - HOSE code generation from GitHub
- surge - Structure generation tool (downloaded as binary at runtime)

## Configuration

**Environment:**
- Environment variables via `.env` file (loaded by docker-compose)
- INCLUDE_OCSR (default "true") - Controls TensorFlow/DECIMER loading
- CMS_INTERNAL_AUTH_TOKEN - Rate limit bypass token (32+ chars recommended)
- ALLOWED_ORIGINS (default "*") - CORS allowed origins (comma-separated)
- GRAFANA_ADMIN_PASSWORD - Required for Grafana container (no default)
- RELEASE_VERSION (default "1.0") - API version identifier
- HOMEPAGE_URL (default "/latest/docs") - Root redirect destination
- WORKERS (default 2) - Uvicorn worker count
- JAVA_HOME - JVM path for CDK toolkit (auto-detected in Docker)

**Build:**
- `requirements.txt` - Python dependencies (pip format)
- `frontend/package.json` + `package-lock.json` - Frontend dependencies
- `.pre-commit-config.yaml` - Pre-commit hook configuration
- `Dockerfile` - Full backend image (includes OCSR/TensorFlow)
- `Dockerfile.lite` - Lightweight backend image (INCLUDE_OCSR=false)
- `docker-compose.yml` - Multi-service orchestration (api, web, prometheus, grafana)

## Platform Requirements

**Development:**
- Python 3.10+ (CI runs 3.10, local dev uses 3.11)
- OpenJDK 11 JRE (for CDK toolkit)
- Node.js 18+ (frontend)
- Docker & docker-compose (for containerized local dev)
- git (for pre-commit and installation of packages from GitHub)

**Production:**
- Docker/container platform (published to Docker Hub: `nfdi4chem/cheminformatics-microservice`)
- Image variants:
  - `latest` / `{version}` - Full image with OCSR
  - `latest-lite` / `{version}-lite` - Lightweight without TensorFlow
  - `frontend` - React SPA served on port 3000
- Monitoring: Prometheus (metrics at `/metrics`) and Grafana (dashboards)
- Minimum 4GB JVM heap (`-Xmx4096M` in CDK wrapper)

---

*Stack analysis: 2026-03-12*
