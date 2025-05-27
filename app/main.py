from __future__ import annotations

import os
from fastapi import FastAPI
from fastapi import status
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import RedirectResponse
from fastapi_versioning import VersionedFastAPI
from prometheus_fastapi_instrumentator import Instrumentator
from typing import Dict
from fastapi_mcp import FastApiMCP

from .routers import chem
from .routers import converters
from .routers import depict
from .routers import tools
from app.exception_handlers import input_exception_handler
from app.exception_handlers import InvalidInputException
from app.schemas import HealthCheck


def create_app_metadata() -> Dict:
    return {
        "title": "Cheminformatics API",
        "description": """This set of essential and valuable microservices is designed to be accessed via API calls to support cheminformatics. Generally, it is designed to work with SMILES-based inputs and could be used to translate between different machine-readable representations, get Natural Product (NP) likeliness scores, visualize chemical structures, and generate descriptors. In addition, the microservices also host an instance of DECIMER (a deep learning model for optical chemical structure recognition).
        """,
        "version": "3.0.2",
        "terms_of_service": "https://docs.api.naturalproducts.net",
        "contact": {
            "name": "Steinbeck Lab",
            "url": "https://cheminf.uni-jena.de",
            "email": "caffeine@listserv.uni-jena.de",
        },
        "license_info": {
            "name": "CC BY 4.0",
            "url": "https://creativecommons.org/licenses/by/4.0",
        },
        "openapi_tags": [
            {
                "name": "Chemical Analysis",
                "description": "Advanced molecular structure analysis and property prediction",
            },
            {
                "name": "Visualization",
                "description": "High-fidelity chemical structure rendering and visualization",
            },
            {
                "name": "Conversion",
                "description": "Sophisticated molecular format transformation tools",
            },
        ],
    }


# Create base FastAPI app
app = FastAPI(**create_app_metadata())

# Add CORS middleware first
origins = ["*"]
app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Initialize FastAPI-MCP if enabled
mcp = None
if os.getenv("ENABLE_MCP", "false").lower() == "true":
    mcp = FastApiMCP(app, name="cheminformatics")
    mcp.mount()

# Add routers
app.include_router(chem.router)
app.include_router(converters.router)
app.include_router(depict.router)
app.include_router(tools.router)

# Import OCSR router if necessary
if os.getenv("INCLUDE_OCSR", "true").lower() == "true":
    from .routers import ocsr

    app.include_router(ocsr.router)

# Setup MCP server if enabled
if mcp is not None:
    mcp.setup_server()

# Add instrumentation
Instrumentator().instrument(app).expose(app)

# Add custom CORS headers middleware
@app.middleware("http")
async def add_cors_headers(request, call_next):
    response = await call_next(request)
    response.headers["Access-Control-Allow-Origin"] = "*"
    response.headers["Access-Control-Allow-Credentials"] = "true"
    response.headers["Access-Control-Allow-Methods"] = "GET, POST, PUT, DELETE, OPTIONS"
    response.headers["Access-Control-Allow-Headers"] = "Content-Type, Authorization"
    return response


# register exception handlers
for sub_app in app.routes:
    if hasattr(sub_app.app, "add_exception_handler"):
        sub_app.app.add_exception_handler(
            InvalidInputException,
            input_exception_handler,
        )


@app.get("/", include_in_schema=False)
async def root():
    return RedirectResponse(url=os.getenv("HOMEPAGE_URL", "/latest/docs"))


@app.get(
    "/health",
    tags=["healthcheck"],
    summary="Perform a Health Check",
    response_description="Return HTTP Status Code 200 (OK)",
    status_code=status.HTTP_200_OK,
    response_model=HealthCheck,
)
def get_health() -> HealthCheck:
    """## Perform a Health Check.

    Endpoint to perform a health check on. This endpoint can primarily be used by Docker
    to ensure a robust container orchestration and management are in place. Other
    services that rely on the proper functioning of the API service will not deploy if this
    endpoint returns any other HTTP status code except 200 (OK).
    Returns:
        HealthCheck: Returns a JSON response with the health status
    """
    return HealthCheck(status="OK")
