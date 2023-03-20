from fastapi import FastAPI
from fastapi.responses import RedirectResponse

# from .config import settings
from .routers import converters, chem, compose, decimer
from fastapi.middleware.cors import CORSMiddleware
from fastapi.openapi.utils import get_openapi

app = FastAPI()

origins = ["*"]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

app.include_router(converters.router)
app.include_router(chem.router)
app.include_router(compose.router)
app.include_router(decimer.router)


@app.get("/", include_in_schema=False)
async def docs_redirect():
    return RedirectResponse(url="/docs")


def custom_openapi():
    if app.openapi_schema:
        return app.openapi_schema
    openapi_schema = get_openapi(
        title="Cheminf Micro Services",
        version="0.3.0",
        description="This set of essential and valuable microservices is designed to be accessed via API calls to support cheminformatics. Generally, it is designed to work with SMILES-based inputs and could be used to translate between different machine-readable representations, get Natural Product (NP) likeliness scores, visualize chemical structures, and generate descriptors. In addition, the microservices also host an instance of STOUT and another instance of DECIMER (two deep learning models for IUPAC name generation and optical chemical structure recognition, respectively).",
        routes=app.routes,
    )
    openapi_schema["info"]["x-logo"] = {
        "url": "https://github.com/Steinbeck-Lab/cheminformatics-python-microservice/raw/main/public/img/logo.png"
    }
    app.openapi_schema = openapi_schema
    return app.openapi_schema


app.openapi = custom_openapi
