from fastapi import FastAPI
from fastapi.responses import RedirectResponse
from fastapi_versioning import VersionedFastAPI

# from .config import settings
from .routers import chem, converters, decimer
from fastapi.middleware.cors import CORSMiddleware

app = FastAPI(
    title="Cheminf Micro Services",
    description="This set of essential and valuable microservices is designed to be accessed via API calls to support cheminformatics. Generally, it is designed to work with SMILES-based inputs and could be used to translate between different machine-readable representations, get Natural Product (NP) likeliness scores, visualize chemical structures, and generate descriptors. In addition, the microservices also host an instance of STOUT and another instance of DECIMER (two deep learning models for IUPAC name generation and optical chemical structure recognition, respectively).",
    terms_of_service="https://github.com/Steinbeck-Lab",
    contact={
        "name": "Steinbeck Lab",
        "url": "https://cheminf.uni-jena.de/",
        "email": "caffeine@listserv.uni-jena.de",
    },
    license_info={
        "name": "CC BY 4.0",
        "url": "https://creativecommons.org/licenses/by/4.0/",
    },
)

origins = ["*"]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

app.include_router(chem.router)
app.include_router(converters.router)
app.include_router(decimer.router)


@app.get("/", include_in_schema=False)
async def docs_redirect():
    return RedirectResponse(url="/docs")


app = VersionedFastAPI(
    app,
    version_format="{major}",
    prefix_format="/v{major}",
    terms_of_service="https://github.com/Steinbeck-Lab",
    contact={
        "name": "Steinbeck Lab",
        "url": "https://cheminf.uni-jena.de/",
        "email": "caffeine@listserv.uni-jena.de",
    },
    license_info={
        "name": "CC BY 4.0",
        "url": "https://creativecommons.org/licenses/by/4.0/",
    },
)
