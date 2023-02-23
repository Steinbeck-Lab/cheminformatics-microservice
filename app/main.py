from fastapi import Depends, FastAPI
from .config import settings
from .routers import converters, chem, compose
from fastapi.middleware.cors import CORSMiddleware

app = FastAPI()

origins = [
    "http://localhost",
    "http://localhost:8080",
    "http://localhost:80",
]

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