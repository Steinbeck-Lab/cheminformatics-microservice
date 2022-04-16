from fastapi import Depends, FastAPI

from .routers import converters, chem

app = FastAPI()


app.include_router(converters.router)
app.include_router(chem.router)

@app.get("/")
async def root():
    return {"message": "Hello World!"}
