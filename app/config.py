from pydantic import BaseSettings, Field, PostgresDsn


class Settings(BaseSettings):
    PGPASSWORD: str = Field(..., env="PGPASSWORD")
    POSTGRES_DB: str = Field(..., env="POSTGRES_DB")
    POSTGRES_USER: str = Field(..., env="POSTGRES_USER")
    POSTGRES_PASSWORD: str = Field(..., env="POSTGRES_PASSWORD")
    POSTGRES_PORT: str = Field(..., env="POSTGRES_PORT")


settings = Settings()
