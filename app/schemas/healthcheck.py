from __future__ import annotations

from pydantic import BaseModel


class HealthCheck(BaseModel):
    """Represents the response model for validating and returning the result of.

    a health check.

    Attributes:
        status (str): The status of the health check. The default is "OK".
    """

    status: str = "OK"
