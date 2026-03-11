import logging
import os

from slowapi import Limiter
from slowapi import _rate_limit_exceeded_handler
from slowapi.errors import RateLimitExceeded
from slowapi.util import get_remote_address
from starlette.requests import Request as StarletteRequest

logger = logging.getLogger(__name__)

INTERNAL_AUTH_TOKEN = os.getenv("CMS_INTERNAL_AUTH_TOKEN")

if not INTERNAL_AUTH_TOKEN:
    logger.warning(
        "CMS_INTERNAL_AUTH_TOKEN is not set — authenticated rate-limit bypass is disabled"
    )
elif len(INTERNAL_AUTH_TOKEN) < 32:
    logger.warning(
        "CMS_INTERNAL_AUTH_TOKEN is shorter than 32 characters — consider using a stronger token"
    )


def _is_authenticated(request: StarletteRequest) -> bool:
    """Return True when the request carries a valid internal auth token."""
    token = request.headers.get("x-internal-auth")
    return bool(token and INTERNAL_AUTH_TOKEN and token == INTERNAL_AUTH_TOKEN)


def custom_rate_limit_key_func(request: StarletteRequest) -> str:
    """Return a rate-limit key for the request.

    Authenticated requests get a unique key per request, which effectively
    exempts them from rate limiting.  Unauthenticated requests are keyed
    by IP address.
    """
    if _is_authenticated(request):
        return f"auth_exempt_{id(request)}"
    return get_remote_address(request)


limiter = Limiter(key_func=custom_rate_limit_key_func)
rate_limit_exceeded_handler = _rate_limit_exceeded_handler
RateLimitExceededExc = RateLimitExceeded
