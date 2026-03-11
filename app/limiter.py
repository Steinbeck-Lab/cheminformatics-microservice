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


def custom_rate_limit_key_func(request: StarletteRequest):
    """Custom key function for rate limiting.

    Authenticated requests share a single bucket so the overall rate limit
    still applies.  Unauthenticated requests are rate-limited per IP.
    """
    token = request.headers.get("x-internal-auth")

    if token and INTERNAL_AUTH_TOKEN and token == INTERNAL_AUTH_TOKEN:
        return "auth_bypass"
    return get_remote_address(request)


limiter = Limiter(key_func=custom_rate_limit_key_func)
rate_limit_exceeded_handler = _rate_limit_exceeded_handler
RateLimitExceededExc = RateLimitExceeded
