import os

from slowapi import Limiter, _rate_limit_exceeded_handler
from slowapi.errors import RateLimitExceeded
from slowapi.util import get_remote_address
from starlette.requests import Request as StarletteRequest

INTERNAL_AUTH_TOKEN = os.getenv("CMS_INTERNAL_AUTH_TOKEN")


def custom_rate_limit_key_func(request: StarletteRequest):
    """
    Custom key function for rate limiting.
    - Authenticated requests get unique keys to effectively bypass rate limiting
    - Other requests are rate limited by IP address
    """
    token = request.headers.get("x-internal-auth")

    # Check for valid authentication token
    if token and INTERNAL_AUTH_TOKEN and token == INTERNAL_AUTH_TOKEN:
        # Return a unique key for each authenticated request
        # This effectively bypasses rate limiting since each request gets its own bucket
        return f"auth_bypass_{id(request)}"
    # For unauthenticated requests, use IP-based rate limiting
    return get_remote_address(request)


limiter = Limiter(key_func=custom_rate_limit_key_func)
rate_limit_exceeded_handler = _rate_limit_exceeded_handler
RateLimitExceededExc = RateLimitExceeded
