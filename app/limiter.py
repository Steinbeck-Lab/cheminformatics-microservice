

import os
from slowapi import Limiter, _rate_limit_exceeded_handler
from slowapi.util import get_remote_address
from slowapi.errors import RateLimitExceeded
from starlette.requests import Request as StarletteRequest

INTERNAL_AUTH_TOKEN = os.getenv("CMS_INTERNAL_AUTH_TOKEN")

def custom_rate_limit_key_func(request: StarletteRequest):
    token = request.headers.get("x-internal-auth")
    if token and INTERNAL_AUTH_TOKEN and token == INTERNAL_AUTH_TOKEN:
        return None  # disables rate limiting for this request
    return get_remote_address(request)

limiter = Limiter(key_func=custom_rate_limit_key_func)
rate_limit_exceeded_handler = _rate_limit_exceeded_handler
RateLimitExceededExc = RateLimitExceeded
