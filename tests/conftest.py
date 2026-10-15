from __future__ import annotations

import os

# DECIMER downloads several GB of models from Zenodo when it is imported, which
# made CI take over an hour. With INCLUDE_OCSR=false the app leaves the OCSR
# router out (see app/main.py), so skip the tests that need DECIMER as well.
if os.getenv("INCLUDE_OCSR", "true").lower() != "true":
    collect_ignore = ["test_decimer.py", "test_ocsr.py", "test_deeplearningtools.py"]
