import os
from pathlib import Path

# Base directory for the application (where the app is running)
BASE_DIR = Path(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Directory for storing uploaded files and their processed output
UPLOAD_DIR = os.path.join(BASE_DIR, "uploads")
PDF_DIR = os.path.join(UPLOAD_DIR, "pdfs")
SEGMENTS_DIR = os.path.join(UPLOAD_DIR, "segments")
IMAGES_DIR = os.path.join(UPLOAD_DIR, "chem_images")

# Create directories if they don't exist
os.makedirs(PDF_DIR, exist_ok=True)
os.makedirs(SEGMENTS_DIR, exist_ok=True)
os.makedirs(IMAGES_DIR, exist_ok=True)

# Maximum upload file size (10 MB)
MAX_UPLOAD_SIZE = 10 * 1024 * 1024  # in bytes
