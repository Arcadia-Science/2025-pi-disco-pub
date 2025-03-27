from dotenv import load_dotenv

import os

load_dotenv()  # load environment variables from .env file

OUTPUT_DIR = "outputs"

# API constants
MAX_RETRIES = 5
BACKOFF_FACTOR = 0.1
STATUS_FORCELIST = [429, 500, 502, 503, 504]

# ESMFold constants
API_URL = os.getenv("API_URL")
BEARER_TOKEN = os.getenv("BEARER_TOKEN")
HEADERS = {
    "Content-Type": "application/json",
    "Authorization": f"Bearer {BEARER_TOKEN}",
}

# AWS credential constants
AWS_ACCESS_KEY_ID = os.getenv("AWS_ACCESS_KEY_ID")
AWS_SECRET_ACCESS_KEY = os.getenv("AWS_SECRET_ACCESS_KEY")
AWS_REGION = os.getenv("AWS_REGION")
