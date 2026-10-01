#!/usr/bin/env bash
# Regenerate .env and, where the registry needs it, refresh the container registry login.
#
# The service holds no credentials of its own -- it reads nothing but the structures it is given
# -- so the only secret involved is the package-index credential, which is already sitting in
# the local pip config.  On a network served by an Artifactory or Nexus mirror the same
# credential authenticates the container registry, which is what REGISTRY_LOGIN_HOST turns on.
#
# On a machine that reaches PyPI and Docker Hub directly, this writes an empty .env and stops;
# compose treats a missing env_file as a hard error, so the file has to exist either way.
set -euo pipefail

cd "$(dirname "$0")"

CONFIG="${CONFIG:-deploy.conf}"
if [ ! -f "$CONFIG" ]; then
    echo "ERROR: $CONFIG is missing. Run: cp deploy.conf.example $CONFIG, then edit it." >&2
    exit 1
fi

# Read values out of deploy.conf rather than sourcing it: the file is also `include`d by make,
# so it holds values like NODE_MAP that contain spaces and would be parsed as commands by sh.
conf() { sed -n "s/^[[:space:]]*$1[[:space:]]*=[[:space:]]*//p" "$CONFIG" | tail -n 1; }

CONTAINER_CLI="$(conf CONTAINER_CLI)"
CONTAINER_CLI="${CONTAINER_CLI:-docker}"
REGISTRY_LOGIN_HOST="$(conf REGISTRY_LOGIN_HOST)"

# Some hosts have python3 but no bare `python` outside an activated venv, and this script is
# meant to be runnable either way.
PY="$(command -v python || command -v python3 || true)"
if [ -z "$PY" ]; then
    echo "ERROR: no python interpreter on PATH." >&2
    exit 1
fi

# An index URL in the pip config means a mirror; its absence means plain PyPI, and the build
# needs nothing passed in.
PIP_INDEX_URL="$("$PY" -m pip config list 2>/dev/null |
    awk -F= '/global.index-url/ {gsub("'"'"'","",$2); print $2; exit}')"

if [ -z "$PIP_INDEX_URL" ]; then
    : > .env
    echo "Wrote empty .env (no package-index mirror configured; the build will use PyPI)."
    exit 0
fi

echo "TEMP_PIP_INDEX_URL=$PIP_INDEX_URL" > .env
chmod 600 .env
echo "Wrote .env"

if [ -z "$REGISTRY_LOGIN_HOST" ]; then
    exit 0
fi

export PIP_INDEX_URL
read -r REGISTRY_USER REGISTRY_PASS <<<"$("$PY" - <<'PY'
import os
from urllib.parse import urlparse, unquote

parsed = urlparse(os.environ["PIP_INDEX_URL"].strip())
if not parsed.username or parsed.password is None:
    raise SystemExit("PIP_INDEX_URL carries no credentials; unset REGISTRY_LOGIN_HOST or "
                     "log in to the registry by hand.")
print(unquote(parsed.username), unquote(parsed.password))
PY
)"
printf '%s' "$REGISTRY_PASS" |
    "$CONTAINER_CLI" login "$REGISTRY_LOGIN_HOST" -u "$REGISTRY_USER" --password-stdin >/dev/null
echo "Registry login refreshed for $REGISTRY_LOGIN_HOST."
