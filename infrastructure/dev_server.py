"""Run the service on localhost with Flask's development server, for work outside the container.

    make serve                       # from the package root
    python infrastructure/dev_server.py

`java` has to be on PATH for the chemical-feature step; prepend the JDK's bin directory if the
only one installed is bundled with something else.

Port and host come from GLOSA_DEV_PORT / GLOSA_DEV_HOST rather than the command line so the
invocation stays flag-free.  The image runs gunicorn against `wsgi:server` instead; this server
is single-process and not meant to be exposed.
"""

from __future__ import annotations

import os

from app import app

if __name__ == "__main__":
    app.run(
        host=os.environ.get("GLOSA_DEV_HOST", "127.0.0.1"),
        port=int(os.environ.get("GLOSA_DEV_PORT", "8057")),
        debug=False,
        threaded=True,
    )
