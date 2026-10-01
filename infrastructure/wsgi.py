"""WSGI entry point for the containerised G-LoSA service.

Kept separate from `app.py` so gunicorn imports a module that does nothing but expose the Flask
object, and so a stale `app.run()` can never start the development server in the image.
"""

from __future__ import annotations

import glosa_runner
from app import app as server

print(f"glosa binary : {glosa_runner.GLOSA_BIN} (present={glosa_runner.GLOSA_BIN.exists()})")
print(f"java runtime : {glosa_runner.JAVA_BIN}")
print(f"feature cls  : {glosa_runner.CLASSES_ACF} (present={glosa_runner.CLASSES_ACF.exists()})")
