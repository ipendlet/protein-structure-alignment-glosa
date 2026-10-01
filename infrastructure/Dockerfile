# syntax=docker/dockerfile:1

# Registry to pull the base images from, as a prefix with a trailing slash.  Empty -- the
# default -- means Docker Hub; compose passes a mirror in from deploy.conf on networks that
# cannot reach Hub.  Declared before the first FROM so both stages can use it.
ARG IMAGE_PREFIX=

# --------------------------------------------------------------------------- Java stage
# G-LoSA's chemical-feature assignment is a Java program and has to run on every upload, so a
# JVM is a runtime dependency, not just a build one.  A whole JDK is ~320 MB for two class
# files, and both only import java.io, java.text and java.util -- all of them in java.base -- so
# jlink cuts a ~50 MB runtime with that single module and the final image carries that instead.
FROM ${IMAGE_PREFIX}eclipse-temurin:21-jdk AS jvm

WORKDIR /build
COPY src/AssignChemicalFeatures.java src/AssignSecondaryStructures.java ./

# Separate output directories, not one javac invocation: both files declare their own top-level
# `Atom` class, so compiling them together fails with "duplicate class: Atom".
RUN mkdir -p /opt/glosa/classes/acf /opt/glosa/classes/ass && \
    javac -nowarn -d /opt/glosa/classes/acf AssignChemicalFeatures.java && \
    javac -nowarn -d /opt/glosa/classes/ass AssignSecondaryStructures.java

RUN jlink --add-modules java.base \
          --strip-debug --no-header-files --no-man-pages \
          --output /opt/jre

# --------------------------------------------------------------------------- runtime
# The full python image rather than -slim, for two reasons.  A build network restricted to a
# mirror usually cannot reach deb.debian.org, and such mirrors rarely carry a Debian repo, so
# apt cannot install anything; and the full image is built on buildpack-deps, so it already
# ships the g++ that the scorer needs.  On -slim there is no compiler and no way to add one.
FROM ${IMAGE_PREFIX}python:3.11

ARG PIP_INDEX_URL
ENV PIP_INDEX_URL="${PIP_INDEX_URL}"
ENV PYTHONUNBUFFERED=1
ENV PYTHONDONTWRITEBYTECODE=1

ENV VIRTUAL_ENV=/opt/venv
RUN python -m venv "$VIRTUAL_ENV" && \
    "$VIRTUAL_ENV/bin/pip" install --upgrade pip setuptools wheel
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

COPY --from=jvm /opt/jre /opt/jre
COPY --from=jvm /opt/glosa/classes /opt/glosa/classes
ENV JAVA_BIN=/opt/jre/bin/java

# Dependencies before source, so editing a template does not re-run the pip install.
COPY infrastructure/requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt

# Build the scorer.  -O2 matters here in a way it rarely does for glue code: the maximum-clique
# search over the product graph is the whole runtime, and the unoptimised binary is several
# times slower on anything larger than a binding site.
COPY src/glosa.cpp /tmp/glosa.cpp
RUN g++ -O2 -o /opt/glosa/glosa /tmp/glosa.cpp && rm /tmp/glosa.cpp

ENV GLOSA_BIN=/opt/glosa/glosa \
    GLOSA_CLASSES_ACF=/opt/glosa/classes/acf \
    GLOSA_CLASSES_ASS=/opt/glosa/classes/ass \
    GLOSA_EXAMPLES=/opt/glosa/example \
    GLOSA_JOBS_DIR=/data/jobs

COPY example /opt/glosa/example

WORKDIR /app
COPY infrastructure/glosa_runner.py infrastructure/app.py infrastructure/wsgi.py \
     infrastructure/smoke_test.py infrastructure/align_cli.py /app/
COPY infrastructure/templates /app/templates
COPY infrastructure/static /app/static

# Job outputs are the only thing written at runtime; compose mounts a volume here so a redeploy
# does not throw away results.
RUN mkdir -p /data/jobs
VOLUME ["/data"]

EXPOSE 8057

# One worker with threads.  Alignments run on a background thread pool inside the process and a
# whole-protein product graph can reach several GB, so a second worker doubles the worst-case
# memory for no gain -- the threads are what keep the UI responsive while a job runs.  The long
# timeout covers the upload, not the alignment, which is already off the request thread.
CMD ["gunicorn", \
     "--bind", "0.0.0.0:8057", \
     "--workers", "1", \
     "--threads", "8", \
     "--timeout", "300", \
     "--access-logfile", "-", \
     "--error-logfile", "-", \
     "wsgi:server"]
