# Build and deploy the G-LoSA alignment service.
#
#   make            → build + up + url   (first-time deploy)
#   make redeploy   → cached rebuild + replace the container (after a code change)
#   make down       → stop and remove it
#
# Everything site-specific -- registry, routing domain, per-machine names, limits -- comes from
# deploy.conf, which is not in the repository.  See deploy.conf.example.

CONFIG ?= deploy.conf

ifeq ($(wildcard $(CONFIG)),)
$(error $(CONFIG) is missing. Run: cp deploy.conf.example $(CONFIG), then edit it)
endif
include $(CONFIG)

SERVICE_NAME            ?= glosa
GLOSA_PORT              ?= 8657
CONTAINER_CLI           ?= docker
COMPOSE_CLI             ?= $(CONTAINER_CLI) compose
TRAEFIK_NETWORK         ?= traefik
TRAEFIK_ENTRYPOINT      ?= web
TRAEFIK_RELAY_CONTAINER ?= traefik_relay

ifeq ($(strip $(ROUTING_DOMAIN)),)
$(error ROUTING_DOMAIN is not set in $(CONFIG))
endif

# Short name this machine is known by in the routing mesh.  `hostname` returns a bare name on
# some hosts and the FQDN on others, so strip any domain before matching.  A host missing from a
# non-empty NODE_MAP resolves to UNMAPPED and is rejected rather than guessed at: the wildcard
# record forwards *.<NODENAME>.<domain> to whichever machine claims that name, so a wrong value
# publishes this deployment under another machine's hostnames.
HOST_SHORT := $(firstword $(subst ., ,$(shell hostname)))

ifeq ($(strip $(NODENAME)),)
ifneq ($(strip $(NODE_MAP)),)
NODENAME := $(or $(patsubst $(HOST_SHORT)=%,%,$(filter $(HOST_SHORT)=%,$(NODE_MAP))),UNMAPPED)
endif
endif

ifeq ($(strip $(PUBLIC_HOSTNAME)),)
PUBLIC_HOSTNAME := $(SERVICE_NAME)$(if $(strip $(NODENAME)),.$(NODENAME)).$(ROUTING_DOMAIN)
endif

URL := http://$(PUBLIC_HOSTNAME)/

# Prefix, not host, so an empty DOCKER_REGISTRY_HOST means Docker Hub rather than a malformed
# image reference with a leading slash.
IMAGE_PREFIX := $(if $(strip $(DOCKER_REGISTRY_HOST)),$(strip $(DOCKER_REGISTRY_HOST))/,)

# docker-compose.yml reads all of these from the environment, so it never carries a site value
# of its own and `docker compose` run by hand gets the same settings make would have used.
export SERVICE_NAME GLOSA_PORT PUBLIC_HOSTNAME IMAGE_PREFIX
export TRAEFIK_NETWORK TRAEFIK_ENTRYPOINT
export GLOSA_TIMEOUT GLOSA_JOB_RETENTION_HOURS GLOSA_MAX_UPLOAD_MB

COMPOSE = $(COMPOSE_CLI) -p $(SERVICE_NAME) -f docker-compose.yml

.PHONY: all nodename url config build up down redeploy restart logs shell verify smoke \
        _require_nodename _require_env _traefik_net

all: build up

nodename:
	@echo "SERVICE_NAME=$(SERVICE_NAME)  NODENAME=$(NODENAME)  GLOSA_PORT=$(GLOSA_PORT)"

config: nodename
	@echo "PUBLIC_HOSTNAME=$(PUBLIC_HOSTNAME)"
	@echo "IMAGE_PREFIX=$(IMAGE_PREFIX)"
	@echo "COMPOSE=$(COMPOSE)"

url:
	@echo ""
	@echo "G-LoSA: $(URL)  (also on host port $(GLOSA_PORT))"

# Refuse to build or deploy on a machine with no NODE_MAP entry, rather than publishing under
# some other node's hostnames.
_require_nodename:
	@if [ "$(NODENAME)" = "UNMAPPED" ]; then \
		echo "ERROR: hostname '$(HOST_SHORT)' has no entry in NODE_MAP in $(CONFIG)."; \
		echo "       Add '$(HOST_SHORT)=<label>' to it; the published hostname is derived"; \
		echo "       from that label."; \
		exit 1; \
	fi

# .env carries the package-index credentials the image build needs, and compose treats a missing
# env_file as a hard error, so check for it here and say what to run.
_require_env:
	@if [ ! -f .env ]; then \
		echo "ERROR: .env is missing. Run ./build-local.sh to generate it."; \
		exit 1; \
	fi

# docker-compose.yml joins the proxy network as external, so it has to exist beforehand.
# Creating it here is idempotent, and keeping it external means `make down` leaves it -- and the
# Traefik listener sharing it -- alone.
_traefik_net:
	@$(CONTAINER_CLI) network inspect $(TRAEFIK_NETWORK) >/dev/null 2>&1 || \
		$(CONTAINER_CLI) network create $(TRAEFIK_NETWORK) >/dev/null
	@echo "shared '$(TRAEFIK_NETWORK)' network present"

# Full from-scratch build.  Slow: it recompiles the scorer, both Java programs and the jlink
# runtime.  For a UI change use `make redeploy`, which reuses those layers.
build: _require_nodename _require_env
	@echo "Building $(SERVICE_NAME) (no-cache)..."
	@$(COMPOSE) build --no-cache

up: _require_nodename _require_env _traefik_net
	@$(COMPOSE) up --detach
	@$(MAKE) --no-print-directory url

down:
	@echo "Stopping and removing $(SERVICE_NAME)..."
	@$(COMPOSE) down

# Fast path after editing app.py, glosa_runner.py or a template: a cached rebuild picks up the
# changed COPY layer only, leaving the compiler and pip layers alone.
#
# --force-recreate is required, not a precaution: podman-compose's `up` leaves an already-running
# container in place even when --build has just produced a new image, so a plain `up --build`
# rebuilds and then silently keeps serving the old code.
redeploy: _require_nodename _require_env _traefik_net
	@$(COMPOSE) up --build --detach --force-recreate
	@$(MAKE) --no-print-directory url

restart:
	@$(COMPOSE) restart

logs:
	@$(COMPOSE) logs -f

shell:
	@$(CONTAINER_CLI) exec -it $(SERVICE_NAME) bash

# Check that the local Traefik listener routes the published hostname to this container.  The
# request goes to 127.0.0.1 with an explicit Host header, which exercises that listener without
# needing to reach the upstream proxy -- some hosts are blocked from reaching it, so a direct
# request to the public URL fails there even when routing is healthy.
verify:
	@if ! $(CONTAINER_CLI) ps --format '{{.Names}}' 2>/dev/null | grep -qx $(TRAEFIK_RELAY_CONTAINER); then \
		echo "  WARNING: no '$(TRAEFIK_RELAY_CONTAINER)' container on this host -- is the listener up?"; \
	fi
	@direct=`curl -s -o /dev/null -w '%{http_code}' --max-time 10 http://127.0.0.1:$(GLOSA_PORT)/ || echo 000`; \
	routed=`curl -s -o /dev/null -w '%{http_code}' --max-time 10 -H "Host: $(PUBLIC_HOSTNAME)" http://127.0.0.1/ || echo 000`; \
	printf '  %-52s HTTP %s\n' "host port $(GLOSA_PORT) (container direct)" "$$direct"; \
	printf '  %-52s HTTP %s\n' "$(PUBLIC_HOSTNAME) (via traefik)" "$$routed"
	@echo "  (expect 200 for both; 404 means no traefik router matched the hostname,"
	@echo "   502 means it matched but the container is not serving, and 000 means"
	@echo "   nothing is listening on that port at all.)"

# Run the bundled binding-site example inside the container.  This is the check that the three
# programs -- javac output, the jlink runtime and the compiled scorer -- actually talk to each
# other; `verify` only proves the web server answers.
smoke:
	@$(CONTAINER_CLI) exec $(SERVICE_NAME) python /app/smoke_test.py
