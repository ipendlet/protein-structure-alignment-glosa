# G-LoSA alignment service

Local structure alignment with [G-LoSA v2.2](https://compbio.lehigh.edu/GLoSA/), wrapped so that
two PDB files are the only thing anyone has to supply. Give it a reference and a mobile
structure and it returns the GA-score, the transformation matrix and a single overlay file with
both bodies in it — as a web service, as a command-line tool, or as a Python call.

```
src/                 G-LoSA sources: the C++ scorer and the two Java feature tools
infrastructure/      the service: Flask app, pipeline, CLI, container build, deployment
example/             the structures shipped with upstream G-LoSA, used by the smoke tests
runs/                output (gitignored)
```

Upstream G-LoSA is three programs that have to be driven in sequence, with intermediate files
passed between them by hand. [`infrastructure/glosa_runner.py`](infrastructure/glosa_runner.py)
hides that behind one `align()` call, and everything else here — the web app, the CLI, the smoke
tests — is a front end to it, so all three drive an identical pipeline.

---

## Build

Needs a C++ compiler, a JDK and Python 3.9+. Nothing else; the service has no database and no
external service to reach.

```bash
make                 # compiles src/glosa.cpp with -O2 and both .java files into classes/
```

`-O2` is not optional in the way it is for glue code. G-LoSA scores by finding a maximum clique
in the product graph of the two structures, that search *is* the runtime, and an unoptimised
binary is several times slower on anything larger than a binding site.

The two Java files are compiled separately into `classes/acf` and `classes/ass` because both
declare their own top-level `Atom` class and compiling them together fails with
`duplicate class: Atom`.

Then install the Python dependencies:

```bash
python3 -m venv .venv && .venv/bin/pip install -r infrastructure/requirements.txt
```

The compiled artefacts land where `glosa_runner.py` looks for them in a source checkout, so
there is nothing to configure. In the container they live under `/opt/glosa` instead, and
`GLOSA_BIN`, `GLOSA_CLASSES_ACF` and `GLOSA_CLASSES_ASS` point at them there — set those if you
put them somewhere else again.

Check the build with the bundled binding-site pair, which takes about a third of a second:

```bash
make example         # expects GA-score 0.896194
```

G-LoSA is deterministic, so a different score means the binary was built differently, not that
the inputs changed.

---

## Run

### From the command line

```bash
python3 infrastructure/align_cli.py \
    -s1 example/1IA1-bs.pdb \
    -s2 example/2BL9-bs.pdb \
    -o runs/1IA1_vs_2BL9 \
    --transfer example/2BL9-lig.pdb
```

| Flag | Meaning |
|------|---------|
| `-s1` | Reference. Stays put; everything else is moved onto it. |
| `-s2` | Mobile. Scored against the reference and superposed onto it. |
| `--transfer` | Carried along by the same matrix but not scored. This is where a ligand goes. |
| `-o` | Output directory. |
| `--timeout` | Seconds before the run is killed. Defaults to `0`, meaning no limit. |
| `--keep-scratch` | Keep `pairs.rst` and `product_graph.rst`, which can reach tens of gigabytes. |
| `--s1cf` / `--s2cf` | Chemical feature files, if you want to override the derived ones. |

Trailing arguments are passed to `glosa` untouched, so `-iter 1 -n 3` and the rest of the
upstream flags still work.

**A PDB is the only thing you need to supply.** Upstream has you run `AssignChemicalFeatures` by
hand and pass the result to `-s1cf`/`-s2cf`, but nothing in that step needs a human decision —
the features follow from the residue and atom names — so the pipeline derives them.

Five files come back:

| File | Contents |
|------|----------|
| `combined.pdb` | **The overlay.** Reference and aligned structure in one file, on different chains. |
| `ali_struct.pdb` | The mobile structure alone, after the transform. |
| `ali_struct_with.pdb` | The carried-along structure after the same transform. |
| `matrix.txt` | The rotation and translation. |
| `convergence.tsv` | The score at each search iteration. |

### From Python

```python
from glosa_runner import align

result = align(structure1="a.pdb", structure2="b.pdb", outdir="runs/a_vs_b")
print(result.ga_score, result.outputs["combined"])
```

### As a web service, locally

```bash
make serve           # http://127.0.0.1:8057
```

Flask's development server, single-process and not meant to be exposed — use the container for
anything anyone else will touch. `make serve DEV_PORT=9000` moves it.

### Runtime, and why any of this matters

Both the runtime and the scratch space grow with the square of the residue count, and the
difference between input sizes is not a detail:

| Pair | Size | Wall time | Scratch |
|------|------|-----------|---------|
| Two binding sites (1IA1 / 2BL9) | ~100 atoms each | 0.3 s | negligible |
| Binding surface vs. whole chain (1NEZ / 2ATP, `-surf`) | ~1000 atoms | ~1 s | small |
| Two whole chains (1NEZ_H / 2ATP_C) | ~5000 atoms | **3 h 06 min** | 55 MB on disk, ~460 MB resident |

That last row is measured, not extrapolated, and it scored 0.897871 — barely different from the
0.896194 the binding sites alone give. Aligning whole chains mostly buys runtime, which is why
the deployed service caps a job and the CLI does not: long pairs belong offline.

---

## Deploy

The container build compiles the scorer and the Java tools itself, `jlink`s a ~50 MB JVM for
them and serves the app under gunicorn. All commands run from `infrastructure/`.

```bash
cd infrastructure
cp deploy.conf.example deploy.conf     # then edit it: see below
bash build-local.sh                    # writes .env; refreshes the registry login if needed
make                                   # build + start (several minutes: it compiles everything)
make smoke                             # run a real alignment inside the container
```

Day to day:

| Command | What it does |
|---------|--------------|
| `make redeploy` | Cached rebuild + replace the container. **The one to use after editing the app or a template.** |
| `make build` | Full `--no-cache` rebuild. Only needed when `requirements.txt` or something in `src/` changes. |
| `make up` / `make down` | Start from the existing image / stop and remove. |
| `make logs` / `make shell` | Follow the gunicorn log / get a shell in the container. |
| `make verify` | Check that both the host port and the published hostname answer 200. |
| `make config` | Print the resolved settings — hostname, registry prefix, compose command. |

[`infrastructure/README.md`](infrastructure/README.md) covers the rest: what is in the image and
why, the local patch to `glosa.cpp` that makes a running job observable, the job board, and
troubleshooting.

### Configuration

Everything site-specific — registry, routing domain, machine names, limits — lives in
`infrastructure/deploy.conf`, which is **not in the repository**. Copy
[`deploy.conf.example`](infrastructure/deploy.conf.example), which documents every key, and
edit. The minimum is one line, because everything else has a working default:

```ini
# infrastructure/deploy.conf
ROUTING_DOMAIN=apps.example.com
```

A fuller one, for a network behind a package and image mirror, with one wildcard record fanning
out to several machines:

```ini
SERVICE_NAME=glosa
GLOSA_PORT=8657

ROUTING_DOMAIN=apps.example.com
NODE_MAP=build-node-01=gpu1 build-node-02=gpu2 workstation-07=ws07

CONTAINER_CLI=podman
DOCKER_REGISTRY_HOST=registry.example.com/docker-virtual
REGISTRY_LOGIN_HOST=registry.example.com

TRAEFIK_NETWORK=traefik
TRAEFIK_ENTRYPOINT=web

GLOSA_TIMEOUT=1800
GLOSA_JOB_RETENTION_HOURS=72
GLOSA_MAX_UPLOAD_MB=32
```

The service is published at `<SERVICE_NAME>.<NODENAME>.<ROUTING_DOMAIN>`, where `NODENAME` is the
label `NODE_MAP` gives the local machine — or at `<SERVICE_NAME>.<ROUTING_DOMAIN>` when
`NODE_MAP` is empty. Setting `PUBLIC_HOSTNAME` directly bypasses the scheme entirely. A machine
that is absent from a non-empty `NODE_MAP` resolves to `UNMAPPED`, and `make build`/`up`/
`redeploy` refuse to run rather than publish this deployment under another machine's hostname.

`make` passes all of it to compose through the environment, so `docker-compose.yml` and the
`Dockerfile` carry no site values of their own. Any key can also be overridden for one
invocation — `make up SERVICE_NAME=glosa-test GLOSA_PORT=8658` runs a second copy alongside the
first.

### Credentials

`build-local.sh` writes `infrastructure/.env`, which is gitignored. On a network that reaches
PyPI and Docker Hub directly the file is empty and that is the end of it. Behind a mirror it
carries the package-index URL, credentials and all, read out of the local pip config; setting
`REGISTRY_LOGIN_HOST` additionally logs the container CLI into the image registry with the same
credential, which is the usual arrangement for Artifactory and Nexus. Re-run it when the token
expires.

The service itself holds no credentials and reaches no external service — it reads nothing but
the structures it is given.

### Reverse proxy

The container publishes Traefik labels and joins an existing external network, which a Traefik
instance on the same machine watches. That instance is long-lived and not managed here; `make
up` only creates the shared network if it is missing.

Deploying without a proxy costs nothing: the labels are inert and the host port still serves, so
`http://<machine>:8657/` works on its own.

---

## Licence and attribution

G-LoSA is by Hui Sun Lee and Wonpil Im (Lehigh University); `src/glosa.cpp` and the two
`src/Assign*.java` files are their code, and `example/` is their example set. `glosa.cpp` carries
one local change, marked `ADDED (not upstream)` at every site, which writes a checkpoint and a
convergence row once per search iteration so that a run in progress is not completely opaque. It
only ever writes — no state is read back and none is modified — so the GA-score is bit-for-bit
what the unpatched program produces, and `make smoke` checks that against a known value on every
build. See [infrastructure/README.md](infrastructure/README.md#the-patch-to-glosacpp).

If you use this, cite the G-LoSA papers:

> Lee HS, Im W. *G-LoSA: An efficient computational tool for local structure-centric biological
> studies and drug design.* Protein Science 25:865–876 (2016).
>
> Lee HS, Im W. *G-LoSA for prediction of protein-ligand binding sites and structures.* Methods
> in Molecular Biology 1611:97–108 (2017).
