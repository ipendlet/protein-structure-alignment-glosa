# G-LoSA alignment service — deployment

Wraps the [G-LoSA v2.2](..) local-structure aligner in a container and publishes it behind a
Traefik listener. Upload two PDB files, get back the GA-score, the transformation matrix and a
single overlay file with both structures in it.

Building and running from a source checkout is in [the top-level README](../README.md); this
file is the deployment, and the reasoning behind how the service behaves.

---

## Up and down

All commands run from this directory.

```bash
cp deploy.conf.example deploy.conf   # once: site settings, not in the repository
bash build-local.sh                  # once per credential window: writes .env, registry login
make                                 # build + up  (first time; several minutes)
make smoke                           # run the bundled example inside the container
```

After that, the day-to-day commands are:

| Command | What it does |
|---------|--------------|
| `make up` | Start the container from the existing image. |
| `make down` | Stop and remove it. The shared proxy network is left alone. |
| `make redeploy` | Cached rebuild + replace the container. **This is the one to use after editing the app or a template.** |
| `make restart` | Restart the running container without touching the image. |
| `make build` | Full `--no-cache` rebuild. Only needed when `requirements.txt` or something in `../src/` changes. |
| `make logs` | Follow the gunicorn log. |
| `make shell` | A shell inside the container. |
| `make verify` | Check that both the host port and the published hostname answer 200. |
| `make smoke` | Run a real alignment inside the container and check the score. |
| `make config` | Print every resolved setting: hostname, registry prefix, compose command. |

### Settings

All of them come from `deploy.conf`, which is gitignored;
[`deploy.conf.example`](deploy.conf.example) documents each key and the default it falls back
to. The ones that come up most:

`SERVICE_NAME` (default `glosa`) is the container name, the Traefik router name and the first
label of the hostname. `GLOSA_PORT` (default `8657`) is the host port, published alongside the
Traefik route so `make verify` can tell "the container is broken" apart from "the routing is
broken". Either can be overridden for one invocation, which is how a second copy runs beside the
first: `make up SERVICE_NAME=glosa-test GLOSA_PORT=8658`.

`make` exports the lot into the environment that compose reads, so neither
`docker-compose.yml` nor the `Dockerfile` holds a site value — which is also why
`docker compose` run by hand without those variables refuses rather than guesses.

---

## Using it

The form has three slots. Only the first two are required.

| Slot | G-LoSA flag | Meaning |
|------|-------------|---------|
| Reference | `-s1` | Stays put. Everything else is moved onto it. |
| Mobile | `-s2` | Scored against the reference and superposed onto it. |
| Carried along | `-s2w` | Not scored, but moved by the same matrix. This is where a ligand goes. |

**A PDB is the only thing you need to supply.** The upstream instructions have you run
`AssignChemicalFeatures` by hand and pass the result to `-s1cf`/`-s2cf`, but nothing in that step
needs a human decision — the features follow from the residue and atom names — so the service
derives them from each structure. The advanced section still accepts feature files if you want to
override that with a hand-edited set; supplying them for the bundled 1IA1/3SRQ pair gives the
same 0.804175 as letting the service derive them.

Five files come back:

| File | Contents |
|------|----------|
| `combined.pdb` | **The overlay.** Reference and aligned structure in one file, on different chains. |
| `ali_struct.pdb` | The mobile structure alone, after the transform. |
| `ali_struct_with.pdb` | The carried-along structure after the same transform. |
| `matrix.txt` | The rotation and translation. |
| `convergence.tsv` | The score at each search iteration, behind the [convergence plot](#the-convergence-plot). |

`combined.pdb` is a single model with the reference on the chains it already had, the aligned
copy moved onto the next free letters and anything carried along after that — typically `A`, `B`
and `C`. A multi-model file was the other option, but most viewers show only the first model
unless told otherwise, which would defeat the point of an overlay. As written, the file opens in
PyMOL or ChimeraX with both bodies visible and colourable by chain, and the results page renders
the same file in the browser.

---

## Watching a job

`/jobs` is a board of every submission, unfinished ones first, updating itself every four
seconds. A job's own page does the same. Both keep running whether or not anyone is watching,
and `/jobs.json` returns the same data for polling from a script.

The useful part is the phase, because **G-LoSA prints nothing at all until it is completely
finished** — the log of the three-hour run below was empty for three hours and then produced its
score. What it does do is leave a predictable trail in its working directory, and the phase can
be read back from which files exist and whether they are still growing:

| Phase | What it means |
|-------|---------------|
| assigning chemical features | Reading the two structures. |
| enumerating residue pairs | `pairs.rst` is being written. |
| building the product graph | `product_graph.rst` is growing. |
| **searching for the maximum clique** | The scratch has gone quiet and the CPU counter is climbing. |
| writing results | `pairs.rst` touched again; the search is over. |

That fourth row is where essentially all the time goes. On the 1NEZ/2ATP pair the product graph
stopped changing about two minutes in and the search then ran for nearly three hours against a
file that never moved again. **A job whose scratch size has stopped changing is not stuck**, and
saying so is the main thing the board is for. The advancing CPU seconds are the proof.

There is deliberately no percentage and no ETA. The search is branch-and-bound over the product
graph, so the work remaining is not knowable in advance and a progress bar would be inventing a
number. The bar on the job page is an indeterminate sweep for that reason. What is shown instead
is elapsed wall time, CPU seconds, resident memory and scratch size — enough to tell a job that
is working from one that is wedged.

### The alignment as it stands

Alongside the phase, a running job reports **the best GA-score found so far** and draws **the
current alignment** in the viewer, redrawn whenever the score improves. Both come from a
checkpoint that [a local patch to `glosa.cpp`](#the-patch-to-glosacpp) writes once per search
iteration.

This turns out to matter more than the phase does. The search runs up to three iterations
(`-iter`), and on the whole-chain 1NEZ/2ATP pair:

| | |
|---|---|
| Best score after ~20 s, during iteration 2 of 3 | **0.8978710747** |
| Final score after the full 3 h 06 min | **0.897871** |

The answer is reached in the first twenty seconds and the remaining three hours confirm it.
Being able to see that is the difference between waiting out a run and stopping it.

### The convergence plot

The job page plots those scores: the running best as a line, every individual attempt as a dot
beneath it, one point per search iteration. It is the picture of the table above — the line
reaching its final height almost immediately and then staying flat while the dots keep
scattering — and it is what tells you a long run has nothing left to find.

The plot is free. The patch was already writing a checkpoint once per iteration, so logging a
row at the same moment costs one `fprintf` on a path that runs **a dozen or so times in a whole
job**, not thousands. The history lives in `convergence.tsv` next to the results and is offered
as a download.

Nothing is plotted in the browser either. The geometry is laid out in Python and emitted as
inline SVG, so the page loads no charting library — one less thing to fail on a network that
blocks CDNs. A running job re-fetches the rendered plot from `/job/<id>/convergence.svg`, and
only when the iteration count has actually moved, so hours of polling cost a handful of fetches.

The x axis is the iteration index rather than elapsed time. The rows are wildly uneven in time —
the clique search can spend an hour between two of them — and against a time axis the entire run
collapses onto the left edge. Each point carries its own timestamp in a tooltip instead.

### Stopping a job

Running jobs have a **Stop** button, on the board and on their own page. Stopping kills the
scorer and keeps the best alignment it had reached: the last checkpoint is converted into a real
`combined.pdb`, `ali_struct.pdb` and `matrix.txt`, the job is recorded as `cancelled` rather than
`failed`, and the score is shown with the iteration it got to. Stopping the pair above after 47
seconds yields 0.8978710747 and a full set of files — the same answer the three-hour run gives.

The matrix written this way uses one more column than upstream's `%15.10f`, which is exactly
wide enough for a value like `-338.4458011053` to fill the field and run into its neighbour. The
extra column guarantees the whitespace the parser splits on.

A stop is a `SIGTERM` to the scorer, which reaches the worker thread as an ordinary failure, so a
marker file is written first to tell "stopped" apart from "crashed". It has to be written before
the kill, not after, or the two race.

A restart is handled rather than ignored. The executor lives in the process, so a redeploy or an
OOM kill abandons whatever was in flight; on startup the service marks those jobs failed, deletes
their scratch and **kills the orphaned scorer**. That last step matters: `glosa` does not die
with the process that spawned it, and being CPU-bound it will otherwise sit on a core for hours
against a job nobody is waiting for. The PID is checked against `/proc` before signalling, since
a recorded PID can be stale by minutes and the kernel recycles numbers.

---

## Runtime, and why there is a timeout

G-LoSA scores by finding a maximum clique in the product graph of the two structures, so both the
runtime and the scratch space grow with the square of the residue count. The difference between
input sizes is not a detail:

| Pair | Size | Wall time | Scratch |
|------|------|-----------|---------|
| Two binding sites (1IA1 / 2BL9) | ~100 atoms each | 0.3 s | negligible |
| Binding surface vs. whole chain (1NEZ / 2ATP, `-surf`) | ~1000 atoms | ~1 s | small |
| Two whole chains (1NEZ_H / 2ATP_C) | ~5000 atoms | **3 h 06 min** | 55 MB on disk, ~460 MB resident |

That last row is measured, not extrapolated, and it scored 0.897871 — barely different from the
0.896194 the binding sites alone give. Aligning whole chains mostly buys runtime.

So the container caps a job at `GLOSA_TIMEOUT` (default 1800 s). On a shared endpoint an
unbounded whole-chain pair is a denial of service rather than a feature. Run those offline
instead, where the limit defaults to off:

```bash
python infrastructure/align_cli.py -s1 a.pdb -s2 b.pdb -o runs/a_vs_b --timeout 0
```

Two consequences of the same cost model are worth knowing. Alignments run on a background thread
and the results page polls, because a long job would otherwise outlive any HTTP timeout; and the
thread pool is one wide, because the memory is what bites first and a second concurrent
whole-protein job is how the container gets OOM-killed rather than how it goes faster.

---

## What is in the image

| File | Role |
|------|------|
| [`Dockerfile`](Dockerfile) | Two stages: compile the Java tools + jlink a runtime, then compile `../src/glosa.cpp` and install the app. |
| [`glosa_runner.py`](glosa_runner.py) | The pipeline: features → alignment → overlay. Shared by the web app and the CLI. |
| [`app.py`](app.py) | Flask app: upload form, job queue, results page. |
| [`wsgi.py`](wsgi.py) | Hands gunicorn the Flask object. |
| [`align_cli.py`](align_cli.py) | The same pipeline from the command line, for runs too long for the web. |
| [`templates/jobs.html`](templates/jobs.html) | The status board. |
| [`smoke_test.py`](smoke_test.py) | Runs the bundled example and checks the score. Driven by `make smoke`. |
| [`dev_server.py`](dev_server.py) | Flask's development server, for work outside the container. |
| [`requirements.txt`](requirements.txt) | Exact pins, taken from the project's own `.venv`. |
| [`docker-compose.yml`](docker-compose.yml) | Container, host port, job volume, Traefik labels. Every site value comes from the environment. |
| [`Makefile`](Makefile) | The targets above. Resolves `deploy.conf` into that environment. |
| [`deploy.conf.example`](deploy.conf.example) | Template for the gitignored `deploy.conf`, with every key documented. |
| [`build-local.sh`](build-local.sh) | Generates `.env`; refreshes the image registry login when one is configured. |
| [`../.dockerignore`](../.dockerignore) | Keeps `runs/`, the locally built binary and `deploy.conf` out of the build context. |

Four decisions in there that are not obvious:

**The image compiles G-LoSA from source rather than copying a binary in.** `src/glosa.cpp` is built
with `-O2`, which matters more here than it does for glue code: the clique search *is* the
runtime, and an unoptimised binary is several times slower on anything larger than a binding
site.

**A jlink runtime, not a JDK.** The chemical-feature assignment is a Java program that has to run
on every upload, so a JVM is a runtime dependency and not just a build one. A full JDK is ~320 MB
for two class files, and both only import `java.io`, `java.text` and `java.util` — all in
`java.base` — so `jlink --add-modules java.base` cuts a ~50 MB runtime and the final image
carries that instead.

**The two Java files are compiled in separate `javac` invocations.** Both declare their own
top-level `Atom` class, so compiling them together fails with `duplicate class: Atom`. Two output
directories, `classes/acf` and `classes/ass`, and `-cp` picks one.

**The full `python:3.11` image, not `-slim`.** A build network restricted to an internal mirror
usually cannot reach `deb.debian.org`, and such mirrors rarely carry a Debian repo, so `apt`
cannot install anything. The full image is built on `buildpack-deps`, so the `g++` that
`glosa.cpp` needs is already there; on `-slim` there is no compiler and no way to add one. The
registry those base images are pulled from is `DOCKER_REGISTRY_HOST` in `deploy.conf`, empty
meaning Docker Hub.

### The patch to `glosa.cpp`

`glosa.cpp` is upstream code with one local change, marked `ADDED (not upstream)` at every site.
A `writeCheckpoint()` function writes `checkpoint.rst` — stage, iteration, current and best
score, and `u_best`/`t_best` — and is called once per iteration of each of the three search
loops. Without it a run in progress is completely opaque, since the program prints nothing and
touches no file until it is finished. The same function appends the row to `convergence.tsv`
that the plot is drawn from.

It only ever writes. No state is read back and none is modified, so the GA-score is bit-for-bit
what the unpatched program produces; `make smoke` checks that against a known value on every
build, which is what makes the patch safe to carry.

Two details worth knowing if you touch it. The scores tracked inside the loops are unnormalised
sums and the printed GA-score is `best_score / N_target`, so the checkpoint applies the same
division — otherwise it would report 37.64 where the run ends up reporting 0.896. And the file
is written to `checkpoint.tmp` and renamed into place, because the reader polls on a timer and
would otherwise sometimes catch a half-written file.

The service's own transform was checked against this: applying the checkpoint matrix to the
mobile structure in Python reproduces `glosa`'s own `ali_struct.pdb` for the same pair to
0.0000 Å.

Every upload is also passed through a normalising step before G-LoSA sees it: only the first
model's `ATOM`/`HETATM` records survive, and a trailing `TER` is appended. G-LoSA's parsers stop
at the first `TER` and require one to be present, so a file without it is silently read as empty
rather than rejected — and uploaded files routinely carry CRLF endings or NMR ensembles that
break the fixed-column parsing.

The container holds no credentials and reaches no external service. Job output lives in a named
volume so a redeploy does not throw away results, and is pruned after
`GLOSA_JOB_RETENTION_HOURS` (default 72).

---

## Prerequisite: the Traefik listener

**None of this is required.** Without a proxy the labels are inert and `http://<machine>:8657/`
serves on its own. What follows is the arrangement the defaults are shaped for.

Routing has two hops. An upstream Traefik forwards a wildcard record —
`*.<NODENAME>.<ROUTING_DOMAIN>` — to the right machine, and a Traefik listener **on that
machine** does the last hop, matching the full hostname against the `traefik.*` labels of local
containers. That listener is a separate, long-lived container, named by
`TRAEFIK_RELAY_CONTAINER`, which this Makefile does not manage; `make up` only creates the shared
network (`TRAEFIK_NETWORK`) if it is missing. Under rootless podman the listener needs
`net.ipv4.ip_unprivileged_port_start` ≤ 80 to bind port 80 at all.

`NODENAME` is the label the wildcard fans out on, and the Makefile resolves it by looking the
local short hostname up in `NODE_MAP`. A host that is absent from a non-empty map resolves to
`UNMAPPED` and `make build`/`up`/`redeploy` refuse to run, because the wildcard forwards to
whichever machine claims a given name — a wrong value would publish this deployment under a
different machine's hostnames. Add an entry to `NODE_MAP` when deploying on a new host. Where the
domain has no such level, leave `NODE_MAP` empty and the hostname loses that component.

One consequence worth knowing: where a machine is firewalled off from the upstream Traefik,
requesting the public URL *from that machine* fails even when routing is perfectly healthy. Use
`make verify`, which goes through the local listener with an explicit `Host` header, or open the
URL from somewhere else.

---

## Troubleshooting

| Symptom | Cause |
|---------|-------|
| `make verify` shows `000` on the host port | The container is not running. `make up`, then `make logs`. |
| `make verify` shows `502` via Traefik | The route matched but gunicorn is not answering — check `make logs`. |
| `make verify` shows `404` via Traefik | No router matched the hostname. Usually `NODENAME` is wrong, or the container is not on the proxy network. Check `make config`. |
| `make verify` shows `000` via Traefik | Nothing is listening on port 80: the local Traefik listener is down. |
| `ERROR: deploy.conf is missing` | `cp deploy.conf.example deploy.conf`, then edit it. |
| `ERROR: .env is missing` | Run `bash build-local.sh`. |
| `make` fails with *"has no entry in NODE_MAP"* | A new machine. Add `<short hostname>=<label>` to `NODE_MAP` in `deploy.conf`. |
| `docker compose` by hand fails on `PUBLIC_HOSTNAME` | It reads the settings from the environment, which only `make` fills in. Use the make targets. |
| A job fails with *"contains no ATOM or HETATM records"* | The upload was not a PDB, or was an `._`-prefixed AppleDouble stub rather than the structure itself. |
| A job fails with *"exceeded the … time limit"* | The pair is too large for the shared endpoint. Run it with `align_cli.py --timeout 0`. |
| The overlay panel is blank but the downloads work | The browser could not fetch 3Dmol.js from either the image or the CDN. The page says so; open `combined.pdb` in PyMOL instead. |
| `make smoke` reports a GA-score mismatch | The scorer was built differently — check that the `g++` line in the Dockerfile still passes `-O2` and that `glosa.cpp` is unmodified. |
