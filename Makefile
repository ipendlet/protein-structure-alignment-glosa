# Build and run G-LoSA from a source checkout.  Deployment lives in infrastructure/Makefile.
#
#   make            → compile the scorer and the Java feature tools
#   make example    → align the bundled binding-site pair, as a check that the three agree
#   make serve      → the web UI on http://127.0.0.1:8057 (development server, not for deployment)
#   make deploy     → hand over to infrastructure/Makefile (container build + publish)

CXX      ?= g++
CXXFLAGS ?= -O2
JAVAC    ?= javac
PYTHON   ?= python3
DEV_PORT ?= 8057

# Where the compiled artefacts land.  These paths are glosa_runner.py's defaults for a source
# checkout, so nothing has to be configured to run from here; the image overrides them with
# GLOSA_BIN / GLOSA_CLASSES_* because it puts the same files under /opt/glosa.
GLOSA_BIN   := glosa
CLASSES_ACF := classes/acf/AssignChemicalFeatures.class
CLASSES_ASS := classes/ass/AssignSecondaryStructures.class

.PHONY: all classes example serve align clean deploy help

all: $(GLOSA_BIN) classes

help:
	@echo "make           compile the scorer and the Java feature tools"
	@echo "make example   align the bundled 1IA1/2BL9 binding sites into runs/"
	@echo "make serve     development web server on port $(DEV_PORT)"
	@echo "make deploy    build and publish the container (see infrastructure/README.md)"
	@echo "make clean     remove compiled artefacts"

# -O2 is not optional here in the way it is for glue code: the maximum-clique search over the
# product graph is the entire runtime, and an unoptimised binary is several times slower on
# anything larger than a binding site.
$(GLOSA_BIN): src/glosa.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

classes: $(CLASSES_ACF) $(CLASSES_ASS)

# Separate javac invocations and separate output directories: both .java files declare their own
# top-level `Atom` class, so compiling them together fails with "duplicate class: Atom".
$(CLASSES_ACF): src/AssignChemicalFeatures.java
	@mkdir -p classes/acf
	$(JAVAC) -nowarn -d classes/acf $<

$(CLASSES_ASS): src/AssignSecondaryStructures.java
	@mkdir -p classes/ass
	$(JAVAC) -nowarn -d classes/ass $<

# The cheap pair -- a few hundred milliseconds -- so this checks the wiring without waiting on a
# whole-protein clique search.
example: all
	$(PYTHON) infrastructure/run_example.py

serve: all
	GLOSA_DEV_PORT=$(DEV_PORT) $(PYTHON) infrastructure/dev_server.py

deploy:
	@$(MAKE) -C infrastructure

clean:
	rm -rf $(GLOSA_BIN) classes
