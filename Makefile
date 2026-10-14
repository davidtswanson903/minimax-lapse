# Every target runs from this directory with nothing outside it.
#   make check      the checks: generated files current, golden outputs held, tests
#   make reproduce  rerun the quick grade (or GRADE=...) and compare with stored results; writes nothing
#   make all        rerun every experiment, regenerate everything, run the tests
#   make run GRADE=g [EXP=x]  run a grade: quick, paper, survey or certification
#   make lean       build the Lean proofs and rewrite the stored report the checks read
# The interpreter: PYTHON if given, else this directory's .venv, else python on the path.
VENV_PYTHON := $(firstword $(wildcard $(CURDIR)/.venv/bin/python $(CURDIR)/.venv/Scripts/python.exe))
PYTHON ?= $(or $(VENV_PYTHON),python)
export PYTHONPATH := src

.PHONY: figures all install results run golden generate check test reproduce paper lean

all: results generate figures test

install:
	$(PYTHON) -m pip install -e ".[dev]"

results:
	$(PYTHON) -m repokit.golden run paper

run:
	$(PYTHON) -m repokit.golden run $(or $(GRADE),paper) $(EXP)

golden:
	$(PYTHON) -m repokit.golden pin $(EXP)

generate:
	$(PYTHON) -m repokit.generate

check:
	$(PYTHON) -m repokit.generate --check
	$(PYTHON) -m pytest -q

test:
	$(PYTHON) -m pytest -q

reproduce:
	$(PYTHON) -m repokit.golden reproduce $(or $(GRADE),quick) $(EXP)

paper: generate figures
	cd paper/draft && pdflatex -interaction=nonstopmode paper.tex && bibtex paper && pdflatex -interaction=nonstopmode paper.tex && pdflatex -interaction=nonstopmode paper.tex

figures:
	$(PYTHON) -m repokit.figures

lean:
	cd lean && lake exe cache get && lake build
	$(PYTHON) -m repokit.lean report
	$(PYTHON) -m repokit.generate
