.PHONY: install fetch run test coverage build

install:
	python3 -m venv .venv
	.venv/bin/python -m pip install -e '.[dev,download]'

fetch:
	.venv/bin/coregene-bayes fetch

run:
	.venv/bin/coregene-bayes run

test:
	.venv/bin/pytest

coverage:
	.venv/bin/pytest --cov=coregene_bayes --cov-report=term-missing

build:
	.venv/bin/python -m build

