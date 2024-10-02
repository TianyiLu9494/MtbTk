.PHONY: clean build all test install

clean:
	rm -rf dist src/MtbTk.egg-info

build:
	python -m build

test:
	python -m pytest

install:
	pip install -e .

all: clean build