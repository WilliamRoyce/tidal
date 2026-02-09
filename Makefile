.PHONY: test lint format typecheck docs clean install all help

help:  ## Show this help message
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'

test:  ## Run test suite
	uv run pytest tests/ -x -q

test-verbose:  ## Run test suite with verbose output
	uv run pytest tests/ -v

test-coverage:  ## Run tests with coverage report
	uv run pytest tests/ --cov=tidal --cov-report=term-missing --cov-report=xml

lint:  ## Check code style with ruff
	uv run ruff check

format:  ## Format code with ruff
	uv run ruff format

format-check:  ## Check if code is formatted
	uv run ruff format --check

typecheck:  ## Run type checking with pyright
	uv run pyright

docs:  ## Build Sphinx documentation
	cd docs && make html

docs-clean:  ## Clean documentation build
	cd docs && make clean

clean:  ## Remove build artifacts and caches
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete
	find . -type d -name "*.egg-info" -exec rm -rf {} +
	rm -rf .pytest_cache coverage.xml .ruff_cache

install:  ## Install dependencies
	uv sync --all-extras

all: lint typecheck test  ## Run all checks (lint, typecheck, test)
