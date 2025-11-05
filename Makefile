.PHONY: help test lint format security clean install dev

help:  ## Show this help message
	@echo 'Usage: make [target]'
	@echo ''
	@echo 'Available targets:'
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-15s\033[0m %s\n", $$1, $$2}'

install:  ## Install package in editable mode
	pip install -e .

dev:  ## Install package with development dependencies
	pip install -e ".[dev]"

test:  ## Run tests
	pytest tests/ -v

test-cov:  ## Run tests with coverage
	pytest tests/ --cov=rue_lib --cov-report=html --cov-report=term

lint:  ## Run linter
	ruff check src/ tests/

lint-fix:  ## Run linter and fix issues
	ruff check --fix src/ tests/

format:  ## Format code
	ruff format src/ tests/

security:  ## Run security checks
	bandit -c pyproject.toml -r src/

precommit:  ## Install pre-commit hooks
	pre-commit install

precommit-run:  ## Run all pre-commit hooks
	pre-commit run --all-files

clean:  ## Clean build artifacts
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info
	rm -rf .pytest_cache/
	rm -rf .ruff_cache/
	rm -rf htmlcov/
	rm -rf .coverage
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete

build:  ## Build package
	python -m build

publish-test:  ## Publish to TestPyPI
	python -m twine upload --repository testpypi dist/*

publish:  ## Publish to PyPI
	python -m twine upload dist/*

# Nix-specific targets (require Nix to be installed)
nix-develop:  ## Enter Nix development environment
	nix develop

nix-test:  ## Run tests using Nix
	nix run .#test

nix-test-cov:  ## Run tests with coverage using Nix
	nix run .#test-cov

nix-build:  ## Build package using Nix
	nix build
