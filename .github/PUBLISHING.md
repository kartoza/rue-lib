# Publishing to PyPI

This guide explains how to publish rue-lib to PyPI (Python Package Index).

## Prerequisites

1. **PyPI Account**: Create accounts on:
   - [TestPyPI](https://test.pypi.org/account/register/) (for testing)
   - [PyPI](https://pypi.org/account/register/) (for production)

2. **API Tokens**: Generate API tokens for both accounts:
   - Go to Account Settings → API Tokens
   - Create a token with "Entire account" scope
   - Save the tokens securely

3. **Configure PyPI credentials**:
   ```bash
   # Create ~/.pypirc file
   cat > ~/.pypirc << 'EOF'
   [distutils]
   index-servers =
       pypi
       testpypi

   [pypi]
   username = __token__
   password = pypi-YOUR-TOKEN-HERE

   [testpypi]
   repository = https://test.pypi.org/legacy/
   username = __token__
   password = pypi-YOUR-TESTPYPI-TOKEN-HERE
   EOF

   chmod 600 ~/.pypirc
   ```

## Build the Package

### Using Nix (Recommended)

```bash
nix develop
make build
```

### Using Python directly

```bash
# Activate virtual environment
source .venv/bin/activate

# Install build tools
pip install build twine

# Build the package
python -m build
```

This creates two files in the `dist/` directory:
- A source distribution (`.tar.gz`)
- A wheel distribution (`.whl`)

## Test on TestPyPI First

Always test your package on TestPyPI before publishing to PyPI:

```bash
# Upload to TestPyPI
python -m twine upload --repository testpypi dist/*

# Or using Make
make publish-test
```

### Install from TestPyPI to test

```bash
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ rue-lib

# Test the installation
python -c "from rue_lib.core import greet; print(greet())"
```

Note: We use `--extra-index-url` because dependencies (like `rich`) are on PyPI, not TestPyPI.

## Publish to PyPI

Once you've verified the package works correctly:

```bash
# Upload to PyPI
python -m twine upload dist/*

# Or using Make
make publish
```

## Version Management

Before publishing:

1. **Update version** in `pyproject.toml`:
   ```toml
   [project]
   version = "0.2.0"  # Update this
   ```

2. **Update version** in `src/rue_lib/__init__.py`:
   ```python
   __version__ = "0.2.0"  # Update this
   ```

3. **Tag the release**:
   ```bash
   git tag v0.2.0
   git push origin v0.2.0
   ```

## Automated Publishing with GitHub Actions

You can automate PyPI publishing using GitHub Actions. Create `.github/workflows/publish.yml`:

```yaml
name: Publish to PyPI

on:
  release:
    types: [published]

jobs:
  publish:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4

      - name: Set up Python
        uses: actions/setup-python@v5
        with:
          python-version: "3.11"

      - name: Install dependencies
        run: |
          python -m pip install --upgrade pip
          pip install build twine

      - name: Build package
        run: python -m build

      - name: Publish to PyPI
        env:
          TWINE_USERNAME: __token__
          TWINE_PASSWORD: ${{ secrets.PYPI_API_TOKEN }}
        run: python -m twine upload dist/*
```

Then add your PyPI API token as a GitHub secret named `PYPI_API_TOKEN`.

## Publishing Checklist

- [ ] Update version in `pyproject.toml` and `__init__.py`
- [ ] Update CHANGELOG (if you have one)
- [ ] Run all tests: `make test`
- [ ] Run linter: `make lint`
- [ ] Run security checks: `make security`
- [ ] Build the package: `make build`
- [ ] Test on TestPyPI: `make publish-test`
- [ ] Install and test from TestPyPI
- [ ] Publish to PyPI: `make publish`
- [ ] Create GitHub release
- [ ] Tag the release in git

## Troubleshooting

### Package already exists

If you get an error that the package already exists, you need to bump the version number.

### Import errors

Make sure all dependencies are correctly specified in `pyproject.toml` and that the package structure is correct.

### Missing files

Check `MANIFEST.in` and ensure all necessary files are included.

## Resources

- [Python Packaging User Guide](https://packaging.python.org/)
- [PyPI Help](https://pypi.org/help/)
- [TestPyPI](https://test.pypi.org/)
- [Twine Documentation](https://twine.readthedocs.io/)
