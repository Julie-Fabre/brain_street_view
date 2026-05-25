# Contributing to Brain Street View

Thanks for your interest in contributing!

## Reporting bugs

Please open an issue on [GitHub](https://github.com/Julie-Fabre/brain_street_view/issues) with:
- A description of the problem
- Steps to reproduce it
- Your Python/MATLAB version and OS

## Requesting features

Open a GitHub issue describing the feature and its use case.

## Contributing code

1. Fork the repository
2. Create a branch (`git checkout -b my-feature`)
3. Make your changes
4. Run the tests: `python -m pytest tests/`
5. Open a pull request

Please keep code concise and match the existing style.

## Running the tests

```bash
pip install -e ".[dev]"
python -m pytest tests/
```

Tests are grouped with markers so you can run only what your environment supports:

- `python -m pytest tests/ -m "not network"` — fast, deterministic tests that need
  neither network access nor the atlas. This is what CI runs on every Python version.
- `python -m pytest tests/ -m network` — integration tests that query the live
  [Allen Brain Atlas API](http://api.brain-map.org/).
- Tests using the `atlas_path` fixture are skipped automatically unless the Allen CCF
  atlas files are available. Set `ALLEN_ATLAS_PATH` to a directory containing them to
  run these locally, and optionally `BSV_SAVE_LOCATION` to reuse a cache of fetched
  experiment data.

## Getting help

Open a [GitHub issue](https://github.com/Julie-Fabre/brain_street_view/issues) or email juliemfabre@gmail.com.
