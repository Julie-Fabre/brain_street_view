# Contributing

Thanks for your interest in contributing!

**Reporting bugs:** Open an issue on [GitHub](https://github.com/Julie-Fabre/brain_street_view/issues) with a description, steps to reproduce, and your Python/MATLAB version and OS.

**Requesting features:** Open a GitHub issue describing the feature and its use case.

**Contributing code:**

1. Fork the repository
2. Create a branch (`git checkout -b my-feature`)
3. Make your changes
4. Run the tests: `python -m pytest tests/`
5. Open a pull request

Please keep code concise and match the existing style.

**Running the tests:** Install dev dependencies with `pip install -e ".[dev]"`, then
run `python -m pytest tests/`. Tests are grouped by marker: `-m "not network"` runs the
fast, deterministic subset that needs no network or atlas (what CI runs on every Python
version); `-m network` runs the live [Allen Brain Atlas API](http://api.brain-map.org/)
integration tests. Tests that need the Allen CCF atlas are skipped unless `ALLEN_ATLAS_PATH`
points to the atlas files.

**Getting help:** Open a [GitHub issue](https://github.com/Julie-Fabre/brain_street_view/issues) or email juliemfabre@gmail.com.
