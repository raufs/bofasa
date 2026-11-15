# Contributions

Contributions via pull-requests are welcome! Please add new feature requests/ideas or notify us of any potential bugs via GitHub issues.

## Running the unit tests

From a source checkout, with the bofasa dependencies installed (e.g. the conda environment from `bofasa_env.yml`):

```bash
pip install pytest
pytest
```

The splitting tests (`tests/test_dog_splitting.py`, `tests/test_protein_splitting.py`) compare against results produced by bofasa v1.3.0, stored with their input trees in `tests/data/dog_trees.tar.gz` and `tests/data/protein_njt.tar.gz` (uncompressed automatically into `tests/data/.extracted/` at test time). If you intentionally change splitting behaviour, regenerate them with `python tests/data/make_test_fixtures.py` (needs `fastme` on PATH) and explain the change in your pull request.

`run_tests.sh` runs bofasa end to end on the four-genome test case (`test_case.tar.gz`) and needs `bofasa setup` to have been run first.
