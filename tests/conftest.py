"""Shared pytest configuration for the bofasa test suite."""
import os
import sys
import warnings

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(TESTS_DIR, "data")
SRC_DIR = os.path.join(os.path.dirname(TESTS_DIR), "src")

# Allow running the tests from a source checkout without installing bofasa.
if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)
# Fixture helpers (tree generation, split-file parsing).
if DATA_DIR not in sys.path:
    sys.path.insert(0, DATA_DIR)

# ete3 3.1.3 emits SyntaxWarnings for invalid escape sequences on Python 3.12+.
warnings.filterwarnings("ignore", category=SyntaxWarning)

# Unpack the compressed fixtures (tests/data/*.tar.gz) into tests/data/.extracted/
# before test modules are collected.
import make_test_fixtures  # noqa: E402

make_test_fixtures.extract()
