"""Shared helpers for reading the golden fixtures.

Every fixture records results as ``{"value": repr(...)}`` or
``{"error": "Type: message"}``, so a comparison is bit-for-bit rather than
"close enough" — silent numerical drift is the main risk of this port.
"""

import json
from pathlib import Path

import pytest

GOLDEN = Path(__file__).parent / "golden"
DATA = Path(__file__).parent / "data"


def load_golden(name):
    with open(GOLDEN / name) as fh:
        return json.load(fh)


def check_result(expected, fn, *args, **kwargs):
    """Assert ``fn(*args)`` reproduces a captured value, or fails the same way."""
    if "error" in expected:
        with pytest.raises(Exception) as excinfo:  # noqa: B017 - the type is in the fixture
            fn(*args, **kwargs)
        assert f"{type(excinfo.value).__name__}: {excinfo.value}" == expected["error"]
    else:
        assert repr(fn(*args, **kwargs)) == expected["value"]
