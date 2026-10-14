"""One Python version: the one pyproject.toml requires is the one continuous integration runs, and an older one is refused."""
import re

import pytest

import repokit
from repokit import paths


def test_ci_runs_the_required_python():
    least = repokit.required_python(paths.ROOT / "pyproject.toml")
    workflow = (paths.ROOT / ".github" / "workflows" / "ci.yml").read_text(encoding="utf-8")
    assert least, "pyproject.toml states no requires-python"
    versions = re.findall(r'python-version:\s*"?([\d.]+)', workflow)
    assert versions and set(versions) == {f"{least[0]}.{least[1]}"}, versions


def test_an_older_python_is_refused(tmp_path):
    pyproject = tmp_path / "pyproject.toml"
    pyproject.write_text('[project]\nrequires-python = ">=3.14"\n', encoding="utf-8")
    with pytest.raises(SystemExit, match="needs Python 3.14 or later"):
        repokit.require_python((3, 10), pyproject)
    repokit.require_python((3, 14), pyproject)
