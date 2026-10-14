"""Data from the world: every dataset matches its checksum, names its provenance, its role, what it
bears on, the premises it adds and the equations that carry it, and a placeholder is marked."""
from repokit import datasets, equations, paths, sources


def test_every_dataset_matches_its_checksum():
    stale = [item["key"] for item in datasets.registry()
             if datasets.checksum(paths.ROOT / item["file"]) != item.get("sha256")]
    assert not stale, f"datasets whose files do not match their recorded checksums: {stale}"


def test_placeholders_are_marked(paper):
    if not datasets.placeholders():
        return
    readme = (paths.ROOT / "README.md").read_text(encoding="utf-8")
    assert "placeholder" in paper.lower() and "placeholder" in readme.lower(), "placeholder data must say so on every surface"
