"""What a reader regenerates from the stored results matches the files here, and each carries its banner."""
from repokit import generate, paths


def test_generated_files_are_fresh():
    assert not generate.stale(), f"run `make generate`: {generate.stale()}"


def test_generated_files_carry_banners():
    out = generate.outputs()
    missing = [paths.rel(p) for p, text in out.items() if p != generate.VALUES_JSON and "GENERATED" not in text.splitlines()[0]]
    assert not missing, f"generated files without a banner: {missing}"
