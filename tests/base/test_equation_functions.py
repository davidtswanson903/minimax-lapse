"""Equations as code: the paper states the mathematics; each equation the code computes is a
function carrying the equation's label and stating it; each such function is checked by
another route; and the table of all three is generated and placed."""
from repokit import equations, paths


def test_every_equation_function_names_an_equation_of_the_paper(paper):
    unknown = sorted({e.label for e in equations.functions()} - set(equations.paper_labels(paper)))
    assert not unknown, f"functions carry labels the paper does not have: {unknown}"


def test_every_computed_equation_has_its_function(paper):
    computed = {e.label for e in equations.functions()}
    missing = [label for label, marked in equations.paper_labels(paper).items() if not marked and label not in computed]
    assert not missing, f"labelled equations with no function, and not marked % not computed: {missing}"


def test_equation_functions_state_their_equation():
    bare = [f"{e.function} ({e.label})" for e in equations.functions() if e.label not in e.doc]
    assert not bare, f"equation functions whose docstring does not state their equation: {bare}"


def test_every_equation_function_is_checked_by_another_route():
    checked = equations.check_functions()
    unchecked = sorted({e.label for e in equations.functions()} - set(checked))
    unstated = [test for tests in checked.values() for test, first in tests if not first.startswith("By ")]
    assert not unchecked, f"equations with no independent check: {unchecked}"
    assert not unstated, f"checks whose first line does not name their route: {unstated}"


def test_the_equation_table_is_placed(paper):
    if not equations.functions():
        return
    assert "table-equations.tex" in paper, "the paper does not place the table of equations"
    assert "Computed by" in (paths.ROOT / "README.md").read_text(encoding="utf-8"), "the README does not place it"


def test_the_paper_loads_what_its_tables_need():
    """The generated tables draw booktabs rules: a paper that places one loads the presentation layer or booktabs."""
    from repokit import latex, sources

    paper = latex.strip_comments(sources.paper())
    if "build/latex/table-" in paper:
        assert "build/latex/presentation.tex" in paper or "{booktabs}" in paper, "the tables need booktabs"
