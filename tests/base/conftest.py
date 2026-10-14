"""Fixtures for the reader's checks. Everything is read once per session."""
import pytest

from repokit import sources


@pytest.fixture(scope="session")
def claims():
    return sources.claims()


@pytest.fixture(scope="session")
def values():
    return sources.values()


@pytest.fixture(scope="session")
def paper():
    return sources.paper()
