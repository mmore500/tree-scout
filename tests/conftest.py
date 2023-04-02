import random

import numpy as np
import pytest


@pytest.fixture(scope="session", autouse=True)
def _session() -> None:
    pass


@pytest.fixture(scope="module", autouse=True)
def _module() -> None:
    pass


@pytest.fixture(autouse=True)
def reseed_random(request: pytest.FixtureRequest) -> None:
    """Reseed random number generator to ensure determinstic test."""
    random.seed(request.node.name)
    np.random.seed(random.randint(1, 2**32))
