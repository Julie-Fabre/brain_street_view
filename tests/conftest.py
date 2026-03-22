import sys
import os
import pytest

# Add parent dir to path so we can import bsv
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

ALLEN_ATLAS_PATH = os.environ.get('ALLEN_ATLAS_PATH', '')
SAVE_LOCATION = os.environ.get('BSV_SAVE_LOCATION', '')

# Known VISp experiment IDs for testing (small set)
TEST_EXPERIMENT_IDS = [100141219, 112423392]


@pytest.fixture
def atlas_path():
    if not os.path.exists(ALLEN_ATLAS_PATH):
        pytest.skip('Allen atlas not found locally')
    return ALLEN_ATLAS_PATH


@pytest.fixture
def save_location():
    return SAVE_LOCATION


@pytest.fixture
def test_experiment_ids():
    return TEST_EXPERIMENT_IDS
