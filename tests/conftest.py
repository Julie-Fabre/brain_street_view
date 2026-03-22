import sys
import os
import pytest

# Add parent dir to path so we can import bsv
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

ALLEN_ATLAS_PATH = '/home/julie/Dropbox/Atlas/allenCCF'
SAVE_LOCATION = '/home/julie/Dropbox/Data/AllenQueries'

# Known VISp experiment IDs for testing (small set)
TEST_EXPERIMENT_IDS = [100141219, 112423392]


@pytest.fixture
def atlas_path():
    return ALLEN_ATLAS_PATH


@pytest.fixture
def save_location():
    return SAVE_LOCATION


@pytest.fixture
def test_experiment_ids():
    return TEST_EXPERIMENT_IDS
