import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.dirname(__file__))
RESULTS = os.path.join(ROOT, "data", "interim", "tests", "results")
os.makedirs(RESULTS, exist_ok=True)
GLOBAL_SEED = 20261030
