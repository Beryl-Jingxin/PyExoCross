"""Run HITRAN API tests from repository root."""

import runpy


if __name__ == "__main__":
    runpy.run_path("test/test_api_hitran.py", run_name="__main__")
