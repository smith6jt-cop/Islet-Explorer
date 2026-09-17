"""Prevent concurrent QuPath writers from losing shared project metadata."""
import importlib.util
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


MODULE = Path(__file__).resolve().parents[1] / "cohort_expansion/run_qupath.py"
spec = importlib.util.spec_from_file_location("cohort_runner", MODULE)
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)

CHILD = """
import importlib.util
from pathlib import Path
import sys
spec = importlib.util.spec_from_file_location('runner', sys.argv[1])
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)
try:
    with runner.project_lock(Path(sys.argv[2])):
        print('acquired')
except RuntimeError:
    raise SystemExit(78)
"""


class ProjectLockTests(unittest.TestCase):
    def attempt(self, root):
        return subprocess.run([sys.executable, "-c", CHILD, str(MODULE), str(root)],
                              capture_output=True, text=True, timeout=10)

    def test_competing_process_is_blocked_until_release(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with runner.project_lock(root):
                self.assertEqual(self.attempt(root).returncode, 78)
            result = self.attempt(root)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(result.stdout.strip(), "acquired")

    def test_different_projects_can_run_independently(self):
        with tempfile.TemporaryDirectory() as directory:
            first, second = Path(directory) / "first", Path(directory) / "second"
            first.mkdir()
            second.mkdir()
            with runner.project_lock(first):
                self.assertEqual(self.attempt(second).returncode, 0)

    def test_exception_releases_lock(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with self.assertRaisesRegex(ValueError, "interrupted"):
                with runner.project_lock(root):
                    raise ValueError("interrupted")
            self.assertEqual(self.attempt(root).returncode, 0)


if __name__ == "__main__":
    unittest.main()
