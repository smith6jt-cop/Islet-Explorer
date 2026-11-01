#!/usr/bin/env bash
# smoke.sh -- run inside the site container (run_in_site.sh --no-seal bash scripts/hpg/smoke.sh): the binds, the
# Python environment, git, QuPath with its OpenCV natives, flock contention and a rename across data/ and output/.
# Exits nonzero at the first failure. Leaves nothing behind.
set -euo pipefail
WS_ROOT=/data/panc_CODEX/Islet-Explorer
WS_QPI=/data/QuPath_Projects/QuPath_project_IsletExplorer
tag=${SLURM_JOB_ID:-$$}
work=$WS_ROOT/output/.hpg_smoke_$tag
mkdir "$work"
trap 'rm -rf "$work" "$WS_ROOT/data/.hpg_smoke_$tag"' EXIT

echo "== binds"
grep -E " ($WS_ROOT|$WS_QPI|/usr) " /proc/self/mountinfo | cut -d" " -f4-6
test -d "$WS_ROOT/scripts/uniform" && test -d "$WS_QPI"
touch "$WS_QPI/.write_test" 2>/dev/null && { rm -f "$WS_QPI/.write_test"; echo "QuPath folder is writable" >&2; exit 1; } || true

echo "== python"
python - <<'EOF'
import sys, os, numpy, pandas, pyarrow, tifffile, cv2, shapely, sklearn, scipy, skimage, anndata, scanpy
print(sys.executable, sys.version.split()[0])
print({m.__name__: m.__version__ for m in (numpy, pandas, pyarrow, tifffile, cv2, shapely, sklearn, scipy, skimage)})
import posix
assert os.stat is posix.stat, 'hook active in --no-seal smoke test'
EOF

echo "== git"
git -C "$WS_ROOT" rev-parse HEAD

echo "== QuPath and OpenCV"
cat > "$work/smoke.groovy" <<'EOF'
import org.bytedeco.opencv.global.opencv_core
println("QuPath " + qupath.lib.common.GeneralTools.getVersion() + ", OpenCV " + opencv_core.getVersionString().getString())
EOF
/home/smith6jt/QuPath/bin/QuPath script "$work/smoke.groovy"

echo "== flock contention"
python - "$work/lock" <<'EOF'
import fcntl, subprocess, sys
path = sys.argv[1]
with open(path, 'w') as held:
    fcntl.flock(held, fcntl.LOCK_EX)
    probe = 'import fcntl,sys\nf=open(sys.argv[1],"w")\ntry:\n fcntl.flock(f,fcntl.LOCK_EX|fcntl.LOCK_NB)\n print("acquired")\nexcept BlockingIOError:\n print("blocked")'
    out = subprocess.run([sys.executable, '-c', probe, path], capture_output=True, text=True).stdout.strip()
    assert out == 'blocked', out
print('flock blocks a second process')
EOF

echo "== rename across data/ and output/"
python - "$WS_ROOT/data/.hpg_smoke_$tag" "$work/renamed" <<'EOF'
import os, sys
src, dst = sys.argv[1:]
open(src, 'w').write('x')
os.replace(src, dst)
assert os.path.isfile(dst) and not os.path.exists(src)
print('os.replace works across data/ and output/')
EOF
echo "smoke OK"
