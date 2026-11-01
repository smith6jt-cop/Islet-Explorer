"""HPG-only start-up hook: on PYTHONPATH (scripts/hpg) inside the site container; inert unless CODEX_SITE_SEAL is set.

See site_identity.py. Python's site.py only prints a warning when sitecustomize raises, so any failure here ends the
process with exit code 78 instead of running without the hook.
"""
import os

if os.environ.get('CODEX_SITE_SEAL'):
    try:
        import site_identity
        site_identity.install_from_environment()
    except BaseException as e:
        import sys
        print(f'CODEX site identity hook failed: {type(e).__name__}: {e}', file=sys.stderr, flush=True)
        os._exit(78)
