#!/usr/bin/env python3
"""Run a Python script under an audit hook that records every path it opens (Python-level opens only; files opened
by C libraries such as pyarrow are not seen). Used to list the files the read-only gates read before the HPG copy.

    python scripts/hpg/audit_opens.py OPENED.txt scripts/hpg/gate_dump.py --out DUMP.json
"""
import os
import runpy
import sys


def main() -> None:
    log, script = sys.argv[1], sys.argv[2]
    sys.argv = sys.argv[2:]
    seen = set()

    def hook(event, args):
        if event == 'open' and isinstance(args[0], (str, bytes, os.PathLike)):
            seen.add(os.path.abspath(os.fsdecode(args[0])))
    sys.addaudithook(hook)
    try:
        runpy.run_path(script, run_name='__main__')
    finally:
        with open(log, 'x') as stream:
            stream.write(''.join(p + '\n' for p in sorted(seen)))


if __name__ == '__main__':
    main()
