"""Make a separate native QuPath review copy from a completed candidate replay."""

import argparse
import json
import os
import shutil
import subprocess
from pathlib import Path


def build(project, dest, evidence, qupath):
    summary = evidence / "transitions.csv"
    if not summary.is_file():
        raise ValueError("Cohort replay must finish before building the review project")
    shutil.copytree(project, dest)
    (dest / "V24_REVIEW_COPY").touch(exist_ok=True)
    (dest / "REVIEW.txt").write_text(
        f"Candidate for owner review. Transition audit and evidence: {evidence.resolve()}\n"
    )
    env = dict(
        os.environ,
        JAVA_TOOL_OPTIONS="-Xmx30g -XX:ActiveProcessorCount=4 -Djava.awt.headless=true",
    )
    with (evidence / "qupath.log").open("w") as log:
        subprocess.run(
            [
                str(qupath),
                "script",
                str(
                    Path(__file__).with_name("qupath") / "apply_review_revision.groovy"
                ),
                "--args",
                str(dest.resolve()),
                "--args",
                str(evidence.resolve()),
            ],
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )
    entries = json.loads((dest / "project.qpproj").read_text())["images"]
    for entry in entries:
        tag = entry["imageName"]
        r = json.loads((evidence / "samples" / (tag + ".json")).read_text())
        if r["image"] != tag or str(r["entry"]) != str(entry["entryID"]):
            raise ValueError(f"{tag}: missing native readback")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("project", "dest", "evidence", "qupath"):
        p.add_argument("--" + name, type=Path, required=True)
    build(**vars(p.parse_args()))
