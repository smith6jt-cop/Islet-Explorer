"""Render native TIFF channels with QuPath cell outlines for stratified UUID review."""

import argparse
import html
import json
import sys
from contextlib import closing
from pathlib import Path
from urllib.parse import unquote, urlparse

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont

sys.path.insert(0, str(Path(__file__).resolve().parent / "typing"))
from redsea_cells import Channels

SENTINEL = "87b35e6a-1f95-46bc-b1d8-85afe222ad09"
MARKERS = (
    "DAPI",
    "CD45",
    "CD3e",
    "CD8a",
    "CD4",
    "INS",
    "CHGA",
    "SST",
    "GCG",
    "CD20",
    "CD68",
    "CD163",
    "LGALS3",
    "SMA",
    "panCK",
    "ECAD",
    "CK19",
    "VIM",
    "CD34",
    "PDPN",
    "CD31",
    "B3TUBB",
)
SIZE = 60


def render(project, dest, typing):
    run_root = typing if (typing / "qupath").is_dir() else typing.parent
    gallery = dest / "gallery"
    gallery.mkdir(exist_ok=True)
    font = ImageFont.truetype("DejaVuSans.ttf", 13)
    pages = []
    for entry in json.loads((project / "project.qpproj").read_text())["images"]:
        name = entry["imageName"]
        geom = json.loads((dest / "samples" / (name + ".json")).read_text())["polygons"]
        a = pd.read_parquet(dest / "audit" / (name + ".parquet"))
        ident = pd.read_parquet(
            dest / "identity" / (name + ".parquet"), columns=["object_id", "label"]
        ).set_index("object_id")
        b = entry["serverBuilder"]
        while "builder" in b:
            b = b["builder"]
        f = Path(unquote(urlparse(b["uri"]).path))
        source = json.loads(
            (
                run_root / "qupath" / "data" / str(entry["entryID"]) / "server.json"
            ).read_text()
        )
        labels = [c["name"] for c in source["metadata"]["channels"]]
        if len(labels) != len(b["metadata"]["channels"]):
            raise ValueError(f"{name}: source and display channel counts differ")
        windows, images = {}, {}
        with closing(Channels(f)) as channels:
            series = channels.tf.series[0]
            if series.axes != "CYX" or series.shape[0] != len(labels):
                raise ValueError(f"{name}: channel dimensions differ")
            for uid, polygon in geom.items():
                pts = np.asarray(polygon)
                center = (pts.min(axis=0) + pts.max(axis=0)) / 2
                x, y = np.clip(
                    center.astype(int) - SIZE // 2,
                    0,
                    [series.shape[2] - SIZE, series.shape[1] - SIZE],
                )
                windows[uid] = (x, y, pts - [x, y])
            for m in MARKERS:
                indices = [
                    i
                    for i, label in enumerate(labels)
                    if label == m
                    or (m == "DAPI" and label.startswith("DAPI") and i == 0)
                ]
                if len(indices) != 1:
                    raise ValueError(
                        f"{name}: ambiguous or missing {m} channel {indices} in {labels}"
                    )
                idx = indices[0]
                plane = channels.native(idx)
                crops = {
                    uid: plane[y : y + SIZE, x : x + SIZE].copy()
                    for uid, (x, y, _) in windows.items()
                }
                del plane
                upper = max(1, np.percentile(np.stack(list(crops.values())), 99.5))
                images[m] = {
                    uid: Image.fromarray(np.uint8(np.clip(crop / upper, 0, 1) * 255))
                    for uid, crop in crops.items()
                }
        records = []
        for uid in sorted(geom, key=lambda u: (u != SENTINEL, u)):
            if uid in a.index:
                row = a.loc[uid]
                title = f"{row.old_label} → {row.destination} · {row.rule or row.old_status or row.disposition}"
            else:
                title = "Control · " + str(ident.at[uid, "label"])
            tile = Image.new("RGB", (len(MARKERS) * 82, 112), "white")
            draw = ImageDraw.Draw(tile)
            for i, m in enumerate(MARKERS):
                patch = (
                    images[m][uid]
                    .resize((82, 82), Image.Resampling.NEAREST)
                    .convert("RGB")
                )
                outline = [tuple(p * 82 / SIZE) for p in windows[uid][2]]
                ImageDraw.Draw(patch).line(
                    outline + [outline[0]], fill=(0, 190, 230), width=1
                )
                tile.paste(patch, (82 * i, 22))
                draw.text((i * 82 + 3, 2), m, font=font, fill="black")
            tile.save(gallery / (uid + ".png"))
            if uid == SENTINEL:
                tile.crop((0, 0, 82 * 7, 107)).resize(
                    (82 * 14, 214), Image.Resampling.NEAREST
                ).save(dest / "sentinel.png")
            records.append(
                f'<p><b>{html.escape(title)}</b> &nbsp; {uid}<br><img loading="lazy" src="{uid}.png"></p>'
            )
        (gallery / (name + ".html")).write_text(
            '<meta charset="utf-8"><h1>'
            + name
            + "</h1><p>Cyan marks the cell boundary. Each channel uses the same intensity range across this donor’s review crops, ending at their pooled 99.5th percentile.</p>"
            + "".join(records)
        )
        pages.append(
            f'<li><a href="{name}.html">{name}</a> ({len(records)} cells)</li>'
        )
        print(name, len(records), flush=True)
    (gallery / "index.html").write_text(
        '<meta charset="utf-8"><h1>V24 native-channel review</h1><ul>'
        + "".join(pages)
        + "</ul>"
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--project", type=Path, required=True)
    p.add_argument("--dest", type=Path, required=True)
    p.add_argument("--typing", type=Path, required=True)
    render(**vars(p.parse_args()))
