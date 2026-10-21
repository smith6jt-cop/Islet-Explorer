"""Make the Panc Floyd header logos from www/logo.png (white background, 1536 x 1024).

Writes www/logo-header-light.png (original colours, transparent background) and www/logo-header-dark.png (for the
dark theme: "Panc" in white and "Floyd" in grey; owner, 2026-10-21). Run with the app venv's Python:
    python scripts/app_assets/make_logo_variants.py
"""
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage

WWW = Path(__file__).resolve().parents[2] / "app" / "shiny_app" / "www"
HEIGHT = 128            # output height in px (shown at 64 px, so sharp on 2x screens)
DARK_TOP = np.array([0xff, 0xff, 0xff])      # "Panc" on the dark theme: white
DARK_BOTTOM = np.array([0xa6, 0xa6, 0xa6])   # "Floyd" on the dark theme: grey

rgb = np.asarray(Image.open(WWW / "logo.png").convert("RGB")).astype(float)
near_white = rgb.min(axis=2) > 235
lab, _ = ndimage.label(near_white)
border = set(np.unique(np.concatenate([lab[0], lab[-1], lab[:, 0], lab[:, -1]]))) - {0}
background = np.isin(lab, list(border))

# Opaque everywhere but the background; anti-aliased pixels next to it get alpha from their distance to white,
# with the white mixed out of their colour.
edge = ~background & ndimage.binary_dilation(background, iterations=3)
alpha = np.where(background, 0.0, 1.0)
a = np.clip((255 - rgb.min(axis=2)) / 255, 0, 1)
alpha[edge] = a[edge]
fg = rgb.copy()
safe = np.maximum(alpha, 1e-6)[..., None]
fg[edge] = np.clip((rgb[edge] - (1 - alpha[edge])[:, None] * 255) / safe[edge], 0, 255)

# The letters: dark components wholly right of the prism (the prism and its beam form one component at the left).
comp, n = ndimage.label(~background)
ys, xs = np.where(~background)
x0 = xs.min()
letters = np.zeros_like(background)
letter_parts = []
for i, sl in enumerate(ndimage.find_objects(comp), start=1):
    if sl is not None and sl[1].start > x0 + 0.35 * (xs.max() - x0):
        letters |= comp == i
        letter_parts.append((i, (sl[0].start + sl[0].stop) / 2))
# The two words are two lines: a letter whose vertical centre is above the middle of the lettering is in "Panc".
ly = np.where(letters)[0]
mid = (ly.min() + ly.max()) / 2
top_word = np.isin(comp, [i for i, cy in letter_parts if cy < mid])
bottom_word = letters & ~top_word

y_lo, y_hi, x_lo, x_hi = ys.min(), ys.max() + 1, xs.min(), xs.max() + 1
for name, recolour in (("logo-header-light.png", False), ("logo-header-dark.png", True)):
    out = fg.copy()
    if recolour:
        out[top_word] = DARK_TOP
        out[bottom_word] = DARK_BOTTOM
    rgba = np.dstack([out, alpha * 255]).round().astype(np.uint8)[y_lo:y_hi, x_lo:x_hi]
    img = Image.fromarray(rgba, "RGBA")
    img = img.resize((round(img.width * HEIGHT / img.height), HEIGHT), Image.LANCZOS)
    img.save(WWW / name, optimize=True)
    print(name, img.size)
