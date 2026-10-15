"""Provisional pixel-support diagnostics; these do not validate islet identity."""
import numpy as np
from scipy import ndimage


def hormone_support(raw, thresholds, pixel_size_um, smoothing_um=.6):
    """Apply explicit marker cutoffs to smoothed raw pixels, never display RGB.

    Exported cell-mean cutoffs are only a starting diagnostic when supplied here;
    they must not be represented as validated pixel-segmentation thresholds.
    """
    raw = np.asarray(raw, dtype=float)
    thresholds = np.asarray(thresholds, dtype=float)
    if (raw.ndim != 3 or raw.shape[-1] != len(thresholds)
            or not np.isfinite(raw).all() or not np.isfinite(thresholds).all()
            or (thresholds < 0).any() or pixel_size_um <= 0 or smoothing_um < 0):
        raise ValueError('Finite raw channels, cutoffs and positive calibration required')
    smoothed = ndimage.gaussian_filter(
        raw, sigma=(smoothing_um/pixel_size_um, smoothing_um/pixel_size_um, 0))
    return (smoothed > thresholds).any(axis=-1)


def signal_distances(support, pixel_size_um):
    support = np.asarray(support, dtype=bool)
    if support.ndim != 2 or pixel_size_um <= 0:
        raise ValueError('Two-dimensional support and positive calibration required')
    if not support.any():
        return np.full(support.shape, np.inf)
    return ndimage.distance_transform_edt(~support, sampling=pixel_size_um)


def mask_support_metrics(mask, support, pixel_size_um, tolerance_um=3):
    mask = np.asarray(mask, dtype=bool)
    support = np.asarray(support, dtype=bool)
    if mask.shape != support.shape or not mask.any() or tolerance_um < 0:
        raise ValueError('A nonempty mask aligned with support is required')
    distances = signal_distances(support, pixel_size_um)
    return dict(signal_pixel_fraction=float(support[mask].mean()),
                area_within_signal_distance_fraction=float((distances[mask] <= tolerance_um).mean()),
                distance_tolerance_um=tolerance_um,
                signal_distance_p90_um=float(np.quantile(distances[mask], .9)) if support.any() else None,
                max_signal_distance_um=float(distances[mask].max()) if support.any() else None)


def signal_footprint(support, pixel_size_um, close_um=2, max_distance_um=3,
                     max_hole_area_um2=80):
    """A bounded signal footprint for comparison, not a complete islet boundary.

    Small gaps/holes can be filled, but every output pixel must stay within the
    specified distance of measured signal. Large dark areas cannot be filled.
    """
    support = np.asarray(support, dtype=bool)
    if pixel_size_um <= 0 or min(close_um, max_distance_um, max_hole_area_um2) < 0:
        raise ValueError('Invalid physical units')
    if not support.any():
        return support.copy()
    radius = int(np.ceil(close_um/pixel_size_um))
    yy, xx = np.mgrid[-radius:radius+1, -radius:radius+1]
    disk = (xx**2+yy**2)*pixel_size_um**2 <= close_um**2
    padded = np.pad(support, radius+1)
    closed = ndimage.binary_closing(padded, structure=disk)
    closed = closed[radius+1:-radius-1, radius+1:-radius-1] | support
    holes = ndimage.binary_fill_holes(closed) & ~closed
    labels, _ = ndimage.label(holes)
    areas = np.bincount(labels.ravel())*pixel_size_um**2
    small = (areas <= max_hole_area_um2)
    small[0] = False
    result = closed | small[labels]
    return result & (signal_distances(support, pixel_size_um) <= max_distance_um)


def keep_seed_near_components(mask, seed_xy, pixel_size_um, radius_um=4):
    """Keep pixel components near operational seeds; do not bridge components."""
    labels, _ = ndimage.label(mask)
    selected = set()
    radius = int(np.ceil(radius_um/pixel_size_um))
    for x, y in np.asarray(seed_xy).reshape(-1, 2):
        x0, x1 = max(0, int(np.floor(x))-radius), min(mask.shape[1], int(np.ceil(x))+radius+1)
        y0, y1 = max(0, int(np.floor(y))-radius), min(mask.shape[0], int(np.ceil(y))+radius+1)
        if x0 >= x1 or y0 >= y1:
            continue
        yy, xx = np.mgrid[y0:y1, x0:x1]
        near = ((xx-x)**2+(yy-y)**2)*pixel_size_um**2 <= radius_um**2
        selected.update(np.unique(labels[y0:y1, x0:x1][near]).tolist())
    selected.discard(0)
    return np.isin(labels, list(selected))
