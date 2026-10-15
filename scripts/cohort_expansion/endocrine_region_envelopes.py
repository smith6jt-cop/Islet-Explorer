"""Experimental whole-region envelopes from hormone foreground, without SAM.

These resolve a geometry question only. Cutoffs, connectivity and biological
identity still require reference annotations before any cohort-wide delivery.
"""
import numpy as np
from scipy import ndimage

from endocrine_signal_support import keep_seed_near_components


def region_envelopes(support, seed_xy, pixel_size_um, closing_radius_um=2,
                     seed_proximity_um=4):
    """Join only short foreground gaps, then fill enclosed component interiors.

    External dark spaces stay external. Filling is separate from foreground
    detection, so vascular/nuclear interiors need not split a whole structure.
    This is an envelope proposal, not evidence that every enclosed pixel is
    endocrine or that the component is a biological islet.
    """
    support = np.asarray(support, dtype=bool)
    if (support.ndim != 2 or not np.isfinite(pixel_size_um) or pixel_size_um <= 0
            or not np.isfinite(closing_radius_um) or closing_radius_um < 0
            or not np.isfinite(seed_proximity_um) or seed_proximity_um < 0):
        raise ValueError('Two-dimensional signal and valid physical distances required')
    if not support.any():
        return np.zeros(support.shape, dtype=np.int32)
    radius = int(np.ceil(closing_radius_um/pixel_size_um))
    yy, xx = np.mgrid[-radius:radius+1, -radius:radius+1]
    disk = (xx**2+yy**2)*pixel_size_um**2 <= closing_radius_um**2
    pad = radius+1
    closed = ndimage.binary_closing(np.pad(support,pad), structure=disk)[pad:-pad,pad:-pad] | support
    closed = keep_seed_near_components(closed,seed_xy,pixel_size_um,seed_proximity_um)
    labels, _ = ndimage.label(closed)
    envelopes = np.zeros_like(support)
    for label, slices in enumerate(ndimage.find_objects(labels),1):
        if slices is None:
            continue
        local = labels[slices] == label
        envelopes[slices] |= ndimage.binary_fill_holes(local)
    return ndimage.label(envelopes)[0]
