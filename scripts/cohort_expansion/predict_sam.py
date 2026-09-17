#!/usr/bin/env python3
"""Run the original SAM VIT_B box protocol synchronously with local model weights."""
import hashlib
import json
import os
from pathlib import Path
import sys
import time

import numpy as np
from PIL import Image
import torch
from mobile_sam import SamPredictor, sam_model_registry
from samapi.utils import mask_to_geometry

root = Path(sys.argv[1])
metadata = json.loads((root / 'sam_prompts.json').read_text())
checkpoint = Path('/home/smith6jt/.samapi/vit_b/sam_vit_b_01ec64.pth')
torch.set_num_threads(32)
device = os.environ.get('ISLET_SAM_DEVICE', 'cpu')
model = sam_model_registry['vit_b'](checkpoint=str(checkpoint)).to(device=device)
predictor = SamPredictor(model)
pixels = np.asarray(Image.open(root / 'sam_render.png').convert('RGB'))
predictor.set_image(pixels)
features = []
started = time.time()
for index, prompt in enumerate(metadata['boxes']):
    masks, quality, _ = predictor.predict(box=np.asarray(prompt['bbox']), multimask_output=False)
    geometry = mask_to_geometry(masks[0], scale=metadata['downsample'])
    if geometry.get('coordinates'):
        features.append({'type': 'Feature', 'geometry': geometry,
                         'properties': {'objectType': 'annotation',
                                        'name': 'Islet', 'sam_quality': float(quality[0]),
                                        'source_prompt_index': prompt['index'],
                                        'measurements': {'SAM: predicted mask quality': float(quality[0]),
                                                         'SAM: source prompt index': prompt['index'],
                                                         'Islet: identity reviewed': 0}}})
    if index % 25 == 0 or index + 1 == len(metadata['boxes']):
        print(f'SAM {index + 1}/{len(metadata["boxes"])} prompts, {time.time() - started:.1f}s', flush=True)
if not features:
    raise RuntimeError('SAM returned no islet geometries')
target = root / 'sam_islets.geojson'
target.write_text(json.dumps({'type': 'FeatureCollection', 'features': features}) + '\n')
provenance = {'model': 'vit_b', 'checkpoint': str(checkpoint),
              'checkpoint_sha256': hashlib.sha256(checkpoint.read_bytes()).hexdigest(),
              'device': device, 'torch_version': torch.__version__, 'downsample': metadata['downsample'],
              'prompt_count': len(metadata['boxes']), 'mask_count': len(features),
              'output_sha256': hashlib.sha256(target.read_bytes()).hexdigest()}
(root / 'sam_provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
