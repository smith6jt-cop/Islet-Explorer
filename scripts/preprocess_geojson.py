#!/usr/bin/env python3
"""
Preprocess massive GeoJSON files for efficient spatial queries in R Shiny.
Extracts only geometries and essential properties, creates spatial index.
"""

import gzip
import json
import pickle
import sys
from pathlib import Path
from typing import Dict, List, Any, Tuple
import argparse

def get_bbox(coordinates: List) -> Dict[str, float]:
    """Calculate bounding box from polygon coordinates."""
    if not coordinates or len(coordinates) == 0:
        return None
    
    # Get first ring of polygon
    ring = coordinates[0]
    if len(ring) < 3:
        return None
    
    xs = [p[0] for p in ring if len(p) >= 2]
    ys = [p[1] for p in ring if len(p) >= 2]
    
    if not xs or not ys:
        return None
    
    return {
        'xmin': min(xs),
        'ymin': min(ys),
        'xmax': max(xs),
        'ymax': max(ys)
    }

def simplify_geojson(input_path: Path, output_path: Path, essential_props: List[str] = None):
    """
    Simplify GeoJSON by stripping heavy properties and creating spatial index.
    
    Args:
        input_path: Path to input .geojson.gz file
        output_path: Path to output .pkl file
        essential_props: List of property names to keep (default: ['classification', 'name'])
    """
    if essential_props is None:
        essential_props = ['classification', 'name', 'Classification', 'Name']
    
    print(f"[PREPROCESS] Processing: {input_path.name}")
    print(f"[PREPROCESS] Loading GeoJSON...")
    
    # Load GeoJSON
    with gzip.open(input_path, 'rt', encoding='utf-8') as f:
        data = json.load(f)
    
    if 'features' not in data:
        raise ValueError("No features found in GeoJSON")
    
    features = data['features']
    n_features = len(features)
    print(f"[PREPROCESS] Found {n_features:,} features")
    
    # Simplify features
    print(f"[PREPROCESS] Simplifying features...")
    simplified = []
    spatial_index = []
    
    for i, feature in enumerate(features):
        if (i + 1) % 5000 == 0:
            progress = 100 * (i + 1) / n_features
            print(f"[PREPROCESS] Processed {i+1:,} / {n_features:,} features ({progress:.1f}%)")
        
        if 'geometry' not in feature or 'coordinates' not in feature['geometry']:
            continue
        
        geometry = feature['geometry']
        if geometry['type'] != 'Polygon':
            continue
        
        coords = geometry['coordinates']
        bbox = get_bbox(coords)
        
        if bbox is None:
            continue
        
        # Extract essential properties only
        props = feature.get('properties', {})
        essential = {}
        for prop in essential_props:
            if prop in props:
                essential[prop] = props[prop]
        
        simplified.append({
            'id': i + 1,
            'geometry': coords,
            'bbox': bbox,
            'properties': essential
        })
        
        spatial_index.append([
            bbox['xmin'], bbox['ymin'], bbox['xmax'], bbox['ymax']
        ])
    
    print(f"[PREPROCESS] Creating output structure...")
    
    # Calculate global bbox
    global_bbox = {
        'xmin': min(idx[0] for idx in spatial_index),
        'ymin': min(idx[1] for idx in spatial_index),
        'xmax': max(idx[2] for idx in spatial_index),
        'ymax': max(idx[3] for idx in spatial_index)
    }
    
    result = {
        'features': simplified,
        'n_features': len(simplified),
        'global_bbox': global_bbox,
        'spatial_index': spatial_index,
        'processed_date': str(Path(__file__).stat().st_mtime)
    }
    
    print(f"[PREPROCESS] Saving to {output_path}...")
    with open(output_path, 'wb') as f:
        pickle.dump(result, f, protocol=pickle.HIGHEST_PROTOCOL)
    
    # Report statistics
    original_size = input_path.stat().st_size / (1024 ** 2)
    new_size = output_path.stat().st_size / (1024 ** 2)
    reduction = 100 * (1 - new_size / original_size)
    
    print(f"[PREPROCESS] Complete! Size: {original_size:.1f} MB -> {new_size:.1f} MB ({reduction:.1f}% reduction)")
    print(f"[PREPROCESS] Retained {len(simplified):,} / {n_features:,} features")
    
    return result

def batch_process(gson_dir: Path, output_dir: Path):
    """Process all GeoJSON files in a directory."""
    output_dir.mkdir(parents=True, exist_ok=True)
    
    geojson_files = list(gson_dir.glob("*.geojson.gz"))
    print(f"[BATCH] Found {len(geojson_files)} files to process")
    
    for gf in geojson_files:
        base_name = gf.stem.replace('.geojson', '')
        output_pkl = output_dir / f"{base_name}_simplified.pkl"
        
        if output_pkl.exists():
            print(f"[BATCH] Skipping (already exists): {base_name}")
            continue
        
        try:
            simplify_geojson(gf, output_pkl)
        except Exception as e:
            print(f"[BATCH] ERROR processing {gf.name}: {e}")
            continue
    
    print("[BATCH] Complete!")

def main():
    parser = argparse.ArgumentParser(description='Preprocess GeoJSON files for efficient spatial queries')
    parser.add_argument('input', help='Input .geojson.gz file or directory')
    parser.add_argument('output', nargs='?', help='Output .pkl file or directory (optional for batch)')
    parser.add_argument('--batch', action='store_true', help='Batch process all files in directory')
    
    args = parser.parse_args()
    
    input_path = Path(args.input)
    
    if args.batch:
        output_dir = Path(args.output) if args.output else input_path.parent / 'gson_cache'
        batch_process(input_path, output_dir)
    else:
        if not args.output:
            output_path = input_path.with_suffix('').with_suffix('.pkl')
        else:
            output_path = Path(args.output)
        
        simplify_geojson(input_path, output_path)

if __name__ == '__main__':
    main()
