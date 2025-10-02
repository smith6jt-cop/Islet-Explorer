# GeoJSON Overlay System

## Overview

The GeoJSON overlay system provides efficient visualization of cell segmentation data (200K+ features per image) over OME-TIFF images in the AVIVATOR viewer. This system was specifically designed to handle massive QuPath-exported GeoJSON files (765MB compressed each) without browser crashes or performance issues.

## The Challenge

Traditional approaches to displaying GeoJSON overlays face several problems:

| Problem | Impact |
|---------|--------|
| Large file sizes | 765 MB compressed, ~3-4 GB uncompressed |
| Many features | 187,000+ cell boundaries per image |
| Browser memory limits | Loading full file causes crashes |
| Network latency | Slow transfer of large files |
| Redundant data | 100+ measurement properties per cell |

## Architecture

The solution uses a **three-tier architecture**: preprocessing, server-side querying, and client-side streaming.

```
┌─────────────────────────────────────────────────────────┐
│  TIER 1: Preprocessing (One-time, offline)             │
│                                                          │
│  Input:  gson/*.geojson.gz (765MB, ~187K features)     │
│           ↓                                              │
│  Python: Extract geometries + essential properties      │
│          Build spatial index (bounding boxes)           │
│           ↓                                              │
│  Output: data/gson_cache/*.pkl (199MB, 74% reduction)  │
└─────────────────────────────────────────────────────────┘
                       ↓
┌─────────────────────────────────────────────────────────┐
│  TIER 2: Server-Side Query (Runtime, R Shiny)          │
│                                                          │
│  Request: case_id, viewport bounds (optional)           │
│           ↓                                              │
│  Load:    Preprocessed .pkl file (fast!)                │
│  Query:   Spatial index for bounding box intersection   │
│  Limit:   Max 5000 features per response                │
│           ↓                                              │
│  Return:  GeoJSON FeatureCollection (JSON)              │
└─────────────────────────────────────────────────────────┘
                       ↓
┌─────────────────────────────────────────────────────────┐
│  TIER 3: Client-Side Rendering (Browser, JavaScript)   │
│                                                          │
│  Receive: Feature subset from server                    │
│  Draw:    HTML5 canvas overlay on viewer                │
│  Toggle:  Show/hide without reloading                   │
│  Update:  Request new features on pan/zoom (future)     │
└─────────────────────────────────────────────────────────┘
```

## Performance Achievements

### File Size Reduction

```python
# Before preprocessing
Original: 765 MB (compressed)
Expanded: ~3-4 GB (in memory)
Features: 187,184

# After preprocessing  
Cached:   199 MB (.pkl format)
Features: 186,940 (99.9% retention)
Reduction: 74%
```

### Runtime Performance

| Operation | Time | Comparison |
|-----------|------|------------|
| Load full GeoJSON | 30-60+ sec | Baseline |
| Load preprocessed cache | <1 sec | **30-60× faster** |
| Query viewport (5K features) | <100 ms | **300-600× faster** |
| Render on canvas | <50 ms | Instant |

### Memory Usage

| Component | Before | After | Savings |
|-----------|--------|-------|---------|
| Server (R) | ~4 GB | ~200 MB | 95% |
| Client (Browser) | ~4 GB | ~50 MB | 98% |
| Network transfer | 765 MB | ~2-5 MB | 99% |

## Implementation Details

### Preprocessing (Python)

**File**: `scripts/preprocess_geojson.py`

```python
def simplify_geojson(input_path, output_path):
    """
    Extract minimal data from GeoJSON:
    - Geometry coordinates only
    - Essential properties (classification, name)
    - Spatial index (bounding boxes)
    """
    # Load full GeoJSON
    with gzip.open(input_path, 'rt') as f:
        data = json.load(f)
    
    # Process each feature
    simplified = []
    for i, feature in enumerate(data['features']):
        coords = feature['geometry']['coordinates']
        bbox = calculate_bbox(coords)
        
        simplified.append({
            'id': i,
            'geometry': coords,
            'bbox': bbox,
            'properties': {
                'classification': feature['properties'].get('classification'),
                'name': feature['properties'].get('name')
            }
        })
    
    # Save as pickle
    with open(output_path, 'wb') as f:
        pickle.dump({
            'features': simplified,
            'spatial_index': [f['bbox'] for f in simplified],
            'n_features': len(simplified)
        }, f)
```

**Why Python?**
- No 2GB string limit (unlike R)
- Fast JSON parsing
- Efficient memory management
- Standard library only (no dependencies)

### Server-Side Query (R)

**File**: `app/shiny_app/app.R`

```r
query_geojson_features <- function(case_id, xmin=NULL, ymin=NULL, 
                                   xmax=NULL, ymax=NULL, max_features=5000) {
  # Load preprocessed cache
  pkl_file <- file.path(gson_cache_dir, paste0(case_id, "_simplified.pkl"))
  data <- load_geojson_cache(pkl_file)  # Uses Python via helper
  
  # Query spatial index
  if (!is.null(xmin)) {
    intersects <- (
      data$spatial_index[, "xmax"] >= xmin &
      data$spatial_index[, "xmin"] <= xmax &
      data$spatial_index[, "ymax"] >= ymin &
      data$spatial_index[, "ymin"] <= ymax
    )
    matching_ids <- which(intersects)
  } else {
    matching_ids <- seq_len(min(max_features, data$n_features))
  }
  
  # Limit and format
  if (length(matching_ids) > max_features) {
    matching_ids <- sample(matching_ids, max_features)
  }
  
  features <- lapply(data$features[matching_ids], format_feature)
  return(list(type="FeatureCollection", features=features))
}
```

### Client-Side Rendering (JavaScript)

**File**: `app/shiny_app/www/geojson_streaming.js`

```javascript
// Load features for current image
Shiny.addCustomMessageHandler('loadGeoJSONOverlay', async function(msg) {
    const params = new URLSearchParams({
        case_id: msg.case_id,
        max_features: 5000
    });
    
    const response = await fetch(`/geojson_query?${params}`);
    const data = await response.json();
    
    currentFeatures = data.features;
    drawFeatures();
});

// Draw on HTML5 canvas
function drawFeatures() {
    const ctx = overlayCanvas.getContext('2d');
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    
    ctx.strokeStyle = '#00FF00';
    ctx.lineWidth = 1.5;
    
    currentFeatures.forEach(feature => {
        feature.geometry.coordinates.forEach(ring => {
            ctx.beginPath();
            ring.forEach((point, i) => {
                if (i === 0) ctx.moveTo(point[0], point[1]);
                else ctx.lineTo(point[0], point[1]);
            });
            ctx.closePath();
            ctx.stroke();
        });
    });
}
```

## Data Flow Example

### User selects image "ND_0112.ome.tiff"

1. **App extracts case ID**: `0112`
2. **Server checks cache**: `data/gson_cache/0112_simplified.pkl` exists
3. **Server sends message**: `loadGeoJSONOverlay{case_id: "0112"}`
4. **Client requests**: `GET /geojson_query?case_id=0112&max_features=5000`
5. **Server loads**: 199 MB pickle file (<1 sec)
6. **Server queries**: Spatial index for first 5000 features
7. **Server returns**: ~2-5 MB JSON (5000 features)
8. **Client receives**: GeoJSON FeatureCollection
9. **Client renders**: 5000 polygons on canvas (<50 ms)

**Total time**: <2 seconds vs 30-60+ seconds for full file

## Spatial Indexing

### Bounding Box Structure

Each feature's bounding box is stored as:

```python
bbox = {
    'xmin': min(coordinates[:, 0]),
    'ymin': min(coordinates[:, 1]),
    'xmax': max(coordinates[:, 0]),
    'ymax': max(coordinates[:, 1])
}
```

### Intersection Test

To find features in a viewport:

```python
def intersects(feature_bbox, viewport_bbox):
    return (
        feature_bbox['xmax'] >= viewport_bbox['xmin'] and
        feature_bbox['xmin'] <= viewport_bbox['xmax'] and
        feature_bbox['ymax'] >= viewport_bbox['ymin'] and
        feature_bbox['ymin'] <= viewport_bbox['ymax']
    )
```

### Complexity

- **Index size**: O(n) where n = number of features
- **Query time**: O(n) linear scan (fast enough for ~187K features)
- **Memory**: 32 bytes per feature (4 × 8-byte floats)

**Future optimization**: R-tree or quad-tree for O(log n) queries

## File Formats

### Input: GeoJSON (.geojson.gz)

```json
{
  "type": "FeatureCollection",
  "features": [
    {
      "type": "Feature",
      "geometry": {
        "type": "Polygon",
        "coordinates": [[[x1,y1], [x2,y2], ...]]
      },
      "properties": {
        "classification": "Immune",
        "name": "Cell_12345",
        "marker1_intensity": 123.45,
        "marker2_intensity": 67.89,
        ... // 100+ more properties
      }
    },
    ... // 187,000 more features
  ]
}
```

### Output: Python Pickle (.pkl)

```python
{
    'features': [
        {
            'id': 1,
            'geometry': [[[x1,y1], [x2,y2], ...]],
            'bbox': {'xmin': ..., 'ymin': ..., 'xmax': ..., 'ymax': ...},
            'properties': {
                'classification': 'Immune',
                'name': 'Cell_12345'
            }
        },
        ...
    ],
    'spatial_index': [
        [xmin1, ymin1, xmax1, ymax1],
        [xmin2, ymin2, xmax2, ymax2],
        ...
    ],
    'n_features': 186940,
    'global_bbox': {'xmin': -55.99, 'ymin': -55.99, ...},
    'processed_date': '2025-10-02T13:22:45'
}
```

## Configuration

### Adjust Feature Limit

Edit `app/shiny_app/www/geojson_streaming.js`:

```javascript
const MAX_FEATURES_PER_QUERY = 5000;  // Reduce if slow
```

**Recommendations**:
- Desktop: 5000-10000
- Mobile: 1000-2000
- High-res displays: 5000
- Low-end devices: 500-1000

### Adjust Overlay Style

```javascript
const OVERLAY_STYLE = {
    strokeColor: '#00FF00',              // Outline color
    fillColor: 'rgba(0, 255, 0, 0.1)',   // Fill color (transparent)
    lineWidth: 1.5,                       // Line thickness
    opacity: 0.8                          // Overall opacity
};
```

## Future Enhancements

### 1. Viewport-Based Streaming (Not Yet Implemented)

Automatically load new features as user pans/zooms:

```javascript
// Detect viewport changes
viewerIframe.addEventListener('viewportchange', (e) => {
    const viewport = e.detail;
    updateGeoJSONViewport(viewport.xmin, viewport.ymin, 
                         viewport.xmax, viewport.ymax);
});
```

### 2. Level-of-Detail (LOD) System

Show simplified geometries when zoomed out:

```python
def simplify_geometry(coords, zoom_level):
    if zoom_level < 0.5:
        return douglas_peucker(coords, epsilon=10)
    elif zoom_level < 1.0:
        return douglas_peucker(coords, epsilon=5)
    else:
        return coords  # Full detail
```

### 3. WebGL Rendering

Use deck.gl for GPU-accelerated rendering of 100K+ features:

```javascript
import {PolygonLayer} from '@deck.gl/layers';

const layer = new PolygonLayer({
    data: features,
    getPolygon: d => d.geometry.coordinates[0],
    stroked: true,
    lineWidthMinPixels: 1
});
```

### 4. Vector Tiles

Pre-render as Mapbox Vector Tiles (.mvt) for even better performance:

```bash
tippecanoe -o cells.mbtiles -l cells --drop-densest-as-needed \\
           --maximum-zoom=16 cells.geojson
```

## Troubleshooting

See [Troubleshooting Guide](../reference/troubleshooting.md) for common issues and solutions.

## Related Documentation

- [Preprocessing Guide](preprocessing.md) - Detailed preprocessing workflow
- [API Reference](api_reference.md) - Function signatures and parameters
- [Data Format](data_format.md) - Input/output data specifications
