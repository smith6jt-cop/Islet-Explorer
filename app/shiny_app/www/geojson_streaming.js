// Efficient streaming GeoJSON overlay for AVIVATOR viewer
// Loads features on-demand based on viewport

(function() {
    'use strict';
    
    console.log('[GEOJSON-STREAM] Streaming overlay handler initialized');
    
    let overlayCanvas = null;
    let overlayCtx = null;
    let currentFeatures = [];
    let overlayEnabled = true;
    let caseId = null;
    let currentViewport = null;
    
    // Configuration
    const MAX_FEATURES_PER_QUERY = 5000;
    const OVERLAY_STYLE = {
        strokeColor: '#00FF00',
        fillColor: 'rgba(0, 255, 0, 0.1)',
        lineWidth: 1.5,
        opacity: 0.8
    };
    
    // Create overlay canvas
    function createOverlayCanvas() {
        if (overlayCanvas) return overlayCanvas;
        
        overlayCanvas = document.createElement('canvas');
        overlayCanvas.id = 'geojson-overlay-canvas';
        overlayCanvas.style.cssText = `
            position: absolute;
            top: 0;
            left: 0;
            pointer-events: none;
            z-index: 10;
            opacity: ${OVERLAY_STYLE.opacity};
        `;
        overlayCtx = overlayCanvas.getContext('2d');
        
        return overlayCanvas;
    }
    
    // Draw features on canvas
    function drawFeatures() {
        if (!overlayCanvas || !overlayCtx || !overlayEnabled || currentFeatures.length === 0) {
            return;
        }
        
        // Clear canvas
        overlayCtx.clearRect(0, 0, overlayCanvas.width, overlayCanvas.height);
        
        // Set drawing style
        overlayCtx.strokeStyle = OVERLAY_STYLE.strokeColor;
        overlayCtx.fillStyle = OVERLAY_STYLE.fillColor;
        overlayCtx.lineWidth = OVERLAY_STYLE.lineWidth;
        
        let drawnCount = 0;
        
        // Draw each feature
        currentFeatures.forEach(feature => {
            if (!feature.geometry || feature.geometry.type !== 'Polygon') return;
            
            const coords = feature.geometry.coordinates;
            if (!coords || coords.length === 0) return;
            
            coords.forEach(ring => {
                if (ring.length < 3) return;
                
                overlayCtx.beginPath();
                ring.forEach((point, i) => {
                    if (point.length >= 2) {
                        const x = point[0];
                        const y = point[1];
                        
                        if (i === 0) {
                            overlayCtx.moveTo(x, y);
                        } else {
                            overlayCtx.lineTo(x, y);
                        }
                    }
                });
                overlayCtx.closePath();
                overlayCtx.fill();
                overlayCtx.stroke();
                drawnCount++;
            });
        });
        
        console.log(`[GEOJSON-STREAM] Drew ${drawnCount} features`);
    }
    
    // Load features for viewport
    async function loadViewportFeatures(viewport) {
        if (!caseId) {
            console.warn('[GEOJSON-STREAM] No case ID set');
            return;
        }
        
        const params = new URLSearchParams({
            case_id: caseId,
            max_features: MAX_FEATURES_PER_QUERY
        });
        
        if (viewport) {
            params.append('xmin', viewport.xmin);
            params.append('ymin', viewport.ymin);
            params.append('xmax', viewport.xmax);
            params.append('ymax', viewport.ymax);
        }
        
        try {
            console.log('[GEOJSON-STREAM] Fetching features for viewport:', viewport);
            const response = await fetch(`/geojson_query?${params.toString()}`);
            
            if (!response.ok) {
                throw new Error(`HTTP ${response.status}: ${response.statusText}`);
            }
            
            const data = await response.json();
            currentFeatures = data.features || [];
            
            console.log(`[GEOJSON-STREAM] Loaded ${currentFeatures.length} / ${data.n_total || 0} features`);
            
            drawFeatures();
            showStatus(`Loaded ${currentFeatures.length} cells`, true);
            
        } catch (error) {
            console.error('[GEOJSON-STREAM] Failed to load features:', error);
            showStatus('Failed to load overlays', false);
        }
    }
    
    // Set up overlay canvas on viewer
    function setupOverlay() {
        const iframe = document.getElementById('avivator_viewer') || 
                      document.querySelector('iframe[src*="avivator"]');
        
        if (!iframe) {
            console.log('[GEOJSON-STREAM] No viewer iframe found, retrying...');
            setTimeout(setupOverlay, 1000);
            return;
        }
        
        const canvas = createOverlayCanvas();
        const rect = iframe.getBoundingClientRect();
        
        canvas.width = rect.width;
        canvas.height = rect.height;
        canvas.style.width = rect.width + 'px';
        canvas.style.height = rect.height + 'px';
        
        iframe.parentNode.style.position = 'relative';
        iframe.parentNode.appendChild(canvas);
        
        console.log('[GEOJSON-STREAM] Overlay canvas attached');
    }
    
    // Handle messages from Shiny
    Shiny.addCustomMessageHandler('loadGeoJSONOverlay', function(message) {
        console.log('[GEOJSON-STREAM] Received load command:', message);
        
        caseId = message.case_id;
        currentViewport = message.viewport || null;
        
        setupOverlay();
        loadViewportFeatures(currentViewport);
    });
    
    Shiny.addCustomMessageHandler('toggleGeoJSONOverlay', function(message) {
        overlayEnabled = message.enabled;
        console.log('[GEOJSON-STREAM] Toggle overlay:', overlayEnabled);
        
        if (overlayCanvas) {
            overlayCanvas.style.display = overlayEnabled ? 'block' : 'none';
        }
        
        showStatus(overlayEnabled ? 'Overlays ON' : 'Overlays OFF', overlayEnabled);
    });
    
    Shiny.addCustomMessageHandler('updateGeoJSONViewport', function(message) {
        currentViewport = message.viewport;
        loadViewportFeatures(currentViewport);
    });
    
    // Show status message
    function showStatus(message, isSuccess) {
        const existing = document.getElementById('overlay-status');
        if (existing) existing.remove();
        
        const status = document.createElement('div');
        status.id = 'overlay-status';
        status.style.cssText = `
            position: fixed;
            top: 60px;
            right: 10px;
            background: ${isSuccess ? 'rgba(0, 150, 0, 0.9)' : 'rgba(150, 0, 0, 0.9)'};
            color: white;
            padding: 8px 12px;
            border-radius: 4px;
            font-size: 12px;
            z-index: 1000;
            box-shadow: 0 2px 6px rgba(0,0,0,0.3);
        `;
        status.textContent = message;
        document.body.appendChild(status);
        
        setTimeout(() => status.remove(), 3000);
    }
    
    window.showOverlayStatus = showStatus;
    
})();
