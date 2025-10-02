// GeoJSON Overlay Handler for Avivator viewer
// Handles client-side overlay rendering to avoid iframe reloading

(function() {
    'use strict';
    
    console.log('[OVERLAY] GeoJSON overlay handler loaded');
    
    let overlayData = null;
    let overlayEnabled = true;
    let overlayCanvas = null;
    
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
            opacity: 0.8;
        `;
        
        return overlayCanvas;
    }
    
    // Draw GeoJSON features on canvas
    function drawOverlay() {
        if (!overlayCanvas || !overlayData || !overlayEnabled) return;
        
        const ctx = overlayCanvas.getContext('2d');
        ctx.clearRect(0, 0, overlayCanvas.width, overlayCanvas.height);
        
        if (!overlayData.features) return;
        
        ctx.strokeStyle = '#00FF00';
        ctx.lineWidth = 2;
        ctx.globalAlpha = 0.8;
        
        overlayData.features.forEach(feature => {
            if (feature.geometry && feature.geometry.type === 'Polygon') {
                feature.geometry.coordinates.forEach(polygon => {
                    if (polygon.length < 3) return;
                    
                    ctx.beginPath();
                    polygon.forEach((point, i) => {
                        if (point.length >= 2) {
                            const x = point[0];
                            const y = point[1];
                            if (i === 0) {
                                ctx.moveTo(x, y);
                            } else {
                                ctx.lineTo(x, y);
                            }
                        }
                    });
                    ctx.closePath();
                    ctx.stroke();
                });
            }
        });
    }
    
    // Load GeoJSON data
    async function loadGeoJSON() {
        const params = new URLSearchParams(window.location.search);
        const geojsonUrl = params.get('geojson_url');
        
        if (!geojsonUrl) {
            console.log('[OVERLAY] No GeoJSON URL found');
            return;
        }
        
        try {
            console.log('[OVERLAY] Loading GeoJSON from:', geojsonUrl);
            const response = await fetch(decodeURIComponent(geojsonUrl));
            
            if (!response.ok) {
                throw new Error(`HTTP ${response.status}: ${response.statusText}`);
            }
            
            overlayData = await response.json();
            console.log('[OVERLAY] Loaded', overlayData.features?.length || 0, 'features');
            
            // Set up overlay canvas
            setupOverlay();
            
        } catch (error) {
            console.error('[OVERLAY] Failed to load GeoJSON:', error);
            showStatus('Failed to load overlays', false);
        }
    }
    
    // Set up overlay on the viewer
    function setupOverlay() {
        // Find the Avivator iframe or container
        const iframe = document.getElementById('avivator_viewer') || 
                      document.querySelector('iframe') ||
                      document.querySelector('[id*="avivator"]');
        
        if (!iframe) {
            console.log('[OVERLAY] No viewer iframe found, retrying...');
            setTimeout(setupOverlay, 1000);
            return;
        }
        
        // Create overlay canvas
        const canvas = createOverlayCanvas();
        
        // Position canvas over iframe
        const rect = iframe.getBoundingClientRect();
        canvas.width = rect.width;
        canvas.height = rect.height;
        canvas.style.width = rect.width + 'px';
        canvas.style.height = rect.height + 'px';
        
        // Position relative to iframe
        iframe.style.position = 'relative';
        iframe.parentNode.insertBefore(canvas, iframe.nextSibling);
        
        // Draw overlay
        drawOverlay();
        showStatus('Islet overlays loaded', true);
    }
    
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
    
    // Handle overlay toggle messages
    window.addEventListener('message', (event) => {
        if (event.data && event.data.type === 'toggleOverlay') {
            overlayEnabled = event.data.enabled;
            console.log('[OVERLAY] Toggle overlay:', overlayEnabled);
            
            if (overlayCanvas) {
                overlayCanvas.style.display = overlayEnabled ? 'block' : 'none';
            }
            
            showStatus(overlayEnabled ? 'Overlays ON' : 'Overlays OFF', overlayEnabled);
        }
    });
    
    // Global function for direct access
    window.showOverlayStatus = showStatus;
    
    // Initialize
    if (document.readyState === 'loading') {
        document.addEventListener('DOMContentLoaded', loadGeoJSON);
    } else {
        setTimeout(loadGeoJSON, 500);
    }
    
})();