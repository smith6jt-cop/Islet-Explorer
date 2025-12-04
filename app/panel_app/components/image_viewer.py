"""
Image viewer component for Islet Explorer Panel App.
Embeds AVIVATOR viewer for OME-TIFF visualization.
"""

import os
import json
import base64
from pathlib import Path
from typing import Optional, List, Dict
import panel as pn

# Channel configuration
CHANNEL_COLORS = {
    "DAPI": "#9CA3AF",   # grey
    "INS": "#EF4444",    # red
    "GCG": "#3B82F6",    # blue
    "SST": "#F59E0B",    # yellow
    "CD31": "#10B981",   # green
    "CD3e": "#8B5CF6",   # purple
    "CD4": "#EC4899",    # pink
    "CD8a": "#06B6D4",   # cyan
}

DEFAULT_CHANNELS = [
    {"name": "DAPI", "color": "#9CA3AF", "visible": True, "index": 0},
    {"name": "INS", "color": "#EF4444", "visible": True, "index": 25},
    {"name": "GCG", "color": "#3B82F6", "visible": True, "index": 19},
    {"name": "SST", "color": "#F59E0B", "visible": True, "index": 13},
]

# Paths
PANEL_APP_DIR = Path(__file__).parent.parent
SHINY_APP_DIR = PANEL_APP_DIR.parent / "shiny_app"
AVIVATOR_DIR = SHINY_APP_DIR / "www" / "avivator"
LOCAL_IMAGES_DIR = SHINY_APP_DIR / "www" / "local_images"


def load_channel_names(channel_file: Optional[str] = None) -> Dict[int, str]:
    """Load channel names from file."""
    if channel_file is None:
        channel_file = SHINY_APP_DIR / "Channel_names"

    channel_names = {}
    try:
        with open(channel_file, "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                # Parse format: "NAME (C#)"
                if "(" in line and ")" in line:
                    name = line.split("(")[0].strip()
                    idx_str = line.split("(")[1].replace(")", "").replace("C", "").strip()
                    try:
                        idx = int(idx_str) - 1  # Convert to 0-indexed
                        channel_names[idx] = name
                    except ValueError:
                        continue
    except Exception:
        pass

    return channel_names


def get_local_images(image_dir: Optional[Path] = None) -> List[str]:
    """Get list of available local OME-TIFF images."""
    if image_dir is None:
        # Check multiple locations
        locations = [
            LOCAL_IMAGES_DIR,
            Path(os.environ.get("LOCAL_IMAGE_ROOT", "")),
        ]
        for loc in locations:
            if loc and loc.exists() and loc.is_dir():
                image_dir = loc
                break

    if image_dir is None or not image_dir.exists():
        return []

    # Find OME-TIFF files
    images = []
    for ext in ["*.ome.tif", "*.ome.tiff", "*.OME.TIFF", "*.OME.TIF"]:
        images.extend(image_dir.glob(ext))

    return sorted([img.name for img in images])


def get_static_dirs() -> Dict[str, str]:
    """Get static directory mappings for Panel server."""
    static_dirs = {}

    # AVIVATOR viewer
    if AVIVATOR_DIR.exists():
        static_dirs["/avivator"] = str(AVIVATOR_DIR)

    # Local images
    if LOCAL_IMAGES_DIR.exists():
        static_dirs["/images"] = str(LOCAL_IMAGES_DIR)

    # Check environment variable for images
    env_image_root = os.environ.get("LOCAL_IMAGE_ROOT")
    if env_image_root and Path(env_image_root).exists():
        static_dirs["/images"] = env_image_root

    return static_dirs


def create_channel_config(channels: List[Dict]) -> str:
    """Create base64-encoded channel configuration for Avivator."""
    config = []
    for ch in channels:
        config.append({
            "name": ch.get("name", f"Channel {ch.get('index', 0)}"),
            "color": ch.get("color", "#FFFFFF"),
            "visible": ch.get("visible", True),
            "contrastLimits": ch.get("contrastLimits", [0, 65535]),
            "domain": ch.get("domain", [0, 65535])
        })
    return base64.b64encode(json.dumps(config).encode()).decode()


class ImageViewerPanel:
    """Panel-based image viewer using Avivator."""

    def __init__(self):
        self.channel_names = load_channel_names()
        self.available_images = get_local_images()
        self._avivator_available = AVIVATOR_DIR.exists()
        self._images_available = len(self.available_images) > 0

        # URL input for external images
        self.url_input = pn.widgets.TextInput(
            name="Image URL",
            placeholder="Enter OME-TIFF URL (http:// or https://)",
            width=400
        )

        # Create widgets
        self.image_selector = pn.widgets.Select(
            name="Local Images",
            options=["-- Select local image --"] + self.available_images,
            value="-- Select local image --",
            width=300,
            disabled=not self._images_available
        )

        self.load_url_button = pn.widgets.Button(
            name="Load URL",
            button_type="primary",
            width=100
        )
        self.load_url_button.on_click(self._on_load_url)

        self.refresh_button = pn.widgets.Button(
            name="Refresh",
            button_type="default",
            width=80
        )
        self.refresh_button.on_click(self._refresh_images)

        # Viewer pane
        self.viewer_pane = pn.pane.HTML(
            self._get_placeholder_html(),
            sizing_mode="stretch_both",
            min_height=700
        )

        # Status indicator
        self.status_pane = pn.pane.Markdown(
            self._get_status_text(),
            styles={"font-size": "12px", "color": "#666"}
        )

        # Watch for image selection changes
        self.image_selector.param.watch(self._on_image_select, "value")

    def _get_status_text(self) -> str:
        """Get status text showing available resources."""
        status = []
        if self._avivator_available:
            status.append("AVIVATOR: Ready")
        else:
            status.append("AVIVATOR: Not found (using external)")

        if self._images_available:
            status.append(f"Local images: {len(self.available_images)}")
        else:
            status.append("Local images: None found")

        return " | ".join(status)

    def _get_placeholder_html(self) -> str:
        """Get placeholder HTML when no image is selected."""
        return """
        <div style="
            display: flex;
            flex-direction: column;
            align-items: center;
            justify-content: center;
            height: 100%;
            min-height: 600px;
            background: linear-gradient(135deg, #1e3a72 0%, #2c5aa0 100%);
            border-radius: 8px;
            color: white;
            font-family: system-ui, -apple-system, sans-serif;
        ">
            <svg width="80" height="80" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5">
                <rect x="3" y="3" width="18" height="18" rx="2"/>
                <circle cx="8.5" cy="8.5" r="1.5"/>
                <path d="M21 15l-5-5L5 21"/>
            </svg>
            <h2 style="margin-top: 20px; margin-bottom: 10px;">AVIVATOR Image Viewer</h2>
            <p style="opacity: 0.8; text-align: center; max-width: 400px; line-height: 1.6;">
                Select a local OME-TIFF image from the dropdown<br>
                or enter a URL to an image hosted elsewhere.
            </p>
            <p style="opacity: 0.6; font-size: 12px; margin-top: 20px;">
                Supports OME-TIFF and OME-Zarr formats
            </p>
        </div>
        """

    def _get_viewer_html(self, image_url: str, use_local_avivator: bool = True) -> str:
        """Get HTML to embed Avivator viewer with image."""

        # Determine which AVIVATOR to use
        if use_local_avivator and self._avivator_available:
            avivator_base = "/avivator/index.html"
        else:
            # Fall back to public AVIVATOR
            avivator_base = "https://avivator.gehlenborglab.org/"

        # Encode the URL properly
        from urllib.parse import quote
        encoded_url = quote(image_url, safe='')

        return f"""
        <iframe
            id="avivator-frame"
            src="{avivator_base}?url={encoded_url}"
            style="
                width: 100%;
                height: 100%;
                min-height: 700px;
                border: none;
                border-radius: 8px;
                background: #000;
            "
            allow="fullscreen"
            sandbox="allow-scripts allow-same-origin allow-popups allow-forms"
        ></iframe>
        """

    def _on_image_select(self, event):
        """Handle local image selection."""
        if event.new and event.new != "-- Select local image --":
            # Use relative URL for local images
            image_url = f"/images/{event.new}"
            self.viewer_pane.object = self._get_viewer_html(image_url, use_local_avivator=True)
            self.url_input.value = ""
        elif not self.url_input.value:
            self.viewer_pane.object = self._get_placeholder_html()

    def _on_load_url(self, event):
        """Handle URL load button click."""
        url = self.url_input.value.strip()
        if url:
            # Reset local image selector
            self.image_selector.value = "-- Select local image --"
            # Load external URL (use external AVIVATOR for CORS)
            self.viewer_pane.object = self._get_viewer_html(url, use_local_avivator=False)

    def _refresh_images(self, event):
        """Refresh available images list."""
        self.available_images = get_local_images()
        self._images_available = len(self.available_images) > 0
        self.image_selector.options = ["-- Select local image --"] + self.available_images
        self.image_selector.disabled = not self._images_available
        self.status_pane.object = self._get_status_text()

    def get_panel(self) -> pn.Column:
        """Return the complete viewer panel."""
        # Header with title
        header = pn.pane.Markdown(
            "## OME-TIFF Viewer",
            styles={"margin": "0 0 10px 0"}
        )

        # Local image controls
        local_controls = pn.Row(
            self.image_selector,
            self.refresh_button,
            align="end"
        )

        # URL input controls
        url_controls = pn.Row(
            self.url_input,
            self.load_url_button,
            align="end"
        )

        # Combined controls
        controls = pn.Column(
            pn.pane.Markdown("**Load from URL:**"),
            url_controls,
            pn.Spacer(height=10),
            pn.pane.Markdown("**Or select local image:**"),
            local_controls,
            pn.Spacer(height=5),
            self.status_pane,
            sizing_mode="stretch_width"
        )

        return pn.Column(
            header,
            controls,
            pn.Spacer(height=10),
            self.viewer_pane,
            sizing_mode="stretch_both",
            styles={
                "background": "#f8f9fa",
                "border-radius": "8px",
                "padding": "15px"
            }
        )


def create_image_viewer() -> pn.Column:
    """Create and return image viewer panel."""
    viewer = ImageViewerPanel()
    return viewer.get_panel()
