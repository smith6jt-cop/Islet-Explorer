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


def load_channel_names(channel_file: Optional[str] = None) -> Dict[int, str]:
    """Load channel names from file."""
    if channel_file is None:
        channel_file = Path(__file__).parent.parent.parent / "shiny_app" / "Channel_names"

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


def get_local_images(image_dir: Optional[str] = None) -> List[str]:
    """Get list of available local OME-TIFF images."""
    if image_dir is None:
        # Check multiple locations
        locations = [
            Path(__file__).parent.parent.parent / "shiny_app" / "www" / "local_images",
            Path(os.environ.get("LOCAL_IMAGE_ROOT", "")),
        ]
        for loc in locations:
            if loc.exists() and loc.is_dir():
                image_dir = loc
                break

    if image_dir is None:
        return []

    image_dir = Path(image_dir)
    if not image_dir.exists():
        return []

    # Find OME-TIFF files
    images = []
    for ext in ["*.ome.tif", "*.ome.tiff", "*.OME.TIFF", "*.OME.TIF"]:
        images.extend(image_dir.glob(ext))

    return sorted([img.name for img in images])


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

    def __init__(
        self,
        avivator_url: str = "/avivator/",
        image_base_url: str = "/images/"
    ):
        self.avivator_url = avivator_url
        self.image_base_url = image_base_url
        self.channel_names = load_channel_names()
        self.available_images = get_local_images()

        # Create widgets
        self.image_selector = pn.widgets.Select(
            name="Select Image",
            options=["-- Select an image --"] + self.available_images,
            value="-- Select an image --",
            width=300
        )

        self.refresh_button = pn.widgets.Button(
            name="Refresh",
            button_type="default",
            width=80
        )
        self.refresh_button.on_click(self._refresh_images)

        # Channel visibility controls
        self.channel_controls = self._create_channel_controls()

        # Viewer pane
        self.viewer_pane = pn.pane.HTML(
            self._get_placeholder_html(),
            sizing_mode="stretch_both",
            min_height=600
        )

        # Watch for image selection changes
        self.image_selector.param.watch(self._on_image_select, "value")

    def _create_channel_controls(self) -> pn.Column:
        """Create channel visibility controls."""
        controls = []

        # Add default channels
        for ch in DEFAULT_CHANNELS:
            name = ch["name"]
            color = ch["color"]

            checkbox = pn.widgets.Checkbox(
                name=name,
                value=ch["visible"],
                styles={"color": color}
            )
            color_picker = pn.widgets.ColorPicker(
                name="",
                value=color,
                width=40
            )
            row = pn.Row(checkbox, color_picker, sizing_mode="stretch_width")
            controls.append(row)

        return pn.Column(
            pn.pane.Markdown("### Channel Controls"),
            *controls,
            width=200
        )

    def _get_placeholder_html(self) -> str:
        """Get placeholder HTML when no image is selected."""
        return """
        <div style="
            display: flex;
            flex-direction: column;
            align-items: center;
            justify-content: center;
            height: 100%;
            min-height: 500px;
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
            <h2 style="margin-top: 20px; margin-bottom: 10px;">Image Viewer</h2>
            <p style="opacity: 0.8; text-align: center; max-width: 400px;">
                Select an OME-TIFF image from the dropdown above to visualize<br>
                multiplex imaging data with channel controls.
            </p>
        </div>
        """

    def _get_viewer_html(self, image_name: str) -> str:
        """Get HTML to embed Avivator viewer with image."""
        image_url = f"{self.image_base_url}{image_name}"

        # Create channel config
        channel_config = create_channel_config(DEFAULT_CHANNELS)

        return f"""
        <iframe
            id="avivator-frame"
            src="{self.avivator_url}?url={image_url}"
            style="
                width: 100%;
                height: 100%;
                min-height: 600px;
                border: none;
                border-radius: 8px;
            "
            allow="fullscreen"
        ></iframe>
        """

    def _on_image_select(self, event):
        """Handle image selection."""
        if event.new and event.new != "-- Select an image --":
            self.viewer_pane.object = self._get_viewer_html(event.new)
        else:
            self.viewer_pane.object = self._get_placeholder_html()

    def _refresh_images(self, event):
        """Refresh available images list."""
        self.available_images = get_local_images()
        self.image_selector.options = ["-- Select an image --"] + self.available_images

    def get_panel(self) -> pn.Column:
        """Return the complete viewer panel."""
        header = pn.Row(
            pn.pane.Markdown("## OME-TIFF Viewer"),
            pn.Spacer(),
            self.image_selector,
            self.refresh_button,
            sizing_mode="stretch_width",
            align="center"
        )

        # Main content with viewer and optional channel controls
        content = pn.Row(
            self.viewer_pane,
            sizing_mode="stretch_both"
        )

        return pn.Column(
            header,
            content,
            sizing_mode="stretch_both",
            styles={
                "background": "#f8f9fa",
                "border-radius": "8px",
                "padding": "15px"
            }
        )


class SimpleImageViewer:
    """Simplified image viewer that embeds external URL or local path."""

    def __init__(self):
        self.url_input = pn.widgets.TextInput(
            name="Image URL",
            placeholder="Enter OME-TIFF URL or select local image...",
            width=400
        )

        self.load_button = pn.widgets.Button(
            name="Load",
            button_type="primary",
            width=80
        )
        self.load_button.on_click(self._load_image)

        self.viewer_pane = pn.pane.HTML(
            self._get_placeholder_html(),
            sizing_mode="stretch_both",
            min_height=600
        )

    def _get_placeholder_html(self) -> str:
        return """
        <div style="
            display: flex;
            align-items: center;
            justify-content: center;
            height: 100%;
            min-height: 500px;
            background: #e5e7eb;
            border-radius: 8px;
            color: #374151;
        ">
            <p>Enter an image URL and click Load to view</p>
        </div>
        """

    def _load_image(self, event):
        url = self.url_input.value
        if url:
            self.viewer_pane.object = f"""
            <iframe
                src="https://avivator.gehlenborglab.org/?url={url}"
                style="width: 100%; height: 100%; min-height: 600px; border: none;"
            ></iframe>
            """

    def get_panel(self) -> pn.Column:
        controls = pn.Row(
            self.url_input,
            self.load_button,
            sizing_mode="stretch_width"
        )

        return pn.Column(
            controls,
            self.viewer_pane,
            sizing_mode="stretch_both"
        )


def create_image_viewer(
    avivator_url: str = "/avivator/",
    image_base_url: str = "/images/"
) -> pn.Column:
    """Create and return image viewer panel."""
    viewer = ImageViewerPanel(
        avivator_url=avivator_url,
        image_base_url=image_base_url
    )
    return viewer.get_panel()
