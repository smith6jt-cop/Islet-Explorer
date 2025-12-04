"""
Image viewer component for Islet Explorer Panel App.
Embeds AVIVATOR viewer for OME-TIFF visualization.
"""

import os
from pathlib import Path
from typing import Optional, List, Dict
import panel as pn

# Paths - use panel_app directories (self-contained)
PANEL_APP_DIR = Path(__file__).parent.parent
ASSETS_DIR = PANEL_APP_DIR / "assets"
LOCAL_IMAGES_DIR = PANEL_APP_DIR / "local_images"

# Public AVIVATOR URL
AVIVATOR_URL = "https://avivator.gehlenborglab.org/"


def get_local_images(image_dir: Optional[Path] = None) -> List[str]:
    """Get list of available local OME-TIFF images."""
    if image_dir is None:
        # Check environment variable FIRST (user override)
        env_root = os.environ.get("LOCAL_IMAGE_ROOT", "")
        if env_root:
            env_path = Path(env_root)
            if env_path.exists() and env_path.is_dir():
                image_dir = env_path
                print(f"[ImageViewer] Using LOCAL_IMAGE_ROOT: {image_dir}")

        # Fall back to default location
        if image_dir is None:
            if LOCAL_IMAGES_DIR.exists() and LOCAL_IMAGES_DIR.is_dir():
                image_dir = LOCAL_IMAGES_DIR
                print(f"[ImageViewer] Using default location: {image_dir}")

    if image_dir is None or not image_dir.exists():
        print(f"[ImageViewer] No image directory found. Set LOCAL_IMAGE_ROOT env var.")
        return []

    # Find OME-TIFF files
    images = []
    for ext in ["*.ome.tif", "*.ome.tiff", "*.OME.TIFF", "*.OME.TIF"]:
        images.extend(image_dir.glob(ext))

    print(f"[ImageViewer] Found {len(images)} images in {image_dir}")
    return sorted([img.name for img in images])


def get_image_dir() -> Optional[Path]:
    """Get the image directory path."""
    env_root = os.environ.get("LOCAL_IMAGE_ROOT", "")
    if env_root:
        env_path = Path(env_root)
        if env_path.exists():
            return env_path
    if LOCAL_IMAGES_DIR.exists():
        return LOCAL_IMAGES_DIR
    return None


class ImageViewerPanel:
    """Panel-based image viewer using AVIVATOR."""

    def __init__(self):
        self.available_images = get_local_images()
        self._images_available = len(self.available_images) > 0
        self._image_dir = get_image_dir()

        # Create widgets
        self.image_selector = pn.widgets.Select(
            name="Select Image",
            options=["-- Select an image --"] + self.available_images,
            value="-- Select an image --",
            width=350,
            disabled=not self._images_available
        )

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
        if self._images_available:
            return f"Found {len(self.available_images)} images in {self._image_dir}"
        else:
            return "No images found. Set LOCAL_IMAGE_ROOT environment variable."

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
            <h2 style="margin-top: 20px; margin-bottom: 10px;">OME-TIFF Viewer</h2>
            <p style="opacity: 0.8; text-align: center; max-width: 400px; line-height: 1.6;">
                Select an OME-TIFF image from the dropdown above<br>
                to visualize multiplex imaging data.
            </p>
        </div>
        """

    def _get_viewer_html(self, image_url: str) -> str:
        """Get HTML to embed Avivator viewer with image."""
        from urllib.parse import quote
        encoded_url = quote(image_url, safe='')

        return f"""
        <iframe
            id="avivator-frame"
            src="{AVIVATOR_URL}?url={encoded_url}"
            style="
                width: 100%;
                height: 100%;
                min-height: 700px;
                border: none;
                border-radius: 8px;
                background: #000;
            "
            allow="fullscreen"
        ></iframe>
        """

    def _on_image_select(self, event):
        """Handle image selection."""
        if event.new and event.new != "-- Select an image --":
            # Construct FULL absolute URL for the image
            # Public AVIVATOR needs complete URL to fetch the image
            try:
                if pn.state.location:
                    protocol = pn.state.location.protocol or "http:"
                    host = pn.state.location.host or "localhost:8080"
                    base_url = f"{protocol}//{host}"
                else:
                    base_url = "http://localhost:8080"
            except:
                base_url = "http://localhost:8080"

            image_url = f"{base_url}/images/{event.new}"
            print(f"[ImageViewer] Selected: {event.new}")
            print(f"[ImageViewer] Full URL: {image_url}")
            self.viewer_pane.object = self._get_viewer_html(image_url)
        else:
            self.viewer_pane.object = self._get_placeholder_html()

    def _refresh_images(self, event):
        """Refresh available images list."""
        self.available_images = get_local_images()
        self._images_available = len(self.available_images) > 0
        self._image_dir = get_image_dir()
        self.image_selector.options = ["-- Select an image --"] + self.available_images
        self.image_selector.disabled = not self._images_available
        self.status_pane.object = self._get_status_text()

    def get_panel(self) -> pn.Column:
        """Return the complete viewer panel."""
        header = pn.pane.Markdown(
            "## OME-TIFF Viewer",
            styles={"margin": "0 0 10px 0"}
        )

        controls = pn.Row(
            self.image_selector,
            self.refresh_button,
            align="center"
        )

        return pn.Column(
            header,
            controls,
            self.status_pane,
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
