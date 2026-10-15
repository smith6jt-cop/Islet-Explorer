// Validate before loading the vendor bundle: its missing-URL behavior loads a
// remote demonstration image. Keep this module independent of browser globals.
export function validatedImageUrl(viewerHref) {
  const viewer = new URL(viewerHref);
  if (!['http:', 'https:'].includes(viewer.protocol) ||
      !/\/avivator\/(?:index\.html)?$/.test(viewer.pathname)) {
    throw new Error('viewer_location_invalid');
  }

  const values = viewer.searchParams.getAll('image_url');
  if (values.length !== 1 || !values[0] || values[0] !== values[0].trim()) {
    throw new Error('viewer_image_required');
  }
  const candidate = values[0];
  // Inspect before URL() normalizes dot segments. Reject encoded separators and
  // residual percent escapes so another decoding layer cannot reveal traversal.
  for (const part of candidate.split('/')) {
    const decoded = decodeURIComponent(part);
    if (decoded === '.' || decoded === '..' || /[\\/%?#\u0000-\u001f\u007f]/.test(decoded)) {
      throw new Error('viewer_image_path_invalid');
    }
  }

  const app = new URL('../', viewer);
  const image = new URL(candidate, app);
  if (image.origin !== viewer.origin || image.username || image.password ||
      image.search || image.hash || !image.pathname.startsWith(app.pathname)) {
    throw new Error('viewer_image_origin_invalid');
  }
  const relative = image.pathname.slice(app.pathname.length);
  if (!/^(?:local_images|images)\/[^/]+$/.test(relative) || !/\.ome\.tiff?$/i.test(relative)) {
    throw new Error('viewer_image_path_invalid');
  }
  return image.href;
}

export function removeXmlElementPrefixes(text, type) {
  if (type !== 'application/xml' || typeof text !== 'string') return text;
  return text.replace(/(<\/?)([A-Za-z_][\w.-]*):([A-Za-z_][\w.-]*)/g,
    (_match, start, _prefix, local) => start + local);
}
