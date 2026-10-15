import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { readFile } from 'node:fs/promises';
import test from 'node:test';
import { removeXmlElementPrefixes, validatedImageUrl } from '../../app/shiny_app/www/avivator/url-policy.mjs';

const viewer = 'https://islet-explorer.rc.ufl.edu/islet-explorer/avivator/index.html';
const image = 'https://islet-explorer.rc.ufl.edu/islet-explorer/local_images/example.ome.tiff';
function viewerWithImage(value, base = viewer) {
  const url = new URL(base);
  url.searchParams.set('image_url', value);
  return url.href;
}

test('accepts only image routes beneath the current app prefix', () => {
  assert.equal(validatedImageUrl(viewerWithImage(image)), image);
  assert.equal(validatedImageUrl(viewerWithImage('/islet-explorer/local_images/example.ome.tiff')), image);
  assert.equal(validatedImageUrl(viewerWithImage('local_images/example.ome.tiff')), image);
  assert.equal(validatedImageUrl(viewerWithImage('images/a%20b.ome.tif')),
    'https://islet-explorer.rc.ufl.edu/islet-explorer/images/a%20b.ome.tif');
  assert.equal(validatedImageUrl(viewerWithImage('images/a.OME.TIF',
    'http://127.0.0.1:3838/avivator/')), 'http://127.0.0.1:3838/images/a.OME.TIF');
  assert.equal(validatedImageUrl(viewerWithImage('images/a.ome.tif',
    'http://localhost:3838/nested/app/avivator/index.html')),
  'http://localhost:3838/nested/app/images/a.ome.tif');
});

test('missing or duplicate image URLs cannot select a vendor demo', () => {
  assert.throws(() => validatedImageUrl(viewer));
  assert.throws(() => validatedImageUrl(viewerWithImage('')));
  assert.throws(() => validatedImageUrl(`${viewerWithImage(image)}&image_url=${encodeURIComponent(image)}`));
});

for (const value of [
  'https://outside.invalid/image.ome.tiff', '//outside.invalid/image.ome.tiff',
  'http://islet-explorer.rc.ufl.edu/islet-explorer/images/a.ome.tif',
  'https://islet-explorer.rc.ufl.edu:444/islet-explorer/images/a.ome.tif',
  'https://user:password@islet-explorer.rc.ufl.edu/islet-explorer/images/a.ome.tif',
  'data:image/tiff;base64,AAAA', 'blob:https://islet-explorer.rc.ufl.edu/example',
  'javascript:alert(1)', '/another-app/images/a.ome.tif',
  '/islet-explorer/local_images_extra/a.ome.tif', 'images/nested/a.ome.tif',
  'images/a.ome.tif?forward=https://outside.invalid', 'images/a.ome.tif#fragment',
  'images/a.xml', 'images/a.ome.tif/extra', ' images/a.ome.tif',
  'images/../local_images/a.ome.tif', 'images/%2e%2e/local_images/a.ome.tif',
  'images/%252e%252e%252flocal_images%252fa.ome.tif', 'images/a%2fb.ome.tif',
  'images/a%5cb.ome.tif', 'images/a\\b.ome.tif', 'images/a%00.ome.tif',
  'images/a\n.ome.tif', 'images/a%.ome.tif', 'images/%3fa.ome.tif',
]) {
  test(`rejects unsafe image URL ${JSON.stringify(value)}`, () => {
    assert.throws(() => validatedImageUrl(viewerWithImage(value)));
  });
}

test('preserves the existing XML element-prefix fix without rewriting other input', () => {
  assert.equal(removeXmlElementPrefixes('<ome:OME><ome:Image ID="Image:0"/></ome:OME>', 'application/xml'),
    '<OME><Image ID="Image:0"/></OME>');
  assert.equal(removeXmlElementPrefixes('<ome:OME/>', 'text/html'), '<ome:OME/>');
});

test('iframe applies its own connection boundary before loading the bootstrap', async () => {
  const html = await readFile(new URL('../../app/shiny_app/www/avivator/index.html', import.meta.url), 'utf8');
  assert.match(html, /http-equiv="Content-Security-Policy"/);
  assert.match(html, /connect-src 'self';/);
  assert.match(html, /script-src 'self' 'wasm-unsafe-eval';/);
  assert.doesNotMatch(html, /(?:^|\s)'unsafe-eval'(?:\s|;)/);
  assert.match(html, /src="\.\/bootstrap\.mjs"/);
  assert.doesNotMatch(html, /src="\.\/assets\/index-/);
  assert.ok(html.indexOf('Content-Security-Policy') < html.indexOf('bootstrap.mjs'));
});

test('the reviewed vendor bytes have not changed beneath the containment layer', async () => {
  const directory = new URL('../../app/shiny_app/www/avivator/', import.meta.url);
  const manifest = await readFile(new URL('VENDOR_SHA256SUMS', directory), 'utf8');
  for (const line of manifest.trim().split('\n')) {
    const [, expected, path] = line.match(/^([0-9a-f]{64})  (assets\/[A-Za-z0-9_.-]+)$/);
    const actual = createHash('sha256').update(await readFile(new URL(path, directory))).digest('hex');
    assert.equal(actual, expected, `${path}: update provenance and rerun browser compatibility review`);
  }
});
