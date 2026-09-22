import { removeXmlElementPrefixes, validatedImageUrl } from './url-policy.mjs';

function showUnavailable() {
  const root = document.getElementById('root');
  root.setAttribute('role', 'status');
  root.textContent = 'Select an available image in Islet Explorer to open the viewer.';
}

function restrictUrlForm() {
  // The vendor panels reference these IDs, but its tab buttons omit them.
  for (const [index, tab] of [...document.querySelectorAll('[role="tablist"] [role="tab"]')].entries()) {
    const panel = document.getElementById(`simple-tabpanel-${index}`);
    if (!panel) continue;
    tab.id = panel.getAttribute('aria-labelledby');
    tab.setAttribute('aria-controls', panel.id);
  }
  // These hidden native inputs are activated by the vendor's visible picker
  // buttons. Remove them from the accessibility tree so they do not create
  // unnamed controls nested inside another interactive element.
  for (const fileInput of document.querySelectorAll('input[type="file"]')) {
    fileInput.hidden = true;
    fileInput.style.display = 'none';
    fileInput.setAttribute('aria-label', 'Open local OME-TIFF images');
  }
  for (const progress of document.querySelectorAll('[role="progressbar"]')) {
    progress.setAttribute('aria-label', 'Loading tissue image');
  }
  for (const canvas of document.querySelectorAll('canvas')) {
    canvas.setAttribute('role', 'img');
    canvas.setAttribute('aria-label', 'Pancreatic tissue image with selectable marker channels');
  }
  for (const selector of document.querySelectorAll('[role="combobox"]')) {
    let channel = selector.parentElement;
    while (channel && !channel.querySelector('input[type="checkbox"]')) channel = channel.parentElement;
    const isLens = Boolean(channel?.textContent.trim().startsWith('Lens:'));
    selector.setAttribute('aria-label', isLens ? 'Lens marker channel' : 'Marker channel');
    const checkbox = channel?.querySelector('input[type="checkbox"]');
    if (checkbox) checkbox.setAttribute('aria-label', isLens ? 'Enable image lens' :
      `Show ${selector.textContent.trim()} marker channel`);
  }
  const input = document.getElementById('ome-input');
  if (!input) return;
  input.disabled = true;
  input.readOnly = true;
  input.setAttribute('aria-label', 'Image selected in Islet Explorer');
  input.title = 'Select an image using the Islet Explorer image selector.';
  const form = input.closest('form');
  if (form) {
    for (const control of form.querySelectorAll('button, input')) control.disabled = true;
  }
}

async function startViewer() {
  let image;
  try {
    image = validatedImageUrl(window.location.href);
  } catch {
    showUnavailable();
    return;
  }

  // Preserve channel_config and normalize the accepted image URL for the vendor.
  const current = new URL(window.location.href);
  current.searchParams.set('image_url', image);
  window.history.replaceState(null, '', current);

  const originalParse = DOMParser.prototype.parseFromString;
  DOMParser.prototype.parseFromString = function (text, type) {
    return originalParse.call(this, removeXmlElementPrefixes(text, type), type);
  };

  // React mounts the URL form asynchronously and can recreate it when the menu
  // changes. Capture prevents submission even before the observer disables it.
  document.addEventListener('submit', event => {
    if (event.target.querySelector('#ome-input')) {
      event.preventDefault();
      event.stopImmediatePropagation();
    }
  }, true);
  const observer = new MutationObserver(restrictUrlForm);
  observer.observe(document.getElementById('root'), { childList: true, subtree: true });

  try {
    await import('./assets/index-DBmmXfNW.js');
    restrictUrlForm();
  } catch {
    observer.disconnect();
    showUnavailable();
  }
}

startViewer();
