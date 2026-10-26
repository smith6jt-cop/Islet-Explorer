/* Islet Explorer color theme switch and white-background PNG export (2026-10-18). See islet-theme.css. */
(function () {
  'use strict';

  function current() {
    return document.documentElement.getAttribute('data-theme') === 'light' ? 'light' : 'dark';
  }

  function label(button) {
    if (!button) return;
    var dark = current() === 'dark';
    button.textContent = dark ? '☀ Light theme' : '☾ Dark theme';
    button.setAttribute('aria-label', dark ? 'Switch to the light theme' : 'Switch to the dark theme');
  }

  // Set the theme, remember it, and tell the server (input$app_theme) so the plots redraw in it.
  window.isletSetTheme = function (theme) {
    theme = theme === 'light' ? 'light' : 'dark';
    document.documentElement.setAttribute('data-theme', theme);
    try { window.localStorage.setItem('islet-theme', theme); } catch (e) { /* private window: not remembered */ }
    var input = document.getElementById('app_theme');
    if (input && input.value !== theme) {
      input.value = theme;
      if (window.jQuery) window.jQuery(input).trigger('change');
    }
    label(document.getElementById('theme_toggle'));
  };

  document.addEventListener('DOMContentLoaded', function () {
    var button = document.getElementById('theme_toggle');
    label(button);
    if (button) {
      button.addEventListener('click', function () {
        window.isletSetTheme(current() === 'dark' ? 'light' : 'dark');
      });
    }
  });

  // ---- PNG export on a white background ----------------------------------------------------------------------------
  // The modebar camera exports what is on screen; on the dark theme that would be a dark figure. Downloads are for
  // publication, so the replacement button (plotly_islet_config) re-draws the figure off screen in light colors.
  var INK = '#222222';
  var DARK_PAPER = '#1f1f1f';  // ISLET_THEME_PALETTE$dark$paper

  function recolourAxis(axis) {
    axis.gridcolor = '#e6e6e6';
    axis.zerolinecolor = '#bdbdbd';
    axis.linecolor = '#444444';
    axis.tickcolor = '#444444';
    axis.tickfont = Object.assign({}, axis.tickfont, { color: INK });
    if (axis.title && typeof axis.title === 'object') {
      axis.title = Object.assign({}, axis.title, { font: Object.assign({}, axis.title.font, { color: INK }) });
    }
  }

  function recolourColorbar(cb) {
    if (!cb) return cb;
    cb = Object.assign({}, cb, { tickfont: Object.assign({}, cb.tickfont, { color: INK }) });
    if (cb.title && typeof cb.title === 'object') {
      cb.title = Object.assign({}, cb.title, { font: Object.assign({}, cb.title.font, { color: INK }) });
    }
    return cb;
  }

  function stripPrivate(obj) {
    if (Array.isArray(obj)) return obj.map(stripPrivate);
    if (obj && typeof obj === 'object') {
      var out = {};
      Object.keys(obj).forEach(function (k) { if (k.charAt(0) !== '_') out[k] = stripPrivate(obj[k]); });
      return out;
    }
    return obj;
  }

  window.isletDownloadLightPng = function (gd, options) {
    var layout = stripPrivate(JSON.parse(JSON.stringify(gd.layout || {})));
    layout.paper_bgcolor = '#ffffff';
    layout.plot_bgcolor = '#ffffff';
    layout.font = Object.assign({}, layout.font, { color: INK });
    if (layout.title && typeof layout.title === 'object') {
      layout.title.font = Object.assign({}, layout.title.font, { color: INK });
    }
    Object.keys(layout).forEach(function (k) { if (/^[xy]axis\d*$/.test(k)) recolourAxis(layout[k]); });
    if (layout.legend) {
      layout.legend.bgcolor = 'rgba(255,255,255,0)';
      layout.legend.font = Object.assign({}, layout.legend.font, { color: INK });
      if (layout.legend.title && typeof layout.legend.title === 'object') {
        layout.legend.title.font = Object.assign({}, layout.legend.title.font, { color: INK });
      }
    }
    (layout.annotations || []).forEach(function (a) { a.font = Object.assign({}, a.font, { color: INK }); });
    var data = (gd.data || []).map(function (trace) {
      var t = Object.assign({}, trace);
      if (t.colorbar) t.colorbar = recolourColorbar(t.colorbar);
      if (t.marker && t.marker.colorbar) t.marker = Object.assign({}, t.marker, { colorbar: recolourColorbar(t.marker.colorbar) });
      // UMAP points are outlined in the screen background (islet_umap_plotly()); outline them in white instead.
      if (t.marker && t.marker.line && t.marker.line.color === DARK_PAPER) {
        t.marker = Object.assign({}, t.marker, { line: Object.assign({}, t.marker.line, { color: '#ffffff' }) });
      }
      return t;
    });
    return window.Plotly.downloadImage({ data: data, layout: layout }, options);
  };
})();
