/* One browser-local hover state per group of linked plots (the Trajectory tab's three, the Plot tab's two).
   Only clicks reach R. */
(function () {
  'use strict';
  if (window.IsletTrajectoryLinks) return;
  const groups = new Map();
  const role = trace => trace && trace.meta && trace.meta.trajectory_role;

  function groupFor(config) {
    if (groups.has(config.group)) return groups.get(config.group);
    const group = {graphs: new Map(), pinned: null, hover: null, active: true,
      frame: null, sequence: 0};
    groups.set(config.group, group);
    // The R Plotly binding clears hover/click inputs even when its event list
    // is empty. Cancel those unused inputs only for these linked widgets.
    window.jQuery(document).on('shiny:inputchanged.' + config.group, function (event) {
      if (!event.name.startsWith('plotly_')) return;
      for (const id of group.graphs.keys()) {
        if (event.name.endsWith('-' + id)) { event.preventDefault(); break; }
      }
    });
    Shiny.addCustomMessageHandler(config.message, function (state) {
      group.pinned = state.key || null;
      group.active = state.active !== false;
      // Selection/clear messages also discard an obsolete hover preview.
      group.hover = null;
      schedule(group);
    });
    return group;
  }

  function schedule(group) {
    if (group.frame !== null) return;
    group.frame = requestAnimationFrame(function () {
      group.frame = null;
      const key = group.active ? (group.hover && group.hover.key) || group.pinned : null;
      group.graphs.forEach(function (record, id) {
        if (!record.gd.isConnected) {
          record.dispose();
          group.graphs.delete(id);
        } else {
          record.requested = key;
          paint(record);
        }
      });
    });
  }

  function indexPoints(record) {
    record.points = new Map();
    (record.gd.data || []).forEach(function (trace) {
      if (role(trace) !== 'islet' || trace.visible === 'legendonly' || trace.visible === false) return;
      const keys = trace.key || [];
      for (let i = 0; i < keys.length; i++) {
        const key = keys[i], x = trace.x[i], y = trace.y[i];
        if (typeof key !== 'string' || !key || !Number.isFinite(x) || !Number.isFinite(y)) continue;
        const sizes = trace.marker && trace.marker.size;
        const size = Array.isArray(sizes) || ArrayBuffer.isView(sizes) ? sizes[i] : sizes;
        record.points.set(key, {x, y, size: Math.max(15, (Number(size) || 4) + 7)});
      }
    });
  }

  function paint(record) {
    if (record.disposed) return;
    const point = record.points.get(record.requested);
    const layout = record.gd._fullLayout, xa = layout && layout.xaxis, ya = layout && layout.yaxis;
    if (!point || !xa || !ya || record.gd.classList.contains('recalculating')) {
      record.overlay.style.display = 'none';
      record.overlay.removeAttribute('data-linked-key');
      return;
    }
    const cx = xa._offset + xa.l2p(point.x), cy = ya._offset + ya.l2p(point.y);
    record.clip.setAttribute('x', xa._offset); record.clip.setAttribute('y', ya._offset);
    record.clip.setAttribute('width', xa._length); record.clip.setAttribute('height', ya._length);
    record.circles.forEach(function (circle) {
      circle.setAttribute('cx', cx); circle.setAttribute('cy', cy);
      circle.setAttribute('r', point.size / 2);
    });
    record.overlay.setAttribute('data-linked-key', record.requested);
    record.overlay.setAttribute('data-point-x', point.x);
    record.overlay.setAttribute('data-point-y', point.y);
    record.overlay.style.display = '';
  }

  function makeOverlay(record, id) {
    const svgNS = 'http://www.w3.org/2000/svg';
    const make = name => document.createElementNS(svgNS, name);
    const svg = make('svg'), defs = make('defs'), clipPath = make('clipPath'), clip = make('rect'), rings = make('g');
    svg.classList.add('trajectory-link-highlight');
    svg.setAttribute('aria-hidden', 'true');
    Object.assign(svg.style, {position: 'absolute', inset: '0', width: '100%', height: '100%', pointerEvents: 'none', zIndex: '3'});
    clipPath.id = id + '-highlight-clip'; clipPath.appendChild(clip); defs.appendChild(clipPath);
    rings.setAttribute('clip-path', 'url(#' + clipPath.id + ')');
    const circles = ['white', '#111827'].map(function (color, i) {
      const circle = make('circle');
      circle.setAttribute('fill', 'none'); circle.setAttribute('stroke', color);
      circle.setAttribute('stroke-width', i === 0 ? '5' : '2');
      rings.appendChild(circle); return circle;
    });
    svg.append(defs, rings);
    // Plotly coordinates start at its inner SVG canvas, not at the widget's
    // outer box (which can include padding for the modebar).
    record.gd.querySelector('.svg-container').appendChild(svg);
    record.overlay = svg; record.clip = clip; record.circles = circles;
  }

  function pointKey(record, point) {
    if (!point || role(point.data) !== 'islet') return null;
    const key = point.data.key && point.data.key[point.pointNumber];
    return typeof key === 'string' && record.points.has(key) ? key : null;
  }

  function attach(el, config) {
    const gd = typeof el.on === 'function' ? el : el.querySelector('.js-plotly-plot');
    if (!gd) return;
    const group = groupFor(config);
    const previous = group.graphs.get(el.id);
    if (previous) previous.dispose();
    if (group.hover && group.hover.source === el.id) group.hover = null;
    const record = {gd, disposed: false, requested: null};
    group.graphs.set(el.id, record);
    indexPoints(record);
    makeOverlay(record, el.id);
    let down = null;

    function hoveredKey() {
      const points = gd._hoverdata || [];
      return points.map(point => pointKey(record, point)).find(Boolean) || null;
    }
    function commit(key) {
      if (!group.active || !record.points.has(key) || record.disposed ||
          gd.classList.contains('recalculating')) return;
      Shiny.setInputValue(config.input, {key, source: el.id, nonce: ++group.sequence}, {priority: 'event'});
    }
    function hover(event) {
      const key = (event.points || []).map(point => pointKey(record, point)).find(Boolean);
      if (!key || !group.active) return;
      group.hover = {key, source: el.id};
      schedule(group);
    }
    function unhover() {
      if (group.hover && group.hover.source === el.id) {
        group.hover = null;
        schedule(group);
      }
    }
    function pointerDown(event) {
      down = event.button === 0 ? {x: event.clientX, y: event.clientY,
        key: hoveredKey(), sent: false, moved: false} : null;
    }
    function pointerMove(event) {
      if (down && (Math.abs(event.clientX - down.x) > 5 || Math.abs(event.clientY - down.y) > 5)) down.moved = true;
    }
    function pointerUp(event) {
      pointerMove(event);
      if (down && down.key && !down.sent && !down.moved) {
        down.sent = true;
        commit(down.key);
      }
    }
    function click(event) {
      const key = (event.points || []).map(point => pointKey(record, point)).find(Boolean);
      // Plotly and the native fallback can report the same physical gesture.
      if (!key || (down && (down.sent || down.moved))) return;
      if (down) down.sent = true;
      commit(key);
    }
    function cancel() { down = null; unhover(); }
    function redraw() { indexPoints(record); schedule(group); }
    function relayout() { schedule(group); }
    const resize = new ResizeObserver(relayout);
    resize.observe(gd);
    gd.on('plotly_hover', hover);
    gd.on('plotly_unhover', unhover);
    gd.on('plotly_click', click);
    gd.on('plotly_afterplot', redraw);
    gd.on('plotly_relayout', relayout);
    gd.addEventListener('pointerdown', pointerDown, true);
    gd.addEventListener('pointermove', pointerMove, true);
    gd.addEventListener('pointerup', pointerUp, true);
    gd.addEventListener('pointercancel', cancel);
    gd.addEventListener('mouseleave', cancel);
    record.dispose = function () {
      record.disposed = true;
      gd.removeListener('plotly_hover', hover);
      gd.removeListener('plotly_unhover', unhover);
      gd.removeListener('plotly_click', click);
      gd.removeListener('plotly_afterplot', redraw);
      gd.removeListener('plotly_relayout', relayout);
      resize.disconnect();
      record.overlay.remove();
      gd.removeEventListener('pointerdown', pointerDown, true);
      gd.removeEventListener('pointermove', pointerMove, true);
      gd.removeEventListener('pointerup', pointerUp, true);
      gd.removeEventListener('pointercancel', cancel);
      gd.removeEventListener('mouseleave', cancel);
    };
    schedule(group);
    // Recover state even if the server's first message preceded widget binding.
    Shiny.setInputValue(config.ready, ++group.sequence, {priority: 'event'});
  }

  window.IsletTrajectoryLinks = {attach};
}());
