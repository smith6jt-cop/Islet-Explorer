/* Card help popovers (card_help() in R/00_globals.R). One open at a time; an outside click or Escape closes it.
   While open, the panel is moved to <body> and placed in window coordinates: a "?" near the left edge or the
   bottom of the window must not push its text off screen, and a card's hover transform (.card:hover) would
   otherwise become the containing block of a fixed panel. Closing puts it back inside its <details>. */
(function () {
  'use strict';
  if (window.IsletCardHelp) return;
  window.IsletCardHelp = true;
  const MARGIN = 16;
  const shown = new Set();

  function place(details) {
    const panel = details._helpPanel, icon = details.querySelector(':scope > summary');
    if (!panel || !icon) return;
    const box = icon.getBoundingClientRect();
    if (!details.isConnected || !box.width) { details.open = false; return; }
    const width = Math.min(420, window.innerWidth - 2 * MARGIN);
    const left = Math.min(Math.max(MARGIN, box.right - width), window.innerWidth - width - MARGIN);
    const below = window.innerHeight - box.bottom - MARGIN;
    const above = box.top - MARGIN;
    Object.assign(panel.style, {position: 'fixed', width: width + 'px', left: left + 'px', right: 'auto'});
    if (below >= 240 || below >= above) {
      Object.assign(panel.style, {top: (box.bottom + 4) + 'px', bottom: 'auto',
        maxHeight: Math.max(160, below - 4) + 'px'});
    } else {
      Object.assign(panel.style, {top: 'auto', bottom: (window.innerHeight - box.top + 4) + 'px',
        maxHeight: Math.max(160, above - 4) + 'px'});
    }
  }

  function show(details) {
    const panel = details._helpPanel || details.querySelector(':scope > .card-help-panel');
    if (!panel) return;
    details._helpPanel = panel;
    panel._helpOwner = details;
    document.body.appendChild(panel);
    shown.add(details);
    place(details);
  }

  function hide(details) {
    const panel = details._helpPanel;
    shown.delete(details);
    if (!panel) return;
    panel.removeAttribute('style');
    // A card re-rendered by Shiny while its help was open leaves no <details> to return to.
    if (details.isConnected) details.appendChild(panel); else panel.remove();
  }

  const closeAll = except => [...shown].forEach(d => { if (d !== except) { d.open = false; hide(d); } });

  // toggle does not bubble; listen in the capture phase.
  document.addEventListener('toggle', function (event) {
    const details = event.target;
    if (!(details instanceof HTMLElement) || !details.matches('details.card-help')) return;
    if (details.open) { closeAll(details); show(details); } else hide(details);
  }, true);
  document.addEventListener('pointerdown', function (event) {
    const target = event.target;
    const panel = target.closest && target.closest('.card-help-panel');
    closeAll((panel && panel._helpOwner) || (target.closest && target.closest('details.card-help')));
  });
  document.addEventListener('keydown', function (event) {
    if (event.key === 'Escape') closeAll(null);
  });
  const replace = () => shown.forEach(place);
  window.addEventListener('resize', replace);
  window.addEventListener('scroll', replace, true);
})();
