/* The page header (logo, tab bar, theme switch in one row since 2026-10-26) is frozen at the top of the scrolling
   main panel, and each tab's control panels (.tab-controls) below it (owner, 2026-10-21). CSS needs the heights to
   stack them: this keeps --islet-header-h and --islet-tabs-h on <html> equal to the rendered heights, which change
   with the theme's logo, the window width and the browser zoom. A tab bar inside the header adds nothing below it,
   so --islet-tabs-h is 0 then; the variable stays for a layout that puts the bar under the header. */
(function () {
  'use strict';
  function measure() {
    const root = document.documentElement;
    const header = document.querySelector('.app-header');
    const tabs = document.getElementById('tabs');
    if (header) root.style.setProperty('--islet-header-h', header.getBoundingClientRect().height + 'px');
    if (tabs) root.style.setProperty('--islet-tabs-h',
      (header && header.contains(tabs) ? 0 : tabs.getBoundingClientRect().height) + 'px');
  }
  function start() {
    measure();
    if (typeof ResizeObserver !== 'function') return;
    const observer = new ResizeObserver(measure);
    ['.app-header', '#tabs'].forEach(function (selector) {
      const el = document.querySelector(selector);
      if (el) observer.observe(el);
    });
  }
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', start);
  else start();
}());
