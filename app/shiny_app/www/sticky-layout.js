/* The page header and tab bar are frozen at the top of the scrolling main panel, and each tab's control panels
   (.tab-controls) below them (owner, 2026-10-21). CSS needs their heights to stack them: this keeps
   --islet-header-h and --islet-tabs-h on <html> equal to the rendered heights, which change with the theme's logo,
   the window width and the browser zoom. */
(function () {
  'use strict';
  function measure() {
    const root = document.documentElement;
    const header = document.querySelector('.app-header');
    const tabs = document.getElementById('tabs');
    if (header) root.style.setProperty('--islet-header-h', header.getBoundingClientRect().height + 'px');
    if (tabs) root.style.setProperty('--islet-tabs-h', tabs.getBoundingClientRect().height + 'px');
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
