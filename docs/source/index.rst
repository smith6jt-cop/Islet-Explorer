Islet-Explorer Documentation
============================

.. image:: https://img.shields.io/badge/R-Shiny-blue
   :alt: R Shiny
   
.. image:: https://img.shields.io/badge/Python-3.6+-green
   :alt: Python

**Islet-Explorer** is an interactive R Shiny application for exploring human pancreatic islet measurements from CODEX multiplexed imaging data. The app provides comprehensive analysis tools for comparing donor groups, islet sizes, and cellular composition with an embedded image viewer for OME-TIFF files.

.. note::
   This documentation covers version 1.0 of Islet-Explorer, including the new efficient GeoJSON overlay system for cell segmentation visualization.

Features
--------

* **Interactive Data Exploration**: Plot and analyze islet measurements by donor group and size
* **Statistical Analysis**: Automated statistical testing with multiple comparison corrections
* **Pseudotime Trajectory Analysis**: Explore cellular differentiation trajectories
* **Integrated Image Viewer**: View OME-TIFF images with Avivator/Viv
* **Cell Segmentation Overlays**: Efficient visualization of 200K+ cell boundaries per image
* **Quality Assurance Tools**: Built-in data validation and diagnostics

Quick Start
-----------

1. **Install R Dependencies**::

     Rscript scripts/install_shiny_deps.R

2. **Preprocess GeoJSON Overlays** (one-time)::

     ./scripts/preprocess_all_geojson.sh

3. **Run the Application**::

     R -q -e "shiny::runApp('app/shiny_app', launch.browser=FALSE)"

4. **Access the App**: Navigate to http://localhost:XXXX in your browser

Table of Contents
-----------------

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   user_guide/installation
   user_guide/getting_started
   user_guide/plot_tab
   user_guide/statistics_tab
   user_guide/trajectory_tab
   user_guide/viewer_tab

.. toctree::
   :maxdepth: 2
   :caption: Technical Documentation

   technical/architecture
   technical/data_format
   technical/geojson_overlay
   technical/preprocessing
   technical/api_reference

.. toctree::
   :maxdepth: 2
   :caption: Development

   development/contributing
   development/testing
   development/deployment

.. toctree::
   :maxdepth: 1
   :caption: Reference

   reference/glossary
   reference/troubleshooting
   reference/faq
   reference/changelog

Citation
--------

If you use Islet-Explorer in your research, please cite:

.. code-block:: bibtex

   @software{islet_explorer_2025,
     title = {Islet-Explorer: Interactive Analysis of Pancreatic Islet CODEX Data},
     author = {Smith Lab},
     year = {2025},
     url = {https://github.com/smith6jt-cop/Islet-Explorer}
   }

Support
-------

* **Issues**: `GitHub Issues <https://github.com/smith6jt-cop/Islet-Explorer/issues>`_
* **Discussions**: `GitHub Discussions <https://github.com/smith6jt-cop/Islet-Explorer/discussions>`_
* **Email**: smith6jt@example.edu

License
-------

This project is licensed under the MIT License - see the LICENSE file for details.

Indices and Tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

