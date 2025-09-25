AIOLOS into PumpKin Tool Documentation
=====================================

Welcome to the AIOLOS into PumpKin Tool documentation. This tool provides a complete pipeline for converting AIOLOS simulation data into PumpKin analysis format and generating comprehensive chemical pathway visualizations.

Overview
--------

The AIOLOS into PumpKin Tool is designed to:

* Process AIOLOS reaction files and simulation data
* Convert data formats for PumpKin pathway analysis
* Generate comprehensive plots and visualizations
* Provide automated pipeline management

Quick Start
-----------

Basic usage::

    python main.py --run-all

Custom configuration::

    python main.py --run-all --num-cells 100 --data-dir ./data/ --output-dir ./results/

Modules
-------

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   conversion_tool
   outputs
   pumpkin
   main

API Reference
=============

Main Pipeline
-------------

.. automodule:: main
   :members:
   :undoc-members:
   :show-inheritance:

Conversion Tool
---------------

.. automodule:: Conversion_Tool.main
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: Conversion_Tool.processing_aiolos_reac_file
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: Conversion_Tool.making_densities_file
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: Conversion_Tool.making_rates_file
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: Conversion_Tool.calculating_rates
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: Conversion_Tool.getting_reactions_from_log
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: Conversion_Tool.reading_rates_from_new_chem_files
   :members:
   :undoc-members:
   :show-inheritance:

Output Generation
-----------------

.. automodule:: Outputs.plotting_main
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: Outputs.species_specific_plots
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: Outputs.pathway_rates_plots
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: Outputs.rates_of_reactions_plots
   :members:
   :undoc-members:
   :show-inheritance:

PumpKin Integration
-------------------

.. automodule:: PumpKin.running_pumpkin
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: PumpKin.running_all_pumpkin
   :members:
   :undoc-members:
   :show-inheritance:

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`