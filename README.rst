.. |date| date::

************************
PPTC PDX Fusion Pipeline
************************

:authors: Komal S Rathi
:contact: Komal Rathi (rathik@email.chop.edu)
:organization: CHOP
:status: Completed
:date: |date|

.. meta::
   :keywords: pdx, mouse, fusion, 2019
   :description: pdx mouse fusion analysis pipeline.

Introduction
============

This repo contains code for:

1. Fusion detection and filtering
2. Create publication quality figures.

Details
=======

update_fusion_output.R: Updates obsolete Hugo gene symbols in raw input data
driver_fusions.R: Driver fusion results from cytogenetics/literature search
venn_diagram.R: Creates venn diagram of driver fusions
summary: Creates combined excel sheet with Fusion detection/filtering pipeline + driver fusions

