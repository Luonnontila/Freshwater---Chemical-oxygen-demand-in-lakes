# Freshwaters - Chemical oxygen demand in lakes
This repository is for the source code of a Luonnontila web service indicator for [chemical oxygen demand in Finnish lakes](https://luonnontila.fi/indikaattorit-elinymparistoittain/sisavedet/kemiallinen-hapenkuluts/). The indicator is used to  for chemical oxygen demand in naturally humic and non-humic Finnish lakes. The indicator describes the general trend in the deposition of organic matter (e.g. humic substances) in Finnish lakes.

## Background

The chemical oxygen demand of inland waters indicates the amount of organic matter in the water bodies and is thus closely related to the cocentration of humic substances. The background concentration of humic substances is strongly influenced by natural factors such as the amount of peatlands in the drainage areas, so the development of the indicator has separately taken into account darker humic and clearer non-humic lakes.

The indicator is produced from an extensive series of water quality measurements of Finnish lakes, the first of which were taken already in 1960. 

## Code information

The code is written in **R** (version 4.2.1) and formatted to a **targets** package pipeline.

### Process description

In short, the process for the indicator production advances from **requesting data from the VESLA-database API** (, administrated by the Finnish Environment Institute) and proceeds through **data filtering and processing steps** to a **fitted generalized additive mixed effects model** taking into account the hierarchical data structure. This 

#### Process in a glimpse

#### Data sources

#### Data processing

#### Indicator calculation

#### Descriptive statistics
