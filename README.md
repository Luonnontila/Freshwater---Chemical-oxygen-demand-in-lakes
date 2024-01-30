# Freshwaters - Chemical oxygen demand in lakes
This repository is for the source code of a Luonnontila web service indicator for [chemical oxygen demand in Finnish lakes](https://luonnontila.fi/indikaattorit-elinymparistoittain/sisavedet/kemiallinen-hapenkuluts/). The indicator describes the general trend in the deposition of organic matter (e.g. humic substances) in Finnish lakes. The background concentration of humic substances is strongly influenced by natural factors such as the amount of peatlands in the drainage areas and thus the indicator differentiates between the trends of the darker humic and clearer non-humic lakes.

The indicator is produced from an extensive series of water quality measurements of Finnish lakes, the first of which were taken already in 1960.

## Structure and files in this repository

The project contains four main folders: **data**, **output**, **script** and **_targets**.  
**data**: Contains data used for successful completion of the pipeline (i.e. a csv-file with a list of waterbodies included in the calculation).  
**output**: Contains outputs prduced in the pipeline (e.g. analysis results).  
**script**: Contains the necessary code determining the functions used in different parts of pipeline.  
**_targets**: Contains objects produced and used by the **targets** package  

Full structure and files in the repository:

```
+-- data
|   +-- Otantalista_SV7_humuspitoisuus_jarvet_vesimuodostumat.csv
+-- output
|   +-- COD_GAM_fit_tidy_ENG.csv
|   +-- COD_GAM_fit_tidy_FIN.csv
|   +-- COD_GAM_fit_tidy_SWE.csv
|   +-- VESLA_COD_status.csv
|   \-- VESLA_COD_trends.csv
+-- README.md
+-- script
|   +-- Functions
|   |   +-- Vesla_functions_COD.r
+-- VESLA_COD.Rproj
+-- _targets
|   +-- meta
|   |   +-- meta
|   |   +-- process
|   |   +-- progress
|   |   +-- tar_temp_11609973cfe
|   |   +-- tar_temp_217c49fa78ab
|   |   +-- tar_temp_23c03dc43f5b
|   |   \-- tar_temp_e141fe8258b
|   +-- objects
|   |   +-- COD_fit_GAM_trend_output_tidy
|   |   +-- COD_fit_trend_output
|   |   +-- COD_fit_trend_output_tidy
|   |   +-- COD_status_output
|   |   +-- COD_trend_new
|   |   +-- COD_trend_old
|   |   +-- COD_trend_output_combined
|   |   +-- COD_trend_output_new
|   |   +-- COD_trend_output_old
|   |   +-- Comparison_COD
|   |   +-- GAM_fit_COD
|   |   +-- outputfile1
|   |   +-- outputfile2
|   |   +-- VESLA_COD_data_combined
|   |   +-- VESLA_COD_minimum_timeframe
|   |   +-- VESLA_query_data
|   |   +-- VESLA_raw_query_data
|   |   +-- VESLA_results
|   |   +-- VESLA_samples
|   |   +-- VESLA_sampling_events
|   |   +-- VESLA_sites_all_results
|   |   +-- VESLA_site_ids
|   |   +-- VESLA_summary_COD_data_Lake
|   |   +-- VESLA_summary_COD_data_Site
|   |   +-- VESLA_summary_COD_data_Site_comparison
|   |   +-- VESLA_summary_COD_data_Site_new
|   |   +-- VESLA_summary_COD_data_Site_old
|   |   \-- VESLA_summary_stats_Lake
|   \-- user
\-- _targets.R
```
## General information on the code

The code is written in **R** (version 4.2.1) and takes advantage of R packages **targets**, **brms**, **readr**, **dplyr**, **tidyr**, **data.table** and **zoo**.

The code for the indicator is primarily in two files: targets.R and Vesla_functions_COD.r. Below find a short description of contents in each of the files.

### _targets.R 

_targets.R -file determines the structure of the pipeline in the format required by the **targets** package. For more information on the logic of targets, see [targets manual](https://books.ropensci.org/targets/).

### Vesla_functions_COD.r
Vesla_functions_COD.r determines the functions used in different steps of the pipeline. The functions can be categorized in four different stages of the pipeline: (1) obtaining (functions starting with `get_VESLA_`), (2) shaping (functions starting with `shape_VESLA_data_`), (3) summarizing (funcitons starting with `summarize_VESLA_data_`), and (4) analyzing (functions starting with `analyze_VESLA_data`) data obtained from the [VESLA 2.0 database](https://rajapinnat.ymparisto.fi/api/vesla/2.0/).

### Process description

In short, the process for the indicator production advances from **requesting data from the VESLA-database API** (, administrated by the Finnish Environment Institute) and proceeds through **data filtering and processing steps** to a **fitted generalized additive mixed effects model** taking into account the hierarchical data structure.

![](test_process.html)

#### Process in a glimpse

#### Data sources

#### Data processing

#### Indicator calculation

#### Descriptive statistics
