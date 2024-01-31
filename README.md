# Freshwaters - Chemical oxygen demand in lakes
This repository is for the source code of a Luonnontila web service indicator for [chemical oxygen demand in Finnish lakes (COD)](https://luonnontila.fi/indikaattorit-elinymparistoittain/sisavedet/kemiallinen-hapenkuluts/). The indicator describes the general trend in the deposition of organic matter (e.g. humic substances) in Finnish lakes. The background concentration of humic substances is strongly influenced by natural factors such as the amount of peatlands in the drainage areas and thus the indicator differentiates between the trends of the darker humic and clearer non-humic lakes.

The indicator is produced from an extensive series of water quality measurements of Finnish lakes, the first of which were taken already in 1960.

## Structure and files in this repository

The project contains the following folders and files:

- **data**: Contains data used for successful completion of the pipeline (i.e. a csv-file with a list of waterbodies included in the calculation).
  - **Otantalista_SV7_humuspitoisuus_jarvet_vesimuodostumat.csv**: Determines the list of waterbodies included in the analysis.
- **output**: Contains outputs produced by running the pipeline.
  - **COD_GAM_fit_tidy_ENG.csv**, **COD_GAM_fit_tidy_FIN.csv**, and **COD_GAM_fit_tidy_SWE.csv**: Indicator based on the fitted statistical model (plotted in the indicator webpage), with formatting matching to language abbreviation in the end of the file (English, Finnish and Swedish language).
  - **VESLA_COD_status.csv**: Describes the status of the indicator in relation to the reference value.
  - **VESLA_COD_trends.csv**: Describes linear trends in the data in two different time periods (1960-1999 and (2001-2023).
- **script**: Contains the necessary code determining the functions used in different parts of pipeline.
  - **Functions**: Contains the function scripts used in the pipeline
    - **Vesla_functions_COD.r**: An R-code containing the script of all functions needed to execute the analysis pipeline.
- **run_targets_pipeline.r**: An R-code file containing commands to run the pipeline and produce outputs descriptive of the pipeline structure.

Hierarchical folder structure and files in the repository:
```
+-- data
|   \-- Otantalista_SV7_humuspitoisuus_jarvet_vesimuodostumat.csv
+-- output
|   +-- COD_GAM_fit_tidy_ENG.csv
|   +-- COD_GAM_fit_tidy_FIN.csv
|   +-- COD_GAM_fit_tidy_SWE.csv
|   +-- VESLA_COD_status.csv
|   \-- VESLA_COD_trends.csv
+-- README.md
+-- run_targets_pipeline.r
+-- script
|   \-- Functions
|   |   \-- Vesla_functions_COD.r
\-- _targets.R
```
## General information on the code

The code is written in **R** (version 4.2.1) and takes advantage of R packages **targets**, **brms**, **readr**, **dplyr**, **tidyr**, **data.table**, **zoo** and **htmltools**.

The code for the indicator is primarily in two files: targets.R and Vesla_functions_COD.r. The former determines the analysis pipeline and the latter the functions used in it.

## Analysis pipeline 

In short, the process for the indicator production advances from **requesting data** from the Finnish Environment Institute administrated [VESLA 2.0 database API](https://rajapinnat.ymparisto.fi/api/vesla/2.0/) and proceeds through **data filtering and processing steps** to a **fitted generalized additive mixed effects model** taking into account the hierarchical data structure. Additional analyses to (1) determine the period of lowest COD across lakes, (2) current COD relative to the period of lowest COD, and (3) linear trend of COD from the year 2001 to present day are also produced in the pipeline. Outputs are then produced from these statistical models.

The structure of the pipeline is in the format required by the **targets** package and it is configured in the `_targets.R` script file. For more information on the logic and formatting of targets, see [targets manual](https://books.ropensci.org/targets/). At the end of the pipeline configuration 

**An interactive directed acyclic graph (DAG) flowchart description of the pipeline can be found [here](https://raw.githack.com/Luonnontila/Freshwater---Chemical-oxygen-demand-in-lakes/main/output/pipeline_network/objects_pipeline.html).** Circles depict different objects produced in the pipeline. Note that the graph can be zoomed (pinched), different objects moved (drag) and different parts highlighted (click on object) to enable easier reading.

### Vesla_functions_COD.r
Vesla_functions_COD.r determines the functions used in different steps of the pipeline. The functions can be categorized in five different stages of the pipeline:  
1. obtaining data (functions starting with `get_VESLA`)
2. restructuring and combining data (functions starting with `shape_VESLA_data`),
3. summarizing data (funcitons starting with `summarize_VESLA_data`), and
4. analyzing data (functions starting with `analyze_VESLA_data`) data obtained from the [VESLA 2.0 database](https://rajapinnat.ymparisto.fi/api/vesla/2.0/)
5. producing outputs from data and analyses (functions starting with `output_table_VESLA`)

**An interactive directed acyclic graph (DAG) flowchart description of how different functions contribute to the process can be found [here](https://raw.githack.com/Luonnontila/Freshwater---Chemical-oxygen-demand-in-lakes/main/output/pipeline_network/network_pipeline.html).** Triangles depict different functions and circles objects produced in the pipeline. Note that the graph can be zoomed (pinch-to-zoom), different objects moved (drag) and different parts highlighted (click on object) to enable easier reading.
