# _targets.R file
library(targets)
source("Script/Functions/Vesla_functions_COD.r", encoding= "ISO-8859-1")
options(tidyverse.quiet = TRUE, download.file.method = "wininet")
tar_option_set(packages = c("dplyr", "brms", "readr", "tidyr", "data.table", "zoo"))
list(
  
  # 1. Obtain data
  
  ## 1.1 Create a query object to send to VESLA
  
  tar_target(
             VESLA_query_data_file,
             "data/Otantalista_SV7_humuspitoisuus_jarvet_vesimuodostumat.csv",
             format="file"
             ),
  tar_target(
             VESLA_raw_query_data,
             read.table(VESLA_query_data_file, sep=",", header=TRUE)
             ),
  tar_target(
             VESLA_query_data,
             unique(VESLA_raw_query_data)
             ),
  
  ## 1.2 Obtain data from http://rajapinnat.ymparisto.fi/api/vesla/2.0/
  
  ## 1.2.1 Get sampling site IDs for the assigned lakes
  tar_target(
             VESLA_site_ids,
             get_VESLA_sampling_sites(VESLA_query_data)#,
             #cue = tar_cue(mode = "always")
             ),
  
  ## 1.2.2 Get all results for the associated sampling sites
  tar_target(
             VESLA_sites_all_results,
             get_VESLA_all_data(VESLA_site_ids$Paikka_Id)
             ),
  
  
  # 2. Shape, restructure, and combine data
  
  ## 2.1 Select appropriate subsets of data
  
  ### 2.1.1 Filter out sampling events prior 1960 
  tar_target(
             VESLA_sampling_events,
             shape_VESLA_data_events(VESLA_sites_all_results, "1960-01-01")
             ),
  
  ### 2.1.2 Select samples collected between 1 and 2 meters depth
  tar_target(
             VESLA_samples,
             shape_VESLA_data_samples(VESLA_sites_all_results, c(1,2))
             ),
  
  ### 2.1.3 Appropriate analysis results
  tar_target(
             VESLA_results,
             shape_VESLA_data_results(VESLA_sites_all_results, c(25,27,3293))
             ),
  
  
  ## 2.2 Combine above pieces together
  tar_target(
             VESLA_COD_data_combined,
             shape_VESLA_data_combine(VESLA_results, VESLA_samples, VESLA_sampling_events, VESLA_site_ids)
             ),
  
  ## 2.3 Create summaries 
  tar_target(
             VESLA_summary_COD_data_Site,
             summarize_VESLA_data(VESLA_COD_data_combined, "Paikka_Id")
             ),
  tar_target(
             VESLA_summary_COD_data_Lake,
             summarize_VESLA_data(VESLA_COD_data_combined, "J_Jarvi_Id")
             ),
  
  ## 2.4 shape and filter data for mean comparisons (i.e. status determination)
  
  ### 2.4.1 Find typical local minimums in time series to allow for status determination
  tar_target(
             VESLA_COD_minimum_timeframe,
             analyze_VESLA_minimum_MA(VESLA_summary_COD_data_Lake, 5)
             ),
  
  ### 2.4.2 Use the obtained minimum to create a comparison data set
  tar_target(
             VESLA_summary_COD_data_Site_comparison,
             shape_VESLA_data_comparison(VESLA_summary_COD_data_Site,
                                         VESLA_COD_minimum_timeframe,
                                         c((max(VESLA_summary_COD_data_Site$na)-4),max(VESLA_summary_COD_data_Site$na)))
             ),
  
  ## 2.5 shape and filter data for trend analyses (i.e. trend determination)
  tar_target(
             VESLA_summary_COD_data_Site_old,
             shape_VESLA_data_time_region(VESLA_summary_COD_data_Site,
                                          c(min(VESLA_summary_COD_data_Site$na),1999))
             ),
  tar_target(
             VESLA_summary_COD_data_Site_new,
             shape_VESLA_data_time_region(VESLA_summary_COD_data_Site,
                                          c(2001, max(VESLA_summary_COD_data_Site$na)))
             ),
  
  
  # 3. Analyze data
  
  ## 3.1 analyze mean differences (i.e. determine status)
  
  tar_target(Comparison_COD,
             analyze_VESLA_status(VESLA_summary_COD_data_Site_comparison)
             ),
  
  
  ## 3.2 analyze (linear) trends
  
  ### 3.2.1 analyze linear trend prior 2000
  tar_target(
             COD_trend_old,
             analyze_VESLA_trend(VESLA_summary_COD_data_Site_old)
             ),
  
  ### 3.2.2 analyze linear trend post 2000
  tar_target(
             COD_trend_new,
             analyze_VESLA_trend(VESLA_summary_COD_data_Site_new)
             ),
  
  
  
  ## 3.3 create a smoothed indicator graph
  
  tar_target(GAM_fit_COD,
             analyze_VESLA_GAM(VESLA_summary_COD_data_Site)),
             #cue = tar_cue(mode = "always")),
  
  
  # 4. Create output
  
  tar_target(
             COD_status_output,
             output_table_VESLA_status_comparison(Comparison_COD, 
                                                  VESLA_summary_COD_data_Site_comparison)
             ),
  
  tar_target(
             COD_trend_output_old,
             output_table_VESLA_trend(COD_trend_old, VESLA_summary_COD_data_Site_old)
             ),
  
  tar_target(
             COD_trend_output_new,
             output_table_VESLA_trend(COD_trend_new, VESLA_summary_COD_data_Site_new)
             ),
  
  tar_target(
             COD_trend_output_combined,
             rbind(#COD_trend_output_all,
                   COD_trend_output_old,
                   COD_trend_output_new)
             ),
  
  
  # 5 Write output to disk
  tar_target(outputfile1,
             write.table(COD_status_output, 
                         file="output/VESLA_COD_status.csv", 
                         row.names = FALSE, 
                         sep=";", 
                         dec=",")
             ),
  tar_target(outputfile2,
             write.table(COD_trend_output_combined, 
                         file="output/VESLA_COD_trends.csv", 
                         row.names = FALSE, 
                         sep=";", 
                         dec=",")
             ),
  
  tar_target(COD_fit_GAM_trend_output_tidy,
             output_table_tidy_VESLA_GAM(GAM_fit_COD)
             )
             
  )

