# Running the targets pipeline for producing an indicator on chemical oxygen demand in Finnish lakes

library(htmltools)
library(targets)

# Command to run pipeline
tar_make()

# Command to save the pipeline of produced objects in an html file 
save_html(html = tar_glimpse(), "output/objects_pipeline.html")

# Command to save a graph depicting the produced objects and the functions used to produce them
save_html(html = tar_visnetwork(), "output/network_pipeline.html")
