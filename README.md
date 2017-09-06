# varbvs-vs-gemma

Experiments comparing varbvs and GEMMA.

## How to build the webpages

Run the following commands in R from the [analysis](analysis)
directory:

```R
library(rmarkdown)
render("demo.sim.plots.Rmd",envir = new.env(),output_dir = "../docs")
```
