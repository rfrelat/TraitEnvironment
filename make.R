# render all files
rmarkdown::render("Tutorials/CleanDataTER.Rmd", output_dir = "docs/")
rmarkdown::render("Tutorials/CWM.Rmd", output_dir = "docs/")
rmarkdown::render("Tutorials/HMSC.Rmd", output_dir = "docs/")
rmarkdown::render("Tutorials/RLQ.Rmd", output_dir = "docs/")
rmarkdown::render("Tutorials/DC-CA_douconca.Rmd", output_dir = "docs/")
rmarkdown::render("Tutorials/DC-CA_old.Rmd", output_dir = "docs/")

# create the index
rmarkdown::render("ReadMe.md", output_dir = "docs/", output_file = "index.html")
