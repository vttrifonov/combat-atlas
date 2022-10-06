Sys.setenv(RETICULATE_PYTHON="./env/bin/python3.7")
Sys.setenv(RENV_CONFIG_SYNCHRONIZED_CHECK = FALSE)
options(rlib_downstream_check = FALSE)
source("renv/activate.R")
if (interactive() && Sys.getenv("RSTUDIO") == "") {
  Sys.setenv(TERM_PROGRAM = "vscode")
  source(file.path(Sys.getenv(
    if (.Platform$OS.type == "windows") "USERPROFILE" else "HOME"
  ), ".vscode-R", "init.R"))
}
