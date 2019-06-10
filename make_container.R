require(containerit)
require(here)

# use a custom startup command
scriptCmd <- CMD_Rscript("other_dependencies.R")

# create Dockerfile for the script
dockerfile_object <- 
  dockerfile(from = here("code_sample_to_containerize.R"), 
             copy = "script_dir",
             maintainer = "agzuurp@unal.edu.co",
             silent = TRUE,
             cmd = scriptCmd)


# Warning messages:
#   1: In FUN(X[[i]], ...) :
#   Failed to identify a source for package GENIE3. Therefore the package cannot be installed in the Docker image.
# 
# 2: In FUN(X[[i]], ...) :
#   Failed to identify a source for package fastGeneMI. Therefore the package cannot be installed in the Docker image.
# 
# 3: In FUN(X[[i]], ...) :
#   Failed to identify a source for package minet. Therefore the package cannot be installed in the Docker image.

print(dockerfile_object)
#containerit::write(dockerfile_object, "dockerfile")
