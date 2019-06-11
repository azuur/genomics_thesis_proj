FROM rocker/verse
LABEL maintainer="agzuurp@unal.edu.co"
#RUN export DEBIAN_FRONTEND=noninteractive
RUN  apt-get -y update \   
  && apt-get install -y git-core \
	libgmp-dev \
	libv8-dev \
#	libcurl4-openssl-dev \
	libgit2-dev \
	libssl-dev \
	libxml2-dev \
	libcurl4-gnutls-dev \
	wget \
	gpg \
	make



# Switch to root for install
#USER root

# Install wget
#RUN apt-get update -y && apt-get install -y \
#	wget \
#	build-essential \
#	--no-install-recommends \
#	&& rm -rf /var/lib/apt/lists/*

# Install glpk from http
# instructions and documentation for glpk: http://www.gnu.org/software/glpk/
WORKDIR /user/local/
RUN wget http://ftp.gnu.org/gnu/glpk/glpk-4.57.tar.gz \
	&& tar -zxvf glpk-4.57.tar.gz

## Verify package contents
# RUN wget http://ftp.gnu.org/gnu/glpk/glpk-4.57.tar.gz.sig \
#	&& gpg --verify glpk-4.57.tar.gz.sig
#	#&& gpg --keyserver keys.gnupg.net --recv-keys 5981E818

WORKDIR /user/local/glpk-4.57
RUN ./configure \
	&& make \
	&& make check \
	&& make install \
	&& make distclean \
	&& ldconfig \
# Cleanup
	&& rm -rf /user/local/glpk-4.57.tar.gz \
	&& apt-get clean

#create a glpk user
#ENV HOME /home/user
#RUN useradd --create-home --home-dir $HOME user \
#    && chmod -R u+rwx $HOME \
#    && chown -R user:user $HOME

# switch back to user
#WORKDIR $HOME



RUN ["install2.r", "abind", "assertthat", "backports", "bibtex", "BiocManager", "bnlearn", "broom", "callr", "cellranger", "cli", "codetools", "colorspace", "corpcor", "crayon", "DEoptimR", "desc", "devtools", "digest", "doParallel", "doRNG", "foreach", "fs", "generics", "glue", "gtable", "haven", "here", "hms", "httr", "igraph", "infotheo", "iterators", "jsonlite", "lars", "lattice", "lazyeval", "lubridate", "magrittr", "memoise", "modelr", "munsell", "nlme", "ParallelPC", "parmigene", "pillar", "pkgbuild", "pkgconfig", "pkgload", "pkgmaker", "plyr", "prettyunits", "processx", "ps", "R6", "Rcpp", "readxl", "registry", "remotes", "rlang", "rngtools", "robustbase", "rprojroot", "rstudioapi", "rvest", "scales", "sessioninfo", "stringi", "testthat", "tidyselect", "usethis", "withr", "xml2", "xtable"]

#, "dplyr", "tidyverse", "ggplot2", "readr", "purrr", "tibble", "tidyr", "stringr", "forcats"

RUN install2.r --error \ 
    -r 'http://cran.rstudio.com' \
    googleAuthR \
    devtools \
    && Rscript -e "BiocManager::install('graph')" \
    && Rscript -e "BiocManager::install('RBGL')" \
    && Rscript -e "BiocManager::install('minet')" \
    && Rscript -e "BiocManager::install('Rglpk')" \
    && Rscript -e "devtools::install_bitbucket('Jonathan-Ish-Horowicz/fastGeneMI', force = T) " \
    && Rscript -e "BiocManager::install(c('minet','graph','RBGL','GENIE3'), ask = FALSE)" \
    ## clean up
    && rm -rf /tmp/downloaded_packages/ /tmp/*.rds




RUN ["installGithub.r", "jpvert/tigress@de70cd256840b08a64f0471c0c63dbb01314b35a"]


WORKDIR /payload/
COPY ["./", "./"]


#CMD ["R"]
#CMD ["R", "--vanilla", "-f", "code_sample_to_containerize.R"]
