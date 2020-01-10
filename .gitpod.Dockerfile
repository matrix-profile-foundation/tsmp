FROM gitpod/workspace-full

LABEL org.label-schema.license="GPL-3.0" \
    org.label-schema.vcs-url="https://github.com/matrix-profile-foundation/tsmp" \
    org.label-schema.vendor="Matrix Profile Foundation" \
    maintainer="Francisco Bischoff <fbischoff@med.up.pt>" \
    description="Docker for testing TSMP"

## Begin Root tasks ##
USER root
### base ###
ENV LC_ALL=en_US.UTF-8
ENV LANG=en_US.UTF-8

RUN echo "deb https://cloud.r-project.org/bin/linux/ubuntu disco-cran35/" >> /etc/apt/sources.list
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN apt-get update
RUN apt-get install -yq \
    r-base \
    && locale-gen en_US.UTF-8 \
    && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/*

COPY --chown=gitpod:gitpod * /workspace/tsmp/

## End Root tasks ##

## Begin User tasks ##
USER gitpod
ENV LC_ALL=en_US.UTF-8
ENV LANG=en_US.UTF-8
ENV HOME=/home/gitpod

RUN Rscript -e 'if(!dir.exists(Sys.getenv("R_LIBS_USER"))) dir.create(Sys.getenv("R_LIBS_USER"), recursive=TRUE)'
RUN Rscript -e 'if (!require("devtools")) install.packages("devtools", lib = Sys.getenv("R_LIBS_USER"))'
RUN Rscript -e 'if (!require("languageserver")) install.packages("languageserver", lib = Sys.getenv("R_LIBS_USER"))'
RUN Rscript -e 'devtools::install_github("jimhester/covr", ref = "master")'
RUN Rscript -e 'setwd("/workspace/tsmp"); devtools::install_deps(dep = TRUE)'
### checks ###
# no root-owned files in the home directory
RUN notOwnedFile=$(find . -not "(" -user gitpod -and -group gitpod ")" -print -quit) \
    && { [ -z "$notOwnedFile" ] \
    || { echo "Error: not all files/dirs in $HOME are owned by 'gitpod' user & group"; exit 1; } }

# RUN Rscript -e 'setwd("/workspace/tsmp"); devtools::check(args = c("--as-cran"), env_vars = NULL)'
RUN echo "/usr/bin/R" >> ~/.bashrc
RUN echo ".First <- function() devtools::load_all()" > ~/.Rprofile
## End user tasks ##

## Give back control to the engine ##
USER root
# https://gitpod.io/#https://github.com/matrix-profile-foundation/tsmp/tree/franzbischoff/gitpod-setup
