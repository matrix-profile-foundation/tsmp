FROM gitpod/workspace-full

LABEL org.label-schema.license="GPL-3.0" \
      org.label-schema.vcs-url="https://github.com/matrix-profile-foundation/tsmp" \
      org.label-schema.vendor="Matrix Profile Foundation" \
      maintainer="Francisco Bischoff <fbischoff@med.up.pt>" \
      description="Docker for testing TSMP"

## Begin Root tasks ##
USER root
### base ###
ENV R_BASE_VERSION=3.6.2
ENV LC_ALL=en_US.UTF-8
ENV LANG=en_US.UTF-8
RUN yes | unminimize \
  && apt-get install -yq \
      r-base \
  && locale-gen en_US.UTF-8 \
  && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/*

RUN Rscript -e 'if (!require("devtools")) install.packages("devtools")'
COPY --chown=gitpod:gitpod * /workspace/tsmp/
RUN Rscript -e 'setwd("/workspace/tsmp"); devtools::install_deps(dep = TRUE)'
# RUN Rscript -e 'devtools::install_github("jimhester/covr")
## End Root tasks ##

## Begin User tasks ##
USER gitpod
ENV R_BASE_VERSION=3.6.2
ENV LC_ALL=en_US.UTF-8
ENV LANG=en_US.UTF-8
ENV HOME=/home/gitpod
### checks ###
# no root-owned files in the home directory
RUN notOwnedFile=$(find . -not "(" -user gitpod -and -group gitpod ")" -print -quit) \
  && { [ -z "$notOwnedFile" ] \
    || { echo "Error: not all files/dirs in $HOME are owned by 'gitpod' user & group"; exit 1; } }
    
RUN Rscript -e 'setwd("/workspace/tsmp"); devtools::check(args = c("--as-cran"), env_vars = NULL)'
CMD ["R"]
## End user tasks ##

## Give back control to the engine ##
USER root
