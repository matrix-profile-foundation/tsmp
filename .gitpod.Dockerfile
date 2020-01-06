FROM gitpod/workspace-full

USER root
### base ###
RUN yes | unminimize \
    && apt-get install -yq \
        r-base \
    && locale-gen en_US.UTF-8 \
    && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/*
RUN Rscript -e 'if (!require("devtools")) install.packages("devtools")'
USER gitpod
ENV LANG=en_US.UTF-8
USER root
