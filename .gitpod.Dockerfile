FROM gitpod/workspace-full

USER root
### base ###
RUN yes | unminimize \
    && apt-get install -yq \
        r-base \
    && locale-gen en_US.UTF-8 \
    && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/*
ENV LANG=en_US.UTF-8

USER gitpod

USER root
