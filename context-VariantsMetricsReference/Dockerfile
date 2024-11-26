# Install Ubuntu 20.04 (focal) with R base and R version management (https://rocker-project.org/images/versioned/r-ver.html)
FROM rocker/r-ver:4.4.0
# set the default repository to CRAN
RUN echo 'options(repos = c(CRAN = "https://cloud.r-project.org"))' >>"${R_HOME}/etc/Rprofile.site"
# sets the working directory for any RUN, CMD, ENTRYPOINT, COPY and ADD instructions
WORKDIR /context-VariantsMetricsReference
# copy R scripts
COPY *.R ./
# install a sufficient Python binary in Ubuntu 20.04 required for the argparse library
RUN apt update
RUN apt install -y software-properties-common build-essential libffi-dev libssl-dev zlib1g-dev libncurses5-dev libncursesw5-dev libreadline-dev libsqlite3-dev libgdbm-dev libdb5.3-dev libbz2-dev libexpat1-dev liblzma-dev libffi-dev libssl-dev
RUN add-apt-repository ppa:deadsnakes/ppa
RUN apt update
RUN apt install -y python3.12 python3.12-venv
# install R libraries
RUN Rscript /context-VariantsMetricsReference/install_packages.R
# set entrypoint
ENTRYPOINT  ["Rscript", "/context-VariantsMetricsReference/VariantsMetricsReference:1.0.R"]
# sets the working directory for any RUN, CMD, ENTRYPOINT, COPY and ADD instructions
WORKDIR /wd
# set Docker image version
ADD VERSION .
# get the usage when launching rather the help
CMD ["/bin/bash"]
