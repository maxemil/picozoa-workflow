From:continuumio/miniconda3:4.6.14
Bootstrap:docker

%labels
    MAINTAINER Max Emil Schön <max-emil.schon@icm.uu.se>
    DESCRIPTION Singularity image containing all requirements for the Picozoa single cell pipeline
    VERSION 0.1

%environment
    PATH=/opt/conda/envs/picozoa/bin:$PATH
    export PATH

%files
    environment.yml /

%post
    apt-get update
    apt-get install -y procps libxtst6
    apt-get clean -y

    /opt/conda/bin/conda env create -f /environment.yml
    /opt/conda/bin/conda clean -a
