FROM rocker/shiny-verse:4.2.2

# Apt packages
RUN apt-get update && apt-get install -y \
    aptitude \
    libcurl4-openssl-dev \
    libxml2-dev \
    libglpk-dev \
    && rm -rf /var/lib/apt/lists/*

# Install Bioconductor and its packages
RUN Rscript -e 'if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")'
RUN Rscript -e 'BiocManager::install()'
RUN Rscript -e 'BiocManager::install("ggtree")'
RUN Rscript -e 'BiocManager::install("InterMineR")'
RUN Rscript -e 'BiocManager::install("limma")'
RUN Rscript -e 'BiocManager::install("msa")'
RUN Rscript -e 'BiocManager::install("biomaRt")'

# Install CRAN packages
RUN Rscript -e 'install.packages("dplyr")'
RUN Rscript -e 'install.packages("stringr")'
RUN Rscript -e 'install.packages("ggplot2")'
RUN Rscript -e 'install.packages("msaR")'
RUN Rscript -e 'install.packages("RColorBrewer")'
RUN Rscript -e 'install.packages("shinycssloaders")'
RUN Rscript -e 'install.packages("scales")'
RUN Rscript -e 'install.packages("ape")'
RUN Rscript -e 'install.packages("reshape2")'
RUN Rscript -e 'install.packages("remotes")'
RUN Rscript -e 'install.packages("bio3d", dependencies=TRUE)'

RUN Rscript -e 'remotes::install_github("nvelden/NGLVieweR")'

# Run rootless
USER shiny

CMD ["Rscript", "/reelgene/app.R"]
