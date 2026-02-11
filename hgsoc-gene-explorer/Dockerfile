# Base image with R + Shiny Server pre-installed
FROM rocker/shiny:4.4.1

# Make a directory in the container
RUN mkdir /home/shiny-app

# Set working directory
WORKDIR /home/shiny-app

# Install system libraries you need (needed for ggplot2 and data.table and dplyr)

RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libpng-dev \
    libxxf86vm1 \
    libnlopt-dev \
    cmake \
    gfortran \
    build-essential \
    && rm -rf /var/lib/apt/lists/*



# Install problem packages first
RUN R -e "install.packages('ggpubr', repos='https://cloud.r-project.org'); sessionInfo()"
RUN R -e "install.packages('BiocManager'); BiocManager::install('ComplexHeatmap'); sessionInfo()"
# Install the necessary R package dependencies
RUN R -e "install.packages(c('shiny', 'shinycssloaders', 'here', 'dplyr', 'data.table', 'ggplot2', 'cowplot', 'ggrepel', 'patchwork', 'colorRamp2'), repos='https://cloud.r-project.org')"

# Copy this app to the Shiny Server directory
COPY . /home/shiny-app/

# Allow shiny-server to read/write logs
RUN chown -R shiny:shiny /home/shiny-app

# Expose Shiny Server port
EXPOSE 3838

# (Optional) disable the default example apps
#RUN rm -rf /srv/shiny-server/index.html \
#           /srv/shiny-server/sample-apps

# Start Shiny Server
CMD ["Rscript", "app.R"]
