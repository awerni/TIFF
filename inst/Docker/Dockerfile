FROM rocker/r-base:4.3.1

ARG RENV_VERSION=1.0.3

ENV dbhost=google
ENV dbport=5432
ENV dbname=bioinfo_hg38
ENV dbuser=reader
ENV dbpass=reader

RUN apt-get update && apt-get install -y \
  libssl-dev \
  libpq-dev \
  libxml2-dev \
  libcurl4-gnutls-dev \
  libfontconfig1-dev \
  libcairo2-dev \
  cmake \
  libnlopt-dev  \
  libharfbuzz-dev libfribidi-dev \ 
  libpng-dev \
  libtiff5-dev \
  libjpeg-dev \
  libfreetype-dev

RUN apt autoremove -y

RUN R -e 'install.packages("remotes")'
RUN R -e 'remotes::install_version("renv", version="'${RENV_VERSION}'")'

RUN mkdir /renvtmp /app
COPY TIFF/renv.lock /renvtmp/renv.lock

# install "heavy deps" separately to avoid strange error: bit out of range 0 FD_SETSIZE on fd_set
RUN cd /renvtmp && R -e 'options(timeout=600); renv::settings$use.cache(FALSE); renv::restore(library = .libPaths()[1], packages = "shiny")'
RUN cd /renvtmp && R -e 'options(timeout=600); renv::settings$use.cache(FALSE); renv::restore(library = .libPaths()[1], packages = "dplyr")'
RUN cd /renvtmp && R -e 'options(timeout=600); renv::settings$use.cache(FALSE); renv::restore(library = .libPaths()[1], packages = "ggtips")'
RUN cd /renvtmp && R -e 'options(timeout=600); renv::settings$use.cache(FALSE); renv::restore(library = .libPaths()[1], packages = "DESeq2")'
RUN cd /renvtmp && R -e 'options(timeout=600); renv::settings$use.cache(FALSE); renv::restore(library = .libPaths()[1], packages = "tidyverse")'
RUN cd /renvtmp && R -e 'options(timeout=600); renv::settings$use.cache(FALSE); renv::restore(library = .libPaths()[1], packages = "shinytest")'
RUN cd /renvtmp && R -e 'options(timeout=600); renv::settings$use.cache(FALSE); renv::restore(library = .libPaths()[1], exclude = "XIFF")'
RUN rm -rf /renvtmp

COPY XIFF /app/XIFF
COPY TIFF /app/TIFF

RUN R -e 'renv::install("/app/XIFF")'
RUN R -e 'renv::install("/app/TIFF")'

#RUN R -e 'shiny::runTests("/app/general")'

RUN rm -rf /app

RUN useradd -m tiffuser
USER tiffuser

EXPOSE 3838
WORKDIR /home/tiffuser

CMD ["Rscript", "-e", "TIFF::run(port=3838, host='0.0.0.0')"]
