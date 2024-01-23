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
  libnlopt-dev \
  cmake \
  curl # renv needs it to avoid FD_SETSIZE errors

RUN apt autoremove -y

RUN R -e 'install.packages("remotes")'
RUN R -e 'remotes::install_version("renv", version="'${RENV_VERSION}'")'

RUN mkdir /renvtmp /install /app
COPY TIFF/renv.lock /renvtmp/renv.lock

RUN cd /renvtmp && R -e 'options(timeout=600); renv::settings$use.cache(FALSE); renv::restore(library = .libPaths()[1], exclude = c("XIFF", "CLIFF"))'
RUN rm -rf /renvtmp

# install phate python lib required by phateR R package
RUN apt-get install -y python3-pip libpython3-dev && \
 python3 -m pip install --no-cache-dir --break-system-packages phate && \
 apt-get remove -y python3-pip && \
 apt autoremove -y

COPY XIFF /install/XIFF
RUN R -e 'renv::install("/install/XIFF")'

COPY CLIFF /install/CLIFF
RUN R -e 'renv::install("/install/CLIFF")'

COPY TIFF /install/TIFF
RUN R -e 'renv::install("/install/TIFF")'

#RUN R -e 'shiny::runTests("/app/general")'

RUN rm -rf /install

RUN chmod -R 655 /app

RUN useradd -m tiffuser
USER tiffuser

WORKDIR /app
EXPOSE 3838

CMD ["Rscript", "-e", "TIFF::run(port=3838, host='0.0.0.0')"]
