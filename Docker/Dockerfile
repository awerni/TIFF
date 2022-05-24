FROM wernitznig/xiff

LABEL maintainer="Andreas Wernitznig"
LABEL email="awernitz@broadinstitute.org"

ENV dbhost=google
ENV dbport=5432
ENV dbname=bioinfo_hg38
ENV dbuser=reader
ENV dbpass=reader

RUN mkdir /app
COPY DESCRIPTION /app/DESCRIPTION
RUN R -e 'remotes::install_deps("/app", upgrade="never", repos=BiocManager::repositories(), dependencies=TRUE)'

COPY . /app
RUN R CMD INSTALL /app
RUN export OPENSSL_CONF="/tmp/openssl.cnf" && R -e 'testthat::test_local("/app", stop_on_failure=TRUE)'
RUN rm -rf /app

# set access and port
COPY Docker/Rprofile.site /usr/lib/R/etc/
EXPOSE 3838

CMD ["Rscript", "-e", "TIFF::run()"]
