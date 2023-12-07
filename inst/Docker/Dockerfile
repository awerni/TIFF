FROM wernitznig/shiny

LABEL maintainer="Andreas Wernitznig"
LABEL email="andreas.wernitznig@pierre-fabre.com"

ENV dbhost=azure
ENV dbport=5432
ENV dbname=bioinfo_hg38
ENV dbuser=reader
ENV dbpass=reader

RUN apt-get update && apt-get install -y \
 cmake

# install additional packages
COPY install_prereqs_tiff.R /install/install_prereqs_tiff.R
RUN Rscript /install/install_prereqs_tiff.R

# set access and port
COPY Rprofile.site /usr/lib/R/etc/

EXPOSE 3838

CMD ["Rscript", "-e", "TIFF::run()"]

