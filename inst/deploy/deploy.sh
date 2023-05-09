#!/usr/bin/sh

cd ..
pwd
cp Docker/Dockerfile /nas1/.scratch/TIFF/
cp Docker/Rprofile.site /nas1/.scratch/TIFF/
cp Docker/install_prereqs_tiff.R /nas1/.scratch/TIFF/
ssh charlotte "cd /nas1/.scratch/TIFF;docker build . --no-cache -t tiff/broadinstitute.org"
#ssh charlotte "cd /nas1/.scratch/TIFF;docker build . -t tiff/broadinstitute.org"
