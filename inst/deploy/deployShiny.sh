#!/usr/bin/sh

cd ..
pwd
cp Docker/Dockerfile_shiny /nas1/.scratch/TIFF/
cp Docker/install_prereqs.R /nas1/.scratch/TIFF/ 
ssh charlotte "cd /nas1/.scratch/TIFF/;docker build . -f Dockerfile_shiny -t wernitznig/shiny"
