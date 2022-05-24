#!/bin/sh

targetDir=/data/bioinf/ShinyApps/tiff

rm -f ${targetDir}/*.R
rm -f ${targetDir}/www/*

cp -a ../shiny/ui.R ../shiny/server.R ../shiny/global.R ../shiny/config.yml ${targetDir}/
cp -a ../shiny/www/TIFF.png ../shiny/www/main.css ${targetDir}/www/
touch ${targetDir}/restart.txt
