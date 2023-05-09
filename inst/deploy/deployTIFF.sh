#!/usr/bin/sh

cd ~/git/ 

export dbname=bioinfo_22Q4.hg38
export dbhost=charlotte.mausloch.wernitznig.aut

R CMD build  --no-build-vignettes TIFF 
ver=$(cat ~/git/TIFF/DESCRIPTION  | egrep "Version:" | sed  -e 's/Version: //g')
name="TIFF_${ver}.tar.gz"
echo $name
mv $name /nas1/.scratch/TIFF/

pwd
#ssh charlotte "rm /var/www/html/R/src/contrib/TIFF_*;mv /nas1/.scratch/TIFF/${name} /var/www/html/R/src/contrib/;/var/www/html/R/update_packages.sh"

rm /var/www/html/R/src/contrib/TIFF_*
mv /nas1/.scratch/TIFF/${name} /var/www/html/R/src/contrib/
/var/www/html/R/update_packages.sh
