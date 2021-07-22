#! /bin/bash

sudo apt-get -y update
sudo apt-get -y install gdebi-core
sudo apt-get -y install r-base r-base-dev
sudo apt-get -y install libcurl4-openssl-dev libssl-dev libxml2-dev
sudo apt-get -y install wget
sudo apt-get -y install jq
wget https://s3.amazonaws.com/rstudio-ide-build/server/bionic/amd64/rstudio-server-1.4.1118-amd64.deb
sudo gdebi -n rstudio-server-1.4.1118-amd64.deb
sudo apt-get -y install python3-distutils
sudo apt-get -y install python3-apt
sudo apt-get -y install python3-pip
pip3 install dsub
# sudo Rscript -e 'install.packages(c("bigQueryR", "bigrquery", "googleCloudStorageR", "data.table", "dplyr", "arrow"))'
sudo Rscript -e 'install.packages(c("bigQueryR", "bigrquery", "googleCloudStorageR", "data.table", "dplyr"))'