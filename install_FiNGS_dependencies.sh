#!/bin/bash

apt-get update
apt-get -y install python3
apt-get -y install python3-pip
pip3 install PyVCF
pip3 install pysam
pip3 install editdistance
pip3 install scipy
pip3 install joblib

apt-get -y install r-base r-base-dev
apt-get -y install libcurl4-openssl-dev
apt-get -y install libxml2-dev
apt-get -y install zlib1g-dev
Rscript --vanilla R/install_packages.R

