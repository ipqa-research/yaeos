#!/bin/bash
apt-get update && apt-get install wget
wget https://raw.githubusercontent.com/ipqa-research/fortran-setup/main/bootstrap_fortran.sh
bash bootstrap_fortran.sh
