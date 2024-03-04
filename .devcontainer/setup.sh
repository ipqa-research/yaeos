#!/bin/bash
sudo apt-get update && sudo apt-get install wget
wget https://raw.githubusercontent.com/ipqa-research/fortran-setup/main/bootstrap_fortran.sh
bash bootstrap_fortran.sh
