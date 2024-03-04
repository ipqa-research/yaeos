#!/bin/bash
apt-get update && apt-get install -y wget
wget https://raw.githubusercontent.com/ipqa-research/fortran-setup/main/bootstrap_fortran.sh
printf "n\nn\ny\ny" | bash bootstrap_fortran.sh
