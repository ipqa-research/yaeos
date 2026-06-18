#!/usr/bin/env bash

set -euo pipefail

TEMPLATE="template.fypp"
OUTDIR="tapenade_inputs"

mkdir -p "${OUTDIR}"

# -----------------------------------------------------------------------------
# PC-SAFT
# -----------------------------------------------------------------------------
model_name="PCSAFT"
attributes="PCSAFT-ATTRIBUTES"
diameter="PCSAFT-DIAMETER"
zetas="PCSAFT-ZETAS"
ar_hard_sphere="Boublik-Mansoori-Carnahan-Starling"
ar_chain="Chapman-PCSAFT"
ar_dispersion="Dispersion-PCSAFT"

fypp \
    -D MODELNAME="'${model_name}'" \
    -D ATTRIBUTES="'${attributes}'" \
    -D DIAMETER="'${diameter}'" \
    -D ZETAS="'${zetas}'" \
    -D AR_HARD_SPHERE="'${ar_hard_sphere}'" \
    -D AR_CHAIN="'${ar_chain}'" \
    -D AR_DISPERSION="'${ar_dispersion}'" \
    "${TEMPLATE}" "${OUTDIR}/${model_name,,}.f90"

sed -i "s/ModelNameSedMe/${model_name}/g" "${OUTDIR}/${model_name,,}.f90"
echo "Generated ${model_name}.f90"

