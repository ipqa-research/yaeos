#!/usr/bin/env bash

set -euo pipefail

TEMPLATE="template.fypp"
OUTDIR="tapenade_inputs"

mkdir -p "${OUTDIR}"

# -----------------------------------------------------------------------------
# PC-SAFT
# -----------------------------------------------------------------------------
model_name="PCSAFT"

fypp \
    -DATTRIBUTES='"attributes/pcsaft_attributes"' \
    -DSETUP='"pcsaft_setup"' \
    -DDIAMETER='"pcsaft_diameter"' \
    "${TEMPLATE}" \
    "${OUTDIR}/${model_name,,}.f90"

sed -i "s/ModelNameSedMe/${model_name}/g" "${OUTDIR}/${model_name,,}.f90"

echo "Generated ${model_name}.f90"

