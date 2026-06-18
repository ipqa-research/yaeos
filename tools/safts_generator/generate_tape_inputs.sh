#!/usr/bin/env bash

set -euo pipefail

TEMPLATE="template.fypp"
TEMPLATE_SETUP="template_setup.fypp"
TEMPLATE_GET_V0="template_get_v0.fypp"
OUTDIR="tapenade_inputs"

mkdir -p "${OUTDIR}"

# -----------------------------------------------------------------------------
# PC-SAFT
# -----------------------------------------------------------------------------
model_name="TPCSAFT"
attributes="PCSAFT-ATTRIBUTES"
diameter="PCSAFT-DIAMETER"
zetas="PCSAFT-ZETAS"
ar_hard_sphere="Boublik-Mansoori-Carnahan-Starling"
ar_chain="Chapman-PCSAFT"
ar_dispersion="Dispersion-PCSAFT"
get_v0="PCSAFT_V0"

fypp \
    -D MODELNAME="'${model_name}'" \
    -D ATTRIBUTES="'${attributes}'" \
    -D DIAMETER="'${diameter}'" \
    -D ZETAS="'${zetas}'" \
    -D AR_HARD_SPHERE="'${ar_hard_sphere}'" \
    -D AR_CHAIN="'${ar_chain}'" \
    -D AR_DISPERSION="'${ar_dispersion}'" \
    "${TEMPLATE}" "${OUTDIR}/${model_name,,}.f90"

fypp \
    -D MODELNAME="'${model_name,,}'" \
    -D ATTRIBUTES="'${attributes}'" \
    "${TEMPLATE_SETUP}" "${OUTDIR}/${model_name,,}_setup.f90"

fypp \
    -D GETV0="'${get_v0}'" \
    "${TEMPLATE_GET_V0}" "${OUTDIR}/${model_name,,}_v0.f90"

sed -i "s/ModelNameSedMe/${model_name}/g" "${OUTDIR}/${model_name,,}.f90"
sed -i "s/ModelNameSedMe/${model_name}/g" "${OUTDIR}/${model_name,,}_setup.f90"
sed -i "s/ModelNameSedMe/${model_name}/g" "${OUTDIR}/${model_name,,}_v0.f90"

bash ../tapenade_diff/gen_tapemodel.sh "${OUTDIR}/${model_name,,}.f90"

sed -i '/^! PUT setup and get_v0 here after tapenade$/{
    i\

    r '"${OUTDIR}/${model_name,,}_setup.f90"'
    a\

    r '"${OUTDIR}/${model_name,,}_v0.f90"'
    d
}' "tapeout/${model_name,,}_diff.f90"

echo "Generated ${model_name}.f90"

