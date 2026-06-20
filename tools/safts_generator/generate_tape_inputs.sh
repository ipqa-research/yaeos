#!/usr/bin/env bash

set -euo pipefail

TEMPLATE="template.fypp"
TEMPLATE_SETUP="template_setup.fypp"
TEMPLATE_GET_V0="template_get_v0.fypp"

OUTDIR="tapenade_inputs"
TAPEOUT="tapeout"

clean=false
nodiff=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --clean)
            clean=true
            shift
            ;;
        --nodiff)
            nodiff=true
            shift
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

if $clean; then
    rm -rf "${OUTDIR}" "${TAPEOUT}"
    exit 0
fi



mkdir -p "${OUTDIR}"

generate_model() {
    local model_name="$1"
    local attributes="$2"
    local diameter="$3"
    local zetas="$4"
    local ar_hard_sphere="$5"
    local ar_chain="$6"
    local ar_dispersion="$7"
    local get_v0="$8"

    local model_lower="${model_name,,}"

    fypp \
        -D MODELNAME="'${model_name}'" \
        -D ATTRIBUTES="'${attributes}'" \
        -D DIAMETER="'${diameter}'" \
        -D ZETAS="'${zetas}'" \
        -D AR_HARD_SPHERE="'${ar_hard_sphere}'" \
        -D AR_CHAIN="'${ar_chain}'" \
        -D AR_DISPERSION="'${ar_dispersion}'" \
        "${TEMPLATE}" "${OUTDIR}/${model_lower}.f90"

    fypp \
        -D MODELNAME="'${model_lower}'" \
        -D ATTRIBUTES="'${attributes}'" \
        "${TEMPLATE_SETUP}" "${OUTDIR}/${model_lower}_setup.f90"

    fypp \
        -D GETV0="'${get_v0}'" \
        "${TEMPLATE_GET_V0}" "${OUTDIR}/${model_lower}_v0.f90"

    sed -i "s/ModelNameSedMe/${model_name}/g" \
        "${OUTDIR}/${model_lower}.f90"

    sed -i "s/ModelNameSedMe/${model_name}/g" \
        "${OUTDIR}/${model_lower}_setup.f90"

    sed -i "s/modelnamesedme/${model_lower}/g" \
        "${OUTDIR}/${model_lower}.f90"

    sed -i "s/modelnamesedme/${model_lower}/g" \
        "${OUTDIR}/${model_lower}_setup.f90"

    sed -i "s/ModelNameSedMe/${model_name}/g" \
        "${OUTDIR}/${model_lower}_v0.f90"

    if ! $nodiff; then
        bash ../tapenade_diff/gen_tapemodel.sh \
            "${OUTDIR}/${model_lower}.f90"

        sed -i '/^! PUT setup and get_v0 here after tapenade$/{
            i\

            r '"${OUTDIR}/${model_lower}_setup.f90"'
            a\

            r '"${OUTDIR}/${model_lower}_v0.f90"'
            d
        }' "tapeout/${model_lower}_diff.f90"

        sed -i \
            "/^[[:space:]]*END[[:space:]]\+TYPE/! s/TYPE ${model_name}/type, extends(ArModelTapenade) :: ${model_name}/g" \
            "tapeout/${model_lower}_diff.f90"

        sed -i \
            "/^[[:space:]]*END[[:space:]]\+TYPE[[:space:]]\+${model_name}$/c\\
    contains\\
        procedure :: ar\\
        procedure :: ar_d\\
        procedure :: ar_b\\
        procedure :: ar_d_b\\
        procedure :: ar_d_d\\
        procedure :: get_v0 => get_v0\\
    END TYPE ${model_name}
    " "tapeout/${model_lower}_diff.f90"

        sed -i \
            "s/USE YAEOS, ONLY : r/use yaeos__constants, only: pr, R/g" \
            "tapeout/${model_lower}_diff.f90"
    fi

    echo "Generated ${model_name}.f90"
}

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

generate_model \
    "$model_name" \
    "$attributes" \
    "$diameter" \
    "$zetas" \
    "$ar_hard_sphere" \
    "$ar_chain" \
    "$ar_dispersion" \
    "$get_v0"
