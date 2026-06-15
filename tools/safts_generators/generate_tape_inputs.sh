#!/usr/bin/env bash

set -euo pipefail

TEMPLATE="template_master.fypp"
OUTDIR="tampenade_inputs"

mkdir -p "${OUTDIR}"

# -----------------------------------------------------------------------------
# PC-SAFT
# -----------------------------------------------------------------------------
fypp \
    -D MODEL_NAME="'PCSAFT'" \
    -D USE_DIAMETER=1 \
    "${TEMPLATE}" \
    "${OUTDIR}/PCSAFT.f90"

echo "Generated PCSAFT.f90"

