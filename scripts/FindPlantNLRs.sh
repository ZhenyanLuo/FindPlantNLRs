#!/usr/bin/env bash

set -euo pipefail

cd /work

cp --update=none /opt/FindPlantNLRs/Snakefile .
cp --update=none /opt/FindPlantNLRs/FindPlantNLRs.config .
cp --update=none /opt/FindPlantNLRs/NLR_Pfams.list .
cp --recursive --update=none /opt/FindPlantNLRs/ref_db .

snakemake --cores all --latency-wait 60 --keep-going