#!/bin/bash
#SBATCH --job-name=transgam-sim
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=72:00:00
#SBATCH --array=1-10
#SBATCH --output=transgam-sim-%A_%a.out

set -euo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV:-transgam-r}"

cd "${SLURM_SUBMIT_DIR}"

MODEL="${MODEL:-weibull}"
N_S="${N_S:-2000}"
N_T="${N_T:-300}"
N_TEST="${N_TEST:-1000}"
N_REP_TOTAL="${N_REP_TOTAL:-100}"
REPS_PER_TASK="${REPS_PER_TASK:-10}"
SETTINGS="${SETTINGS:-all}"
TQ="${TQ:-{0.25,0.50,0.75}}"
METHODS="${METHODS:-stat,transfer,gam_gam_basis,oracle}"
OUT_ROOT="${OUT_ROOT:-${SLURM_SUBMIT_DIR}/results/${MODEL}_${N_T}}"

REP_START=$(( (SLURM_ARRAY_TASK_ID - 1) * REPS_PER_TASK + 1 ))
REP_END=$(( SLURM_ARRAY_TASK_ID * REPS_PER_TASK ))

if [ "${REP_START}" -gt "${N_REP_TOTAL}" ]; then
  echo "REP_START=${REP_START} exceeds N_REP_TOTAL=${N_REP_TOTAL}; nothing to run."
  exit 0
fi

if [ "${REP_END}" -gt "${N_REP_TOTAL}" ]; then
  REP_END="${N_REP_TOTAL}"
fi

CHUNK_TAG="rep_${REP_START}_${REP_END}"
mkdir -p "${OUT_ROOT}/chunks"

echo "Running MODEL=${MODEL}, N_T=${N_T}, ${CHUNK_TAG}, METHODS=${METHODS}"

Rscript scripts/main_simulation_server.R \
  --model "${MODEL}" \
  --methods "${METHODS}" \
  --settings "${SETTINGS}" \
  --n-rep "${N_REP_TOTAL}" \
  --rep-start "${REP_START}" \
  --rep-end "${REP_END}" \
  --n-s "${N_S}" \
  --n-t "${N_T}" \
  --n-test "${N_TEST}" \
  --tq "${TQ}" \
  --out-root "${OUT_ROOT}/chunks/${CHUNK_TAG}" \
  --cores "${SLURM_CPUS_PER_TASK}"
