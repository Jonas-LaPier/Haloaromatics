#!/bin/bash
# Generic Gaussian 16 job-array runner.
# Do not submit directly; use:   python3 hx.py submit <stage>
# (hx.py sets --array, --cpus-per-task, --mem, --time, --job-name and the manifest.)
#
# Usage: sbatch [options] slurm/g16_array.sbatch <manifest.txt>
#   manifest: one .gjf path per line (relative to the project root)

#SBATCH --job-name=g16
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --output=slurm_logs/%x.%A_%a.out
#SBATCH --error=slurm_logs/%x.%A_%a.err

set -u
MANIFEST="$1"
cd "$SLURM_SUBMIT_DIR"

GJF=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MANIFEST")
if [ -z "$GJF" ] || [ ! -f "$GJF" ]; then
    echo "ERROR: no input for task $SLURM_ARRAY_TASK_ID in $MANIFEST"; exit 1
fi

NAME=$(basename "$GJF" .gjf)
STAGE_DIR=$(dirname "$(dirname "$GJF")")
LOG="$STAGE_DIR/logs/$NAME.log"
mkdir -p "$STAGE_DIR/logs" "$STAGE_DIR/chks"

# node-local scratch for Gaussian's .rwf/.int files
export GAUSS_SCRDIR="${L_SCRATCH:-${SCRATCH:-/tmp}}/g16_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$GAUSS_SCRDIR"

echo "task=$SLURM_ARRAY_TASK_ID job=$SLURM_JOB_ID node=$(hostname) input=$GJF start=$(date)"
module load chemistry gaussian
g16 < "$GJF" > "$LOG"
RC=$?
echo "g16 exit=$RC end=$(date)"

if [ "${FORMCHK:-0}" = "1" ] && [ -f "$STAGE_DIR/chks/$NAME.chk" ]; then
    formchk "$STAGE_DIR/chks/$NAME.chk" "$STAGE_DIR/chks/$NAME.fchk"
fi
rm -rf "$GAUSS_SCRDIR"
exit $RC
