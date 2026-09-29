#!/bin/bash
# Generic Gaussian 16 job-array runner.
# Do not submit directly; use:   python3 hx.py submit <stage>
# (hx.py sets --array, --cpus-per-task, --mem, --time, --job-name and the manifest.)
#
# Usage: sbatch [options] slurm/g16_array.sbatch <manifest.txt>
#   manifest: one .gjf path per line (relative to the project root)
#   As an array job, task N runs line N. Submitted without --array (hx.py submit --serial),
#   it runs every line in turn, for batches of short jobs such as AM1.

#SBATCH --job-name=g16
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --output=slurm_logs/%x.%A_%a.out
#SBATCH --error=slurm_logs/%x.%A_%a.err

set -u
MANIFEST="$1"
cd "$SLURM_SUBMIT_DIR"

run_one() {
    local GJF="$1" TAG="$2"
    if [ -z "$GJF" ] || [ ! -f "$GJF" ]; then
        echo "ERROR: no input for task $TAG in $MANIFEST"; return 1
    fi
    local NAME STAGE_DIR LOG RC
    NAME=$(basename "$GJF" .gjf)
    STAGE_DIR=$(dirname "$(dirname "$GJF")")
    LOG="$STAGE_DIR/logs/$NAME.log"
    mkdir -p "$STAGE_DIR/logs" "$STAGE_DIR/chks"

    # node-local scratch for Gaussian's .rwf/.int files
    export GAUSS_SCRDIR="${L_SCRATCH:-${SCRATCH:-/tmp}}/g16_${SLURM_JOB_ID}_${TAG}"
    mkdir -p "$GAUSS_SCRDIR"

    echo "task=$TAG job=$SLURM_JOB_ID node=$(hostname) input=$GJF start=$(date)"
    g16 < "$GJF" > "$LOG"
    RC=$?
    echo "g16 exit=$RC end=$(date)"

    if [ "${FORMCHK:-0}" = "1" ] && [ -f "$STAGE_DIR/chks/$NAME.chk" ]; then
        formchk "$STAGE_DIR/chks/$NAME.chk" "$STAGE_DIR/chks/$NAME.fchk"
    fi
    rm -rf "$GAUSS_SCRDIR"
    return $RC
}

module load chemistry gaussian
if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
    run_one "$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MANIFEST")" "$SLURM_ARRAY_TASK_ID"
    RC=$?
else
    RC=0; N=0
    while IFS= read -r GJF || [ -n "$GJF" ]; do
        [ -z "$GJF" ] && continue
        N=$((N + 1))
        run_one "$GJF" "$N" || RC=1
    done < "$MANIFEST"
fi
exit $RC
