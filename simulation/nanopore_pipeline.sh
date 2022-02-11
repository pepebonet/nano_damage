#!/bin/bash

# Configuration variables ------------------------------------------------------

# Maximal number of jobs to execute at the same time
MAX_JOBS=5
# Maximal number of jobs per second
MAX_JOBS_PER_SECOND=10
# Seconds to wait before checking for output file
LATENCY_WAIT=3700

# Check preconditions ----------------------------------------------------------

# Enforce existence of lsf_log directory
mkdir -p lsf_log

# Enforce existence of TMPDIR
export TMPDIR=/home/jbonet/tmp
mkdir -p ${TMPDIR}

# Create one log directory per Snakemake run -----------------------------------


test -z "${LSB_JOBID}" && LSB_JOBID=$(date +%Y-%m-%d_%H-%M)
LOGDIR=lsf_log/${LSB_JOBID}
mkdir -p ${LOGDIR}


# Activate bash cmd printing, debug info ---------------------------------------

set -e

# Activate environment
source activate tombo

echo "Run snakemake"
# Kick off Snakemake -----------------------------------------------------------

lsf="sbatch --mem-per-cpu={cluster.h_vmem} -c {cluster.threads} "
lsf="$lsf -o ${PWD}/${LOGDIR}/out_%j.txt -e ${PWD}/${LOGDIR}/err_%j.txt"

snakemake \
    -s ./Snakefile \
    --cluster-config ./config_lsf.yaml \
    --cluster "$lsf" \
    --latency-wait ${LATENCY_WAIT} \
    -j ${MAX_JOBS} \
    -k \
    -p \
    $*

# Print date after finishing, for good measure ---------------------------------
echo "All done. Have a nice day."
