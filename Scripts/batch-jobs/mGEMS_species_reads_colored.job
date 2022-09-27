#!/bin/bash
#SBATCH --job-name=wgs-msweep-species-reads-colored
#SBATCH --output=/wrk/users/temaklin/wgs-msweep/logs/wgs-msweep-species-reads-colored.out
#SBATCH -e /wrk/users/temaklin/wgs-msweep/logs/wgs-msweep-species-reads-colored.err
#SBATCH -t 07:59:59
#SBATCH -M ukko
#SBATCH -p short
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu=4096
#SBATCH --array=11-1679

set -euxo pipefail

export LC_ALL="C"
export LANG="C"

module load CMake GCC

basedir=$WRKDIR/wgs-msweep

accession=$(sed -n "${SLURM_ARRAY_TASK_ID} p" $basedir/out/wgs_reads_accessions.txt)

nthreads=${SLURM_CPUS_PER_TASK}

out=$basedir/out/mSWEEP/$accession
grouping=$basedir/reference/species_msweep_indicators.txt

in=$basedir/out/mSWEEP/$accession
outdir=$basedir/out/mGEMS/$accession

## Check if $outdir is empty and clear contents if not
if [ -n "$(find $outdir -prune -empty -type d 2>/dev/null)" ];
then ## Create $outdir if it doesn't exist
    if [ -d "$outdir" ];
    then ## $outdir already exists; do nothing
	echo $outdir" exists" > /dev/null
    else
	mkdir $outdir
    fi
else ## empty contents of $outdir
    rm -rf $outdir/*
fi

fwdaln=$basedir/out/themisto/$accession/$accession""_1.aln.xz
revaln=$basedir/out/themisto/$accession/$accession""_2.aln.xz
probs=$in/$accession""_probs.csv.xz
abundances=$in/$accession""_abundances.txt
index=$basedir/reference/colored-index-isolates

$USERAPPL/mGEMS/bin/mGEMS-amd bin --min-abundance 0.000001 --themisto-alns $fwdaln,$revaln -o $outdir --probs $probs -a $abundances --index $index -i $grouping

if [ -n "$(find $outdir -prune -empty -type d 2>/dev/null)" ];
then ## no files found in outdir
    echo "no bins created" > /dev/null
else ## compress the output bins
    xz --extreme --threads $nthreads $outdir/*.bin
fi
