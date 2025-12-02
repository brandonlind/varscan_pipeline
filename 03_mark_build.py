"""Create and sbatch picard and samtools command files.

### purpose
# use picard to mark/remove duplicates, build bam index for GATK
###

### usage
# 03_mark_build.py /path/to/sortfile /path/to/pooldir/
###
"""

import sys, os, balance_queue, subprocess, shutil
from os import path as op
from coadaptree import makedir, get_email_info, pklload

thisfile, pooldir, samp = sys.argv
parentdir = op.dirname(pooldir)
bash_variables = op.join(parentdir, 'bash_variables')
sortfiles = pklload(op.join(pooldir, '%s_sortfiles.pkl' % samp))
joined = " -I ".join(sortfiles)

# MarkDuplicates
dupdir = op.join(pooldir, '03_dedup_rg_filtered_indexed_sorted_bamfiles')
tmpdir = op.join(dupdir, 'tmp')
pool = op.basename(pooldir)
dupfile = op.join(dupdir, "%s_rd.bam" % samp)
statfile = dupfile.replace(".bam", "_stats.txt")
dupstat = op.join(dupdir, "%s_rd_dupstat.txt" % samp)

# create sh file
email_text = get_email_info(op.dirname(pooldir), '03')
text = f'''#!/bin/bash
#SBATCH --mem=50G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --job-name={pool}-{samp}-mark
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH -o %x_%j.out
{email_text}

hostname
date

echo MARKDUPS
module load picard/3.1.1
export _JAVA_OPTIONS="-Xms256m -Xmx27g"
java -Djava.io.tmpdir={tmpdir} -jar $PICARD MarkDuplicates  \
-I {joined} \
-O {dupfile} \
--MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 1000 \
-M {dupstat} \
-REMOVE_DUPLICATES true

date

echo BAMINDEX

cd {dupdir}

java -jar $PICARD BuildBamIndex -I {dupfile}
module unload picard

date

echo SAMTOOLS_STATS
module load samtools/1.19.2
samtools stat -@ 8 {dupfile} > {statfile}
module unload samtools

date

source {bash_variables}

python $HOME/pipeline/04_realignTargetCreator.py {pooldir} {samp} {dupfile}

date
'''

# create shdir and file
shdir = op.join(pooldir, 'shfiles/03_mark_build_shfiles')
for d in [shdir, dupdir, tmpdir]:
    makedir(d)
file = op.join(shdir, '%(pool)s-%(samp)s-mark.sh' % locals())
with open(file, 'w') as o:
    o.write("%s" % text)

# sbatch file
os.chdir(shdir)
print('shdir = ', shdir)
subprocess.call([shutil.which('sbatch'), file])
