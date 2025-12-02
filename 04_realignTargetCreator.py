"""Create and sbatch gatk realignTargetCreator command files.

### purpose
# use the GATK to create target intervals for realignment around indels
###

### usage
# python 04_realignTargetCreator.py /path/to/pooldir/ sampID
###
"""

import os, sys, subprocess, shutil
from os import path as op
from coadaptree import makedir, pklload, get_email_info

thisfile, pooldir, samp, dupfile = sys.argv

# RealignerTargetCreator
aligndir = op.join(pooldir, '04_realign')
listfile = op.join(aligndir, f'{samp}_realingment_targets.list')

# get ref
parentdir = op.dirname(pooldir)
pool = op.basename(pooldir)
ref = pklload(op.join(parentdir, 'poolref.pkl'))[pool]
bash_variables = op.join(parentdir, 'bash_variables')

email_text = get_email_info(parentdir, '04')
text = f'''#!/bin/bash
#SBATCH --mem=50G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --job-name={pool}-{samp}-realign
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH -o %x_%j.out
{email_text}

hostname
date

echo REALIGNER
source $HOME/conda_init.sh
conda activate gatk3_8
export _JAVA_OPTIONS="-Xms256m -Xmx48g"
gatk3 -T RealignerTargetCreator -R {ref} --num_threads 32 -I {dupfile} -o {listfile}
date

# next step
source {bash_variables}
python $HOME/pipeline/05_indelRealign.py {pooldir} {samp} {dupfile} {ref}

date
'''

# create shdir and shfile
shdir = op.join(pooldir, 'shfiles/04_realignTarget_shfiles')
for d in [aligndir, shdir]:
    makedir(d)
file = op.join(shdir, f'{pool}-{samp}-realign.sh')
with open(file, 'w') as o:
    o.write("%s" % text)

# sbatch file
os.chdir(shdir)
print('shdir =', shdir)
subprocess.call([shutil.which('sbatch'), file])
