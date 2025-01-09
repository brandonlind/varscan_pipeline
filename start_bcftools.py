"""Create and sbatch bcftools command files if all realigned bamfiles have been created.

If a single sample per pool_name:
    - set min freq to 0
If multiple samps per pool_name:
    - set min freq to 1/(total_ploidy_across_samples_in_a_pool)

# usage
# python start_bcftools.py parentdir pool
#

# fix
# check_seff does not handle slurm timeouts
"""


import sys, os, time, random, subprocess, shutil
from os import path as op
from datetime import datetime as dt
from coadaptree import makedir, fs, pklload, get_email_info
from balance_queue import getsq


def gettimestamp(f):
    """Get last time modified."""
    return time.ctime(op.getmtime(f))


def getmostrecent(files):
    """Determine the most recent file from a list of files."""
    if not isinstance(files, list):
        files = [files]
    if len(files) > 1:
        whichout = files[0]
        dt1 = dt.strptime(gettimestamp(whichout)[4:], "%b %d %H:%M:%S %Y")
        for o in files[1:]:
            dt2 = dt.strptime(gettimestamp(o)[4:], "%b %d %H:%M:%S %Y")
            if dt2 > dt1:
                whichout = o
                dt1 = dt2
        return whichout
    elif len(files) == 1:
        return files[0]
    else:
        # if len(files) == 0
        return None


def getfiles(samps, shdir, grep):
    """Determine if all realign bam jobs have been created and sbatched.

    Positional arguments:
    samps - list of sample names (length is number of expected shfiles)
    shdir - directory where .sh and .out files are
    grep - program name - keyword used to find correct files

    Returns:
    files - dictionary where key = sh file, val = most recent outfile
    """
    found = [sh for sh in fs(shdir) if sh.endswith(".sh") and grep in sh]
    outs = [out for out in fs(shdir) if out.endswith('.out') and grep in out]
    if len(found) != len(samps):
        print('not all shfiles have been created, exiting %s' % sys.argv[0])
        exit()
    files = dict((f, getmostrecent([out for out in outs if op.basename(f).replace(".sh", "") in out]))
                 for f in found)
    if None in files.values():
        print('not all shfiles have been sbatched, exiting %s' % sys.argv[0])
        exit()
    return files


def check_seff(outs):
    """Execute slurm seff command on each outfile's slurm_job_id to ensure it ran without error.
    Exit otherwise.
    """
    print('checking seff')
    jobid = os.environ['SLURM_JOB_ID']
    for i,f in enumerate(outs):
        pid = f.split("_")[-1].replace(".out", "")
        if not pid == jobid:
            seff, seffcount = '', 0
            while isinstance(seff, list) is False:
                # sometimes slurm sucks
                seff = subprocess.check_output([shutil.which('seff'), pid]).decode('utf-8').split('\n')
                if seffcount == 10:
                    print('slurm is screwing something up with seff, exiting %s' % sys.argv[0])
                    exit()
                time.sleep(1)
                seffcount += 1
            state = [x.lower() for x in seff if 'State' in x][0]
            if 'exit code 0' not in state:
                status = 'died' if 'running' not in state else 'is running'
                print('cannot proceed with %s' % sys.argv[0])
                print('job %s (%s) for %s' % (status, state, f))
                print('exiting %s' % sys.argv[0])
                exit()
        if (i+1) % 10 == 0:
            print('\t%s/%s' % (i+1, len(outs)))


def checkpids(outs, queue):
    """If any of the other bcftools jobs are pending or running, exit."""
    print('checking pids')
    pids = [q[0] for q in queue]
    jobid = os.environ['SLURM_JOB_ID']
    for out in outs:
        pid = out.split("_")[-1].replace(".out", "")
        if pid in pids and pid != jobid:  # if the job is running, but it's not this job
            print('the following file is still in the queue - exiting %s' % sys.argv[0],
                  '\n', '\t%(out)s' % locals())
            exit()


def check_queue(outs, pooldir):
    """Get jobs from the queue, except those that are closing (assumes jobs haven't failed)."""
    print('checking queue')
    sq = getsq(grepping=['bedfile', op.basename(pooldir)])
    if len(sq) > 0:
        checkpids(outs, sq)
    # no need for an else statement here, if len(sq) == 0: no need to check the pids


def get_bamfiles(samps, pooldir):
    """Using a list of sample names, find the realigned bamfiels.

    Return:
    files - dictionary with key = samp_name, val = /path/to/bamfile
    """
    print('getting bamfiles')
    found = fs(op.join(pooldir, '04_realign'))
    files = dict((samp, f.replace(".bai", ".bam"))
                 for samp in samps
                 for f in found
#                  if samp in f 
                 if op.basename(f).startswith("%s_" % samp)
                 and f.endswith('.bai'))
    if not len(files) == len(samps):
        print('len(files) != len(samps)')
        print('files = ', files)
        print('\nsamps = ', samps)
        print('\nfound = ', found)
        exit()
    return files


def checkfiles(pooldir):
    """Call get_bamfiles."""
    # get the list of file names
    print('checking files')
    pool = op.basename(pooldir)
    samps = pklload(op.join(op.dirname(pooldir), 'poolsamps.pkl'))[pool]
    shdir = op.join(pooldir, 'shfiles/05_indelRealign_shfiles')
    files = getfiles(samps, shdir, 'indelRealign')
    check_queue(files.values(), pooldir)  # make sure job isn't in the queue (running or pending)
    #check_seff(files.values())  # make sure the jobs didn't die
    return get_bamfiles(samps, pooldir)


def create_reservation(pooldir, exitneeded=False):
    """Create a file so that other realign jobs can't start bcftools too."""
    print('creating reservation')
    shdir = makedir(op.join(pooldir, 'shfiles/bcftools'))
    file = op.join(shdir, '%s_bcftools_reservation.sh' % pool)
    jobid = os.environ['SLURM_JOB_ID']
    if not op.exists(file):
        with open(file, 'w') as o:
            o.write("%s" % jobid)
    else:
        exitneeded = True
    time.sleep(random.random()*15)
    with open(file, 'r') as o:
        fjobid = o.read().split()[0]
    if not fjobid == jobid or exitneeded is True:
        # just in case two jobs try at nearly the same time
        print('\tanother job has already created bcftools_reservation.sh for %s' % pool)
        exit()
    return shdir


def get_prereqs(bedfile, parentdir, pool, program):
    """Get object names."""
    num = bedfile.split("_")[-1].split(".bed")[0]
    ref = pklload(op.join(parentdir, 'poolref.pkl'))[pool]
    pooldir = op.join(parentdir, pool)
    outdir = makedir(op.join(pooldir, program))
    vcf = op.join(outdir, f'{pool}_{program}_bedfile_{num}.vcf')
    return (num, ref, vcf)


def get_small_bam_cmds(bamfiles, bednum, bedfile):
    """Get samtools commands to reduce a bamfile to intervals in the bedfile."""
    smallbams = []
    #cmds = '''module load java\nmodule load samtools/1.9\n'''
    cmds = []
    for bam in bamfiles:
        pool = op.basename(bam).split("_realigned")[0]
        smallbam = f'$SLURM_TMPDIR/{pool}_realigned_{bednum}.bam'
        cmd = f'''samtools view -b -L {bedfile} {bam} > {smallbam}'''
        # cmd = f'''samtools view -b -L {bedfile} {bam} > {smallbam}\n'''
        # cmds = cmds + cmd
        cmds.append(cmd)
        smallbams.append(smallbam)

    pool = op.basename(bam).split('_')[0]
    cmd_dir = makedir(f'{parentdir}/{pool}/cmd_files')
    cmd_file = f'{cmd_dir}/{bednum}_view_cmds.sh'
    with open(cmd_file, 'w') as o:
        o.write('\n'.join(cmds))
                
    # return (smallbams, cmds)
    return (smallbams, cmd_file)

# export BCFTOOLS_PLUGINS='/home/lindb/src/bcftools-1.11/plugins'
def get_bcftools_cmd(bamfiles, bedfile, bednum, vcf, ref, pooldir, program):
    # smallbams, smallcmds = get_small_bam_cmds(bamfiles, bednum, bedfile)
    smallbams, cmd_file = get_small_bam_cmds(bamfiles, bednum, bedfile)
    smallbams = ' '.join(smallbams)
    sampfile = op.join(pooldir, 'samples_file.txt')
    # determine MAF
    mafpkl = op.join(op.dirname(pooldir), 'maf.pkl')
    if op.exists(mafpkl):
        maf = float(pklload(mafpkl))
    else:
        maf = 0.0

    #masked_vcf = op.basename(vcf).replace('.vcf', '_masked.vcf')
    #/home/lindb/src/bcftools-1.11/bcftools +setGT $SLURM_TMPDIR/{op.basename(vcf)} -o $SLURM_TMPDIR/{masked_vcf} -- -t q -i 'FORMAT/DP<5 || FORMAT/GQ<20' -n .
    #/home/lindb/src/bcftools-1.11/bcftools filter -e 'F_MISSING > 0.20 || MAF <= 0' $SLURM_TMPDIR/{masked_vcf} > {vcf}
    
    cmd = f'''module load parallel
module load samtools/1.9
cat {cmd_file} | parallel -j {threads} --progress --eta
module unload samtools/1.9

/home/lindb/src/bcftools-1.11/bcftools mpileup --min-MQ 20 --min-BQ 20 -B -f {ref} {smallbams} -a "DP,AD" | \
/home/lindb/src/bcftools-1.11/bcftools call -G - -Ov -mv -f GQ,GP --samples-file {sampfile} > $SLURM_TMPDIR/{op.basename(vcf)}
/home/lindb/src/bcftools-1.11/bcftools filter -e 'F_MISSING > 0.40 || MAF <= 0' $SLURM_TMPDIR/{op.basename(vcf)} > {vcf}
'''
    # final vcf
    outdir = makedir(op.join(pooldir, program))
    finalvcf = op.join(outdir, op.basename(vcf))  # TODO: I think this is redundant, leaving since it's worked before

    return (cmd, finalvcf)


def make_adaptree_sh(bamfiles, bedfile, shdir, pool, pooldir, program, parentdir):
    """Create sh file for bcftools command."""
    bednum, ref, vcf = get_prereqs(bedfile, parentdir, pool, program)
    
    cmd, finalvcf = get_bcftools_cmd(bamfiles, bedfile, bednum, vcf, ref, pooldir, program)

    outtable = finalvcf.replace('.vcf', '.txt')
    filt_outdir = op.dirname(finalvcf)

    bash_variables = op.join(parentdir, 'bash_variables')
    text = f'''#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={threads}
#SBATCH --job-name={pool}-{program}_bedfile_{bednum}
#SBATCH --time='1-00:00:00'
#SBATCH --mem=4000M
#SBATCH --output={pool}-{program}_bedfile_{bednum}_%j.out

# run bcftools (v1.11), filter
module load StdEnv/2018.3
{cmd}

# gzip outfiles to save space
module load nixpkgs/16.09  gcc/7.3.0 htslib/1.9
cd $(dirname {finalvcf})
bgzip -f {finalvcf} --threads {threads}

# index and convert to table
module load StdEnv/2023
module load gatk/4.4.0.0
gatk IndexFeatureFile --input {finalvcf}.gz
gatk VariantsToTable --variant {finalvcf}.gz -F CHROM -F POS -F REF -F ALT -F AF -F QUAL -F TYPE -F FILTER -F ADP -F WT -F HET -F HOM -F NC \
-GF GT -GF GQ -GF SDP -GF DP -GF PL -GF PVAL -GF AD -GF RD -GF GP\
-O {outtable}
module unload gatk/4.4.0.0

# mask genotypes with low genotype quality or depth, refilter loci for <=40% missing data
source $HOME/activate_py3124.sh
python $HOME/pipeline/filter_bcftools.py {outtable} {filt_outdir} {threads}

# if any other bcftools jobs are hanging due to priority, change the account
source {bash_variables}
python $HOME/pipeline/balance_queue.py {program} {parentdir}

'''
    file = op.join(shdir, f'{pool}-{program}_bedfile_{bednum}.sh')
    with open(file, 'w') as o:
        o.write("%s" % text)
    return file, finalvcf

def sbatch(file):
    """Sbatch file."""
    os.chdir(op.dirname(file))
    pid = subprocess.check_output([shutil.which('sbatch'), file]).decode('utf-8').replace("\n", "").split()[-1]
    print("sbatched %s" % file)
    #time.sleep(10)
    return pid


def get_bedfiles(parentdir, pool):
    """Get a list of paths to all of the bed files for ref.fa."""
    ref = pklload(op.join(parentdir, 'poolref.pkl'))[pool]
    beddir = op.join(op.dirname(ref), 'bedfiles_%s' % op.basename(ref).split(".fa")[0])
    return [f for f in fs(beddir) if f.endswith('.bed')]

def create_sh(bamfiles, shdir, pool, pooldir, program, parentdir):
    """Create and sbatch shfiles, record pid to use as dependency for combine job."""
    bedfiles = get_bedfiles(parentdir, pool)
    pids = []
    finalvcfs = []
    for bedfile in bedfiles:
#         file = make_sh(bamfiles, bedfile, shdir, pool, pooldir, program, parentdir)
        file,finalvcf = make_adaptree_sh(bamfiles, bedfile, shdir, pool, pooldir, program, parentdir)
        pids.append(sbatch(file))
        finalvcfs.append(finalvcf + '.gz')
    return pids, finalvcfs

def create_combine(pids, parentdir, pool, program, shdir, finalvcfs):
    """Create command file to combine bcftools jobs once they're finished.

    Positional arguments:
    pids = list of slurm job id dependencies (the jobs that need to finish first)
    ...
    """
    pooldir = op.join(parentdir, pool)
    email_text = get_email_info(parentdir, 'final')
    catout = op.join(op.dirname(finalvcfs[0]), f'{pool}-{program}_all_bedfiles.vcf.gz')
    dependencies = "#SBATCH --dependency=afterok:" + ",".join(pids)
    bash_variables = op.join(parentdir, 'bash_variables')
    joined = " ".join(finalvcfs)
    text = f'''#!/bin/bash
#SBATCH --job-name={pool}-combine-{program}
#SBATCH --time=12:00:00
#SBATCH --mem=20000M
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={threads}
#SBATCH --output={pool}-combine-{program}_%j.out
{dependencies}
{email_text}


# source {bash_variables}

# python $HOME/pipeline/combine_varscan.py {pooldir} {program} {pool}


/home/lindb/src/bcftools-1.11/bcftools concat {joined} -O z -o {catout} --threads {threads}
    
'''
    combfile = op.join(shdir, f'{pool}-combine-{program}.sh')
    with open(combfile, 'w') as o:
        o.write("%s" % text)
    sbatch(combfile)
    print(f'sbatched {program} combinefile with dependencies: ' + ','.join(pids))

    pass

def main(parentdir, pool):
    """Start <program> if it's appropriate to do so."""

    # check to see if all bam files have been created; if not: exit()
    bamfiles = checkfiles(op.join(parentdir, pool))

    # create reservation so other files don't try and write files.sh, exit() if needed
    shdir = create_reservation(op.join(parentdir, pool))

    # create .sh files
    for program in ['bcftools']:
        print('starting %s commands' % program)
        # create .sh file and submit to scheduler
        pids,finalvcfs = create_sh(bamfiles.values(),
                                   shdir,
                                   pool,
                                   op.join(parentdir, pool),
                                   program,
                                   parentdir)

        # create .sh file to combine bcftools parallels using jobIDs as dependencies
        #create_combine(pids, parentdir, pool, program, shdir, finalvcfs)


if __name__ == "__main__":
    # args
    thisfile, parentdir, pool = sys.argv

    threads = 48

    main(parentdir, pool)
