#!/usr/bin/env python
"""
This program accepts a directory of proteins sequences (one seq per file). 
It sets up a script to run hhblits, hhpred.
It then submits the scripts.
Then aggregates the results.
"""

# Imports --------------------------------------------------------------------------------------------------
import argparse
from pathlib import Path
import subprocess
import random
import pandas as pd
from io import StringIO
import time

#import sys
#from anvio_fastalib import FastaOutput
#from Bio.Seq import Seq

# Args -----------------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser(
        description='')
parser.add_argument('-i',
                    '--input',
                    help="Path/to/dir/",
                    required=True)

parser.add_argument('-o',
                    '--output',
                    help="name of the output directory",
                    required=True)

parser.add_argument('-w',
                    '--work',
                    help="0 = Dry-run ; 1 = Execute. Default = 0",
                    required=False)

parser.add_argument('-nob',
                    '--no_build',
                    action='store_true',
                    help="Don't build MSA with -n iterations through -db",
                    required=False)

parser.add_argument('-s',
                    '--suffix',
                    help="Input files have these suffixes, default = 'faa'",
                    required=False)

parser.add_argument('-msa',
                    '--multiple_seq_align',
                    action='store_true',
                    help="Input is a multiple sequence alignment without consensus",
                    required=False)

parser.add_argument('-xcon',
                    '--consensus',
                    action='store_true',
                    help="Input is a multiple sequence alignment with consensus",
                    required=False)

parser.add_argument('-db',
                    '--database',
                    help="HHsuite-formatted database. Default = /home/benlersm/data/hhsuite/UniRef30_2020_02",
                    required=False)

parser.add_argument('-pdb',
                    '--pdb_database',
                    help="HHsuite-formatted database. Default = /home/benlersm/data/hhsuite/pdb70",
                    required=False)

parser.add_argument('-n',
                    '--iterations',
                    help="Number of iterations. Default = '3'",
                    required=False)

parser.add_argument('-p',
                    '--parse',
                    help="Try to parse the output into a CSV table",
                    required=False)

#parser.add_argument('-l',
#                    '--list',
#                    help="List of accessions to search",
#                    required=False)

args = parser.parse_args()

# Constants ------------------------------------------------------------------------------------------------
if args.database:
    hhblits_uniprot_db = args.database + '*'
else:
    hhblits_uniprot_db = Path('/home/benlersm/data/hhsuite/UniRef30_2020_02')
    
if args.pdb_database:
    hhblits_pdb_db = args.pdb_database + '*'
else:
    hhblits_pdb_db = Path('/home/benlersm/data/hhsuite/pdb70')
    
if args.iterations:
    iterations = args.iterations
else:
    iterations = 3

if args.suffix:
    suffix = args.suffix
else:
    suffix = 'faa'

if args.work:
    work = int(args.work)
else:
    work = 0
assert(work == 0 or work == 1), f"Work must be '0' or '1' and you put {args.work}"

hhblits_CDD_db = Path('/home/benlersm/data/hhsuite/CDD')

if args.no_build:
    assert(args.multiple_seq_align), f"""
    You chose not to build an MSA (option --no_build) but did not specify
    that the input files in {args.input} are MSAs already. 
    Add option -msa to your command
    """
    
if args.consensus:
    assert(args.multiple_seq_align), f"""
    You said the input files in {args.input} have consensus sequences, 
    but did not specify they are MSAs.
    Add option -msa to your command
    """
# Functions ------------------------------------------------------------------------------------------------

def directory_to_chunks(directory, extension, chunk):
    """
    Returns a list containing chunk-sized lists of files /
    of the given extension in directory
    
    E.g., directory_to_chunks('/path/to/dir', 'txt', 20)
    """
    assert Path(directory).is_dir(), f'{directory} isnt a directory'
    files = list(Path(directory).glob('*.' + extension))
    assert len(files) > 0, f"Couldnt find any *.{extension} files in {directory}"
    
    file_list = []
    for i in range(0, len(files), chunk):
        file_list.append((files[i:i + chunk]))
    return file_list

def make_submit_script(file):
    """
    Accepts a writable object, prints the commands to configure a batch script
    Loads HHsuite, copies the databases to lscratch
    """
    
    print("#! /bin/bash", file = f)
    #print("#SBATCH --partition=multinode", file = f)
    #print("#SBATCH --constraint=x2695", file = f)
    #print("#SBATCH --ntasks=64", file = f)
    #print("#SBATCH --ntasks-per-core=1", file = f)
    #print("#SBATCH --exclusive", file = f)
    #print("#SBATCH --qos=turbo", file = f)
    print("#SBATCH --cpus-per-task=32", file = f)
    print("#SBATCH --mem=20g", file = f)
    print("#SBATCH --gres=lscratch:300", file = f)
    print("#SBATCH --time=06:00:00", file = f)
    print("module load hhsuite || exit 1", file = f)
    if not args.no_build:
        print(f"cp {hhblits_uniprot_db}* /lscratch/$SLURM_JOB_ID", file = f)
    print(f"cp {hhblits_pdb_db}* /lscratch/$SLURM_JOB_ID", file = f)
    print(f"cp {hhblits_CDD_db}* /lscratch/$SLURM_JOB_ID", file = f)
    
def check_jobs(jobs, sleep):
    """
    accepts a CSV-list of jobs. Checks status every SLEEP seconds,
    Exits when no jobs are running/pending
    """
    job_str = ",".join(jobs)
    running_jobs = 1
    pending_jobs = 1
    while int(running_jobs) > 0 or int(pending_jobs) > 0:
        p2 = subprocess.run(f'jobhist {job_str}',
                        check = True, 
                        shell = True,
                        universal_newlines = True,
                        stdout=subprocess.PIPE)
        pending_jobs = p2.stdout.count("PENDING")
        running_jobs = p2.stdout.count("RUNNING")
        completed_jobs = p2.stdout.count("COMPLETED")
        failed_jobs = p2.stdout.count("FAILED")

        print(f'{pending_jobs} pending jobs, {running_jobs} running, {completed_jobs} completed, {failed_jobs} failed')

        time.sleep(sleep)
    print(f"""
    Done. {len(jobs)} were submitted. {pending_jobs} are pending, {running_jobs} are running,
    {completed_jobs} completed successfully and {failed_jobs} failed
    """)
    
    
def parse_hhr_output(hhrfile, num_hits=10):
    """
    This function parses the horrible output of HHblits/HHpred (hhr format)
    
    It returns a pandas dataframe of the top 10 hits, or the number defined by num_hits
    """
    
    #Define the space-delimited field widths for Pandas
    colspecs = [(0, 3),  #RANK
            (4, 12),  #ACCESSION
            (12, 35), #DESCRIPTION
            (36, 40), #PROB
            (41, 48), #EVAL
           (50, 55), #PVAL
            (56, 64), #SCORE
            (64, 68), #SS
            (69, 75), #cOLS
            (76, 85), #QUERY HMM
            (86, 100)] #TEMPLATE HMM
    
    #Define the column names
    names = ['hitrank',
             'accession', 
             'desc',
             'prob',
            'eval',
            'pval',
            'score',
            'ss',
            'cols',
            'query_hmm',
            'template_hmm']
    
    df_list = []
    with open(hhrfile) as f:
        lines = f.readlines()
        query = lines[0].strip().split()[1]
        
        #lines 9-19 contain the table of top 10 hits
        for i in range(9, 19):
            line = lines[i]
            df = pd.read_fwf(StringIO(line), 
                                 headers = None, 
                                 colspecs = colspecs,
                                 names = names
                                )
            df_list.append(df)
    df4 = pd.concat(df_list)
    
    #split the accession column into two
    df4[["accession", "chain"]] = df4["accession"].str.split("_", expand = True)
    
    #set a column containing the query name
    df4["query"] = query
    
    #reorder the columns
    df5 = df4[["query", 
               "hitrank", 
               "accession", 
               "chain", 
               "desc", 
               "prob",
               "eval",
               "pval",
               "score",
               "ss",
               "cols",
               "query_hmm",
               "template_hmm"]]
    
    #select only X number of hts
    df6 = df5.query('hitrank <= @num_hits')
    
    return df6


def add_pdb_metadata_to_hhr_df(hhr_df):
    """
    This function adds PDB metadata to an HHR-formatted dataframe,
    returning another dataframe.
    
    The PDB metadata file "entries.idx" was downloaded on 05/28/20 from here:
    ftp://ftp.wwpdb.org/pub/pdb/derived_data/index/entries.idx
    
    This page pointed me to it:
    https://www.rcsb.org/pages/general/summaries
    
    I removed the first rows containing the column names so that Pandas is happy.
    """
    metadata_df = pd.read_csv('/home/benlersm/data/hhsuite/entries.idx', 
            sep = '\t', 
            header = None,  
            #index_col = 0,
            names = ["IDCODE", 
                     "HEADER", 
                     "DATE",
                     "COMPOUND",
                     "SOURCE",
                     "AUTHORLIST",
                     "RESOLUTION", 
                     "EXPERIMENT"
                    ]
           )

    merge = pd.merge(hhr_df, 
                     metadata_df, 
                     how = 'left', 
                     left_on = 'accession', 
                     right_on = 'IDCODE'
                    )
    
    merge2 = merge.drop(columns = ['IDCODE', 
                                   'HEADER', 
                                   'DATE', 
                                   'AUTHORLIST', 
                                   'RESOLUTION', 
                                   'EXPERIMENT'])
    return merge2

def combine_hhr_results(file_list):
    df_list = []
    for f in file_list:
        df = parse_hhr_output(f, 1)
        df_list.append(df)
    df2 = pd.concat(df_list)
    df3 = add_pdb_metadata_to_hhr_df(df2)

    return df3

#-----------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    #Make a list of n-sized file lists. This will be the number of queries per machine.
    l1 = directory_to_chunks(args.input, suffix, 10)

    #Make numbers for the bash scripts
    x = random.randint(1, 1000)
    i = 1
    
    #set paths
    uniprot_db_local = '/lscratch/$SLURM_JOB_ID/' + hhblits_uniprot_db.stem
    pdb_db_local = '/lscratch/$SLURM_JOB_ID/' + hhblits_pdb_db.stem
    cdd_db_local =  '/lscratch/$SLURM_JOB_ID/' + hhblits_CDD_db.stem
    
    #Set CPUs to the num_threads given in SBATCH script
    cpu = "$SLURM_CPUS_PER_TASK"
    #cpu = "$SLURM_NTASKS"

    #Make the hhblits bash scripts
    for file_list in l1:
        with open('submit_hhblits_' + str(x) + '_' + str(i) + '.sh', 'a') as f:
            make_submit_script(f)
            for file in file_list:
                a3m = file.with_suffix('.a3m')
                hhpred = file.with_suffix('.hhpred')
                
                #If we want to build an MSA vs db
                if not args.no_build: 
                    #And the input is an MSA
                    if args.multiple_seq_align:
                        M_param = '49'
                        #Does the MSA have a consensus?
                        if args.consensus:
                            M_param = 'first'
                        print(f'hhblits -i {file} -M {M_param} -d {uniprot_db_local} -o /dev/null -oa3m {a3m} -n {iterations} -norealign -cpu {cpu}', file = f)
                        
                    #The input is not an MSA    
                    else:
                    #To do: assert single fasta seq as input
                    #if not a3m.exists():
                        print(f'hhblits -i {file} -d {uniprot_db_local} -o /dev/null -oa3m {a3m} -n {iterations} -norealign -cpu {cpu}', file = f)
                
                #Run HHPred
                if not hhpred.exists():
                    if args.multiple_seq_align:
                        M_param = '49'
                        if args.consensus:
                            M_param = 'first'
                        print(f'hhblits -i {file} -M {M_param} -d {pdb_db_local} -d {cdd_db_local} -o {hhpred} -n 1 -norealign -cpu {cpu}', file = f)
                    else:
                        print(f'hhblits -i {a3m} -d {pdb_db_local} -d {cdd_db_local} -o {hhpred} -n 1 -norealign -cpu {cpu}', file = f)

    #Make sure each bash script has an hhblits command, otherwise remove.
    p = Path('.')
    bash_script_list = list(p.glob('submit_hhblits_' + str(x) + '*' + '.sh'))
    for script in bash_script_list:
        if "hhblits" not in str(open(script).readlines()):
            script.unlink()
            
    #Check there is something to do
    bash_script_list2 = list(p.glob('submit_hhblits_' + str(x) + '*' + '.sh'))
    assert len(bash_script_list2) > 0, f"""
    The input directory {args.input} has  *.hhpred outputs for every file it found.
    That means there is nothing to do. Good news? If not, delete *.hhpred files to 
    force a re-run
    """
    
    #Submit the scripts. #0 is dry-run; 1 is execute
    jobids = []
    if args.work == 1:
        for script in bash_script_list2:
            #Submit the jobs, recording jobIDs
            p1 = subprocess.run(f"sbatch {script}", 
                               shell = True, 
                               check = True,
                               universal_newlines = True,
                               stdout=subprocess.PIPE
                               )

            jobids.append(p1.stdout.strip())

            time.sleep(1)

            check_jobs(jobids, 300)
            
            #Cleanup the job submission scripts
            for script in bash_script_list:
                script.unlink()
        
    #Parse the output
    if args.parse:
        file_list = list(Path(args.input).glob('*.hhpred'))
        df = combine_hhr_results(file_list)
        df.to_csv(args.output, index = False)
    
    
    
