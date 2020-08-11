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

parser.add_argument('-b',
                    '--build',
                    action='store_true',
                    help="Build an MSA with -n iterations through -db",
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
                    help="""
                    /Path/to/HHsuite-formatted-database.
                    E.g., /databases/UniRef30_2020_02
                    NOT /databases/UnirRef30_2020_02_a3m.ffdata
                    Default = /home/benlersm/data/hhsuite/UniRef30_2020_02
                    """,
                    required=False)

parser.add_argument('-pdb',
                    '--pdb_database',
                    help="HHsuite-formatted database. Default = /home/benlersm/data/hhsuite/pdb70",
                    required=False)

parser.add_argument('-n',
                    '--iterations',
                    help="Number of iterations. Default = '3'",
                    required=False)

parser.add_argument('-no_ssd',
                    '--no_ssd',
                    action='store_true',
                    help="Do not transfer the HHsuite databases to /lscratch",
                    required=False)

parser.add_argument('-bt',
                    '--batch_size',
                    help="Number of queries to batch together into a single script. Default = '10'",
                    required=False)

parser.add_argument('-t',
                    '--threads',
                    help="Number of threads per machine. Default = '32'",
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

#Paths to HHsuite formatted DBs
hhblits_uniprot_db = Path('/home/benlersm/data/hhsuite/UniRef30_2020_02')
hhblits_pdb_db = Path('/home/benlersm/data/hhsuite/pdb70')
hhblits_cdd_db = Path('/home/benlersm/data/hhsuite/CDD')

#Number of iterations to build an MSA
iterations = 3

#Input file suffixes
suffix = 'faa'

#Effort - 0 is dry-run, 1 is submit
work = 0

#Number of queries to batch together
batch = 10

#Number of threads per machine
threads = 32
# Inputs  ----------------------------------------------------------------------------------------------

if args.database:
    hhblits_uniprot_db = args.database
    #To do: Check the input is /path/to/database/db 
    #Not /path/to/database
    #and check if _a3m.ffdata exists
    
if args.pdb_database:
    hhblits_pdb_db = args.pdb_database
    
if args.iterations:
    iterations = int(args.iterations)
    assert(type(iterations) == int), f"Iterations must be a number, and you put {args.iterations}"

if args.suffix:
    suffix = args.suffix

if args.work:
    work = int(args.work)
    assert(work == 0 or work == 1), f"Work must be '0' or '1' and you put {args.work}"
    
if args.batch_size:
    batch = int(args.batch_size)
    
if args.threads:
    threads = int(args.threads)

assert Path(args.input).is_dir(), f"You provided '{args.input}' as an input, but it isn't a directory :/"

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

def make_submit_script(file, threads, *dbs):
    """
    Accepts a writable object, prints the commands to configure a batch script
    Loads HHsuite, copies the databases to lscratch
    """
    
    print(f"#! /bin/bash", file = f)
    print(f"#SBATCH --cpus-per-task={threads}", file = f)
    print("#SBATCH --mem=20g", file = f)
    if args.no_ssd:
        #don't need lscratch space
        pass
    else:
        #CDD = ~2g, UniRef ~177g, PDB ~60g
        print("#SBATCH --gres=lscratch:300", file = f)
    print("#SBATCH --time=06:00:00", file = f)
    print("module load hhsuite || exit 1", file = f)

    if args.no_ssd:
        #we dont need to transfer the DBs to /lscratch
        pass
    
    #Otherwise, transfer
    else:
        print(f"cp {hhblits_pdb_db}* /lscratch/$SLURM_JOB_ID", file = f)
        print(f"cp {hhblits_cdd_db}* /lscratch/$SLURM_JOB_ID", file = f)
        
        #if we need to build an MSA, transfer the uniprot_db
        if args.build:
            print(f"cp {hhblits_uniprot_db}* /lscratch/$SLURM_JOB_ID", file = f)
    
    #For MPI runs
    #print("#SBATCH --partition=multinode", file = f)
    #print("#SBATCH --constraint=x2695", file = f)
    #print("#SBATCH --ntasks=64", file = f)
    #print("#SBATCH --ntasks-per-core=1", file = f)
    #print("#SBATCH --exclusive", file = f)
    #print("#SBATCH --qos=turbo", file = f)
    
def check_jobs(jobs, sleep):
    """
    accepts a python list of jobs. Checks status every SLEEP seconds,
    Exits when no jobs are either running or pending
    """
    job_str = ",".join(jobs)
    
    #Check every sleep seconds
    time.sleep(sleep)
    
    p1 = subprocess.run(f'jobhist {job_str}',
                        check = True, 
                        shell = True,
                        universal_newlines = True,
                        stdout=subprocess.PIPE)
    pending_jobs = p1.stdout.count("PENDING")
    running_jobs = p1.stdout.count("RUNNING")
    completed_jobs = p1.stdout.count("COMPLETED")
    failed_jobs = p1.stdout.count("FAILED")
        
    print(f'{pending_jobs} pending jobs, {running_jobs} running, {completed_jobs} completed, {failed_jobs} failed')
        
    while int(running_jobs) > 0 or int(pending_jobs) > 0:
        return check_jobs(jobs, sleep)
    
    print(f"""
    Done. {len(jobs)} were submitted. {pending_jobs} are pending, {running_jobs} are running,
    {completed_jobs} completed successfully and {failed_jobs} failed
    """)
    
    
def parse_hhr_output(hhrfile, num_hits=10):
    """
    This function parses the horrible output of HHblits/HHpred (hhr format)
    
    It returns a pandas dataframe of the top 10 hits, or the number defined by num_hits
    """
    
    df_list = []
    with open(hhrfile, encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()
        query = lines[0].strip().split()[1]
        
        #lines 9-509 contain the table of top 500 hits
        for i in range(9, 509):
            line = lines[i]
            
            #Analyze the lines only if they start with a digit
            try:
                if line.split()[0].isdigit():
                    hit_raw, prob, evalue, pvalue, score, ss, cols, q_hmm, t_hmm, tlen = line.rsplit(maxsplit = 9)
                    
                    #Most times template_hmm and tlen are separated by a space;
                    # 92-310 (449)
                    #But if the template_hmm is long, they are not, and it looks like this:
                    # 802-1021(1151)
                    #Fucking hell. If so, skip that line
                    
                    try:
                        #will raise a ValueError if prob doesnt equal a float b/c of the above issue
                        prob2 = float(prob)
                        
                        #parse out the hit and accession
                        hitrank = (hit_raw.split(maxsplit = 2)[0])
                        accession = (hit_raw.split(maxsplit = 2)[1])
                        desc = (hit_raw.split(maxsplit = 2)[2])
                        
                        #make a pandas dataframe
                        d = {'hitrank' : hitrank,
                             'accession' : accession,
                             'desc': desc,
                             'prob' : prob2,
                             'evalue' : evalue,
                             'pvalue' : pvalue,
                             'score' : score,
                             'ss' : ss,
                             'cols' : cols,
                             'query_hmm' : q_hmm,
                             'template_hmm' : t_hmm,
                             'template_length' : tlen}
                        df = pd.DataFrame(data = d, index = [0])
                        
                        #add the df's to a list
                        df_list.append(df)
        
                    #pass over the line if prob2 isn't  a float
                    except ValueError:
                        pass
                    
            #If there are fewer than 500 hits, the subsequent line will be empty
            #and will raise an IndexError
            except IndexError:
                break
    
    #Combine the results together.
    df2 = pd.concat(df_list)
    
    #set a column containing the query name
    df2["query"] = query
    df2["file"] = Path(hhrfile).stem
    
    #reorganize the columns
    df3 = df2[["file",
              "query",
              "accession",
              "desc",
              "prob",
              "evalue", 
              "pvalue", 
              "score", 
              "ss", 
              "cols", 
              "query_hmm", 
              "template_hmm", 
              "template_length"]]
    
    #select only X number of hts
    df4 = df3.query('index <= @num_hits')
    
    return df4


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
    df2 = pd.concat(df_list, ignore_index = True)
    #df3 = add_pdb_metadata_to_hhr_df(df2)

    return df2


def make_hhblits_command(query, output, n, threads, swarmfile, *dbs, output_a3m = True):
    """
    Accepts a python pathlib object as a query, prints the hhblits command to swarmfile
    """
    
    #Unpack the databases
    db_list = [str(x) for x in dbs]
    databases = " -d ".join(db_list)
        
    #Figure out what the input is and set the -M flag appropriately
    if args.multiple_seq_align:
        M_param = '-M 49'
        #Does the MSA have a consensus?
        if args.consensus:
            M_param = '-M first'
        #Is the MSA already in a3m format?
        if query.suffix == ".a3m":
            M_param = ""
    
    #If the query is an a3m MSA but the -msa flag wasn't set
    elif query.suffix == ".a3m":
            M_param = ""
    else:
        M_param = ""
    
    #Set the output msa
    oa3m = ""
    if output_a3m == True:
        a3m = query.with_suffix('.a3m')
        assert a3m.is_file() == False, f"""
        You're trying to build an a3m-formatted MSA, but an a3m-formatted MSA already exists in the directory. Try disabling -build or delete the existing a3m file
        """
        oa3m = f"-oa3m {a3m}"
    
    #Print the command
    print(f'hhblits -i {query} {M_param} -d {databases} -o {output} {oa3m} -n {n} -norealign -cpu {threads}',
          file = swarmfile)


#-----------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    #Make a list of batch-sized file lists. This will be the number of queries per machine.
    l1 = directory_to_chunks(args.input, suffix, batch)
    
    #set paths to HHsuite formatted databases
   
    
    #Set CPUs to the num_threads given in SBATCH script
    cpu = "$SLURM_CPUS_PER_TASK"
    
    #If using MPI
    #cpu = "$SLURM_NTASKS"
    
    #Make numbers for the bash scripts
    x = random.randint(1, 1000)
    i = 1

    #Make the hhblits bash scripts
    for file_list in l1:
        with open('submit_hhblits_' + str(x) + '_' + str(i) + '.sh', 'a') as f:
            make_submit_script(f, threads)
            
            if args.no_ssd:
                #The paths to the HHsuite databases are the same as default
                hhblits_uniprot_db_local = hhblits_uniprot_db
                hhblits_pdb_db_local = hhblits_pdb_db
                hhblits_cdd_db_local =  hhblits_cdd_db
            else:
                #The paths were moved to lscratch, update
                hhblits_uniprot_db_local = Path('/lscratch/$SLURM_JOB_ID/' + hhblits_uniprot_db.stem)
                hhblits_pdb_db_local = Path('/lscratch/$SLURM_JOB_ID/' + hhblits_pdb_db.stem)
                hhblits_cdd_db_local =  Path('/lscratch/$SLURM_JOB_ID/' + hhblits_cdd_db.stem)
            
            for file in file_list:
                
                if args.build:
                    #We need to build an MSA
                    make_hhblits_command(file, 
                                         '/dev/null', 
                                         iterations, 
                                         cpu, 
                                         f, 
                                         hhblits_uniprot_db_local, 
                                         output_a3m = True)
                    
                    #set file equal to the constructed a3m-formatted MSA
                    file = file.with_suffix('.a3m')
                    
                make_hhblits_command(file, 
                                     file.with_suffix('.hhr'), 
                                     1, 
                                     cpu, 
                                     f, 
                                     hhblits_pdb_db_local,
                                     hhblits_cdd_db_local,
                                     output_a3m = False)
                    
        #Increment i so that the multiple scripts get made
        i += 1
    
    #check if an hhblits command is present, otherwise delete
    #for script in bash_script_list:
    #    if "hhblits" not in str(open(script).readlines()):
    #        script.unlink()
    
    
    #Submit the scripts. #0 is dry-run; 1 is execute
    if work == 1:
        
        bash_script_list = list(Path('.').glob('submit_hhblits_' + str(x) + '*' + '.sh'))
        jobids = []
        
        for script in bash_script_list:
            #Submit the jobs, recording jobIDs
            p1 = subprocess.run(f"sbatch {script}", 
                               shell = True, 
                               check = True,
                               universal_newlines = True,
                               stdout = subprocess.PIPE
                               )
            print(f"submitted {script} with job id {p1.stdout.strip()}")
            jobids.append(p1.stdout.strip())

            time.sleep(1)

        check_jobs(jobids, 60)
        
        #Parse the output
        if args.parse:
            file_list = list(Path(args.input).glob('*.hhr'))
            df = combine_hhr_results(file_list)
            df.to_csv(args.output, index = False)
            
        #Cleanup the job submission scripts
        for script in bash_script_list:
            script.unlink()
        for jobid in jobids:
            Path('slurm-' + jobid + '.out').unlink()

    
    
