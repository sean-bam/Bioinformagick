def colabfold_logfile_to_table(logfile, output):
    
    with open(output, 'w') as o:
        with open(logfile) as f:
            for line in f:
                if 'Query' in line:
                    jobname = line.strip().split()[4]

                if 'took' in line:
                    model = line.strip().split()[2]

                    time_raw = line.strip().split()[4]
                    time = time_raw.split(".")[0]

                    pLDDT = line.strip().split()[7]
                    print(jobname,model,time,pLDDT, sep = ",", file = o)
                    
def colabfold_logtable_to_df(table):
    
    df = pd.read_csv(table, 
                    header = None, 
                    names = ['profile', 
                             'model', 
                             'time', 
                             'score']
                    )
    return df

def get_top_ranked_models(af_output_dir):
    """
    Searches through a directory containing *.pdb files produced from alphafold/colabfold
    Selects the top ranked file based on the presence of "rank_1" in the filename
    Returns a list
    """
    files = []
    for file in Path(af_output_dir).rglob('*.pdb'):
        if "rank_1" in str(file):
            files.append(file)
    return files