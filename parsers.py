import pandas as pd
import drivers
from pathlib import Path
import subprocess
from Bio.SeqIO.FastaIO import SimpleFastaParser
from io import StringIO


def file_exists(file_path):
    if not file_path:
        raise RuntimeError("No input file is declared...")
    if not Path(file_path).exists():
        raise RuntimeError(f"No such file: {file_path}")
    return True
    
def multifasta_2_singlefasta(fasta, output_dir, word=1, ext='.fna'):
    """
    """
    
    o = Path(output_dir)
    assert o.is_dir() == True, f"{output_dir} is not a directory"
    
    with open(fasta) as f:
        for title, seq in SimpleFastaParser(f):
            index = 1 - word
            filename = title.split()[index]
            outfile = o / Path(filename).with_suffix(ext)
            assert not outfile.exists(), f"An output named {outfile} already exists!"
            with open(outfile, 'w') as output:
                print(f">{title}", file = output)
                print(seq, file = output)

def cogpsicog_nt2df(cog_psicog_output):
    """
   This function parses the output of cog_psicognitor run with the flags `-nt` and `-qpath`
   It removes complex hits (e.g., query_coords = "15526-15711=22225-22650")
   and returns a dataframe
    """
    df = pd.read_csv(cog_psicog_output, 
                     index_col = False, 
                     names = ["query", 
                              "db_name", 
                              "query2", 
                              "len", 
                              "query_coords", 
                              "aln_len", 
                              "subject", 
                              "subject2",
                              "hit_code",
                              "subject_len",
                              "subject_coords"]
                     ).drop(columns = ['query2', 'subject2'])
    
    #Remove columns where the query coordinates are complex
    df2 = df[~df.query_coords.str.contains('=')]
    df3 = df2[~df2.subject_coords.str.contains('=')]

    df3[["qstart", "qend"]] = df3["query_coords"].str.split('-', expand =True)
    df3[["sstart", "ssend"]] = df3["subject_coords"].str.split('-', expand =True)
    df4 = df3.drop(columns = ['query_coords', 'subject_coords'])
    
    return df4

def cogpsicog_to_df(cog_psicog_output):
    """
   This function parses the output of cog_psicognitor run with default parameters
    """
    df = (pd.read_csv(cog_psicog_output, 
                     index_col = False, 
                     names = ["query", 
                              "db_name", 
                              "query2", 
                              "len", 
                              "query_coords", 
                              "aln_len", 
                              "profile", 
                              "subject2",
                              "hit_code"]
                              #"subject_len",
                              #"subject_coords"]
                     )
              .drop(columns = ['query2', 'subject2'])
         )
    
    return df

def import_profile_descriptions(pandas_df, csv):
    """
    !!! This function is behaving and I'm not sure why... 
    
    This function accepts a csv-delimited files with the following columns:
    profile	source	name	description
    
    And joins the table with a pandas dataframe on the column 'profile'
    
    If the dataframe already has the columns source/name/description, the values
    are replaced.
    """
    #Make sure df has a column labelled profile
    assert('profile' in pandas_df.columns.tolist()), f"Could not find a column labelled 'profile' in {df}"
    
    #read in the profile table, check headers
    df_profile = pd.read_csv(csv)
    expected_headers = ['profile', 'source', 'name', 'description']
    for header in expected_headers:
        assert(header in df_profile.columns.tolist()), f"could not find a column labelled {header}"

    #figure out what to do: merge, or update values?
    if 'source' not in pandas_df.columns.tolist() and 'name' not in pandas_df.columns.tolist() and 'description' not in pandas_df.columns.tolist():
        #Combine the two tables with a merge
        df2 = pd.merge(pandas_df, 
                        df_profile, 
                        how = 'left', 
                        on = 'profile',
                       )
    else:
        #update the values
        pandas_df.update(df_profile, overwrite = False)
        df2 = pandas_df.copy()

    return df2

def cogpsicog_filter(df, marker_list):
    """
    Given a cog_psicognitor dataframe, this function grabs rows with the first marker gene in the list.
    Then, the best hit is taken (first by hit_code, then by length)
    Then, it moves onto the next marker in the list
    Finally, it concatenates the results into a single, sorted dataframe
    """
    
    df.astype({'aln_len' : 'int64'})
    
    df_list = []
    
    for marker in marker_list:
        #Sort by nucleotide accession then hit_code (low to high) then prot_length (high to low)
        df_marker = (df.query('profile_superfamily == @marker and aln_len > 99')
                                .sort_values(['query', 
                                            'hit_code', 
                                            'aln_len'
                                            ], 
                                            ascending = [True, 
                                                         True,
                                                         False
                                                        ]
                                           )
                               .drop_duplicates(subset = ['query', 'profile_superfamily'])
                    )
        #list all the dataframes    
        df_list.append(df_marker)
        
        df_out = (pd.concat(df_list).sort_values(['query']))
    
    return df_out

def df2blastdbcmd(dataframe, querycol, qstartcol, qendcol, blastdb, output):
    """  
    The idea is to give this function a dataframe and the column names to run blastdbcmd.
    
    Blastdbcmd does not accept qstart > qend
    
    This is almost working. For some reason, some of the rows still have a qstart > qend.
    """
    
    #Subset the dataframe to three columns
    df = dataframe[[querycol, qstartcol, qendcol]]
    
    #Split the dataframe, where qstart < qend
    df1 = df.where(df.iloc[:,1] < df.iloc[:,2])
    
    #Where qstart is > qend, swap the names of the columns
    df2 = df.where(df.iloc[:,1] > df.iloc[:,2]).rename(columns={qstartcol : qendcol, qendcol : qstartcol})
    
    #Combine the dataframes back together
    
    (pd.concat([df1, df2], sort = False)
       .dropna()
       .to_csv('tmp12345.csv', index = False)
    )
    
    #Overwrite the output file
    if Path(output).is_file():
        Path(output).unlink()
    
    with Path(output).open('w') as outfile:
        
        subprocess.run(['blastdbcmd', 
                        '-entry_batch', 
                        'tmp12345.csv', 
                         '-db', 
                        blastdb, 
                        ], 
                        stdout = outfile,
                        text = True,
                        check = True
                    )

def df2fa2frag(dataframe, querycol, qstartcol, qendcol, fastafile, output):
    """
    This function accepts a dataframe and the names of the columns to use to extract subsequences with `fa2frag`.
    Example with 
    a dataframe (df1) 
    containing a fasta header in the column 'query'
    containing start coordinates in the column `qstart'
    containing end coordinates in the column `qend`
    
    df2fa2frag(df1, 'query', 'qstart', 'qend', 'input.fna', 'output.fna')
    """
    
    df = dataframe[[querycol, qstartcol, qendcol]]
    
    df.to_csv('tmp12345.csv', index = False, header = False)
    
    with open('tmp6789.csv', 'w') as fragslist:
        with open('tmp12345.csv') as infile:
            for line in infile:
                query, start, stop = line.split(sep = ",")
                print(query + ":"+start + ".." + stop, end = '', file = fragslist)
    
    
    with open(output, 'w') as outfile:
        subprocess.run(['fa2frag', 
                        fastafile, 
                        '-f=tmp6789.csv', 
                        '-w=1',  
                        ], 
                        stdout = outfile,
                        text = True,
                        check = True
                    )
        
    Path('tmp12345.csv').unlink()
    Path('tmp6789.csv').unlink()
    

def prodigal2vcontact_g2g(prodigalfasta, output):
    """
    This function accepts a prodigal-formatted fasta file and outputs a 
    three column table for vcontact: protein_id, contig_id, keyword
    """
    
    #Parse the input file to a temporary file
    parse_prodigal(prodigalfasta, 'tmp12345.csv')
    
    #Make sure the output doesn't exist
    if Path(output).is_file():
        Path(output).unlink()
    
    with open(output, 'a') as outfile:
        #Write the expected header
        print("protein_id,contig_id,keywords", file = outfile)
        
        #Read in the temporary file, write the output
        with open('tmp12345.csv') as infile:
            for line in infile:
                contig,orf,protein_id,start,stop,sign,info = line.split(sep = ",")
                print(protein_id, contig, "N/A", sep = ",", file = outfile)
    Path('tmp12345.csv').unlink()


def parse_prodigal(prodigal_fasta, output):
    """
    This function parses the default headers of prodigal and writes a CSV table
    """
    
    #Make sure the output doesn't exist
    if Path(output).is_file():
        Path(output).unlink()
        
    with open(prodigal_fasta) as fastafile:
            for line in fastafile:
                if line.startswith('>'):
                    newline = line.strip('>')
                    protein_id, start, stop, strand, info = newline.split(' # ')
                    try:
                        contig, orf = protein_id.split('_')
                    except ValueError:
                        contig_str1, contig_str2, orf = protein_id.split('_')
                        contig = contig_str1 + "_" + contig_str2
                    if int(strand) == 1:
                        sign = '+'
                    else:
                        sign = '-'
                    with open(output, 'a') as outfile:
                        print(contig,orf,protein_id,start,stop,sign,info, sep=",", file = outfile, end = '')
                        
def prodigal_fasta2df(prodigal_fasta):
    """
    This function takes the fasta output from prodigal and returns a pandas dataframe
    """
    data = []
    with open(prodigal_fasta) as fastafile:
            for line in fastafile:
                if line.startswith('>'):
                    newline = line.strip('>')
                    protein_id, start, stop, strand, info = newline.split(' # ')
                    contig, orf = protein_id.split('_')
                    if int(strand) == 1:
                        sign = '+'
                    else:
                        sign = '-'
                    data.append([contig,orf,protein_id,start,stop,sign,info])
            df = pd.DataFrame(data, columns=['contig', 'orf', 'protein_id', 'start', 'stop', 'sign', 'info'])
            
            return df
        
def prodigal_gff2df(gff):
    """
    This function takes the GFF output from prodigal and returns a pandas dataframe
    """
    
    data = []
    with open(gff) as f:
        for line in f:
            if line.startswith('##gff-version'):
                pass
            elif line.startswith('# Sequence'):
                split_header = line.split(";")
                contig_length_raw = split_header[1]
                contig_length = contig_length_raw.split('=')[1]
            elif line.startswith('# Model'):
                pass
            else:
                contig, prodigal_version, CDS, start, stop, confidence, strand, number, info = line.split()
                data.append([contig, contig_length, start, stop, strand])
                
    df = pd.DataFrame(data, columns=['contig', 
                                     'contig_length', 
                                     'start', 
                                     'stop', 
                                     'strand'
                                    ]
                     )
    
    df2 = df.astype({'contig_length' : 'int64', 
                     'start' : 'int64', 
                     'stop' : 'int64'
                    }
                   )
    
    return df2
                        
def vcontact2cytoscape(vcontact_input, output):
    """
    This function accepts the output by vcontact and formats the input
    """
    df_vc = pd.read_csv('vcontact_tree/mash_reps/vContact_Output/genome_by_genome_overview.csv').drop(columns = ['Unnamed: 0'])
    
    df_phage = pd.read_csv('pevzner_candidate_phage_annotated.csv')
    
    #Count how many times each mash representative sequence occurs
    df_phage["mash_count"] = (df.groupby(by = ['mash_cluster_rep'])['mash_cluster_rep']
                                .transform('count')
                             )
    
    df_phage2 = (df_phage[['db_name', 
                           'mash_cluster_rep', 
                           'Order', 
                           'Family', 
                           'Subfamily', 
                           'Genus',
                           'mash_count'
                          ]
                         ]
                 .rename(columns = {'mash_cluster_rep' : "Genome"})
                 .replace({'db_name' : {'code15' : 'pevzner'}})
                 .drop_duplicates()
                )
    
    
    df_merge = (pd.merge(df_vc, df_phage2, on = ['Genome'], how = 'left')
                   .fillna(value = {'db_name' : 'Refseq'}))
    
    df_merge.to_csv(output)
    
def bamfiles_to_coveragematrix(directory, output):
    
    file_exists(directory)
                           
    file_list = list(Path(directory).glob('*.bam'))
    
    assert len(file_list) > 0, f"There are no .bam files in {directory}"

    df_list = []
    for file in file_list:
        name = file.stem
                           
        subprocess.run(f'samtools idxstats {file} > samtools_tmp',
                       shell = True,
                       check = True)
        df = (pd.read_table('samtools_tmp', 
                           header = None, 
                           names = ['length', name, 'unmapped_reads']
                           ).drop(columns=['unmapped_reads', 'length']))
        
        df_list.append(df)
        
    Path('samtools_tmp').unlink()
   
    coverage_matrix_df = pd.concat(df_list, sort = False, axis = 1)
    
def bamfiles_to_coveragematrix2(directory, output):
    
    file_exists(directory)
                           
    file_list = list(Path(directory).glob('*.bam'))
    
    assert len(file_list) > 0, f"There are no .bam files in {directory}"

    df_list = []
    for file in file_list:
        name = file.stem
                           
        subprocess.run(f'samtools idxstats {file} > samtools_tmp',
                       shell = True,
                       check = True)
        df = (pd.read_table('samtools_tmp', 
                           header = None, 
                           names = ['length', name, 'unmapped_reads']
                           ).drop(columns=['unmapped_reads', 'length']))
        
        df_list.append(df)
        
    Path('samtools_tmp').unlink()
   
    coverage_matrix_df = pd.concat(df_list, sort = False, axis = 1)
    #coverage_matrix_df.to_csv(output, index_label = 'contig')
    
    return coverage_matrix_df

    
    
def taxonomy2tree_rename(protein2taxonomy_output, output):
    """
    This function expects the entrez output from the script "protein2taxonomy3.sh"; i.e, 10 columns of:
    ProteinID,TaxID,ScientificName,Kingdom,Phylum,Class,Order,Family,Subfamily,Genus
    
    And returns a two-column table for the script tree_rename:
    ProteinID    ProteinID@Order_Family_Subfamily_Genus
    """
    df = pd.read_csv(protein2taxonomy_output).dropna()
    df['name'] = df['Accession'] + "@" + df['Family'] + "_" + df['Subfamily'] + "_" + df['Genus']
    df.to_csv(output, sep='\t', columns = ['Accession', 'name'], index = False, header = False)
    
def df2island(df, output, color_dict = {}):
    """
    This function expects a pandas dataframe with the following columns:
    'contig', 'orf', 'protein_id', 'start', 'stop', 'sign', 'name', 'annotation_short', 'category']
    """
        
    #Make coordinates column. Sometimes Start/Stop is an object, sometimes int64. So, set explicitly
    df1 = df.astype({'start' : 'str', 'stop' : 'str'})
    df1['coordinates'] = df1['start'] + '..' + df1['stop']
    
    #make protein_length column. Same problem as above.
    df2 = df1.astype({'start' : 'int64', 'stop' : 'int64'})
    df2['protein_length'] = ((df2['stop'] - df2['start']) / 3).astype('int')
    
    #make name columns
    df2.loc[:,"name_fam"] = df.name
    df2.loc[:,"name_name"] = df.annotation_short
    
    #If I don't have an annotation_short but do have a name, use the name in both columns
    df2["name_name"] = df2.annotation_short.fillna(df2.name)
    
    #make color column
    if color_dict:
        df2["color"] = [color_dict.get(x) for x in df2.category]
        df2["color"] = df2.color.fillna("#A9A9A9")
    else:
        df2["color"] = "#A9A9A9"
    
    #make mystery column
    df2["mystery"] = "X"
    
    #reorder the columns
    df3 = df2.reindex(columns=['contig', 'coordinates', 'sign', 'protein_id', 'mystery', 'orf', 'name_fam', 'protein_length', 'name_name', 'color'])
    
    contig_name = df.iloc[0].contig
    
    #print the file
    df3.to_csv('island.tmp', index = False, header = False, sep = "\t")
    
    #Write the header
    #if not Path(output).is_dir():
        #Path(output).mkdir(parents=True, exist_ok=True)
        
    with open(output, 'w') as outfile:
        outfile.writelines(["===\t" + contig_name + "\tnum1\t\t" + contig_name + "\t\t\t\t\t\n"])
        
        #Add the data
        with open('island.tmp') as infile:
            for line in infile:
                outfile.write(line)
    
    #remove the temporary file
    Path('island.tmp').unlink()
    

                        
def get_cdd_metadata_as_df():
    
    df = pd.read_csv('/panfs/pan1.be-md.ncbi.nlm.nih.gov/phage_hunting/databases/cddid.tbl', 
                names = ['UID', 
                         'profile', 
                         'name', 
                         'description', 
                         'length'],
               sep = "\t")
    
    return df

def blast_outfmt6csv_to_df(outfmt6_csv):
    df = pd.read_csv(outfmt6_csv, 
                    sep = '\t',
                   header = None,
                   names = ['query',
                           'subject',
                           'pident', 
                            'length',
                            'mismatch',
                            'gapopen',
                            'qstart',
                            'qend',
                            'sstart',
                            'send',
                            'evalue',
                            'bitscore'
                           ]
                   )
    
    return df

def blast_outfmt6_to_df(blastoutfmt6):
    df = pd.read_csv(blastoutfmt6,
                 names=['query',
                        'subject',
                        'pident', 
                        'length',
                        'mismatch',
                        'gapopen',
                        'qstart',
                        'qend',
                        'sstart',
                        'send',
                        'evalue',
                        'bitscore'
                        ],
                 header = None,
                 sep = '\t')
    return df

def blast_outfmt6_subprocess_to_df(outfmt6_subprocess):
    """
    This function converts a blast outfmt6 output in a subprocess object
    to a pandas dataframe.
    """
    
    #Conver the subprocess output into a list, one entry for each hit
    blastoutfmt6_list = outfmt6_subprocess.stdout.strip().split('\n')

    #make an empty list of pandas datframes
    df_list = []

    for hit in blastoutfmt6_list:
        #Convert each hit into a list of lists for Pandas
        hit_list = [hit.split()]

        df = pd.DataFrame(hit_list, columns=['query',
                                          'subject',
                                          'pident', 
                                           'length',
                                           'mismatch',
                                           'gapopen',
                                           'qstart',
                                           'qend',
                                           'sstart',
                                           'send',
                                           'evalue',
                                           'bitscore'
                                          ]
                  )

        df_list.append(df)
    
    df2 = pd.concat(df_list)

    return df2

def parse_hhr_output_fixedwidth(hhrfile, num_hits=10):
    """
    This function parses the horrible output of HHblits/HHpred (hhr format)
    
    It returns a pandas dataframe of the top 10 hits, or the number defined by num_hits
    """
    names = ['hitrank',
             'hit', 
             'prob',
             'eval',
             'pval',
             'score',
             'ss',
             'cols',
             'query_hmm',
             'template_hmm',
             'template_length']
    
    colspecs = [(0, 3),  #RANK
            (4, 12),  #ACCESSION
            (12, 35), #DESCRIPTION
            (35, 41), #PROB
            (41, 48), #EVAL
           (50, 55), #PVAL
            (56, 64), #SCORE
            (64, 68), #SS
            (69, 75), #cOLS
            (76, 85), #QUERY HMM
            (86, 100)] #TEMPLATE HMM
    
    with open(hhrfile, encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()
        query = lines[0].strip().split()[1]
        
        #lines 9-509 contain the table of top 500 hits
        with open('tmp.hhr.table', 'w') as o:
            for i in range(9, 509):
                line = lines[i]
                
                #Print the first 500 lines, only if they start with a digit
                try:
                    if line.split()[0].isdigit():
                        print(line, end = '', file = o)
                        
                #If there are fewer than 500 hits, the subsequent line will be empty
                #and will raise an IndexError
                except IndexError:
                    break
    
    df = (pd.read_fwf('tmp.hhr.table', header = None, names = names, colspecs = colspecs,)
            #.dropna()
            .drop(columns = "hitrank")
            .reset_index()
         )

    #Make a new column for the accession
    #df["accession"] = df.hit.str.split(' ', expand = True)[0]
    df["query"] = query
    df["file"] = Path(hhrfile).stem
    
    #reorganize the columns
    df2 = df[["file",
              "query",
              #"accession", 
              "hit", 
              "prob",
              "eval", 
              "pval", 
              "score", 
              "ss", 
              "cols", 
              "query_hmm", 
              "template_hmm", 
              "template_length"]]
    
    #select only X number of hts
    df3 = df2.query('index <= @num_hits')
    
    #remove the tmp file
    #Path('tmp.hhr.table').unlink()
    return df3

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
                        
                        #parse out the hit. accession, and description
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
    
    #Combine the results together
    df2 = (pd.concat(df_list, ignore_index = True)
             #.drop(columns = "hitrank")
             #.reset_index()
             #.drop(columns = 'index')
          )
    
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
    metadata_df = pd.read_csv('/panfs/pan1.be-md.ncbi.nlm.nih.gov/hhsuite_db/entries.idx', 
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

def combine_hhr_results(file_list, num_hits):
    """
    Accepts a list of files in the HHblits `hhr` output format, and a number of top hits to return
    Returns a pandas datafarme
    """
    df_list = []
    for f in file_list:
        df = parse_hhr_output(f, int(num_hits))
        df_list.append(df)
    df2 = pd.concat(df_list, ignore_index = True)
    #df3 = add_pdb_metadata_to_hhr_df(df2)

    return df2

def parse_yuris_run_mmclust_output_to_dict(run_mmclust_output):
    
    """ 
    This function takes a cluster file that has looks like this:
    4    Prot1 Prot2
    15   Prot4 Prot5 Prot6...
    
    The first sequence is assumed to be the representative sequence.
    
    It returns a dictionary with the representative as value and subsequent sequences as keys:
    {"Prot1 : Prot1", 
    "Prot2 : Prot1", 
    Prot4: Prot4, 
    Prot5: Prot4, } 
    etc..
    """
    
    cluster_dict = {}
    with open(run_mmclust_output) as f:
        for line in f:
            size, cluster_str = line.strip().split('\t')
            rep = cluster_str.split()[0]
            for seq in cluster_str.split():
                cluster_dict[seq] = rep                
    return cluster_dict

def parse_yuris_run_mmclust_to_df(run_mmclust_output, cluster_name='cls_'):
    """
    This function takes a cluster file that has looks like this:
    4    Prot1 Prot2
    15   Prot4 Prot5 Prot6...
    
    The first sequence is assumed to be the representative sequence.
    
    It returns a dataframe with the columns
    protein representative cluster_size, cluster_name 
    """
    #Get the results from Yuri's script
    mm_dict = parse_yuris_run_mmclust_output_to_dict(run_mmclust_output)
    
    #convert to a dataframe
    df = (pd.DataFrame.from_dict(mm_dict, 
                            orient = 'index', 
                            columns = ['representative'])
        .reset_index()
        .rename(columns = {'index' : 'protein'}, 
                copy = False)
     )
    
    #Get the size of the cluster
    df["cluster_size"] = df.groupby('representative')['representative'].transform('count')
    
    #Name each cluster
    grouped = df.groupby('representative')
    i=1
    df_list = []
    for name, group in grouped:
        group["cluster_name"] = str(cluster_name) + str(i)
        df_list.append(group)
        i += 1
    df2 = pd.concat(df_list)
    
    return df2

def parse_yuris_prof_align_to_df(prof_align_cls):
    """
    to do
    """
    cluster_dict = {}
    with open(prof_align_cls) as f:
        for line in f:
            size, node, members = line.strip().split(maxsplit = 2)
            for seq in members.split():
                cluster_dict[seq] = node

    #convert to a dataframe
    df = (pd.DataFrame.from_dict(cluster_dict, 
                            orient = 'index', 
                            columns = ['node'])
        .reset_index()
        .rename(columns = {'index' : 'protein'}, 
                copy = False)
     )

    #Get the size of the cluster
    df["node_size"] = df.groupby('node')['protein'].transform('count')
    return df

def parse_yuris_run_psiprofile(run_psiprofile_output):
    profile_name = Path(run_psiprofile_output).stem
    df = pd.read_csv(run_psiprofile_output,
                 names=['query',
                        'subject',
                        'qlen', 
                        'length',
                        'qstart',
                        'qend',
                        'sstart',
                        'send',
                        'evalue',
                        'bitscore'
                        ],
                 header = None,
                 comment = '#',
                 sep = '\t')
    
    #change the query to the MSA input file name
    df.loc[:,'query'] = profile_name
    
    #pick the best hit 
    df2 = (df.sort_values(by = ['subject', 'evalue'], ascending = [True, True])
             .drop_duplicates(subset = ['subject'])
             .reset_index(drop = True)
          )
    
    return df2

def combine_psiprofile_outputs(file_list):
    """
    Accepts a list of files corresponding to yuri's run_psiprofile
    Parses each file, combines them all and returns a dataframe
    """
    df_list = []
    for file in file_list:
        df = parse_yuris_run_psiprofile(file)
        df_list.append(df)
    
    df_psi = (pd.concat(df_list, ignore_index = True)
                .sort_values(by = ['subject', 'evalue'], ascending = [True, True])
             )
    
    return dF_psi

def entrez_protein_to_taxonomy_to_tree_rename(protein2taxonomy_output, output):
    """
    This function expects 10 columns of:
    ProteinID,TaxID,ScientificName,Kingdom,Phylum,Class,Order,Family,Subfamily,Genus
    
    And returns a two-column table for the script tree_rename:
    ProteinID    ProteinID@Order_Family_Subfamily_Genus
    """
    df = pd.read_csv(protein2taxonomy_output).dropna()
    df['name'] = df['Accession'] + "@" + df['Order'] + "_" + df['Family'] + "_" + df['Subfamily'] + "_" + df['Genus']
    
    df.to_csv(output, 
              sep='\t', 
              columns = ['Accession', 'name'], 
              index = False, 
              header = False)
    
def drep_to_table(drep_dir, output, suffix=".fna"):
    """
    This function accepts an output directory generated by dRep
    https://drep.readthedocs.io/en/latest/index.html
    
    It converts the outputs into a three column CSV table of
    genome, cluster, representative
    
    And also removes the file extension specified by suffix
    """
    # Genomes and cluster designations
    cdb_df = pd.read_csv(drep_dir + '/data_tables/Cdb.csv')
    cdb_df.genome = cdb_df.genome.str.replace(suffix, "")

    # Winning genomes
    wdb_df = pd.read_csv(drep_dir + '/data_tables/Wdb.csv')
    wdb_df.genome = wdb_df.genome.str.replace(suffix, "")

    #Merge the two tables, keeping only the three columns we want
    df1 = (pd.merge(cdb_df, 
                    wdb_df, 
                    how = 'left', 
                    left_on = 'secondary_cluster', 
                    right_on = 'cluster')
           .rename(columns = {"genome_y" : "representative",
                             "genome_x" : "contig"}, 
                   copy = False)
           .loc[:,['contig',
                   'cluster',
                   'representative'
                  ]
               ]
          )
    
    #Add a column for the size of the cluster
    df1["cluster_size"] = df1.groupby('cluster')['contig'].transform('count')

    df1.to_csv(output, index = False)
    
def sr_to_df(srfile):
    p = Path(srfile)
    
    df_dict = {}
    df_list = []
    with open(p) as f:
        for line in f:
            header, seq = line.strip().split()
            df_dict[header] = p.stem
            
        df = (pd.DataFrame.from_dict(df_dict, orient = 'index')
                .reset_index()
                .rename(columns = {'index' : 'accession', 0 : 'file'}, copy = False)
             )

    return df

def chunks(lst, n):
    """
    Yield successive n-sized chunks from lst.
    credit: https://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks
    """
    for i in range(0, len(lst), n):
        yield lst[i:i + n]
        

def entrez_2_island(accession_list, output_dir, chr_start=0, chr_stop=100000000):
    """
    Accepts a list of entrez nuccore accessions.
    Runs drivers.entrez_nuccore_to_protein.
    Outputs an island file named 'output_dir/accession.island'
    """
    p = Path(output_dir)
    assert p.is_dir() and p.exists(), f"{output_dir} is not an existing directory"
    
    drivers.entrez_nuccore_to_protein(accession_list, 'tmp_feature_table.ft', chr_start, chr_stop)
    
    df = pd.read_csv('tmp_feature_table.ft', 
                 sep = '\t') 
    #             header = None, 
    #             names = ['contig', 
    #                      'raw_coords', 
    #                      'protein', 
    #                      'annotation', 
    #                      'strand']
    #            )
    assert not df.empty, f"one of the accessions did not return any proteins"

    #Make a column named strand
    df.loc[df.query('Raw_coords.str.contains("complement")').index, "sign"] = "-"
    df.sign.fillna(value = "+", inplace = True)
    
    #fix the raw_coords column
    df.drop(df[df.Raw_coords.str.contains("join")].index, inplace = True)
    df.drop(df[df.Raw_coords.str.contains("|", regex = False)].index, inplace = True)
    df["Raw_coords"] = (df["Raw_coords"].str.replace("complement", "")
                                      .str.replace("(", "")
                                      .str.replace(")", "")
                                      .str.replace(">", "")
                                      .str.replace("<", "")
                                      .str.replace("..", "-", regex = False))


    df[["start", "stop"]] = df["Raw_coords"].str.split("-", expand = True)

    
    #Format for df2island
    #columns needed = 'contig', 'orf', 'protein_id', 'start', 'stop', 'sign', 'name', 'annotation_short', 'category'
    df['contig'] = df['Accession']
    df['orf'] = df["Protein_id"]
    df['annotation_short'] = df["Product"]
    df['name'] = df["Product"]
    df['category'] = ""
    df2 = df.drop(columns = ['Raw_coords', 'Accession'])#[['Accession', 'Protein_id', 'start', 'stop', 'sign', 'Product']]
    
    #Make an island file for each contig
    grouped = df2.groupby('contig')
    for name, group in grouped:
        output = Path(name).with_suffix('.island')
        df2island(group, p / output)
    
    #remove the feature table
    Path('tmp_feature_table.ft').unlink()