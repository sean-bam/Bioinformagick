"""
Classes to parse HHsuite outputs
"""

class HHsuiteHHR(object):
    """
    Parameters
    ==========
    HHR: str
        Path to the HHR output of HHsuite
    """
    def __init__(self, hhsuite_hhr):

        self.hhsuite_hhr = hhsuite_hhr

    def parse_hhr_output_fixedwidth(self.hhsuite_hhr, num_hits=10):
        """
        !!Deprecated!!
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

        with open(self.hhsuite_hhr, encoding='utf-8', errors='ignore') as f:
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

    def parse_hhr_output(self.hhsuite_hhr):
        """
        This function parses the horrible output of HHblits/HHpred (hhr format)
        Note: some alignments are discarded b/c of inconsistent delimiting in the hhr output

        It returns a pandas dataframe
        """
        #initialize an empty dataframe
        columns = [#'hitrank', 
               'accession', 
               'desc', 
               'prob', 
               'evalue', 
               'pvalue', 
               'score', 
               'ss', 
               'cols', 
               'query_hmm', 
               'template_hmm',
               'template_length'
              ]

        df = pd.DataFrame(columns = columns)

        with open(self.hhsuite_hhr, encoding='utf-8', errors='ignore') as f:
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
                            d = {#'hitrank' : hitrank,
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
                            #df = pd.DataFrame(data = d, index = [0])
                            df = df.append(d, ignore_index = True)

                        #pass over the line if prob2 isn't  a float
                        except ValueError:
                            pass

                #If there are fewer than 500 hits, the subsequent line will be empty
                #and will raise an IndexError
                except IndexError:
                    break

        #I now have a dataframe with hits, or an empty dataframe
        if not df.empty:

            #set a column containing the query name
            df["query"] = query
            df["file"] = Path(self.hhsuite_hhr).stem

            #reorganize the columns, set dtypes
            df = df[["file",
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
                      "template_length"]].astype({"prob" : "float", 
                                                "evalue" : "float",
                                                "cols" : "float"})

        return df

    def get_best_hhpred_hit(self.hhsuite_hhr, min_prob = 90, min_aln_length = 50):
        """
        Filters HHR output for hits > 90 prob and aln_length > 50
        Selects the best non-PDB hit, if present
        Otherwise, returns best PDB hit or empty dataframme
        """
        df = parse_hhr_output(self.hhsuite_hhr)

        #filter the hits
        df2 = df.query('prob >= @min_prob and cols >= @min_aln_length')

        #get the best non PDB hit
        df3 = df2.query('~accession.str.contains("_")').head(1)

        #if only pdb hits
        if df3.empty:
            df3 = df2.head(1)

        #if no hits pass the criteria, this is an empty dataframe
        return df3


    def add_pdb_metadata_to_hhr_df(hhr_df):
        """
        !!Deprecated!!
        
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