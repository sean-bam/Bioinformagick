"""
Classes to parse HHsuite outputs
"""

class HMMerTbl(object):
    """
    Parameters
    ==========
    hmmsearch_tblout: str
        Path to the HMMer table output of HMMer
    """
    def __init__(self, hmmsearch_tblout):

        self.hmmsearch_tblout = hmmsearch_tblout


    def parse_hmmer_tblout(self.hmmsearch_tblout):
        """
        Converts the --tblout output of hmmsearch into a dataframe
        """

        names = ["target_name",
             "t_accession",
             "query_name",
             "q_accession",
             "e-value",
             "score",
             "bias",
             "dom_E-value",
             "dom_score",
             "dom_bias",
             "exp",
             "reg",
             "clu",
             "ov",
             "env",
             "dom",
             "rep",
             "inc",
             "t_description"]

        df = pd.read_csv(self.hmmsearch_tblout, 
                        delim_whitespace = True, 
                        comment = "#",
                        names = names
                       )
        return df

    
class HMMerDomtbl(object):
    """
    Parameters
    ==========
    hmmsearch_domtblout: str
        Path to the HMMer domain table output of HMMer
    """
    def __init__(self, hmmsearch_domtblout):

        self.hmmsearch_domtblout = hmmsearch_domtblout

    def parse_hmmer_domtblout(self.hmmsearch_domtblout):
        """
        Converts the --tblout output of hmmsearch into a dataframe
        """

        names = ["target_name",
                 "t_accession",
                 "t_len",
                 "query_name",
                 "q_accession",
                 "q_len",
                 "e-value",
                 "score",
                 "bias",
                 "dom_num",
                 "dom_count",
                 "c-e_value",
                 "i-e_value",
                 "dom_score",
                 "dom_bias",
                 "hmm_from",
                 "hmm_to",
                 "ali_from",
                 "ali_to",
                 "env_from",
                 "env_to",
                 "acc",
                 "t_description"]

        df = pd.read_csv(self.hmmsearch_domtblout, 
                        delim_whitespace = True, 
                        comment = "#",
                        names = names
                       )

        #calculate query/subject coverages
        df["t_cov"] = (df["env_to"] - df["env_from"]) / df["t_len"]
        df["q_cov"] = (df["hmm_to"] - df["hmm_from"]) / df["q_len"]

        return df