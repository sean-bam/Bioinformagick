"""
Classes to parse the output of blast
"""

import pandas as pd

class BlastTab(object):
    """
    Parameters
    ==========
    blasttab: str
        Path to the tabular output of blast
    """
    def __init__(self, blasttab):

        self.blasttab = blasttab
        
    def blast_outfmt6_to_df(self):
        """
        
        """
        df = pd.read_csv(self.blasttab,
                         names=['qaccver',
                                'saccver',
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
                         sep = '\t'
                        )
        
        return df
    
class BlastSubprocess(object):
    
    """
    work in progress
    """
    
    def __init__(self, blasttab):
        pass
    
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