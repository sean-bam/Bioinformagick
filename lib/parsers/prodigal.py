"""
Class for parsing prodigal outputs
"""
import pandas as pd

class ProdigalGFF(object):
    """
    Parameters
    ==========
    prodigal_gff: str
        Path to the gff output of prodigal
    """
    def __init__(self, prodigal_gff):

        self.prodigal_gff = prodigal_gff
    
    def prodigal_gff2df(self):
        """
        This function takes the GFF output from prodigal and returns a pandas dataframe
        """

        data = []
        with open(self.prodigal_gff) as f:
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
                    data.append([contig, start, stop, strand])

        df = pd.DataFrame(data, 
                          columns=['contig', 
                                   'start', 
                                   'stop', 
                                   'strand'
                                    ]
                         )

        df2 = df.astype(
                        {'contig' : 'str', 
                         'start' : 'int64', 
                         'stop' : 'int64',
                         'strand' : 'str'
                        }
                       )

        return df2

class ProdigalFasta(object):
    
    def __init__(self, prodigal_gff):

        self.prodigal_gff = prodigal_fasta
    
    def prodigal_fasta2df(self):
    """
    This function takes the fasta output from prodigal and returns a pandas dataframe
    """
        data = []
        with open(self.prodigal_fasta) as fastafile:
                for line in fastafile:
                    if line.startswith('>'):
                        newline = line.strip('>')
                        protein_id, start, stop, strand, info = newline.split(' # ')
                        contig, orf = protein_id.split('_')
                        
                        #if int(strand) == 1:
                        #    sign = '+'
                        #else:
                        #    sign = '-'

                        data.append([contig,
                                     orf,
                                     protein_id,
                                     start,
                                     stop,
                                     strand,
                                     info
                                    ]
                                   )

                df = pd.DataFrame(data, 
                                  columns=['contig', 
                                           'orf', 
                                           'protein_id', 
                                           'start', 
                                           'stop', 
                                           'strand', 
                                           'info'
                                          ]
                                 )
                
                df2 = df.astype(
                        {'contig' : 'str', 
                         'start' : 'int64', 
                         'stop' : 'int64',
                         'strand' : 'str'
                        }
                       )

                return df
            
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
