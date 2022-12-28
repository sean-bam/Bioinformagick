"""
This program implements MSA filtering as described in 
Esterman et al (2021) https://doi.org/10.1093/ve/veab015

In summary:

1. Calculate HH94 score of a seq in an MSA
2. Calculate the Q score of an AA (x) in an MSA column as Qx = (HH94 score) * (blosum62 score for X against all AAs)
 - For each position in the MSA, an aligned AA is given a frequency score that is incremented by the HH94 score of the sequence. If the position is a gap, all AAs are possible, so each AA frequency is incremented by their default frequency in swiss prot*HH94 score. What we are left with is a vector of 20 amino acids and their frequencies for each column, where the frequencies of aligned AAs are typically much higher than their background swiss prot frequency.
3. The consensus AA for that position is the AA with the highest Qx score
4. Calculate the expectation of the score of an MSA column against a randomly selected AA (R) as Qr = (vector of relative frequencies of AAs) * Qscore
5. Calculate the homogeneity of an MSA column as H = (consensus - expectation) / (blosum score - expectation)
6. Keep/discard consensus AA, depending on threshold

In other words:

Basically, the "best" amino acid is the one with the highest score (BLOSUM62 against the alignment column), calculated with sequence weights. Homogeneity prorates this score between the maximum (strictly homogeneous) and the non-pathological minimum (random assortment of amino acids). Consensus amino acid is registered when the homogeneity is above the threshold, otherwise it's set to "X".

Some more background:

All of this involves sequence weighting. Sequence weighting is typically performed in MSAs to downweight highly redundant sequences and upweight diverse sequences. There are tree-based and distance-based weighting schemes. According to Henikoff & Henikoff 1994 paper, for both schemes:

"the weight assigned to a sequence is a measure of the distance between the sequence and a root or generalized sequence. Each distance is based on the entire sequence in question. However, the seq weight are typically applied to PSSMs in which each position is considered independently".

Or, according to https://doi.org/10.1186/s12859-021-04183-8:
"the score of a sequence is the average of the scores of each position of the sequence, the score of a position being 1/rd, with r the number of different characters at the considered alignment column and d the number of times the character of the considered sequence and position appears in the considered alignment column. The idea of this weighting scheme is to give equal weight to all characters observed at one alignment column, dividing this weight equally among those sequences sharing that character at that position."

What they fail to mention from the Henikoff paper is the part about the consensus:
"for every position of the alignment, the PSSM entry for each residue was the sequence- weighted observed frequency of that residue divided by its expected frequency tabulated from SWISS-PROT".
"""

# Imports --------------------------------------------------------------------------------------------------
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
from Bio.Align import substitution_matrices
from Bio.SeqIO.FastaIO import SimpleFastaParser

# Args -----------------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser(
        description='')
parser.add_argument('-i',
                    '--input',
                    help="Path/to/input",
                    required=False)

parser.add_argument('-o',
                    '--output',
                    help="name of the output file",
                    required=False)

#args = parser.parse_args()

# Constants ------------------------------------------------------------------------------------------------

# Functions -----------

def get_tables_for_msa_filtering():
    
    #set a background frequency of AAs
    background_aa_freq = {
                            "A" : 0.07422,
                            "C" : 0.02469,
                            "D" : 0.05363,
                            "E" : 0.05431,
                            "F" : 0.04742,
                            "G" : 0.07415,
                            "H" : 0.02621,
                            "I" : 0.06792,
                            "K" : 0.05816,
                            "L" : 0.09891,
                            "M" : 0.02499,
                            "N" : 0.04465,
                            "P" : 0.03854,
                            "Q" : 0.03426,
                            "R" : 0.05161,
                            "S" : 0.05723,
                            "T" : 0.05089,
                            "V" : 0.07292,
                            "W" : 0.01303,
                            "Y" : 0.03228,
                            "X" : 0.00001,
                         }
        
    #get a blosum62 matrix
    blosum62 = substitution_matrices.load("BLOSUM62")

    #set a expectation of AAs
    #expectation_aa = {}
    #for aa1 in background_aa_freq.keys():
    #    expectation_aa[aa1] = 0
    #    for aa2 in background_aa_freq.keys():
    #        expectation_aa[aa1] = expectation_aa[aa1] + blosum62[aa1,aa2]*background_aa_freq[aa2]
    
    expectation_aa = aa_freq_to_score(background_aa_freq)

    return background_aa_freq, blosum62, expectation_aa

def aa_freq_to_score(col_aa_frq):
    """
    expects a dictionary of AA : frequency
    returns a dictionary of AA : score
    
    score is sum of blosum62 scores of an index AA versus all other AAs
    """
    #set all AA weights as zero
    weights = {}
    for aa in col_aa_frq.keys():
        weights[aa] = 0
        
    #get a blosum62 matrix
    blosum62 = substitution_matrices.load("BLOSUM62")
    
    #assign weights to each AA
    for aa1 in col_aa_frq.keys():
        for aa2 in col_aa_frq.keys():
            weights[aa1] = weights[aa1] + (blosum62[aa1,aa2] * col_aa_frq[aa2])
            #print(f"AA1 is {aa1}, AA2 is {aa2}, frequency of AA2 is {col_aa_frq[aa2]}, blosum62 of {aa1}{aa2} is {blosum62[aa1,aa2]}, updating weight of AA1 to {weights[aa1]}")
    
    return weights


def aln_to_df(aln, msa_format='fasta'):
    """
    #deprecated
    convertsa  biopython alignment object into a pandas dataframe
    """
    
    headers = [rec.id for rec in aln]
    seqs = [list(rec) for rec in aln]
    df = pd.DataFrame.from_records(seqs, index = headers).astype('str')
    return df

def fa_to_df(fa, is_msa=False):
    """
    Converts a fasta-formatted file to a pandas dataframe
    """
    headers = []
    seqs = []
    with open(fa) as in_handle:
        for title, seq in SimpleFastaParser(in_handle):
            headers.append(title)
            seqs.append(list(seq.upper()))

        
    df = pd.DataFrame.from_records(seqs, index = headers)
    
    #check to make sure the MSA is aligned properly
    if is_msa:
        assert df.isna().sum().sum() == 0, f"All sequences must be the same length"
   
    return df

def df_to_fa(df, output):
    """
    converts a pandas dataframe of sequences to a fasta-formatted file
    """
    
    #combine the columns into a single string
    series = df.T.apply("".join)
    
    #write the output
    with open(output, 'w') as o:
        for header, seq in series.items():
            print(f">{header}", seq, sep = '\n', file = o)
            
def get_gap_cols_from_msa(df, grcut=1):
    """
    -grcut = no more than r fraction of gaps
    
    returns a list of columns that are more than r fraction of gaps
    """
    #df = msa_to_df(msa)
    num_seqs, num_columns = df.shape
    EPSILON = 1e-6

    #calculate the minimum number of non-NA values per column
    max_gaps_in_column = round((num_seqs * grcut) + EPSILON, 3)

    #convert gaps to NA, count how many NAs per column
    gap_count_series = df.replace("-", np.nan).isnull().sum()

    #get columns where number of NAs is >= max_gaps
    gap_col_list = gap_count_series.loc[lambda s: s >= max_gaps_in_column].index.tolist()
    #df2 = df.loc[:, good_cols]
    
    ##RO
    
    return gap_col_list

def calc_HH94(df, grswe=0.51):
    """
    - sequences have an initial weight of 1/nseq
    - Calculates the score of a position in an MSA as 1/rd
        r = # of dif characters (excluding gaps) in a column
        d = number of times the character of the considered sequence 
        and position appears in the considered alignment column
    
    - columns in MSA that are >grswe fraction of gap positions are ignored
    
    - the scores of a sequence's aligned positions are added to its weight
    
    Returns a series object containing the sequence name and HH94 weight
    """
    #df = msa_to_df(msa)
    
    #convert gaps to NA so they are not counted in the weights 
    df2 = df.replace("-", np.nan)
    
    #get a list of gap columns that are more than grswe fraction of gaps
    gap_cols = get_gap_cols_from_msa(df, grswe)
    
    #initial sequence weight is just 1 / nseq
    initial_weight = 1 / len(df)

    
    colaa_weights_nest_dict = {}
    for column in df2.columns:
        if column not in gap_cols:
        
            #get a count of the frequency of each AA in the column
            #by default, NAs (gaps) are not counted
            col_series = df2[[column]].value_counts()

            #calculate the AA weights in the column
            #the length of the series is the number of unique AAs
            colaa_weights = round(1 / (col_series * len(col_series)), 3)

            #store the column-specific results in a dict of dict
            colaa_weights_dict = colaa_weights.to_dict()
            colaa_weights_nest_dict[column] = colaa_weights_dict

    #replace the AAs in the MSA with their column-specific weights
    df3 = df2.replace(colaa_weights_nest_dict)
    
    #for each sequence, sum all the values and add the initial wieght
    hh94_series = df3.sum(axis = 1).round(3) + initial_weight
    
    #normalize
    hh94_mean = hh94_series.mean()
    if hh94_mean < 0:
        hh94_mean = 1
    hh94_norm = round(hh94_series / hh94_series.mean(), 4)
    
    #return a dictionary containing the sequence and its normalized weight
    return hh94_norm.to_dict()
    
def get_col_aa_frequencies(column, hh94_scores, background_aa_freq):
    """
    Accepts a dictionary or pandas series containing seq : aligned char
    
    calculates the column=specific AA weights, given HH94 scores of sequences
    
    returns a dictionary of AA : weight
    """
    
    #set all AA frequencies as zero
    aafreq = {}
    for aa in background_aa_freq.keys():
        aafreq[aa] = 0
    
    #set gap score as zero
    ngap = 0
    
    for seq,char in column.items():
        if char == "-":
            
            #increment the background frequencies of all AAs
            for aa, score in aafreq.items():
                aafreq[aa] = aafreq[aa] + background_aa_freq[aa]*hh94_scores[seq]
            
            #increment the gap frequency 
            ngap += hh94_scores[seq]
        else:
            try:
                aafreq[char] = aafreq[char] + hh94_scores[seq]
            except KeyError:
                pass
    
    #add in ngaps to the dict
    aafreq["-"] = ngap
    
    return aafreq

def get_best_aa_and_score(aascores):
    """
    Expects a dictionary of AA : score
    Selects the AA with the highest score
    
    Returns AA and score
    """
    
    #get top AA and score, as a tuple
    bestaa_tuple = max(aascores.items(), key = lambda k : k[1])
    
    return bestaa_tuple[0],bestaa_tuple[1]
    
def calc_homogeneity(aa, score, nseqs, expectaascores, blosum62):
    """
    calculates the homogeneity of a column in an MSA
    homogeneity ranges from 0-1
    """
    
    homogeneity = (score/nseqs - expectaascores[aa]) / (blosum62[aa, aa] - expectaascores[aa])
    return round(homogeneity, 4)

def assign_consensus(bestaa, grcon, homo, hocon, ngaps, nseqs, ambiguousAA="X"):
    """
    
    """
    #Default consensus is gap
    caa = "-"

    if ngaps <= nseqs * grcon:
        if homo > hocon:
            caa = bestaa
        else:
            caa = ambiguousAA

    return caa

def filter_msa(msa, msa_out, grcut=1.1, hocut=-100, gcon=0.499, hcon=0.333, conplus=False, grswe=0.51):
    """
    MSA filtering + trimming, as described in 
    Esterman et. al (2021) https://doi.org/10.1093/ve/veab015
    """

    bad_cols = []
    
    if grcut <= 1:
        #print(f'getting a list of positions that are > {grcut * 100}% gaps')
        
        df = fa_to_df(msa, is_msa=True)
        
        #get a list of gap positions
        gapcols = get_gap_cols_from_msa(df, grcut)
        
        bad_cols += gapcols
        
    if hocut >= 0 or conplus:
        
        #then we need to iterate through each position in the MSA
        #print(f'iterating through the alignment')
        
        #consensus seq is empty
        conseq = []
        
        #get the necessary tables
        background_aa_freq,blosum62,expectedAAscores = get_tables_for_msa_filtering()
        
        #convert the msa to a dataframe 
        df = fa_to_df(msa, is_msa=True)
        nseqs, ncols = df.shape
        
        #get the HH94 scores of the sequences
        hh94_scores = calc_HH94(df, grswe)
        
        for col in range(0,ncols):
            
            #get the column
            column = df[col]

            #calculate the AA frequencies
            col_aa_frq = get_col_aa_frequencies(column, hh94_scores, background_aa_freq)

            #the last item in the dictionary is the gap frequencies
            #record it and remove from the dictionary,
            #otherwise the next fxn throws an error
            ngaps = col_aa_frq.popitem()[1]

            #convert the frequencies to scores
            col_aa_scores = aa_freq_to_score(col_aa_frq)

            #get the top scoring AA
            bestaa, bestaascore = get_best_aa_and_score(col_aa_scores)
                
            #assess the homogeneity of the AA, if the bestscoring AA is above expectation
            #if bestaascore/nseqs > expectedAAscores[bestaa]:
            homo = calc_homogeneity(bestaa, bestaascore, nseqs, expectedAAscores, blosum62)
            
            #record columns that are below the homogeneity threshold
            if homo < hocut:
                bad_cols.append(col)

            #assign the consensus
            con = assign_consensus(bestaa, gcon, homo, hcon, ngaps, nseqs)
            conseq.append(con)
        
        if conplus:
            
            #make a dataframe from the consensus
            df_c = pd.DataFrame.from_records([conseq], index = ["CONSENSUS"])
            
            #update the original MSA datafarme
            df = pd.concat([df_c,df])
        
        #remove gap/nonhomogenous columns
        if len(bad_cols) > 0:
            df = df.drop(columns = bad_cols)
            
        #write the output
        with open(msa_out, 'w') as f:
            df_to_fa(df, msa_out)
        
        #return df
        
    else:
        print('doing nothing, printing alignment')
        df = fa_to_df(msa, is_msa=True)
        with open(msa_out, 'w') as f:
            df_to_fa(df, msa_out)
        
