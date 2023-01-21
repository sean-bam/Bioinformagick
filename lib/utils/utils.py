"""
A smorgasbord of random useful fxns.
"""

import string
#import random

def df2gff(df, output):
    """
    takes a dataframe with the following columns:
    "accession", "type", "start", "stop", "strand", "id", "product"
    strand can be +/- or 1/-1. Otherwise, set to ".".
    
    prints a table with gff columns
    1. seqid
    2. source
    3. type
    4. start
    5. stop
    6. score
    7. strand
    8. phase
    9. attributes
    """
    
    df = (df.loc[:,["accession",
                    "type",
                   "start",
                   "stop",
                   "strand",
                   "id",
                   "product"
                  ]
               ]
            .reset_index(drop = True)
            .astype({"start" : "int64",
                   "stop" : "int64",
                     "strand" : "str"
                  }
                   )
         )
    
    #fix the strand columns
    df['strand'] = df.strand.replace("1", "+")
    df['strand'] = df.strand.replace("+1", "+")
    df["strand"] = df.strand.replace("-1", "-")
    df["strand"] = df.strand.replace("nan", ".")
    
    #fill empty "products" and remove any comma's
    #to be compatible with GFF file format specs
    df["product"] = df["product"].fillna(".")
    df["product"] = df["product"].str.replace(",", "+")
    
    #make the gff columns
    df["source"] = "DB"
    df["score"] = "."
    df["phase"] = 0
    df["attributes"] = "ID=" + df["id"] + ";" + "Name=" + df["product"]
    
    #rearrange the columns
    df = df.loc[:, ["accession",
                    "source",
                    "type",
                    "start",
                    "stop",
                    "score",
                    "strand",
                    "phase",
                    "attributes"
                   ]
               ]
    
    
    #print the output
    df.to_csv(output, 
              header = False, 
              index = False, 
              sep = '\t')
    #return df
    
def df2gbk(df, output, flip=False):
    """
    takes a dataframe with the following columns:
    "accession", "length", "start", "stop", "strand", "type", "product"
    strand can be +/- or 1/-1. Otherwise, set to None.
    """

    df = (df.loc[:,["accession", 
                   "length",
                   "start",
                   "stop",
                   "strand",
                   "type",
                   "product"
                  ]
               ]
            .reset_index(drop = True)
            .astype({"start" : "int64",
                   "stop" : "int64",
                     "strand" : "str"
                  }
                   )
         )
    
    #fix the strand columns
    df['strand'] = df.strand.replace("+", "+1")
    df['strand'] = df.strand.replace("1", "+1")
    df["strand"] = df.strand.replace("-", "-1")
    
    #fill empty "products"
    df["product"] = df["product"].fillna("hypothetical")
    
    accession = df.at[0,"accession"]
    length = df.at[0,"length"]
    
    #initilize a seq record
    unk_dna = UnknownSeq(length, character="N")
    my_sequence_record = SeqRecord(unk_dna, id = accession)
    my_sequence_record.annotations["molecule_type"] = "DNA"

    #add annotations
    flip = False
    for index, row in df.iterrows():
        start = row[2]
        stop = row[3]
        strand = row[4]
        feature_type = row[5]
        product = row[6]
        
        #this is the only way I deal with strand apprpriately,
        #after spending way too much time on this
        if strand == "+1" or strand == "-1":
            strand = int(strand)
        else:
            strand = None
        
        #make all non-CDS features "CDS" features, fow now
        if feature_type != "CDS":
            feature_type = "CDS"
            
        #put in a Gene
        my_feature = SeqFeature(FeatureLocation(start, stop), 
                                type="gene", 
                                strand=strand,
                                qualifiers = {"locus_tag" : product}
                               )
        my_sequence_record.features.append(my_feature)
        
        #put in the feautres
        my_feature = SeqFeature(FeatureLocation(start, stop), 
                                type=feature_type, 
                                strand=strand,
                                qualifiers = {"product" : product,
                                              "codon_start" : 1,
                                              }
                               )
        print(my_feature.translate(cds=True))
        my_sequence_record.features.append(my_feature)
        
    if flip:
        #reverse complement, keeping the ID's the same
        my_sequence_record.reverse_complement(id=True,
                                             name=True,
                                             description=True,
                                             annotations=True)
    #eqIO.write(my_sequence_record, output, "genbank")
    return my_sequence_record

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
    
def pdb_id_generator(size=4, chars=string.ascii_uppercase + string.digits):
    """
    Generates a random 4 character string
    from: 
    https://stackoverflow.com/questions/2257441/random-string-generation-with-upper-case-letters-and-digits
    """
    return ''.join(random.choice(chars) for _ in range(size))

def slice_pdb(pdbid, pdbfile, start, end, output):
    """
    Accepts a PDBfile and start/end coordinates
    Outputs a new file with just the structure between start-end
    
    End can be larger than length of the molecule
    """
    from Bio.PDB.PDBParser import PDBParser
    import Bio.PDB.Dice as BPD
    
    parser = PDBParser()
    structure = parser.get_structure(pdbid, pdbfile)
    BPD.extract(structure, 'A', start, end, output)
    
def get_prot_coordinates_from_gbk(genbank_file, accessions, output):
    """
    Given an iterable of protein accessions, gets the nucleotide coordinates from a genbank file
    """
    with open(output, 'a') as f:
        for seq_record in SeqIO.parse(genbank_file, "genbank"):
            contig = seq_record.id.split("|")[0]
            
            for feature in seq_record.features:
                prot_id_list = feature.qualifiers.get('protein_id')
                if prot_id_list:
                    prot_id = prot_id_list[0]
                    if prot_id in accessions:
                        print(contig,
                              prot_id, 
                              feature.location.start.position, 
                              feature.location.end.position, 
                              file = f)
                        
def gbk_to_prots(genbank, output):
    """
    Writes fasta formatted protein seqs to output
    """
    with open(output, 'w') as f:
        for seq_record in SeqIO.parse(genbank, "gb"):
            for feature in seq_record.features:
                if feature.type == "CDS":
                    print(">", feature.qualifiers.get('protein_id')[0], sep = "", file = f)
                    print(feature.qualifiers.get('translation')[0], file = f)
                    
def concatenate_files(file_list, output):
    
    output = Path(output)
    
    with output.open('wb') as wfd:
        for infile in file_list:
            with open(infile, 'rb') as fd:
                shutil.copyfileobj(fd, wfd)

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
                
