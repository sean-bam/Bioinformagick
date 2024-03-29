import parsers
import pandas as pd
from pathlib import Path
import subprocess
import shutil

def entrez_protein_to_taxonomy(protein_accession_list, output):
    """
    This function expects a list of protein accession and returns 10 columns of:
    ProteinID,TaxID,ScientificName,Kingdom,Phylum,Class,Order,Family,Subfamily,Genus
    """
    assert type(protein_accession_list) == list or type(protein_accession_list) == set, f"{protein_accession_list} is not a list/set"
    assert len(protein_accession_list) > 0, f"{protein_accession_list} is empty!"
    
    if not Path(output).parent.is_dir():
        Path(output).parent.mkdir(parents = True)
    
    with open(output, 'w') as f:
        print(f"Accession", 
              "TaxID", 
              "ScientificName", 
              "Kingdom", 
              "Phylum", 
              "Class", 
              "Order", 
              "Family", 
              "Subfamily", 
              "Genus",
              sep = "~",
              file = f)
        
        for accession in protein_accession_list:
            p1 = subprocess.run(f'elink -db protein -target taxonomy -name protein_taxonomy -id {accession} | efetch -format xml | xtract -pattern Taxon -tab "~" -first TaxId ScientificName -group Taxon -KING "(-)" -PHYL "(-)" -CLSS "(-)" -ORDR "(-)" -FMLY "(-)" -SFMLY "(-)" -GNUS "(-)" -block "*/Taxon" -match "Rank:kingdom" -KING ScientificName -block "*/Taxon" -match "Rank:phylum" -PHYL ScientificName -block "*/Taxon" -match "Rank:class" -CLSS ScientificName -block "*/Taxon" -match "Rank:order" -ORDR ScientificName -block "*/Taxon" -match "Rank:family" -FMLY ScientificName -block "*/Taxon" -match "Rank:subfamily" -SFMLY ScientificName -block "*/Taxon" -match "Rank:genus" -GNUS ScientificName -group Taxon -tab "~" -element "&KING" "&PHYL" "&CLSS" "&ORDR" "&FMLY" "&SFMLY" "&GNUS"',
                               shell = True,
                                text = True,
                                capture_output = True,
                               check = True)
            print(accession, 
                  p1.stdout, 
                  sep = "~", 
                  #end = '', 
                  file = f)

def entrez_assembly_to_taxonomy(assembly_accession_list, output):
    """
    This function expects a list of protein accession and returns 10 columns of:
    ProteinID,TaxID,ScientificName,Kingdom,Phylum,Class,Order,Family,Subfamily,Genus
    """
    assert type(assembly_accession_list) == list or type(protein_accession_list) == set, f"{protein_accession_list} is not a list/set"
    assert len(assembly_accession_list) > 0, f"{assembly_accession_list} is empty!"
    
    if not Path(output).parent.is_dir():
        Path(output).parent.mkdir(parents = True)
    
    with open(output, 'w') as f:
        print(f"Accession", 
              "TaxID", 
              "ScientificName", 
              "Kingdom", 
              "Phylum", 
              "Class", 
              "Order", 
              "Family", 
              "Subfamily", 
              "Genus",
              sep = "~",
              file = f)
        
        for accession in assembly_accession_list:
            p1 = subprocess.run(f'elink -db assembly -target taxonomy -name assembly_taxonomy -id {accession} | efetch -format xml | xtract -pattern Taxon -tab "~" -first TaxId ScientificName -group Taxon -KING "(-)" -PHYL "(-)" -CLSS "(-)" -ORDR "(-)" -FMLY "(-)" -SFMLY "(-)" -GNUS "(-)" -block "*/Taxon" -match "Rank:kingdom" -KING ScientificName -block "*/Taxon" -match "Rank:phylum" -PHYL ScientificName -block "*/Taxon" -match "Rank:class" -CLSS ScientificName -block "*/Taxon" -match "Rank:order" -ORDR ScientificName -block "*/Taxon" -match "Rank:family" -FMLY ScientificName -block "*/Taxon" -match "Rank:subfamily" -SFMLY ScientificName -block "*/Taxon" -match "Rank:genus" -GNUS ScientificName -group Taxon -tab "~" -element "&KING" "&PHYL" "&CLSS" "&ORDR" "&FMLY" "&SFMLY" "&GNUS"',
                               shell = True,
                                text = True,
                                capture_output = True,
                               check = True)
            print(accession, 
                  p1.stdout, 
                  sep = "~", 
                  #end = '', 
                  file = f)            
            
def entrez_nuccore_to_taxonomy(accession_list, output):
    """
    Accepts a nucleotide ID
    Runs an entrez search and returns a string of:
    ProteinID,TaxID,ScientificName,Kingdom,Phylum,Class,Order,Family,Subfamily,Genus
    
    Empty values are marked with a '-'
    
    Notes: 
     - "Realms" don't have keywords yet, so I am not pulling those
     - This script is useful, because empty values are marked with a "-", 
     whereas pulling a taxa string from other places doesn't do that, 
     so its harder to parse
     -I can't parallelize this on the command line with 'join-into-groups-of', 
     because the links between the query and result get broken.
     -I Use "~" as a delimter, b/c some idiots use "," in the names of their species
    """
    assert type(accession_list) == list or type(accession_list) == set, f"{accession_list} is not a list/set"
    assert len(accession_list) > 0, f"{accession_list} is empty!"
    
    if not Path(output).parent.is_dir():
        Path(output).parent.mkdir(parents = True)
        
    assert not Path(output).exists(), f"The file {output} exists already and we don't want to overwrite it"
    
    #split the accessions into chunks of n
    chunker_gen = parsers.chunks(accession_list, 200)
    
    #keep a count of how many accessions don't return a result
    missed = 0
    taxids = []
    with open('tmp_accession_to_taxid_123.csv', 'w') as t:       
        for lst in chunker_gen:
            
            #get a table of accession->taxid
            command = f'efetch -format docsum -db nuccore -id {lst} | xtract -pattern DocumentSummary -element AccessionVersion TaxId '
            p1 = subprocess.run(command,
                                capture_output = True,
                                shell = True,
                                check = True,
                                text = True
                               )
            print(p1.stdout, file = t)
        
        #get a list of taxids
    df = pd.read_csv('tmp_accession_to_taxid_123.csv', sep = '\t', header = None, names = ["accession", "taxid"])
    taxids = df.taxid.unique().tolist()
        
    with open('tmp_taxid_to_lineage_123.csv', 'w') as f:
        
        #split the taxids into a reasonable amount
        chunker_gen = parsers.chunks(taxids, 200)
        
        for lst in chunker_gen:
        
            command = f'efetch -format xml -db taxonomy -id {lst} | xtract -pattern Taxon -tab "~" -first TaxId ScientificName -group Taxon -KING "(-)" -PHYL "(-)" -CLSS "(-)" -ORDR "(-)" -FMLY "(-)" -SFMLY "(-)" -GNUS "(-)" -block "*/Taxon" -match "Rank:kingdom" -KING ScientificName -block "*/Taxon" -match "Rank:phylum" -PHYL ScientificName -block "*/Taxon" -match "Rank:class" -CLSS ScientificName -block "*/Taxon" -match "Rank:order" -ORDR ScientificName -block "*/Taxon" -match "Rank:family" -FMLY ScientificName -block "*/Taxon" -match "Rank:subfamily" -SFMLY ScientificName -block "*/Taxon" -match "Rank:genus" -GNUS ScientificName -group Taxon -tab "~" -element "&KING" "&PHYL" "&CLSS" "&ORDR" "&FMLY" "&SFMLY" "&GNUS"'
        
            p1 = subprocess.run(command,
                                capture_output = True,
                                shell = True,
                                check = True,
                                text = True
                               )
     
    
            print(p1.stdout, 
                sep = "~", 
                #end = '', 
                file = f)
        
    #print the output
 
    cols = ["taxid", 
      "ScientificName", 
      "Kingdom", 
      "Phylum", 
      "Class", 
      "Order", 
      "Family", 
      "Subfamily", 
      "Genus",
       ]
        
    df2 = pd.read_csv('tmp_taxid_to_lineage_123.csv', header = None, names = cols, sep = "~")

    df3 = pd.merge(df, df2, how = 'left', on = 'taxid')
    df3.to_csv(output, index = False)
    
    #remove tmp files
    Path('tmp_accession_to_taxid_123.csv').unlink()
    Path('tmp_taxid_to_lineage_123.csv').unlink()

            
def entrez_nuccore_to_sra(protein_accession_list, output):
    """
    Given a nucleotide accession, (e.g., `OMBN01000154.1`), retrieve the largest SRA run accession.
    This sometimes works...
    """
    assert type(protein_accession_list) == list or type(protein_accession_list) == set, f"{protein_accession_list} is not a list/set"
    assert len(protein_accession_list) > 0, f"{protein_accession_list} is empty!"
    
    if not Path(output).parent.is_dir():
        Path(output).parent.mkdir(parents = True)
    
    with open(output, 'w') as f:
        print(f"Accession",
               "SRA_run",
               "SRA_total_spots",
               sep = ',',
               file = f)
        
        for accession in protein_accession_list:
                p1 = subprocess.run(f'esearch -db nuccore -query {accession} | elink -target biosample | elink -target sra | efetch -format docsum | xtract -pattern DocumentSummary -block Runs -tab "," -element Run@acc Run@total_spots | sort -k 2 -t, | head -n1',
                                    shell = True,
                                    text = True,
                                    capture_output = True,
                                    check = True)
                print(accession, 
                      p1.stdout, 
                      sep = ",",) 
                      #end = '', 
                      #file = f)
                        

    
def entrez_biosample_to_sra(protein_accession_list, output):
    """
    Given a biosample accession, (e.g., `OMBN01000154.1`), outputs a table of:
    Accession,BioProject,Run,bases,LibraryLayout,InsertSize
    """
    assert type(protein_accession_list) == list or type(protein_accession_list) == set, f"{protein_accession_list} is not a list/set"
    assert len(protein_accession_list) > 0, f"{protein_accession_list} is empty!"
    
    if not Path(output).parent.is_dir():
        Path(output).parent.mkdir(parents = True)
    
    with open(output, 'w') as f:
        print(f"Accession",
               "BioProject",
               "Run",
               "Bases",
               "LibraryLayout",
               "InsertSize",
               sep = '\t',
               file = f)
              
        for accession in protein_accession_list:
                p1 = subprocess.run(f'esearch -db biosample -query "{accession} AND biosample sra[filter]"  | elink -target sra | efetch -format runinfo -mode xml | xtract -pattern Row -element BioProject,Run,bases,LibraryLayout,InsertSize',
                                    
                                    shell = True,
                                    text = True,
                                    capture_output = True,
                                    check = True)
                
                for line in p1.stdout.strip().split('\n'):
                    print(accession, line, sep = '\t', file = f)
                    
def entrez_nuccore_to_protein(accession_list, output, chr_start=0, chr_stop=100000000000):
    """
    Returns a tab-separated table of 
    Nucleotide_accssion, coordinates, protein_accession, product
    
    Optional chromosome start/stop filtering
    """
    
    if type(accession_list) != list:
        accession_list = [accession_list]
    #assert type(accession_list) == list, f"{accession_list} is not a list"
    #assert len(accession_list) > 0, f"{accession_list} is empty!"
    
    if not Path(output).parent.is_dir():
        Path(output).parent.mkdir(parents = True)
    
    with open(output, 'w') as f:
        print(f"Accession",
               "Raw_coords",
               "Protein_id",
               "Product",
               sep = '\t',
               file = f)
        
        chunker_gen = parsers.chunks(accession_list, 200)
        for lst in chunker_gen:
            query = ','.join(lst)
            
                
            p1 = subprocess.run(f"efetch -db nuccore -id {query} -format gbc -chr_start {chr_start} -chr_stop {chr_stop} | xtract -insd CDS INSDFeature_location protein_id product",
                                    shell = True,
                                    text = True,
                                    capture_output = True,
                                    check = True)
                
            print(p1.stdout, file = f)
                        

def entrez_assembly_to_file(accession_list, output, filetype = "genbank", database = "GenBank"):
    """
    Given an assembly accession, retrieves the given filetype via wget
    """
    
    assert type(accession_list) == list, f"{accession_list} is not a list"
    assert len(accession_list) > 0, f"{accession_list} is empty!"
    
    if not Path(output).parent.is_dir():
        Path(output).parent.mkdir(parents = True)
        
    assert database == "GenBank" or database == "RefSeq", f"Database must be 'RefSeq' or 'GenBank'"
    assert filetype == "genbank" or  filetype == "nucleotide" or filetype == 'protein', f"""
    filetype must be 'genbank' 'nucleotide' or 'protein'
    """

    chunker_gen = parsers.chunks(accession_list, 200)
    for lst in chunker_gen:
        query = ','.join(accession_list)
        p1 = subprocess.run(f"efetch -db assembly -id {query} -format docsum | xtract -pattern DocumentSummary -element FtpPath_{database}",
                           shell = True,
                           capture_output = True,
                           check = True,
                            text = True,
                          )
        #flag queries with no result
        if len(p1.stdout.split()) < len(accession_list):
            print("The number of entrez results is less than the number of search queries!")
            
        i = 0
        for ftp in p1.stdout.split():
            assembly_base = ftp.split("/").pop()
            
            if filetype == "genbank":
                url = ftp + "/" + assembly_base + "_genomic.gbff.gz"
            elif filetype == "nucleotide":
                url = ftp + "/" + assembly_base + "_genomic.fna.gz"
            elif filetype == "protein":
                url = ftp + "/" + assembly_base + "_protein.faa.gz"
                
            #retry 1 once, don't overwrite, wait 1 sec
            subprocess.run(f'wget {url} -t 1 -w 1 -nc -P {output}',
                          shell = True,
                          check = True)
            i += 1
    print(f"downloaded {i} files to {output}")
    
def entrez_cdd_superfamily_to_members(accessions):
    """
    given an iterable of superfamily accessions (e.g., cl00656),
    gets all the cdd profiles belonging to that supercluster
    returns a dataframe
    note: singleton clusters return no hit
    """
    superfamily_dict = {}
    for accession in accessions:
        p1 = subprocess.run(f'elink -db cdd -target cdd -id {accession} -name cdd_cdd_superfamily_2 | efetch -format docsum | xtract -pattern DocumentSummary -element Accession',
                      check = True,
                      text = True,
                      shell = True,
                      capture_output = True)

        for profile in p1.stdout.strip().split("\n"):
            superfamily_dict.update({profile : accession})
        
    df = pd.DataFrame.from_dict(superfamily_dict, orient = 'index', columns = ['superfamily'])
    df2 = df.reset_index().rename(columns = {'index' : 'profile'})
    
    return df2

