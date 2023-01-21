from BCBio import GFF
from Bio import SeqIO, SeqFeature
from Bio.SeqRecord import SeqRecord

def get_proteins_from_fasta_with_gff(fasta,gff,output,start=0,stop=1000000000000):
    """
    Accepts a fasta file and gff file
    extracts the proteins to a new file, 
    optionally filtering by location
    """
    
    #import the gff annotations to a seqrecord
    seq_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    seq_record_gen = GFF.parse(gff, base_dict=seq_dict)
    
    
    prot_records = []
    for annotated_record in seq_record_gen:
        #print(annotated_record.id)
        for feature in annotated_record.features:
            if feature.location.start.position > start and feature.location.end.position < stop: 
                try:
                    prot = feature.translate(annotated_record, 
                                  cds = True, 
                                  table = 11, 
                                  to_stop = True)
                    
                #If the last ORF is partial w/o a stop codon (runs off the edge)
                except Bio.Data.CodonTable.TranslationError:
                    #print(f'translation error on {feature.id}')
                    pass
                orf_id = f'{annotated_record.id}__{feature.location.start.position}_{feature.location.end.position}_{feature.location.strand}'
                #orf_id = f'{annotated_record.id}__{feature.id}'
                prot_records.append(SeqRecord(seq=prot.seq,id=orf_id,description=""))
    SeqIO.write(prot_records, output, "fasta")