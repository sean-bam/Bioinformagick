# Database setup

DALI requires PDB files be converted into an internal format. There are a few ways to do this. Ultimately, I chose option 3, which using both option 1 and option 2 outputs.

Basically, make a mirror of the PDB database using the --rsync option of import.pl (option 1). This will try to import every file in the PDB archive, but will crash. Then, get a list of cluster representatives from the PDB website (option 2). Use this list to run import.pl for the PDB archive.

In the future, try just getting a copy of the PDB archive without running import.pl --rsync.

## Option 1

The first option is to run the impoort.pl script with --rsync, which downloads the entire PDB. For each PDB entry, it gets a *.ent.gz* file and then converts every single entry to the internal data format. I'll try this option first:

import.pl --rsync --pdbmirrordir dbs/pdb/ --dat dbs/dali/dali_pdb/DAT/ --clean

However, I get this error:
```
Reading DSSP file
At line 288 of file ../src/util.f
Fortran runtime error: Bad integer for item 1 in list input
Error in puu: /usr/local/DALI/5.0/bin/puu 5a1vU 32240.tmp 5a1v.dssp
```
There is something wrong with one of the PDB chains that is breaking the automatic import. Once import.pl breaks, it creates a `dali.lock` file that prevents subsequent imports. The second problem is that it appears that 250000 files is the max that can be put into a folder on our system, and there are ~450000 PDB entries.

## Option 2
Rather than doing the whole PDB via import.pl --rsync, just do a subset.

PDB regularly makes a subset clustered with MMSeqs2 at various thresholds, available [here](https://www.rcsb.org/docs/programmatic-access/file-download-services#:~:text=4HHB.A/download-,Sequence%20clusters%20data,-Results%20of%20the). I can download the PDB files using the PDB script `batch_download.sh`, available [here](https://www.rcsb.org/docs/programmatic-access/batch-downloads-with-shell-script). The problem with this approach is that some pdb entries are too large to download the `.pdb.gz` file and would take ~1 day to complete

`batch_download.sh` requires a file with a comma-separated list of IDs. So, I converted the clustering table, where I assume the first entry is the representative, to this format.

>./batch_download.sh -f pdb70/bc-70.csv -o pdb70/ -p

Again, some files crash `import.pl` and create a dali.lock file, so i have to import them individually. :/

## Option 3
Option 2 is too slow. PDB provides an rsync script that mirrors the archive, available [here](https://www.rcsb.org/docs/programmatic-access/file-download-services#:~:text=the%20ftp%20protocol.-,Automated%20download%20of%20data,-The%20URLs%20in). My guess is this is what `import.pl --rysnc` uses.

So, use the `*.ent` files generated in option 1, and the list of representatives in option 2, to import the files.

Again, If I use `import.pl --pdblist pdb70_subset.list`, where `pdb70_subset.list` is a list of file paths to `*.ent.gz` files, I run into the same issue where a bad file crashes the program. So, I wrote a python script to import them individually, and delete the `dali.lock` file if it appears.

Import.pl extracts all chains from a PDB file, but only some of them are actually representatives. I wrote a script to remove the non-representative DAT files. 

# Running Dali

An example command:

```bash
dali.pl --cd1 T43NA --dat1 output/exp1/models10/dat/ --dat2 dbs/dali/dali_pdb/pdb70/ --clean --oneway --hierarchical --repset pdb70.list --db pdb70.list --outfmt "summary,alignments"
```

- cd1 : Name of the (fake) PDB ID of the .DAT file
 - dat1 : location of the .DAT file. If using a PDB file as a query instead, this is where a temporary .DAT file is made (defaults to ./DAT, throws an error if this dir isn't present
 - dat2 : location of the database .DAT files
 - repset : only consider a subset of files
 - db : I'm not sure if this is needed or not when using repset.

Other params
-pdbid1 : name of the output file (defaults to mol1A)

# Structure metadata

It would be highly desirable to have domains mapped onto known/predicted structures. There are a couple of strategies I can think of:

1. The Pfam FTP offers a mapping of Pfam domains --> known PDB structures [here](http://ftp.ebi.ac.uk/pub/databases/Pfam/mappings/pdb_pfam_mapping.txt). There are ~75700 mappings. The drawback is this is only pfam domains and known structures.
2. The Pfam FTP has a file of for RosettaFold. However, it is just a3m-formatted MSAs of each pfam entry.
3. On the 'AlphaFold structures' tab of a given Pfam, there is a list of proteins "that match this family and have AlphaFold structures". This is nice, because it includes predicted structures in the AlphaFold database. However, there is no mapping of the coordinates of the pfam domain onto the predicted structure and I don't know how to download this data. I think I would first have to get a table of Pfam-->Uniprot ProteinID mappings, then get ProteinID mappings --> AlphaFold structure. To do the former, there is a [Spark-SQL query interface](https://sparql.uniprot.org/sparql), but I don't know how to do construct the query, so I contacted the helpdesk.
4. Entrez offers a [linkname](https://www.ncbi.nlm.nih.gov/entrez/query/static/entrezlinks.html) to go from CDD-->known PDB structure. However, I can't find the coordinates of the domain on the structure, and it only includes known PDB structures.
