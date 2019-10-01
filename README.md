# Bioinformagick
Scripts for bioinformagicking

## <a name='TOC'>Table of Contents</a>

1. [BASH](#bash)
2. [Random Linux stuff](#Linux)
3. [E-utilities](#E-utilities)
4. [Biowulf](#Biowulf)
4. [CBB-Dev](#Biowulf)

...

# <a name="bash">Bash</a>

## Awk
Output lines where field 2 is less than 99
>awk '$2 > 98' fileIN

Remove blank lines
>awk ’NF > 0’ file

Find things in field 6 that match the string "Nov". For matching lines, add field 5 to variable `sum`. Print after reading all lines. 
>awk ’$6 == "Nov" { sum += $5 } END { print sum }’ file

Find things in Field 1 that have `J`. Replacing `~` with `!~` would do the opposite (things that don't have `J`). 
>awk ’$1 ~ /J/’ inventory-shipped file

remove duplicate lines
awk '!seen[$0]++' filename.txt

Search for a thing in tab-delimited column 2, return column 3
awk -F "\t" '$2==Refseq' {print $3}

Extract field 2 with tab separator
>awk ’BEGIN { FS = "\t" } ; { print $2 }’ file
or
>awk -F"\t" '{print $2}' file

Change output to tab delimited 
>awk 'BEGIN { OFS="\t"} {print $1,$2}'

Print lines where field $5 is empty
awk ’$5 == ""’ file

# <a name="linux">Linux</a>

Read in variables as an array
```
#!/bin/bash
while IFS=$'\t' read -r -a myArray
do
        echo "${myArray[0]}" "${myArray[1]}" "${myArray[2]}"
done < "fileIn.txt"
```
move a list of files
>cat list.txt | xargs mv -t /path/to/folder

>move files to another computer
```
sftp benler@anthill 
sftp> get filename
sftp> exit 
```
display file size of subfolders
>du -h --max-depth=1 | sort -hr

compress a folder using tar
>tar -zcvf archive.tar.gz directory/

Check if two directories are equivalent
>diff -r -q /path/to/folder1 /path/to/folder2

>sort a file on field 1 only
sort fileIN.txt -k 1,1 > fileOUT.txt

compare two sorted files and output unique lines
>comm file1 file2

add an alias to your .bashrc
>alias command_shorthand='/path/to/script.sh -options'

change colors of your terminal by modifying your .bashrc file
```
go to https://geoff.greer.fm/lscolors
paste the linux LS_colors output to your .bashrc file. example
export LSCOLORS=
```

count number of CPUS on a machine
>lscpu | egrep '^Thread|^Core|^Socket|^CPU\('

Split a fasta sequence into subsequences of length 60, each with their own header (note, only works for single sequence fasta files)
>grep -v '^>' chr1.fa | tr -d '\n' | fold -w 60 | nl -n rz -s '
' | sed 's/^/>fragment_/;N'

output a list of files with numbers 1-10, in a column
>echo file{1..10}.fasta | xargs -n1

display a particular line in a long file using Less
>less +13257894 -N bacteroides.faa

remove lines 1-13 in a file using sed, and create a backup
>sed -i.backup -e '1,13d' fileIN.fasta

Add "WP_" to the beginning of every new line
>sed -i -e 's/^/WP_/' fileIN.txt

# <a name="E-utilities">E-utilities)</a>

Get the coordinates of a domain annotated in the `Protein` database record
>efetch -db protein -id WP_025069012.1 -format gpc | xtract -insd Region region_name INSDInterval_from INSDInterval_to

# <a name="Biowulf">Biowulf</a>
#### Modules
list all modules
>module avail 

list only default modules
>module -d avail

Search for modules, case-insensitive
>module spider searchterm

Set up a default set of modules:
```
#Make a folder named 'mymodules' and a file called `ngs` that looks like this 
#%Module 1.0
module load bowtie/2-2.2.5
module load samtools/0.1.18

#Add that folder to the environment list
>module use --append /home/user/mymodules

#Load all the modules in the file `ngs`
>module load ngs
```
#### Job submission
Notes from 09/18/walkin:
1. Swarm files run each line in your file.swarm **independently**.
2. Sbatch is equivalent to qsub. In your sbatch command, you can specificy the # of CPUs, amount of memory, and which nodes to allocate to.
3. sinteractive is like the normal command line. Use this to debug/make sure everything works.
4. To figure out which nodes to allocate to, use `freen`. This lists the partition (i.e, queue) and the number of free nodes. For each node, it has the number of CPUs (i.e., threads) and memory. So, for the *norm* queue, I can set up to `--cpus-per-task= 56` and `--mem 246Gb`. I was told I should set CPUs ~16-32 and mem ~30 gb for genomics from the walkin guy.

#### data
The hi-performance drive is /data/
The high-volume / low-performance is /scratch
The good-performance/good-volume is /lscratch

# <a name="CBB-Dev">CBB-Dev</a>
Search Facilities
>ncbi-facilities-search

List all faciliities
>help_config

Change version on the fly
>facswitch python 3.7

Start a jupyter notebook
>jupyter lab --no-browser --port=8888 --ip=0.0.0.0

Submit a batch job with 4 day time limit, at 4gb RAM, no email
>qsub -v SGE_FACILITIES -P unified -l h_rt=259200,h_vmem=4G -n ./myjob.sh

Kill all jobs
>qdel -f -u benlersm

How to submit multithreaded jobs. The qsub command below reserves 4 x 4G = 16G of reserved and free memory (m_mem_free and mem_free = 4G X 4 slots), so only hosts with at least 16G of RAM would be considered by the SGE scheduler to run this job. It sets a virtual memory limit of 24G ( h_vmem=6G x 4 slots).  
>qsub -P unified -l h_vmem=6G,mem_free=4G,m_mem_free=4G -pe multicore 4 -R y ./fourway-job.sh
