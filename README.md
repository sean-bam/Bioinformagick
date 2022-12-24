# Bioinformagick
Scripts for bioinformagicking

## <a name='TOC'>Table of Contents</a>

1. [BASH](#bash)
2. [Random Linux stuff](#Linux)
3. [E-utilities](#E-utilities)
4. [GCP](https://sean-bam.github.io/Bioinformagick/GCP)
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
>awk ’$5 == ""’ file

Get the best hit of blast result
>awk '!a[$1]++'"

# <a name="linux">Linux</a>

Rename a fasta header using the filename
var=$(echo "file.fa" | cut -d'.' --complement -f2-).mito
sed -i "1s/.*/>$var/" file.fa

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
