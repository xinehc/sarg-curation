# sarg-curation
SARG+ construction

## Prerequisite
```bash
mkdir -p protein
for folder in env_nr nr refseq_protein
do
    curl ftp://ftp.ncbi.nlm.nih.gov/blast/db/ \
        | grep '[^ ]*.gz$' -o \
        | grep ^$folder \
        | xargs -P 32 -I {} wget -qN --show-progress https://ftp.ncbi.nlm.nih.gov/blast/db/{} -P protein/$folder
    ## decompress all files
    for file in protein/$folder/*.tar.gz; do tar -xvf $file -C protein/$folder; done
    rm -rf protein/$folder/*.tar.gz
done

blastdbcmd -db protein/env_nr/env_nr -entry all > env_nr.full.fa
blastdbcmd -db protein/nr/nr -entry all > nr.full.fa
blastdbcmd -db protein/refseq_protein/refseq_protein -entry all > refseq_protein.full.fa
diamond makedb --in refseq_protein.full.fa --db refseq_protein
rm -rf protein
```