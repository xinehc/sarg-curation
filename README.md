# SARG+
<a href="https://doi.org/10.5281/zenodo.13121333"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.13121333.svg" alt="DOI"></a>

## Introduction

SARG+ is a manually curated database of Antibiotic Resistance Genes (ARGs), designed to enhance read-based environmental surveillance at species-level resolution. It extends existing databases ([NDARO](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/), [CARD](https://card.mcmaster.ca/), and [SARG](https://smile.hku.hk/ARGs/Indexing)) by incorporating a comprehensive collection of protein sequences from RefSeq that are annotated through the same evidence (BlastRules or Hidden Markov Models provided by the NCBI Prokaryotic Genome Annotation Pipeline, [PGAP](https://github.com/ncbi/pgap)) as experimentally validated ARGs. This expansion addresses the limitations of existing databases, which often include only a single or a few representative sequences per ARG, and allows for the use of more stringent cutoffs while maintaining sensitivity.

**SARG+** ([sarg.fa](https://github.com/xinehc/sarg-curation/blob/master/sarg.fa)) consists of two components:

1. **SARG+ reference** ([sarg_ref.fa](https://github.com/xinehc/sarg-curation/blob/master/sarg_ref.fa)): Contains only experimentally validated sequences.
2. **SARG+ extension** ([sarg_ext.fa](https://github.com/xinehc/sarg-curation/blob/master/sarg_ext.fa)): Contains computationally derived homologs of the reference sequences.

## Prerequisites

### Install Dependencies

Create a new conda environment with the necessary dependencies:

```bash
conda create -n sarg-curation -c bioconda -c conda-forge blast diamond mmseqs2 seqkit
conda activate sarg-curation
```

Install additional Python modules and Jupyter for running the notebooks:

```bash
conda install jupyter regex json5 wget tqdm biopython pandas
```

### Download NCBI Databases
> [!TIP]
> Extracting sequences and annotation evidence can be time-consuming. To take advantage of multiple CPU cores, run `blastdbcmd` and `diamond makedb` in parallel.

Download `nr` and `env_nr` from [NCBI FTP](https://ftp.ncbi.nlm.nih.gov/blast/db/) and extract sequences:

```bash
for db in nr env_nr
do
    curl ftp://ftp.ncbi.nlm.nih.gov/blast/db/ \
        | grep -o '[^ ]*.gz$' \
        | grep ^$db \
        | xargs -P 48 -I {} wget -qN --show-progress https://ftp.ncbi.nlm.nih.gov/blast/db/{} -P tmp/protein/$db

    for file in tmp/protein/$db/*.tar.gz; do tar -xvf $file -C tmp/protein/$db; done
done

blastdbcmd -db tmp/protein/nr/nr -entry all > tmp/nr.fa
blastdbcmd -db tmp/protein/env_nr/env_nr -entry all > tmp/env_nr.fa
rm -rf tmp/protein
```

Download `refseq_protein` from [NCBI RefSeq](https://ftp.ncbi.nlm.nih.gov/refseq/release/) and extract annotation evidence:

```bash
for kingdom in archaea bacteria
do
    curl ftp://ftp.ncbi.nlm.nih.gov/refseq/release/${kingdom}/ \
        | grep '[^ ]*wp_protein.*.gz$' -o \
        | xargs -P 48 -I {} wget -qN --show-progress https://ftp.ncbi.nlm.nih.gov/refseq/release/${kingdom}/{} -P tmp/refseq
done

cat tmp/refseq/*.faa.gz > tmp/refseq_protein.faa.gz
diamond makedb --in tmp/refseq_protein.faa.gz --db tmp/refseq_protein
rm -rf tmp/refseq_protein.faa.gz

python -c "
import gzip
import glob
import pandas as pd

from Bio import SeqIO
from tqdm.contrib.concurrent import process_map

def parser(file):
    lines = []
    with gzip.open(file, 'rt') as f:
        for record in SeqIO.parse(f, 'genbank'):
            header = record.description.rsplit(' [')[0].split('MULTISPECIES: ')[-1]
            if 'structured_comment' in record.annotations:
                annotation = record.annotations['structured_comment']['Evidence-For-Name-Assignment']
                identifier, evidence = annotation.get('Source Identifier', 'NA'), annotation.get('Evidence Accession', 'NA')
            else:
                identifier, evidence = 'NA', 'NA'

            gene = [x for x in record.features if x.type == 'gene']
            symbol = gene[0].qualifiers['gene'][0] if gene else 'NA'

            lines.append([record.id, header, evidence, symbol, identifier])
    pd.DataFrame(lines).to_csv(file.replace('.gpff.gz', '.tsv'), sep='\t', header=None, index=False)

r = process_map(parser, glob.glob('tmp/refseq/*.gpff.gz'), max_workers=48, chunksize=1)
"
```

### Download NDARO and CARD

NDARO need to be downloaded from https://www.ncbi.nlm.nih.gov/pathogens/refgene/ (`refgenes.tsv` and `protein.faa`). CARD can be obtained from https://card.mcmaster.ca/download (`aro_index.tsv` and `protein_fasta_protein_homolog_model.fasta`). These files need to be manually unzipped and placed in the [reference](https://github.com/xinehc/sarg-curation/tree/master/reference) folder.

## Run

We provide a series of Jupyter notebooks for step-wise construction of SARG+:

1. `a0-parse-refs.ipynb`
   - Parses NDARO and CARD metadata and sequences to create a raw reference. Curates the reference according to `sarg.json`, producing the initial SARG+ reference database.
2. `a1-standardize-headers.ipynb`
   - Standardizes the headers of SARG+ reference sequences according to `nr` and `env_nr`.
4. `b0-parse-evidence.ipynb`
   - Finds sequences annotated through the same evidence (BlastRules and Hidden Markov Models) as SARG+ reference sequences.
5. `b1-remove-dups.ipynb`
   - Removes duplicated and cross-mapped sequences by clustering.

> [!NOTE]
> Scripts `z0` and `z1` are used for testing and are not required for building SARG+.

## Output

- **Summary File:** `misc/summary.tsv` provides a summary of SARG+ reference sequences, including their types, subtypes, and sources.
- **Evidence File:** `misc/evidence.tsv` lists all evidence used for creating SARG+ extension.
- **Counts:** `misc/sarg.txt`, `misc/sarg_ref.txt`, and `misc/sarg_ext.txt` display the counts of different ARGs for SARG+, SARG+ reference, and SARG+ extension, respectively.
- **Sequences:**
  - `sarg.fa`: combined reference and extension sequences.
  - `sarg_ref.fa`: the reference component of SARG+.
  - `sarg_ext.fa`: the extension component of SARG+.

## Contribution

New ARG sequences can be integrated into SARG+ by editing [sarg.json](https://github.com/xinehc/sarg-curation/blob/master/sarg.json) and [reference/reference.fasta](https://github.com/xinehc/sarg-curation/blob/master/reference/reference.fasta).

- `sarg.json`: Specifies the ARG type (class/family) for each new gene, relevant literature, rationales, and links to sequences.
- `reference.fasta`: Contains the protein sequences of these new genes.

For example:

```json5
// aminonucleoside
"aminonucleoside": {
    "added": {
        // "Molecular analysis of the pac gene encoding a puromycin N-acetyl transferase from Streptomyces alboniger"
        // Detoxification of puromycin.
        // https://www.uniprot.org/citations/2676728
        "sp|P13249|PUAC_STRAD": "pac", // Puromycin N-acetyltransferase OS=Streptomyces alboniger OX=132473 GN=pac PE=1 SV=2 | GNAT family N-acetyltransferase

        // "The pur8 gene from the pur cluster of Streptomyces alboniger encodes a highly hydrophobic polypeptide which confers resistance to puromycin"
        // May be involved in active puromycin efflux energized by a proton-dependent electrochemical gradient. In addition, it could be implicated in secreting N-acetylpuromycin, the last intermediate of the puromycin biosynthesis pathway, to the environment.
        // https://www.uniprot.org/citations/7916693
        "sp|P42670|PUR8_STRAD": "pur8", // Puromycin resistance protein pur8 OS=Streptomyces alboniger OX=132473 GN=pur8 PE=3 SV=1 | MFS transporter

        // "The ard1 gene from Streptomyces capreolus encodes a polypeptide of the ABC-transporters superfamily which confers resistance to the aminonucleoside antibiotic A201A"
        // The gene ard1 induced antibiotic resistance that was highly specific for A201A.
        // https://www.ncbi.nlm.nih.gov/nuccore/X84374
        "CAA59109.1": "ard1", // Ard1 protein [Saccharothrix mutabilis subsp. capreolus] | Ard1 protein

        // "The aminonucleoside antibiotic A201A is inactivated by a phosphotransferase activity from Streptomyces capreolus NRRL 3817, the producing organism. Isolation and molecular characterization of the relevant encoding gene and its DNA flanking regions"
        // A novel resistance determinant (ard2) to the aminonucleoside antibiotic A201A was cloned from Streptomyces capreolus NRRL 3817, the producing organism, and expressed in Streptomyces lividans.
        // https://www.ncbi.nlm.nih.gov/nuccore/X84374
        "CAD62197.1": "ard2", // Ard2 protein [Saccharothrix mutabilis subsp. capreolus] | Ard2 protein
    }
}
```

If you identify any suspicious entries or wish to contribute sequences to SARG+, please consider creating a pull request by editing `sarg.json` and `reference/reference.fasta`, or opening an issue.

## FAQ

### Which ARGs are included/excluded?

- All ARGs associated with point mutations in essential genes (primarily antibiotic targets) are excluded. Examples include mutations in *gyrA*, *parC*, and *rpoB*.

- Regulators (e.g., activators, repressors) are excluded since they do not confer direct resistance. Exceptions include *tipA* and *albAB*, which act as self-regulated sequesters. Putative accessory genes such as *vanZ* are also removed.

- Fused genes are removed since they can create ambiguities when being aligned using reads (note that SARG+ is designed for read-based profiling). For instance, *catB/aac(6')-I* ([WP_071593228.1](https://www.ncbi.nlm.nih.gov/protein/WP_071593228.1)) is a fusion of *catB* ([WP_264840997.1](https://www.ncbi.nlm.nih.gov/protein/WP_264840997.1)) and *aac(6')-I* ([WP_033917551.1](https://www.ncbi.nlm.nih.gov/protein/WP_033917551.1)).

- Highly similar ARGs are grouped into subtype clusters to reduce the chance of false identifications. For example, the two alleles of *bla*<sub>OXA</sub>, *bla*<sub>OXA-1</sub> ([WP_001334766.1](https://www.ncbi.nlm.nih.gov/protein/WP_001334766.1)) and *bla*<sub>OXA-1042</sub> ([WP_231869587.1](https://www.ncbi.nlm.nih.gov/protein/WP_231869587.1)), differ by a single amino acid. Such subtle differences can be difficult to resolve using reads. By default, clustering uses thresholds of 95% sequence identity and 95% query/subject coverage.

### Why some genes have uncommon names?

ARG names are standardized to ensure unambiguous identification by type and subtype. Examples include:

- *mdtP* refers to an RND efflux pump in *Escherichia coli* and an MFS transporter in *Bacillus subtilis*. To avoid confusion, the *Bacillus* variant is renamed *mdt(P)*, aligning it with *mdt(A)* (another MFS transporter).
- Genes with labeling errors in RefSeq are corrected (e.g., *efrCD* is misspelled as *erfCD*).
- *qacA* and *qacB* are renamed to *qacA/B* due to their high sequence similarity.

For naming consistency, certain ARG names are revised:

- *tnrB2* and *tnrB3* are renamed to *tnrB-1* and *tnrB-2* to reflect their roles as ABC transporter components.
- *cap21* denotes *orf21* in a *Streptomyces griseus* biosynthetic gene cluster ([AB476988](https://www.ncbi.nlm.nih.gov/nuccore/AB476988)), which lacks an established name.


## Citation
Chen, X., Yin, X., Xu, X., & Zhang, T. (2025). Species-resolved profiling of antibiotic resistance genes in complex metagenomes through long-read overlapping with Argo. *Nature Communications*, 16(1), 1744. https://doi.org/10.1038/s41467-025-57088-y

Alcock, B. P., Huynh, W., Chalil, R., Smith, K. W., Raphenya, A. R., Wlodarski, M. A., ... & McArthur, A. G. (2023). CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database. *Nucleic acids research*, 51(D1), D690-D699. https://doi.org/10.1093/nar/gkac920

Feldgarden, M., Brover, V., Gonzalez-Escalona, N., Frye, J. G., Haendiges, J., Haft, D. H., ... & Klimke, W. (2021). AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. *Scientific reports*, 11(1), 12728. https://doi.org/10.1038/s41598-021-91456-0
