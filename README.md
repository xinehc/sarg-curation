# SARG+

## Introduction

SARG+ is a manually curated database of Antibiotic Resistance Genes (ARGs), designed to enhance read-based environmental surveillance at species-level resolution. It extends existing databases ([NDARO](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/), [CARD](https://card.mcmaster.ca/), and [SARG](https://smile.hku.hk/ARGs/Indexing)) by incorporating a comprehensive collection of protein sequences from RefSeq that are annotated through the same evidence sources (BlastRules or Hidden Markov Models provided by the NCBI Prokaryotic Genome Annotation Pipeline, [PGAP](https://github.com/ncbi/pgap)) as experimentally validated ARGs. This expansion addresses the limitations of existing databases, which often include only a single or a few representative sequences per ARG, and allows for the use of more stringent cutoffs while maintaining sensitivity.

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
conda install jupyter regex json5 biopython wget
```

### Download NCBI Databases

Download the `nr`, `env_nr`, and `refseq_protein` databases from [NCBI FTP](https://ftp.ncbi.nlm.nih.gov/blast/db/) and extract sequences:

> [!NOTE]
> Extracting sequences and creating `diamond` databases can be time-consuming. Consider running `blastdbcmd` and `diamond makedb` in parallel to utilize multiple cores.

```bash
mkdir -p tmp/protein
for folder in nr env_nr refseq_protein; do
    mkdir -p tmp/protein/$folder
    curl ftp://ftp.ncbi.nlm.nih.gov/blast/db/ \
        | grep -o '[^ ]*.gz$' \
        | grep ^$folder \
        | xargs -P 32 -I {} wget -qN --show-progress https://ftp.ncbi.nlm.nih.gov/blast/db/{} -P tmp/protein/$folder

    for file in tmp/protein/$folder/*.tar.gz; do tar -xvf $file -C tmp/protein/$folder; done
    rm -rf tmp/protein/$folder/*.tar.gz
done

blastdbcmd -db tmp/protein/nr/nr -entry all > tmp/nr.fa
blastdbcmd -db tmp/protein/env_nr/env_nr -entry all > tmp/env_nr.fa
blastdbcmd -db tmp/protein/refseq_protein/refseq_protein -entry all > tmp/refseq_protein.fa
diamond makedb --in tmp/refseq_protein.fa --db tmp/refseq_protein
rm -rf tmp/protein
```

### Download NDARO and CARD

NDARO can be downloaded manually from https://www.ncbi.nlm.nih.gov/pathogens/refgene/ (click `Download` for both the metadata `refgenes.tsv` and reference protein sequences `refgene_catalog.zip`). CARD can be obtained from https://card.mcmaster.ca/download. These files need to be unzipped and placed in the [reference](https://github.com/xinehc/sarg-curation/tree/master/reference) folder.


## Run

We provide a series of Jupyter notebooks for step-wise construction of SARG+:

1. **`a0-parse-refs.ipynb`**
   - Parses NDARO and CARD metadata and sequences to create a raw reference. Curates the reference according to `sarg.json`, producing the initial SARG+ reference database.
2. **`a1-standardize-headers.ipynb`**
   - Standardizes headers of the SARG+ reference according to `nr`.
3. **`b1-get-remarks.ipynb`**
   - Retrieves annotation evidence of SARG+ sequences.
4. **`b2-parse-remarks.ipynb`**
   - Finds sequences annotated through the same evidence sources (BlastRules and Hidden Markov Models) as the reference sequences.
5. **`b3-remove-dups.ipynb`**
   - Removes duplicated and cross-mapped sequences by clustering.

> [!NOTE]
> Scripts `z0` and `z1` are used for testing and are not required for building SARG+.

## Output

- **Summary File:** `misc/summary.tsv` provides a summary of SARG+ reference sequences, including their types, subtypes, and sources.
- **Evidence File:** `misc/evidence.tsv` lists all evidence used for creating SARG+ extension.
- **Counts:** `misc/sarg.txt`, `misc/sarg_ref.txt`, and `misc/sarg_ext.txt` display the counts of different ARGs for SARG+, SARG+ reference, and SARG+ extension, respectively.
- **FASTA Files:**
  - `sarg.fa`: The combined reference and extension database.
  - `sarg_ref.fa`: The reference component of SARG+.
  - `sarg_ext.fa`: The extension component of SARG+.

## Contribution

New ARG sequences can be integrated into SARG+ by editing [`sarg.json`](https://github.com/xinehc/sarg-curation/blob/master/sarg.json) and [`reference/reference.fasta`](https://github.com/xinehc/sarg-curation/blob/master/reference/reference.fasta).

- **`sarg.json`**: Specifies the ARG type (class/family) for each new gene, relevant literature, rationales, and links to sequences.
- **`reference.fasta`**: Contains the protein sequences of these new genes.

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

If you spot any suspicious cases or wish to add sequences to SARG+, please consider creating a pull request by editing `sarg.json` and `reference/reference.fasta`, or opening an issue. We will update SARG+ regularly.

## FAQ

### Does SARG+ cover point mutations?

No, we exclude all ARGs related to point mutations in essential genes (mainly antibiotic targets). Examples include mutations in ***gyrA***, ***parC***, and ***rpoB***.

### Why doesn't SARG+ include detailed gene numbers like ***bla***<sub>**OXA-1**</sub>?

We group highly similar ARG subtypes (genes) into clusters to reduce the chance of false identifications. For instance, ***bla***<sub>**OXA-1**</sub> and ***bla***<sub>**OXA-1024**</sub> differ by a single amino acid, and this subtle difference can be difficult to detect using reads. By default, we apply 95% identity and 95% query/subject cover as cutoffs for subtype clustering.

### Why does SARG+ exclude fused genes?

SARG+ is designed for read-based ARG profiling. Fused genes can create ambiguities when being aligned using reads. For example, [WP_071593228.1](https://www.ncbi.nlm.nih.gov/protein/WP_071593228.1) (***catB/aac(6')-I***) likely results from the fusion of [WP_264840997.1](https://www.ncbi.nlm.nih.gov/protein/WP_264840997.1) (***catB***) and [WP_033917551.1](https://www.ncbi.nlm.nih.gov/protein/WP_033917551.1) (***aac(6')-I***). Reads, especially short ones, may not reliably distinguish between these genes, potentially leading to false identifications.

### Which other genes are omitted?

We omit regulators (e.g., activators, repressors) since they do not confer direct antibiotic resistance (with the exceptions of ***tipA*** and ***albAB***, which function as self-regulated sequesters). We also remove genes that are putatively accessory, such as ***vanZ***.

### Why some genes have uncommon names?

We standardize gene names to ensure all of them are identifiable through a combination of ARG type and subtype. For example:

- ***mdtP*** refers to both a component of RND efflux pumps in *Escherichia coli* and an MFS transporter in *Bacillus subtilis*. To avoid confusion, we rename the *Bacillus* version to ***mdt(P)***, aligning it with ***mdt(A)*** (also an MFS transporter).
- Some genes are mislabeled by RefSeq; for example, ***efrCD*** is misspelled as ***erfCD***. We correct such cases.
- ***qacA*** and ***qacB*** are renamed to ***qacA/B*** due to their high sequence similarity.

Maintaining naming consistency is a key goal of SARG+, so we modify some gene names accordingly. For instance:

- ***tnrB2*** and ***tnrB3*** are renamed to ***tnrB-1*** and ***tnrB-2*** to reflect their two-component nature as ABC transporters.
- ***cap21*** refers to ***orf21*** of a biosynthetic gene cluster ([AB476988](https://www.ncbi.nlm.nih.gov/nuccore/AB476988)) in *Streptomyces griseus*, which lacks a proper name.