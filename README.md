# Downloading and parsing GenBank files from Python

## Installation
1. Fork git repo into local machine (click on fork) and clone, or simply clone main branch with
```
git clone https://github.com/Robaina/GenBankpy.git
```
2. CD to project directory and set conda environment if not already set:
```
conda env create -n ncbi -f environment.yml
```

3. Activate environment:
```
conda activate ncbi
```


```python
# conda activate ncbi
from genbankpy.parser import GenBankFastaWriter, GBK

"""
This package requires:

pip install ncbi-acc-download
"""

# First we need to define the NCBI entry ids to download the data
entry_ids = [
    'AE001863.1',
    'AF000579.1',
    'AF242489.1', 
    'AP003593.1', 
    'NC_000911.1',
    'NC_007288.1'
]
gbkwriter = GenBankFastaWriter.fromAccessionIDs(entry_ids=entry_ids)
# gbkwriter = GenBankFastaWriter.fromGBKdirectory('gbk_data')
```

    Downloading GenBank files
    Downloading entry: NC_007288.1 (6 / 6) (5 / 6)


```python
# Write fasta containing all peptide sequences of these two organisms
gbkwriter.writeSequencesInFasta(
    gene_keywords={'product': ['any']},
    output_fasta='results/allPeptides.faa', 
    sequence='protein',
    entry_ids=['AE001863.1', 'AP003593.1']
)

# Write fasta containing all nucleotide sequences of these two organisms
gbkwriter.writeSequencesInFasta(
    gene_keywords={'product': ['any']},
    output_fasta='results/allNucleotides.fasta', 
    sequence='nucleotide',
    entry_ids=['AE001863.1', 'AP003593.1']
)

# Write fasta containing nucleotide sequences of the two organisms corresponding to Urease alpha
gbkwriter.writeSequencesInFasta(
    gene_keywords={'product': ['urease', 'alpha']},
    output_fasta='results/ureC.fasta', 
    sequence='nucleotide'
)

# Write fasta containing peptide sequences of the two organisms corresponding to Urease alpha
gbkwriter.writeSequencesInFasta(
    gene_keywords={'product': ['urease', 'alpha']},
    output_fasta='results/ureC.faa', 
    sequence='protein',
    entry_ids=['AE001863.1', 'AP003593.1']
)

# Write fasta containing nucleotide sequences of all five corresponding to 16S
gbkwriter.writeSequencesInFasta(
    gene_keywords={'product': ['16S']},
    output_fasta='results/16s.fasta', 
    sequence='nucleotide',
    entry_ids=None
)
```

# Initializing from list of species names


```python
sp_list = ['Emiliania huxleyi']

gbkwriter = GenBankFastaWriter.fromSpecies(species_list=sp_list,
                                           only_representatives=True,
                                           data_dir='emiliania_data')


gbkwriter.writeSequencesInFasta(
    gene_keywords={'product': ['any']},
    output_fasta='results/allPeptidesEmiliania.faa', 
    sequence='protein'
)
```

# Parsing GenBank files


```python
gbk = GBK('gbk_data/AE001863.1.gbk')
```


```python
gbk.cds.get_by_gene_id('DRA0303')
```




    [SeqFeature(FeatureLocation(ExactPosition(113558), ExactPosition(113924), strand=-1), type='CDS')]


