# Downloading and parsing GenBank files


```python
import pandas as pd
from genbankpy.parser import GenBankFastaWriter, GBK

# entry_ids = pd.read_csv('tests/entry_ids.txt', sep='\t', header=None).iloc[:, 1].values.tolist()
entry_ids = [
    'AE001863.1',
    'AF000579.1',
    'AF242489.1', 
    'AP003593.1', 
    'NC_000911.1'
]
gbkwriter = GenBankFastaWriter.fromAccessionIDs(entry_ids=entry_ids)
# gbkwriter = GenBankFastaWriter.fromGBKdirectory('gbk_data')
```

    Downloading GenBank files



```python
gbkwriter.writeSequencesInFasta(
    gene_keywords={'product': ['urease', 'alpha']},
    output_fasta='results/ureC.fasta', 
    sequence='nucleotide'
)

gbkwriter.writeSequencesInFasta(
    gene_keywords={'product': ['urease', 'alpha']},
    output_fasta='results/ureC.faa', 
    sequence='protein',
    entry_ids=['AE001863.1', 'AP003593.1']
)
```


```python
gbkwriter.writeSequencesInFasta(
    gene_keywords={'product': ['16S']},
    output_fasta='results/16s.fasta', 
    sequence='nucleotide'
)
```

# Parsing GenBank files


```python
gbk = GBK('gbk_data/AE001863.1.gbk')
```


```python
gbk.cds.get_by_keywords({'product': ['urease']})
```




    SeqFeature(FeatureLocation(ExactPosition(125128), ExactPosition(125965), strand=-1), type='CDS')




```python
gbk.cds.get_by_gene_id('DRA0303')
```




    SeqFeature(FeatureLocation(ExactPosition(113558), ExactPosition(113924), strand=-1), type='CDS')


