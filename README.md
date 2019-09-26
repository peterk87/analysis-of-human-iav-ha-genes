# Download, parse, extract Ontario human IAV from NCBI Influenza DB

- Download all nucleotide sequences in the NCBI Influenza DB 
- Parse Human HA gene sequences from Ontario beloning to H1 or H3
- Group sequences by subtype and flu season
- Perform multiple sequence alignments (MSA)



```python
import re
from pathlib import Path

import pandas as pd
import numpy as np

from Bio import SeqIO
from Bio import Entrez
```

Email is required to use NCBI Entrez API


```python
Entrez.email = 'peter.kruczkiewicz@canada.ca'
```

Downloaded all Influenza nucleotide sequences from the NCBI FTP site for the NCBI Influenza DB (2019-09-26T08:51:00+5):

```bash
wget ftp://ftp.ncbi.nih.gov/genomes/INFLUENZA/influenza.fna.gz
```

Output FASTA headers for H1 or H3 HA genes (`(segment 4|\(HA\)|hemagglutinin)`) Influenza A from Ontario (`(Ontario|Toronto|Canada-ON)`) to `ontario-H1-H3-HA-gene-seg4-or-HA.txt`

```bash
zcat influenza.fna.gz | grep -P "^>.*A.*/(Ontario|Toronto|Canada-ON)/.*\(H(1|3)N\w\).*(segment 4|\(HA\)|hemagglutinin)" > ontario-H1-H3-HA-gene-seg4-or-HA.txt
```

Extract GI/NCBI ID from FASTA headers with regex


```python
header = '>gi|52078172|gb|AY619969|Influenza A virus (A/swine/Ontario/K01477/01(H3N3)) hemagglutinin (HA) gene, complete cds'
```

Check that the number after `>gi|` can be parsed


```python
REGEX_GI = re.compile(r'>gi\|(\d+)')
m = REGEX_GI.match(header)
if m:
    print(m.group(1))
```

    52078172


Parse all GIs


```python
gis = []
with open('ontario-H1-H3-HA-gene-seg4-or-HA.txt') as f:
    for l in f:
        m = REGEX_GI.match(l)
        if m:
            gis.append(m.group(1))
```

Should have 892 unique GIs


```python
assert len(set(gis)) == len(gis)
len(gis)
```




    892



Fetch Genbank entries for each of the 892 GI/NCBI IDs 


```python
with Entrez.efetch(db='nucleotide',
                   id=gis,
                   rettype='gb',
                   retmode='text') as efetch_handle:
    entrez_gb_recs = [x for x in SeqIO.parse(efetch_handle, 'genbank')]
```

Should have retrieved 892 Genbank records


```python
len(entrez_gb_recs)
```




    892



Peek into one of the Genbank entries


```python
gb_rec = entrez_gb_recs[0]
```


```python
gb_rec
```




    SeqRecord(seq=Seq('ATGAAGACCATTATTGTTCTGAGTTGTTTTTTCTGTCTGGCTTTCAGCCAAAAT...TGA', IUPACAmbiguousDNA()), id='AY619969.1', name='AY619969', description='Influenza A virus (A/swine/Ontario/K01477/01(H3N3)) hemagglutinin (HA) gene, complete cds', dbxrefs=[])




```python
gb_rec.name
```




    'AY619969'




```python
gb_rec.description
```




    'Influenza A virus (A/swine/Ontario/K01477/01(H3N3)) hemagglutinin (HA) gene, complete cds'




```python
gb_rec.annotations
```




    {'molecule_type': 'RNA',
     'topology': 'linear',
     'data_file_division': 'VRL',
     'date': '28-DEC-2004',
     'accessions': ['AY619969'],
     'sequence_version': 1,
     'keywords': [''],
     'source': 'Influenza A virus (A/swine/Ontario/K01477/01(H3N3))',
     'organism': 'Influenza A virus (A/swine/Ontario/K01477/01(H3N3))',
     'taxonomy': ['Viruses',
      'ssRNA viruses',
      'ssRNA negative-strand viruses',
      'Orthomyxoviridae',
      'Influenzavirus A'],
     'references': [Reference(title='Characterization of avian H3N3 and H1N1 influenza A viruses isolated from pigs in Canada', ...),
      Reference(title='Direct Submission', ...)]}




```python
gb_rec.features
```




    [SeqFeature(FeatureLocation(ExactPosition(0), ExactPosition(1701), strand=1), type='source'),
     SeqFeature(FeatureLocation(ExactPosition(0), ExactPosition(1701), strand=1), type='gene'),
     SeqFeature(FeatureLocation(ExactPosition(0), ExactPosition(1701), strand=1), type='CDS')]




```python
gb_rec.features[0].qualifiers
```




    OrderedDict([('organism',
                  ['Influenza A virus (A/swine/Ontario/K01477/01(H3N3))']),
                 ('mol_type', ['genomic RNA']),
                 ('strain', ['A/swine/Ontario/K01477/01']),
                 ('serotype', ['H3N3']),
                 ('db_xref', ['taxon:292589'])])



Check that all Genbank records have a source feature as the first sequence feature


```python
all([x.features[0].type == 'source' for x in entrez_gb_recs])
```




    True



Parse metadata from Genbank source sequence feature qualifiers


```python
def parse_gb_md(rec):
    source_feature = rec.features[0]
    assert source_feature.type == 'source', rec
    return {k:v[0] for k,v in rec.features[0].qualifiers.items()}
```


```python
def genbank_md(rec): 
    out = parse_gb_md(rec)
    out['accession'] = rec.id
    return out
```

Convert list of metadata dicts to Pandas DataFrame


```python
df_gb_md_892 = pd.DataFrame([genbank_md(x) for x in entrez_gb_recs])
```

Parse/coerce `collection_date` values into standard DateTime format (`pd.Timestamp`)


```python
dates = pd.to_datetime(df_gb_md_892.collection_date, errors='coerce')
years = [x.year if not pd.isnull(x) else None for x in dates]
df_gb_md_892['collection_year'] = years
df_gb_md_892.collection_date = [str(x).split()[0] if not pd.isnull(x) else None for x in dates]
```

Set values for collection year, month and day in DataFrame


```python
df_gb_md_892['collection_month'] = dates.dt.month
df_gb_md_892['collection_day'] = dates.dt.day
```


```python
dts = pd.to_datetime(df_gb_md_892.collection_date)
```

### Compute flu season

Come up with a function to output flu season value where the collection date is used to derive the flu season based on the following:

- if month is August (8) or less, then flu season is `{year - 1}-{year}`
- otherwise, flu season is `{year}-{year + 1}`

For example,
- a collection date of 2018-05-01 should have a flu season of 2017-2018
- a collection date of 2016-10-01 should have a flu season of 2016-2017


```python
dt = dts.loc[891]
```


```python
dt
```




    Timestamp('2018-05-01 00:00:00')




```python
dt.month
```




    5




```python
def flu_season(dt):
    m = dt.month
    y = dt.year
    return f'{y}-{y+1}' if m > 8 else f'{y-1}-{y}'
```


```python
flu_season(dt)
```




    '2017-2018'



Add flu season for each non-null collection date in the dataframe


```python
flu_seasons = [flu_season(dt) if not pd.isnull(dt) else None for dt in dts]
```


```python
df_gb_md_892['flu_season'] = flu_seasons
```

Sort by collection date in descending order (most recent to least)


```python
df_gb_md_892.sort_values('collection_date', inplace=True, ascending=False)
```

Peek at DataFrame of sequence metadata


```python
df_gb_md_892
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>organism</th>
      <th>mol_type</th>
      <th>strain</th>
      <th>serotype</th>
      <th>db_xref</th>
      <th>accession</th>
      <th>country</th>
      <th>host</th>
      <th>segment</th>
      <th>collection_date</th>
      <th>isolation_source</th>
      <th>isolate</th>
      <th>note</th>
      <th>lab_host</th>
      <th>lat_lon</th>
      <th>collection_year</th>
      <th>collection_month</th>
      <th>collection_day</th>
      <th>flu_season</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>891</td>
      <td>Influenza A virus</td>
      <td>viral cRNA</td>
      <td>A/swine/Ontario/SD0298/2018</td>
      <td>H3N2</td>
      <td>taxon:11320</td>
      <td>MK462790.1</td>
      <td>Canada: Ontario</td>
      <td>swine</td>
      <td>4</td>
      <td>2018-05-01</td>
      <td>MDCK cells</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>2018.0</td>
      <td>5.0</td>
      <td>1.0</td>
      <td>2017-2018</td>
    </tr>
    <tr>
      <td>744</td>
      <td>Influenza A virus</td>
      <td>viral cRNA</td>
      <td>A/Ontario/026/2018</td>
      <td>H3N2</td>
      <td>taxon:11320</td>
      <td>MG889769.1</td>
      <td>Canada: Ontario</td>
      <td>Homo sapiens</td>
      <td>4</td>
      <td>2018-01-11</td>
      <td>nasopharyngeal swab</td>
      <td>026</td>
      <td>original specimen</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>2018.0</td>
      <td>1.0</td>
      <td>11.0</td>
      <td>2017-2018</td>
    </tr>
    <tr>
      <td>752</td>
      <td>Influenza A virus</td>
      <td>viral cRNA</td>
      <td>A/Ontario/034/2018</td>
      <td>H3N2</td>
      <td>taxon:11320</td>
      <td>MG889777.1</td>
      <td>Canada: Ontario</td>
      <td>Homo sapiens</td>
      <td>4</td>
      <td>2018-01-11</td>
      <td>nasopharyngeal swab</td>
      <td>034</td>
      <td>original specimen</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>2018.0</td>
      <td>1.0</td>
      <td>11.0</td>
      <td>2017-2018</td>
    </tr>
    <tr>
      <td>751</td>
      <td>Influenza A virus</td>
      <td>viral cRNA</td>
      <td>A/Ontario/033/2018</td>
      <td>H3N2</td>
      <td>taxon:11320</td>
      <td>MG889776.1</td>
      <td>Canada: Ontario</td>
      <td>Homo sapiens</td>
      <td>4</td>
      <td>2018-01-11</td>
      <td>nasopharyngeal swab</td>
      <td>033</td>
      <td>original specimen</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>2018.0</td>
      <td>1.0</td>
      <td>11.0</td>
      <td>2017-2018</td>
    </tr>
    <tr>
      <td>754</td>
      <td>Influenza A virus</td>
      <td>viral cRNA</td>
      <td>A/Ontario/038/2018</td>
      <td>H3N2</td>
      <td>taxon:11320</td>
      <td>MG889779.1</td>
      <td>Canada: Ontario</td>
      <td>Homo sapiens</td>
      <td>4</td>
      <td>2018-01-10</td>
      <td>nasopharyngeal swab</td>
      <td>038</td>
      <td>original specimen</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>2018.0</td>
      <td>1.0</td>
      <td>10.0</td>
      <td>2017-2018</td>
    </tr>
    <tr>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <td>11</td>
      <td>Influenza A virus (A/Ontario/RV1273/2005(H3N2))</td>
      <td>genomic DNA</td>
      <td>A/Ontario/RV1273/2005</td>
      <td>H3N2</td>
      <td>taxon:381529</td>
      <td>DQ469962.1</td>
      <td>Canada</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>None</td>
    </tr>
    <tr>
      <td>12</td>
      <td>Influenza A virus (A/swine/Ontario/33853/2005(...</td>
      <td>genomic DNA</td>
      <td>A/swine/Ontario/33853/2005</td>
      <td>H3N2</td>
      <td>taxon:381533</td>
      <td>DQ469994.1</td>
      <td>Canada</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>None</td>
    </tr>
    <tr>
      <td>13</td>
      <td>Influenza A virus (A/turkey/Ontario/31232/2005...</td>
      <td>genomic DNA</td>
      <td>A/turkey/Ontario/31232/2005</td>
      <td>H3N2</td>
      <td>taxon:381534</td>
      <td>DQ470002.1</td>
      <td>Canada</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>None</td>
    </tr>
    <tr>
      <td>14</td>
      <td>Influenza A virus (A/swine/Ontario/00130/97(H3...</td>
      <td>genomic RNA</td>
      <td>A/swine/Ontario/00130/97</td>
      <td>NaN</td>
      <td>taxon:133777</td>
      <td>AF251395.2</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>None</td>
    </tr>
    <tr>
      <td>21</td>
      <td>Influenza A virus (A/Ontario/1252/2007(H3N2))</td>
      <td>viral cRNA</td>
      <td>A/Ontario/1252/2007</td>
      <td>H3N2</td>
      <td>taxon:496427</td>
      <td>EU399751.1</td>
      <td>Canada</td>
      <td>NaN</td>
      <td>4</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>None</td>
    </tr>
  </tbody>
</table>
<p>892 rows × 19 columns</p>
</div>



### Save entire dataframe to CSV file


```python
df_gb_md_892.to_csv('2019-09-26-McLaughlin-NCBI-Influenza-DB-H3-or-H1-Ontario-metadata-from-genbank.csv', index=False)
```

### Filter for human-derived samples

Frequencies of distinct host values


```python
df_gb_md_892.host.value_counts()
```




    Homo sapiens                       698
    swine                               27
    turkey                              10
    human; gender M; age 25              6
    Swine                                6
                                      ... 
    Homo sapiens; gender M; age 84Y      1
    Homo sapiens; gender M; age 9Y       1
    human; gender M; age 0               1
    Homo sapiens; gender F; age 28Y      1
    human; gender M; age 55              1
    Name: host, Length: 109, dtype: int64



IAV with host info not provided in GenBank file source feature


```python
df_gb_md_892[pd.isnull(df_gb_md_892.host)]
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>organism</th>
      <th>mol_type</th>
      <th>strain</th>
      <th>serotype</th>
      <th>db_xref</th>
      <th>accession</th>
      <th>country</th>
      <th>host</th>
      <th>segment</th>
      <th>collection_date</th>
      <th>isolation_source</th>
      <th>isolate</th>
      <th>note</th>
      <th>lab_host</th>
      <th>lat_lon</th>
      <th>collection_year</th>
      <th>collection_month</th>
      <th>collection_day</th>
      <th>flu_season</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>0</td>
      <td>Influenza A virus (A/swine/Ontario/K01477/01(H...</td>
      <td>genomic RNA</td>
      <td>A/swine/Ontario/K01477/01</td>
      <td>H3N3</td>
      <td>taxon:292589</td>
      <td>AY619969.1</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>None</td>
    </tr>
    <tr>
      <td>1</td>
      <td>Influenza A virus (A/swine/Ontario/42729A/01(H...</td>
      <td>genomic RNA</td>
      <td>A/swine/Ontario/42729A/01</td>
      <td>H3N3</td>
      <td>taxon:292590</td>
      <td>AY619977.1</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>None</td>
    </tr>
    <tr>
      <td>2</td>
      <td>Influenza A virus (A/swine/Ontario/Biovet1/05(...</td>
      <td>genomic RNA</td>
      <td>A/swine/Ontario/Biovet1/05</td>
      <td>H3N2</td>
      <td>taxon:354557</td>
      <td>DQ241762.1</td>
      <td>Canada</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>None</td>
    </tr>
    <tr>
      <td>3</td>
      <td>Influenza A virus (A/swine/Ontario/57561/03(H1...</td>
      <td>genomic RNA</td>
      <td>A/swine/Ontario/57561/03</td>
      <td>H1N1</td>
      <td>taxon:358575</td>
      <td>DQ280195.1</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>None</td>
    </tr>
    <tr>
      <td>4</td>
      <td>Influenza A virus (A/swine/Ontario/55383/04(H1...</td>
      <td>genomic RNA</td>
      <td>A/swine/Ontario/55383/04</td>
      <td>H1N2</td>
      <td>taxon:358580</td>
      <td>DQ280212.1</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>None</td>
    </tr>
    <tr>
      <td>5</td>
      <td>Influenza A virus (A/swine/Ontario/53518/03(H1...</td>
      <td>genomic RNA</td>
      <td>A/swine/Ontario/53518/03</td>
      <td>H1N1</td>
      <td>taxon:358577</td>
      <td>DQ280219.1</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>None</td>
    </tr>
    <tr>
      <td>6</td>
      <td>Influenza A virus (A/swine/Ontario/52156/03(H1...</td>
      <td>genomic RNA</td>
      <td>A/swine/Ontario/52156/03</td>
      <td>H1N2</td>
      <td>taxon:358581</td>
      <td>DQ280227.1</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>None</td>
    </tr>
    <tr>
      <td>7</td>
      <td>Influenza A virus (A/swine/Ontario/48235/04(H1...</td>
      <td>genomic RNA</td>
      <td>A/swine/Ontario/48235/04</td>
      <td>H1N2</td>
      <td>taxon:358582</td>
      <td>DQ280236.1</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>None</td>
    </tr>
    <tr>
      <td>8</td>
      <td>Influenza A virus (A/swine/Ontario/23866/04(H1...</td>
      <td>genomic RNA</td>
      <td>A/swine/Ontario/23866/04</td>
      <td>H1N1</td>
      <td>taxon:358578</td>
      <td>DQ280243.1</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>None</td>
    </tr>
    <tr>
      <td>9</td>
      <td>Influenza A virus (A/swine/Ontario/11112/2004(...</td>
      <td>genomic RNA</td>
      <td>A/swine/Ontario/11112/04</td>
      <td>H1N1</td>
      <td>taxon:358579</td>
      <td>DQ280250.1</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>None</td>
    </tr>
    <tr>
      <td>10</td>
      <td>Influenza A virus (A/duck/Ontario/05/00 (H3N2))</td>
      <td>genomic RNA</td>
      <td>(A/duck/Ontario/05/00 (H3N2)</td>
      <td>NaN</td>
      <td>taxon:273361</td>
      <td>AJ697864.1</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>None</td>
    </tr>
    <tr>
      <td>11</td>
      <td>Influenza A virus (A/Ontario/RV1273/2005(H3N2))</td>
      <td>genomic DNA</td>
      <td>A/Ontario/RV1273/2005</td>
      <td>H3N2</td>
      <td>taxon:381529</td>
      <td>DQ469962.1</td>
      <td>Canada</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>None</td>
    </tr>
    <tr>
      <td>12</td>
      <td>Influenza A virus (A/swine/Ontario/33853/2005(...</td>
      <td>genomic DNA</td>
      <td>A/swine/Ontario/33853/2005</td>
      <td>H3N2</td>
      <td>taxon:381533</td>
      <td>DQ469994.1</td>
      <td>Canada</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>None</td>
    </tr>
    <tr>
      <td>13</td>
      <td>Influenza A virus (A/turkey/Ontario/31232/2005...</td>
      <td>genomic DNA</td>
      <td>A/turkey/Ontario/31232/2005</td>
      <td>H3N2</td>
      <td>taxon:381534</td>
      <td>DQ470002.1</td>
      <td>Canada</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>None</td>
    </tr>
    <tr>
      <td>14</td>
      <td>Influenza A virus (A/swine/Ontario/00130/97(H3...</td>
      <td>genomic RNA</td>
      <td>A/swine/Ontario/00130/97</td>
      <td>NaN</td>
      <td>taxon:133777</td>
      <td>AF251395.2</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>None</td>
    </tr>
    <tr>
      <td>21</td>
      <td>Influenza A virus (A/Ontario/1252/2007(H3N2))</td>
      <td>viral cRNA</td>
      <td>A/Ontario/1252/2007</td>
      <td>H3N2</td>
      <td>taxon:496427</td>
      <td>EU399751.1</td>
      <td>Canada</td>
      <td>NaN</td>
      <td>4</td>
      <td>None</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>None</td>
    </tr>
  </tbody>
</table>
</div>




```python
human_regex_pattern = r'.*([Hh]uman|[Hh]omo).*'
```


```python
df_gb_md_892.host.str.match(human_regex_pattern).sum()
```




    828



What does the metadata for the non-human regex pattern matching sequences look like?


```python
for i,r in df_gb_md_892[~df_gb_md_892.host.str.match(human_regex_pattern, na=False)].iterrows():
    print(f'{r.strain: <40} {r.host: <20} {r.collection_date}')
```

    A/swine/Ontario/SD0298/2018              swine                2018-05-01
    A/swine/Ontario/DM_21/2017               swine                2017-04-26
    A/swine/Ontario/DM_11/2017               swine                2017-01-16
    A/turkey/Ontario/FAV-006-4/2016          turkey               2016-04-07
    A/turkey/Ontario/FAV-006-10/2016         turkey               2016-04-07
    A/turkey/Ontario/FAV-005-2/2016          turkey               2016-04-05
    A/swine/Ontario/G3/2014                  swine                2014-11-07
    A/swine/Ontario/G10/2014                 swine                2014-09-04
    A/swine/Ontario/G12/2014                 swine                2014-09-03
    A/swine/Ontario/G13/2014                 swine                2014-09-03
    A/swine/Ontario/G14/2014                 swine                2014-08-12
    A/swine/Ontario/G16/2014                 swine                2014-07-24
    A/swine/Ontario/G11/2014                 swine                2014-07-15
    A/swine/Ontario/G15/2014                 swine                2014-06-18
    A/swine/Ontario/118-38/2012              swine                2012-12-06
    A/swine/Ontario/120-55/2012              swine                2012-11-28
    A/swine/Ontario/16-24/2012               swine                2012-11-16
    A/swine/Ontario/115-2/2012               swine                2012-10-20
    A/swine/Ontario/84/2012                  swine                2012-10-10
    A/swine/Ontario/68/2012                  swine                2012-10-10
    A/swine/Ontario/114-13/2012              swine                2012-10-09
    A/swine/Ontario/13-1/2012                swine                2012-10-05
    A/swine/Ontario/204-76/2012              swine                2012-08-22
    A/swine/Ontario/46/2012                  swine                2012-07-05
    A/swine/Ontario/62/2012                  swine                2012-07-03
    A/swine/Ontario/107-22/2012              swine                2012-04-15
    A/swine/Ontario/105-56/2012              swine                2012-02-25
    A/swine/Ontario/104-25/2012              swine                2012-01-13
    A/swine/Ontario/103-18/2011              swine                2011-11-23
    A/swine/Ontario/11-105317/2011           swine                2011-11-23
    A/turkey/Ontario/FAV-10/2011             turkey               2011-07-21
    A/turkey/Ontario/FAV-3/2011              turkey               2011-07-06
    A/turkey/Ontario/FAV-9/2011              turkey               2011-02-24
    A/turkey/Ontario/FAV117-1C/2009          turkey               2009-12-07
    A/turkey/Ontario/FAV114-17/2009          turkey               2009-11-03
    A/turkey/Ontario/FAV110-4/2009           turkey               2009-11-03
    A/turkey/Ontario/FAV110/2009             turkey               2009-10-23
    A/mallard/Ontario/15873/2005             mallard; gender F; age Hatch year 2005-08-29
    A/ring-necked duck/Ontario/15530/2005    ring-necked Duck; gender M; age Hatch year 2005-08-26
    A/equine/Ontario/40754/2002              horse                2002-01-01
    A/equine/Ontario/51610/2001              horse                2001-01-01
    A/equine/Ontario/V-50/1997               horse                1997-01-01
    A/swine/Ontario/7/1981                   Swine                1981-01-01
    A/swine/Ontario/6/1981                   Swine                1981-01-01
    A/swine/Ontario/3/1981                   Swine                1981-01-01
    A/swine/Ontario/2/1981                   Swine                1981-01-01
    A/swine/Ontario/1/1981                   Swine                1981-01-01
    A/swine/Ontario/4/1981                   Swine                1981-01-01
    A/swine/Ontario/K01477/01                nan                  None
    A/swine/Ontario/42729A/01                nan                  None
    A/swine/Ontario/Biovet1/05               nan                  None
    A/swine/Ontario/57561/03                 nan                  None
    A/swine/Ontario/55383/04                 nan                  None
    A/swine/Ontario/53518/03                 nan                  None
    A/swine/Ontario/52156/03                 nan                  None
    A/swine/Ontario/48235/04                 nan                  None
    A/swine/Ontario/23866/04                 nan                  None
    A/swine/Ontario/11112/04                 nan                  None
    (A/duck/Ontario/05/00 (H3N2)             nan                  None
    A/Ontario/RV1273/2005                    nan                  None
    A/swine/Ontario/33853/2005               nan                  None
    A/turkey/Ontario/31232/2005              nan                  None
    A/swine/Ontario/00130/97                 nan                  None
    A/Ontario/1252/2007                      nan                  None



```python
df_human = df_gb_md_892[df_gb_md_892.host.str.match(human_regex_pattern, na=False)]
```


```python
df_human.host.value_counts().to_dict()
```




    {'Homo sapiens': 698,
     'human; gender M; age 25': 6,
     'Homo sapiens; gender F; age 22': 4,
     'Homo sapiens; gender M; age 21Y': 3,
     'Homo sapiens; gender M; age 26Y': 3,
     'human; gender F; age 39': 2,
     'human; gender F; age 25': 2,
     'Homo sapiens; gender M; age 83Y': 2,
     'human; gender F; age 14': 2,
     'human; gender M; age 17': 2,
     'Homo sapiens; gender M; age 79Y': 2,
     'Homo sapiens; gender M; age 95Y': 2,
     'human; gender F; age 33': 2,
     'Homo sapiens; gender M; age 86Y': 2,
     'human; gender F; age 5': 2,
     'Homo sapiens; gender M; age 73Y': 2,
     'Homo sapiens; gender M; age 76Y': 2,
     'human; gender F; age 24': 2,
     'Homo sapiens; gender F; age 91Y': 2,
     'Homo sapiens; gender M; age 2Y': 2,
     'Homo sapiens; gender M; age 4Y': 2,
     'Homo sapiens; gender M; age 19Y': 1,
     'human; gender M; age 5': 1,
     'Homo sapiens; gender M; age 32Y': 1,
     'Homo sapiens; gender M; age 11mo': 1,
     'human; gender M; age 14': 1,
     'human; gender F; age 88': 1,
     'Homo sapiens; gender F; age 21Y': 1,
     'human; gender M; age 72': 1,
     'Homo sapiens; gender F; age 96Y': 1,
     'Homo sapiens; DOB 01-Mar-51; M': 1,
     'Homo sapiens; DOB 05-Feb-02; M': 1,
     'Homo sapiens; gender F; age 83Y': 1,
     'Homo sapiens; gender F; age 27Y': 1,
     'Homo sapiens; gender F; age 76Y': 1,
     'Homo sapiens; gender M; age 27Y': 1,
     'Homo sapiens; gender M; age 89Y': 1,
     'Homo sapiens; gender M; age 29': 1,
     'human; gender M; age 70': 1,
     'Homo sapiens; gender F; age 18': 1,
     'Homo sapiens; gender M; age 16Y': 1,
     'Homo sapiens; DOB 25-Jul-96; M': 1,
     'human; gender F; age 1': 1,
     'Homo sapiens; gender F; age 85Y': 1,
     'human; gender M; age 50': 1,
     'Homo sapiens; gender F; age 51Y': 1,
     'human; gender M; age 42': 1,
     'Homo sapiens; gender F; age 24': 1,
     'human; gender M; age 20': 1,
     'Homo sapiens; gender M; age 78Y': 1,
     'Homo sapiens; gender M; age 37Y': 1,
     'Homo sapiens; gender F; age 7mo': 1,
     'Homo sapiens; gender F; age 26Y': 1,
     'human; gender M; age 56': 1,
     'Homo sapiens; gender F; age 73Y': 1,
     'Homo sapiens; gender F; age 90Y': 1,
     'human; gender M; age 13': 1,
     'Homo sapiens; gender female; age 18': 1,
     'human; gender F; age 2M': 1,
     'human; gender F; age 41': 1,
     'human; gender F; age 44': 1,
     'human; gender M; age 9': 1,
     'Homo sapiens; gender F; age 54Y': 1,
     'human; gender M; age 15': 1,
     'Homo sapiens; gender F; age 88Y': 1,
     'Homo sapiens; gender M; age 2M': 1,
     'human; gender M; age 6': 1,
     'Homo sapiens; DOB 02-Dec-06; M': 1,
     'Homo sapiens; gender F; age 15Y': 1,
     'Homo sapiens; gender F; age 17Y': 1,
     'human; gender M; age 1': 1,
     'Homo sapiens; gender M; age 5Y': 1,
     'Homo sapiens; gender F; age 7Y': 1,
     'human; gender M; age 7': 1,
     'Homo sapiens; gender F; age 22Y': 1,
     'Homo sapiens; gender F; age 3Y': 1,
     'Homo sapiens; DOB 13-Oct-82; F': 1,
     'Homo sapiens; gender F; age 1Y': 1,
     'Homo sapiens; gender M; age 70Y': 1,
     'Homo sapiens; gender F; age 20Y': 1,
     'Homo sapiens; gender F; age 81Y': 1,
     'Homo sapiens; gender M; age 45Y': 1,
     'human; gender F; age 49': 1,
     'Homo sapiens; gender M; age 4M': 1,
     'Homo sapiens; gender M; age 62Y': 1,
     'Homo sapiens; gender M; age 71Y': 1,
     'human; gender M; age 46': 1,
     'human; gender F; age 31': 1,
     'Homo sapiens; gender F; age 51': 1,
     'human; gender F; age 53': 1,
     'Homo sapiens; gender M; age 28': 1,
     'human; gender M; age 55': 1,
     'Homo sapiens; gender F; age 74Y': 1,
     'Homo sapiens; DOB 14-Dec-62; F': 1,
     'Homo sapiens; gender M; age 23Y': 1,
     'Homo sapiens; gender F; age 32Y': 1,
     'Homo sapiens; gender M; age 8Y': 1,
     'Homo sapiens; gender M; age 35Y': 1,
     'Homo sapiens; gender M; age 84Y': 1,
     'Homo sapiens; gender M; age 9Y': 1,
     'human; gender M; age 0': 1,
     'Homo sapiens; gender F; age 28Y': 1,
     'Homo sapiens; gender M; age 18Y': 1}




```python
df_human.to_csv('2019-09-26-McLaughlin-NCBI-Influenza-DB-H3-or-H1-Ontario-Human-metadata-from-genbank.csv', index=False)
```

## Write all Genbank records to a file


```python
SeqIO.write(entrez_gb_recs, '2019-09-26-McLaughlin-NCBI-Influenza-DB-H3-or-H1-Ontario-metadata-from-genbank.gb', 'genbank')
```




    892



## Partition human sequences by flu season


```python
id_to_gb = {x.id:x for x in entrez_gb_recs}
```


```python
df_human.flu_season.value_counts()
```




    2015-2016    210
    2016-2017    163
    2014-2015    135
    2010-2011     67
    2012-2013     57
    2009-2010     54
    2017-2018     50
    2008-2009     40
    2011-2012     30
    2013-2014     16
    2007-2008      6
    Name: flu_season, dtype: int64



Group by `flu_season` and `serotype`


```python
g = df_human.groupby(['flu_season','serotype'])
```

Create lists of accessions grouped by `flu_season` and `serotype`


```python
grouped_acc = g.accession.apply(list)
```


```python
grouped_acc
```




    flu_season  serotype
    2007-2008   H1N1        [FJ800811.1, FJ800819.1, FJ800810.1, FJ800818....
    2008-2009   H1N1        [HQ239567.1, CY060534.2, CY060502.2, CY060526....
    2009-2010   H1N1        [CY081093.2, HQ239460.2, CY060726.2, HQ239459....
                H3N2                     [JQ658889.1, JQ658888.1, CY054550.1]
    2010-2011   H1N1        [CY081101.2, CY081069.2, CY081085.2, CY081077....
                H3N2        [CY111003.1, CY111002.1, CY111001.1, CY111000....
    2011-2012   H1N1        [JX875001.1, KF551079.1, KF551095.1, KF551093....
                H3N2        [KF551077.1, KF551078.1, KF551076.1, KF551075....
    2012-2013   H1N1        [KF886379.1, KF886366.1, KF886376.1, KF886371....
                H3N2        [KJ734749.1, KF886354.1, KJ734748.1, KJ734747....
    2013-2014   H1N1        [KP864399.1, KP864397.1, KP864398.1, KP864396....
                H3N2                     [KP864423.1, KP864424.1, KP864425.1]
    2014-2015   H3N2        [KU729355.1, KU729354.1, KU729353.1, KU729458....
    2015-2016   H1N1        [MF195582.1, MF195581.1, MF195580.1, MF195579....
                H3N2        [MF195739.1, MF195738.1, MF195737.1, MF195736....
    2016-2017   H1N1                                             [MH216446.1]
                H3N2        [MH216444.1, MH216443.1, MH216442.1, MH216440....
    2017-2018   H3N2        [MG889769.1, MG889777.1, MG889776.1, MG889779....
    Name: accession, dtype: object



### Write partitioned human HA sequences to GenBank and FASTA files


```python
!mkdir human-HA-sequences-by-flu-season
```


```python
for (flu_season, serotype), accessions in grouped_acc.items():
    print(f'{serotype}| Flu season: {flu_season} (N={len(accessions)})')
    gbs = [id_to_gb[acc] for acc in accessions]
    filename = f'human-HA-sequences-by-flu-season/HA-{serotype}-flu_season-{flu_season}'
    SeqIO.write(gbs, f'{filename}.gb', 'genbank')
    SeqIO.write(gbs, f'{filename}.fa', 'fasta')
```

    H1N1| Flu season: 2007-2008 (N=6)
    H1N1| Flu season: 2008-2009 (N=40)
    H1N1| Flu season: 2009-2010 (N=51)
    H3N2| Flu season: 2009-2010 (N=3)
    H1N1| Flu season: 2010-2011 (N=5)
    H3N2| Flu season: 2010-2011 (N=62)
    H1N1| Flu season: 2011-2012 (N=18)
    H3N2| Flu season: 2011-2012 (N=12)
    H1N1| Flu season: 2012-2013 (N=18)
    H3N2| Flu season: 2012-2013 (N=39)
    H1N1| Flu season: 2013-2014 (N=13)
    H3N2| Flu season: 2013-2014 (N=3)
    H3N2| Flu season: 2014-2015 (N=135)
    H1N1| Flu season: 2015-2016 (N=198)
    H3N2| Flu season: 2015-2016 (N=12)
    H1N1| Flu season: 2016-2017 (N=1)
    H3N2| Flu season: 2016-2017 (N=162)
    H3N2| Flu season: 2017-2018 (N=50)


Peek into output directory at `human-HA-sequences-by-flu-season/`


```python
!ls -lh human-HA-sequences-by-flu-season/
```

    total 5.4M
    -rw-r--r-- 1 pkruczkiewicz grp_pkruczkiewicz 5.7K Sep 26 11:43 HA-H1N1-flu_season-2007-2008.fa
    -rw-r--r-- 1 pkruczkiewicz grp_pkruczkiewicz  21K Sep 26 11:43 HA-H1N1-flu_season-2007-2008.gb
    -rw-r--r-- 1 pkruczkiewicz grp_pkruczkiewicz  67K Sep 26 11:43 HA-H1N1-flu_season-2008-2009.fa
    -rw-r--r-- 1 pkruczkiewicz grp_pkruczkiewicz 206K Sep 26 11:43 HA-H1N1-flu_season-2008-2009.gb
    -rw-r--r-- 1 pkruczkiewicz grp_pkruczkiewicz  91K Sep 26 11:43 HA-H1N1-flu_season-2009-2010.fa
    -rw-r--r-- 1 pkruczkiewicz grp_pkruczkiewicz 279K Sep 26 11:43 HA-H1N1-flu_season-2009-2010.gb
    -rw-r--r-- 1 pkruczkiewicz grp_pkruczkiewicz 9.2K Sep 26 11:43 HA-H1N1-flu_season-2010-2011.fa
    -rw-r--r-- 1 pkruczkiewicz grp_pkruczkiewicz  27K Sep 26 11:43 HA-H1N1-flu_season-2010-2011.gb
    -rw-r--r-- 1 pkruczkiewicz grp_pkruczkiewicz  33K Sep 26 11:43 HA-H1N1-flu_season-2011-2012.fa
    -rw-r--r-- 1 pkruczkiewicz grp_pkruczkiewicz 103K Sep 26 11:43 HA-H1N1-flu_season-2011-2012.gb
    -rw-r--r-- 1 pkruczkiewicz grp_pkruczkiewicz  33K Sep 26 11:43 HA-H1N1-flu_season-2012-2013.fa
    -rw-r--r-- 1 pkruczkiewicz grp_pkruczkiewicz 101K Sep 26 11:43 HA-H1N1-flu_season-2012-2013.gb
    -rw-r--r-- 1 pkruczkiewicz grp_pkruczkiewicz  24K Sep 26 11:43 HA-H1N1-flu_season-2013-2014.fa
    -rw-r--r-- 1 pkruczkiewicz grp_pkruczkiewicz  75K Sep 26 11:43 HA-H1N1-flu_season-2013-2014.gb
    -rw-r--r-- 1 pkruczkiewicz grp_pkruczkiewicz 206K Sep 26 11:43 HA-H1N1-flu_season-2015-2016.fa
    -rw-r--r-- 1 pkruczkiewicz grp_pkruczkiewicz 851K Sep 26 11:43 HA-H1N1-flu_season-2015-2016.gb
    -rw-r--r-- 1 pkruczkiewicz grp_pkruczkiewicz 1.1K Sep 26 11:43 HA-H1N1-flu_season-2016-2017.fa
    -rw-r--r-- 1 pkruczkiewicz grp_pkruczkiewicz 3.4K Sep 26 11:43 HA-H1N1-flu_season-2016-2017.gb
    -rw-r--r-- 1 pkruczkiewicz grp_pkruczkiewicz 5.2K Sep 26 11:43 HA-H3N2-flu_season-2009-2010.fa
    -rw-r--r-- 1 pkruczkiewicz grp_pkruczkiewicz  15K Sep 26 11:43 HA-H3N2-flu_season-2009-2010.gb
    -rw-r--r-- 1 pkruczkiewicz grp_pkruczkiewicz 112K Sep 26 11:43 HA-H3N2-flu_season-2010-2011.fa
    -rw-r--r-- 1 pkruczkiewicz grp_pkruczkiewicz 313K Sep 26 11:43 HA-H3N2-flu_season-2010-2011.gb
    -rw-r--r-- 1 pkruczkiewicz grp_pkruczkiewicz  22K Sep 26 11:43 HA-H3N2-flu_season-2011-2012.fa
    -rw-r--r-- 1 pkruczkiewicz grp_pkruczkiewicz  68K Sep 26 11:43 HA-H3N2-flu_season-2011-2012.gb
    -rw-r--r-- 1 pkruczkiewicz grp_pkruczkiewicz  71K Sep 26 11:43 HA-H3N2-flu_season-2012-2013.fa
    -rw-r--r-- 1 pkruczkiewicz grp_pkruczkiewicz 213K Sep 26 11:43 HA-H3N2-flu_season-2012-2013.gb
    -rw-r--r-- 1 pkruczkiewicz grp_pkruczkiewicz 5.0K Sep 26 11:43 HA-H3N2-flu_season-2013-2014.fa
    -rw-r--r-- 1 pkruczkiewicz grp_pkruczkiewicz  17K Sep 26 11:43 HA-H3N2-flu_season-2013-2014.gb
    -rw-r--r-- 1 pkruczkiewicz grp_pkruczkiewicz 198K Sep 26 11:43 HA-H3N2-flu_season-2014-2015.fa
    -rw-r--r-- 1 pkruczkiewicz grp_pkruczkiewicz 702K Sep 26 11:43 HA-H3N2-flu_season-2014-2015.gb
    -rw-r--r-- 1 pkruczkiewicz grp_pkruczkiewicz  22K Sep 26 11:43 HA-H3N2-flu_season-2015-2016.fa
    -rw-r--r-- 1 pkruczkiewicz grp_pkruczkiewicz  69K Sep 26 11:43 HA-H3N2-flu_season-2015-2016.gb
    -rw-r--r-- 1 pkruczkiewicz grp_pkruczkiewicz 290K Sep 26 11:43 HA-H3N2-flu_season-2016-2017.fa
    -rw-r--r-- 1 pkruczkiewicz grp_pkruczkiewicz 810K Sep 26 11:43 HA-H3N2-flu_season-2016-2017.gb
    -rw-r--r-- 1 pkruczkiewicz grp_pkruczkiewicz  90K Sep 26 11:43 HA-H3N2-flu_season-2017-2018.fa
    -rw-r--r-- 1 pkruczkiewicz grp_pkruczkiewicz 259K Sep 26 11:43 HA-H3N2-flu_season-2017-2018.gb


### MAFFT multiple sequence alignment (MSA) of partitioned human HA genes

Perform MSA with MAFFT (v7.407) of each set of sequences with the L-INS-i strategy for high accuracy. 


```python
!mafft --version
```

    v7.407 (2018/Jul/23)



```python
!mafft --help
```

    ------------------------------------------------------------------------------
      MAFFT v7.407 (2018/Jul/23)
      https://mafft.cbrc.jp/alignment/software/
      MBE 30:772-780 (2013), NAR 30:3059-3066 (2002)
    ------------------------------------------------------------------------------
    High speed:
      % mafft in > out
      % mafft --retree 1 in > out (fast)
    
    High accuracy (for <~200 sequences x <~2,000 aa/nt):
      % mafft --maxiterate 1000 --localpair  in > out (% linsi in > out is also ok)
      % mafft --maxiterate 1000 --genafpair  in > out (% einsi in > out)
      % mafft --maxiterate 1000 --globalpair in > out (% ginsi in > out)
    
    If unsure which option to use:
      % mafft --auto in > out
    
    --op # :         Gap opening penalty, default: 1.53
    --ep # :         Offset (works like gap extension penalty), default: 0.0
    --maxiterate # : Maximum number of iterative refinement, default: 0
    --clustalout :   Output: clustal format, default: fasta
    --reorder :      Outorder: aligned, default: input order
    --quiet :        Do not report progress
    --thread # :     Number of threads (if unsure, --thread -1)



```python
!mkdir human-HA-sequences-by-flu-season/mafft-msa
```

Use GNU Parallel to run parallel instances of `mafft` MSA


```python
!parallel -v mafft --thread -1 --maxiterate 1000 --localpair {} ">" human-HA-sequences-by-flu-season/mafft-msa/{/} ::: human-HA-sequences-by-flu-season/*.fa
```

    mafft --thread -1 --maxiterate 1000 --localpair human-HA-sequences-by-flu-season/HA-H1N1-flu_season-2016-2017.fa > human-HA-sequences-by-flu-season/mafft-msa/HA-H1N1-flu_season-2016-2017.fa
    OS = linux
    The number of physical cores =  28
    outputhat23=16
    treein = 0
    compacttree = 0
    Warning: Only 1 sequence found.
    minimumweight = 0.000010
    autosubalignment = 0.000000
    nthread = 8
    randomseed = 0
    blosum 62 / kimura 200
    poffset = 0
    niter = 16
    sueff_global = 0.100000
    nadd = 16
    Warning: Only 1 sequence found.
    
    Strategy:
     L-INS-i (Probably most accurate, very slow)
     Iterative refinement method (<16) with LOCAL pairwise alignment information
    
    If unsure which option to use, try 'mafft --auto input > output'.
    For more information, see 'mafft --help', 'mafft --man' and the mafft page.
    
    The default gap scoring scheme has been changed in version 7.110 (2013 Oct).
    It tends to insert more gaps into gap-rich regions than previous versions.
    To disable this change, add the --leavegappyregion option.
    
    mafft --thread -1 --maxiterate 1000 --localpair human-HA-sequences-by-flu-season/HA-H3N2-flu_season-2013-2014.fa > human-HA-sequences-by-flu-season/mafft-msa/HA-H3N2-flu_season-2013-2014.fa
    OS = linux
    The number of physical cores =  28
    outputhat23=16
    treein = 0
    compacttree = 0
    stacksize: 8192 kb
    generating a scoring matrix for nucleotide (dist=200) ... done
    All-to-all alignment.
    tbfast-pair (nuc) Version 7.407
    alg=L, model=DNA200 (2), 2.00 (6.00), -0.10 (-0.30), noshift, amax=0.0
    28 thread(s)
    
    outputhat23=16
    Loading 'hat3.seed' ... 
    done.
    Writing hat3 for iterative refinement
    generating a scoring matrix for nucleotide (dist=200) ... done
    Gap Penalty = -1.53, +0.00, +0.00
    tbutree = 1, compacttree = 0
    Constructing a UPGMA tree ... 
        0 / 3
    done.
    
    Progressive alignment ... 
    STEP     2 /2 (thread    1) 
    done.
    tbfast (nuc) Version 7.407
    alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
    16 thread(s)
    
    minimumweight = 0.000010
    autosubalignment = 0.000000
    nthread = 8
    randomseed = 0
    blosum 62 / kimura 200
    poffset = 0
    niter = 16
    sueff_global = 0.100000
    nadd = 16
    Loading 'hat3' ... done.
    generating a scoring matrix for nucleotide (dist=200) ... done
    
        0 / 3
    Segment   1/  1    1-1661
    001-0002-1 (thread    4) identical     
    Converged.
    done
    dvtditr (nuc) Version 7.407
    alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
    8 thread(s)
    
    
    Strategy:
     L-INS-i (Probably most accurate, very slow)
     Iterative refinement method (<16) with LOCAL pairwise alignment information
    
    If unsure which option to use, try 'mafft --auto input > output'.
    For more information, see 'mafft --help', 'mafft --man' and the mafft page.
    
    The default gap scoring scheme has been changed in version 7.110 (2013 Oct).
    It tends to insert more gaps into gap-rich regions than previous versions.
    To disable this change, add the --leavegappyregion option.
    
    mafft --thread -1 --maxiterate 1000 --localpair human-HA-sequences-by-flu-season/HA-H3N2-flu_season-2009-2010.fa > human-HA-sequences-by-flu-season/mafft-msa/HA-H3N2-flu_season-2009-2010.fa
    OS = linux
    The number of physical cores =  28
    outputhat23=16
    treein = 0
    compacttree = 0
    stacksize: 8192 kb
    generating a scoring matrix for nucleotide (dist=200) ... done
    All-to-all alignment.
    tbfast-pair (nuc) Version 7.407
    alg=L, model=DNA200 (2), 2.00 (6.00), -0.10 (-0.30), noshift, amax=0.0
    28 thread(s)
    
    outputhat23=16
    Loading 'hat3.seed' ... 
    done.
    Writing hat3 for iterative refinement
    generating a scoring matrix for nucleotide (dist=200) ... done
    Gap Penalty = -1.53, +0.00, +0.00
    tbutree = 1, compacttree = 0
    Constructing a UPGMA tree ... 
        0 / 3
    done.
    
    Progressive alignment ... 
    STEP     2 /2 (thread    1) 
    done.
    tbfast (nuc) Version 7.407
    alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
    16 thread(s)
    
    minimumweight = 0.000010
    autosubalignment = 0.000000
    nthread = 8
    randomseed = 0
    blosum 62 / kimura 200
    poffset = 0
    niter = 16
    sueff_global = 0.100000
    nadd = 16
    Loading 'hat3' ... done.
    generating a scoring matrix for nucleotide (dist=200) ... done
    
        0 / 3
    Segment   1/  1    1-1702
    001-0002-0 (thread    2) identical     
    Converged.
    done
    dvtditr (nuc) Version 7.407
    alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
    8 thread(s)
    
    
    Strategy:
     L-INS-i (Probably most accurate, very slow)
     Iterative refinement method (<16) with LOCAL pairwise alignment information
    
    If unsure which option to use, try 'mafft --auto input > output'.
    For more information, see 'mafft --help', 'mafft --man' and the mafft page.
    
    The default gap scoring scheme has been changed in version 7.110 (2013 Oct).
    It tends to insert more gaps into gap-rich regions than previous versions.
    To disable this change, add the --leavegappyregion option.
    
    mafft --thread -1 --maxiterate 1000 --localpair human-HA-sequences-by-flu-season/HA-H1N1-flu_season-2007-2008.fa > human-HA-sequences-by-flu-season/mafft-msa/HA-H1N1-flu_season-2007-2008.fa
    OS = linux
    The number of physical cores =  28
    outputhat23=16
    treein = 0
    compacttree = 0
    stacksize: 8192 kb
    generating a scoring matrix for nucleotide (dist=200) ... done
    All-to-all alignment.
    tbfast-pair (nuc) Version 7.407
    alg=L, model=DNA200 (2), 2.00 (6.00), -0.10 (-0.30), noshift, amax=0.0
    28 thread(s)
    
    outputhat23=16
    Loading 'hat3.seed' ... 
    done.
    Writing hat3 for iterative refinement
    generating a scoring matrix for nucleotide (dist=200) ... done
    Gap Penalty = -1.53, +0.00, +0.00
    tbutree = 1, compacttree = 0
    Constructing a UPGMA tree ... 
        0 / 6
    done.
    
    Progressive alignment ... 
    STEP     5 /5 (thread    4) 
    done.
    tbfast (nuc) Version 7.407
    alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
    16 thread(s)
    
    minimumweight = 0.000010
    autosubalignment = 0.000000
    nthread = 8
    randomseed = 0
    blosum 62 / kimura 200
    poffset = 0
    niter = 16
    sueff_global = 0.100000
    nadd = 16
    Loading 'hat3' ... done.
    generating a scoring matrix for nucleotide (dist=200) ... done
    
        0 / 6
    Segment   1/  1    1- 854
    001-0008-1 (thread    8) identical     
    Converged.
    done
    dvtditr (nuc) Version 7.407
    alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
    8 thread(s)
    
    
    Strategy:
     L-INS-i (Probably most accurate, very slow)
     Iterative refinement method (<16) with LOCAL pairwise alignment information
    
    If unsure which option to use, try 'mafft --auto input > output'.
    For more information, see 'mafft --help', 'mafft --man' and the mafft page.
    
    The default gap scoring scheme has been changed in version 7.110 (2013 Oct).
    It tends to insert more gaps into gap-rich regions than previous versions.
    To disable this change, add the --leavegappyregion option.
    
    mafft --thread -1 --maxiterate 1000 --localpair human-HA-sequences-by-flu-season/HA-H1N1-flu_season-2010-2011.fa > human-HA-sequences-by-flu-season/mafft-msa/HA-H1N1-flu_season-2010-2011.fa
    OS = linux
    The number of physical cores =  28
    outputhat23=16
    treein = 0
    compacttree = 0
    stacksize: 8192 kb
    generating a scoring matrix for nucleotide (dist=200) ... done
    All-to-all alignment.
    tbfast-pair (nuc) Version 7.407
    alg=L, model=DNA200 (2), 2.00 (6.00), -0.10 (-0.30), noshift, amax=0.0
    28 thread(s)
    
    outputhat23=16
    Loading 'hat3.seed' ... 
    done.
    Writing hat3 for iterative refinement
    generating a scoring matrix for nucleotide (dist=200) ... done
    Gap Penalty = -1.53, +0.00, +0.00
    tbutree = 1, compacttree = 0
    Constructing a UPGMA tree ... 
        0 / 5
    done.
    
    Progressive alignment ... 
    STEP     4 /4 (thread    3) 
    done.
    tbfast (nuc) Version 7.407
    alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
    16 thread(s)
    
    minimumweight = 0.000010
    autosubalignment = 0.000000
    nthread = 8
    randomseed = 0
    blosum 62 / kimura 200
    poffset = 0
    niter = 16
    sueff_global = 0.100000
    nadd = 16
    Loading 'hat3' ... done.
    generating a scoring matrix for nucleotide (dist=200) ... done
    
        0 / 5
    Segment   1/  1    1-1745
    001-0006-0 (thread    4) identical     
    Converged.
    done
    dvtditr (nuc) Version 7.407
    alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
    8 thread(s)
    
    
    Strategy:
     L-INS-i (Probably most accurate, very slow)
     Iterative refinement method (<16) with LOCAL pairwise alignment information
    
    If unsure which option to use, try 'mafft --auto input > output'.
    For more information, see 'mafft --help', 'mafft --man' and the mafft page.
    
    The default gap scoring scheme has been changed in version 7.110 (2013 Oct).
    It tends to insert more gaps into gap-rich regions than previous versions.
    To disable this change, add the --leavegappyregion option.
    
    mafft --thread -1 --maxiterate 1000 --localpair human-HA-sequences-by-flu-season/HA-H1N1-flu_season-2013-2014.fa > human-HA-sequences-by-flu-season/mafft-msa/HA-H1N1-flu_season-2013-2014.fa
    OS = linux
    The number of physical cores =  28
    outputhat23=16
    treein = 0
    compacttree = 0
    stacksize: 8192 kb
    generating a scoring matrix for nucleotide (dist=200) ... done
    All-to-all alignment.
    tbfast-pair (nuc) Version 7.407
    alg=L, model=DNA200 (2), 2.00 (6.00), -0.10 (-0.30), noshift, amax=0.0
    28 thread(s)
    
    outputhat23=16
    Loading 'hat3.seed' ... 
    done.
    Writing hat3 for iterative refinement
    generating a scoring matrix for nucleotide (dist=200) ... done
    Gap Penalty = -1.53, +0.00, +0.00
    tbutree = 1, compacttree = 0
    Constructing a UPGMA tree ... 
       10 / 13
    done.
    
    Progressive alignment ... 
    STEP    12 /12 (thread   11) 
    done.
    tbfast (nuc) Version 7.407
    alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
    16 thread(s)
    
    minimumweight = 0.000010
    autosubalignment = 0.000000
    nthread = 8
    randomseed = 0
    blosum 62 / kimura 200
    poffset = 0
    niter = 16
    sueff_global = 0.100000
    nadd = 16
    Loading 'hat3' ... done.
    generating a scoring matrix for nucleotide (dist=200) ... done
    
       10 / 13
    Segment   1/  1    1-1764
    001-0022-1 (thread    6) identical     
    Converged.
    done
    dvtditr (nuc) Version 7.407
    alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
    8 thread(s)
    
    
    Strategy:
     L-INS-i (Probably most accurate, very slow)
     Iterative refinement method (<16) with LOCAL pairwise alignment information
    
    If unsure which option to use, try 'mafft --auto input > output'.
    For more information, see 'mafft --help', 'mafft --man' and the mafft page.
    
    The default gap scoring scheme has been changed in version 7.110 (2013 Oct).
    It tends to insert more gaps into gap-rich regions than previous versions.
    To disable this change, add the --leavegappyregion option.
    
    mafft --thread -1 --maxiterate 1000 --localpair human-HA-sequences-by-flu-season/HA-H3N2-flu_season-2011-2012.fa > human-HA-sequences-by-flu-season/mafft-msa/HA-H3N2-flu_season-2011-2012.fa
    OS = linux
    The number of physical cores =  28
    outputhat23=16
    treein = 0
    compacttree = 0
    stacksize: 8192 kb
    generating a scoring matrix for nucleotide (dist=200) ... done
    All-to-all alignment.
    tbfast-pair (nuc) Version 7.407
    alg=L, model=DNA200 (2), 2.00 (6.00), -0.10 (-0.30), noshift, amax=0.0
    28 thread(s)
    
    outputhat23=16
    Loading 'hat3.seed' ... 
    done.
    Writing hat3 for iterative refinement
    generating a scoring matrix for nucleotide (dist=200) ... done
    Gap Penalty = -1.53, +0.00, +0.00
    tbutree = 1, compacttree = 0
    Constructing a UPGMA tree ... 
       10 / 12
    done.
    
    Progressive alignment ... 
    STEP    11 /11 (thread   10) 
    done.
    tbfast (nuc) Version 7.407
    alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
    16 thread(s)
    
    minimumweight = 0.000010
    autosubalignment = 0.000000
    nthread = 8
    randomseed = 0
    blosum 62 / kimura 200
    poffset = 0
    niter = 16
    sueff_global = 0.100000
    nadd = 16
    Loading 'hat3' ... done.
    generating a scoring matrix for nucleotide (dist=200) ... done
    
       10 / 12
    Segment   1/  1    1-1702
    001-0020-1 (thread    5) identical     
    Converged.
    done
    dvtditr (nuc) Version 7.407
    alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
    8 thread(s)
    
    
    Strategy:
     L-INS-i (Probably most accurate, very slow)
     Iterative refinement method (<16) with LOCAL pairwise alignment information
    
    If unsure which option to use, try 'mafft --auto input > output'.
    For more information, see 'mafft --help', 'mafft --man' and the mafft page.
    
    The default gap scoring scheme has been changed in version 7.110 (2013 Oct).
    It tends to insert more gaps into gap-rich regions than previous versions.
    To disable this change, add the --leavegappyregion option.
    
    mafft --thread -1 --maxiterate 1000 --localpair human-HA-sequences-by-flu-season/HA-H3N2-flu_season-2015-2016.fa > human-HA-sequences-by-flu-season/mafft-msa/HA-H3N2-flu_season-2015-2016.fa
    OS = linux
    The number of physical cores =  28
    outputhat23=16
    treein = 0
    compacttree = 0
    stacksize: 8192 kb
    generating a scoring matrix for nucleotide (dist=200) ... done
    All-to-all alignment.
    tbfast-pair (nuc) Version 7.407
    alg=L, model=DNA200 (2), 2.00 (6.00), -0.10 (-0.30), noshift, amax=0.0
    28 thread(s)
    
    outputhat23=16
    Loading 'hat3.seed' ... 
    done.
    Writing hat3 for iterative refinement
    generating a scoring matrix for nucleotide (dist=200) ... done
    Gap Penalty = -1.53, +0.00, +0.00
    tbutree = 1, compacttree = 0
    Constructing a UPGMA tree ... 
       10 / 12
    done.
    
    Progressive alignment ... 
    STEP    11 /11 (thread   10) 
    done.
    tbfast (nuc) Version 7.407
    alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
    16 thread(s)
    
    minimumweight = 0.000010
    autosubalignment = 0.000000
    nthread = 8
    randomseed = 0
    blosum 62 / kimura 200
    poffset = 0
    niter = 16
    sueff_global = 0.100000
    nadd = 16
    Loading 'hat3' ... done.
    generating a scoring matrix for nucleotide (dist=200) ... done
    
       10 / 12
    Segment   1/  1    1-1702
    001-0020-0 (thread    3) identical     
    Converged.
    done
    dvtditr (nuc) Version 7.407
    alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
    8 thread(s)
    
    
    Strategy:
     L-INS-i (Probably most accurate, very slow)
     Iterative refinement method (<16) with LOCAL pairwise alignment information
    
    If unsure which option to use, try 'mafft --auto input > output'.
    For more information, see 'mafft --help', 'mafft --man' and the mafft page.
    
    The default gap scoring scheme has been changed in version 7.110 (2013 Oct).
    It tends to insert more gaps into gap-rich regions than previous versions.
    To disable this change, add the --leavegappyregion option.
    
    mafft --thread -1 --maxiterate 1000 --localpair human-HA-sequences-by-flu-season/HA-H1N1-flu_season-2012-2013.fa > human-HA-sequences-by-flu-season/mafft-msa/HA-H1N1-flu_season-2012-2013.fa
    OS = linux
    The number of physical cores =  28
    outputhat23=16
    treein = 0
    compacttree = 0
    stacksize: 8192 kb
    generating a scoring matrix for nucleotide (dist=200) ... done
    All-to-all alignment.
    tbfast-pair (nuc) Version 7.407
    alg=L, model=DNA200 (2), 2.00 (6.00), -0.10 (-0.30), noshift, amax=0.0
    28 thread(s)
    
    outputhat23=16
    Loading 'hat3.seed' ... 
    done.
    Writing hat3 for iterative refinement
    generating a scoring matrix for nucleotide (dist=200) ... done
    Gap Penalty = -1.53, +0.00, +0.00
    tbutree = 1, compacttree = 0
    Constructing a UPGMA tree ... 
       10 / 18
    done.
    
    Progressive alignment ... 
    STEP    17 /17 (thread    0) 
    done.
    tbfast (nuc) Version 7.407
    alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
    16 thread(s)
    
    minimumweight = 0.000010
    autosubalignment = 0.000000
    nthread = 8
    randomseed = 0
    blosum 62 / kimura 200
    poffset = 0
    niter = 16
    sueff_global = 0.100000
    nadd = 16
    Loading 'hat3' ... done.
    generating a scoring matrix for nucleotide (dist=200) ... done
    
       10 / 18
    Segment   1/  1    1-1702
    001-0032-0 (thread    6) identical     
    Converged.
    done
    dvtditr (nuc) Version 7.407
    alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
    8 thread(s)
    
    
    Strategy:
     L-INS-i (Probably most accurate, very slow)
     Iterative refinement method (<16) with LOCAL pairwise alignment information
    
    If unsure which option to use, try 'mafft --auto input > output'.
    For more information, see 'mafft --help', 'mafft --man' and the mafft page.
    
    The default gap scoring scheme has been changed in version 7.110 (2013 Oct).
    It tends to insert more gaps into gap-rich regions than previous versions.
    To disable this change, add the --leavegappyregion option.
    
    mafft --thread -1 --maxiterate 1000 --localpair human-HA-sequences-by-flu-season/HA-H1N1-flu_season-2011-2012.fa > human-HA-sequences-by-flu-season/mafft-msa/HA-H1N1-flu_season-2011-2012.fa
    OS = linux
    The number of physical cores =  28
    outputhat23=16
    treein = 0
    compacttree = 0
    stacksize: 8192 kb
    generating a scoring matrix for nucleotide (dist=200) ... done
    All-to-all alignment.
    tbfast-pair (nuc) Version 7.407
    alg=L, model=DNA200 (2), 2.00 (6.00), -0.10 (-0.30), noshift, amax=0.0
    28 thread(s)
    
    outputhat23=16
    Loading 'hat3.seed' ... 
    done.
    Writing hat3 for iterative refinement
    generating a scoring matrix for nucleotide (dist=200) ... done
    Gap Penalty = -1.53, +0.00, +0.00
    tbutree = 1, compacttree = 0
    Constructing a UPGMA tree ... 
       10 / 18
    done.
    
    Progressive alignment ... 
    STEP    17 /17 (thread    3) 
    done.
    tbfast (nuc) Version 7.407
    alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
    16 thread(s)
    
    minimumweight = 0.000010
    autosubalignment = 0.000000
    nthread = 8
    randomseed = 0
    blosum 62 / kimura 200
    poffset = 0
    niter = 16
    sueff_global = 0.100000
    nadd = 16
    Loading 'hat3' ... done.
    generating a scoring matrix for nucleotide (dist=200) ... done
    
       10 / 18
    Segment   1/  1    1-1702
    001-0032-0 (thread    8) identical     
    Converged.
    done
    dvtditr (nuc) Version 7.407
    alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
    8 thread(s)
    
    
    Strategy:
     L-INS-i (Probably most accurate, very slow)
     Iterative refinement method (<16) with LOCAL pairwise alignment information
    
    If unsure which option to use, try 'mafft --auto input > output'.
    For more information, see 'mafft --help', 'mafft --man' and the mafft page.
    
    The default gap scoring scheme has been changed in version 7.110 (2013 Oct).
    It tends to insert more gaps into gap-rich regions than previous versions.
    To disable this change, add the --leavegappyregion option.
    
    mafft --thread -1 --maxiterate 1000 --localpair human-HA-sequences-by-flu-season/HA-H3N2-flu_season-2012-2013.fa > human-HA-sequences-by-flu-season/mafft-msa/HA-H3N2-flu_season-2012-2013.fa
    OS = linux
    The number of physical cores =  28
    outputhat23=16
    treein = 0
    compacttree = 0
    stacksize: 8192 kb
    generating a scoring matrix for nucleotide (dist=200) ... done
    All-to-all alignment.
    tbfast-pair (nuc) Version 7.407
    alg=L, model=DNA200 (2), 2.00 (6.00), -0.10 (-0.30), noshift, amax=0.0
    28 thread(s)
    
    outputhat23=16
    Loading 'hat3.seed' ... 
    done.
    Writing hat3 for iterative refinement
    generating a scoring matrix for nucleotide (dist=200) ... done
    Gap Penalty = -1.53, +0.00, +0.00
    tbutree = 1, compacttree = 0
    Constructing a UPGMA tree ... 
       30 / 39
    done.
    
    Progressive alignment ... 
    STEP    38 /38 (thread    9) 
    done.
    tbfast (nuc) Version 7.407
    alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
    16 thread(s)
    
    minimumweight = 0.000010
    autosubalignment = 0.000000
    nthread = 8
    randomseed = 0
    blosum 62 / kimura 200
    poffset = 0
    niter = 16
    sueff_global = 0.100000
    nadd = 16
    Loading 'hat3' ... done.
    generating a scoring matrix for nucleotide (dist=200) ... done
    
       30 / 39
    Segment   1/  1    1-1734
    001-0074-1 (thread    1) identical     
    Converged.
    done
    dvtditr (nuc) Version 7.407
    alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
    8 thread(s)
    
    
    Strategy:
     L-INS-i (Probably most accurate, very slow)
     Iterative refinement method (<16) with LOCAL pairwise alignment information
    
    If unsure which option to use, try 'mafft --auto input > output'.
    For more information, see 'mafft --help', 'mafft --man' and the mafft page.
    
    The default gap scoring scheme has been changed in version 7.110 (2013 Oct).
    It tends to insert more gaps into gap-rich regions than previous versions.
    To disable this change, add the --leavegappyregion option.
    
    mafft --thread -1 --maxiterate 1000 --localpair human-HA-sequences-by-flu-season/HA-H1N1-flu_season-2008-2009.fa > human-HA-sequences-by-flu-season/mafft-msa/HA-H1N1-flu_season-2008-2009.fa
    OS = linux
    The number of physical cores =  28
    outputhat23=16
    treein = 0
    compacttree = 0
    stacksize: 8192 kb
    generating a scoring matrix for nucleotide (dist=200) ... done
    All-to-all alignment.
    tbfast-pair (nuc) Version 7.407
    alg=L, model=DNA200 (2), 2.00 (6.00), -0.10 (-0.30), noshift, amax=0.0
    28 thread(s)
    
    outputhat23=16
    Loading 'hat3.seed' ... 
    done.
    Writing hat3 for iterative refinement
    generating a scoring matrix for nucleotide (dist=200) ... done
    Gap Penalty = -1.53, +0.00, +0.00
    tbutree = 1, compacttree = 0
    Constructing a UPGMA tree ... 
       30 / 40
    done.
    
    Progressive alignment ... 
    STEP    39 /39 (thread   11) 
    done.
    tbfast (nuc) Version 7.407
    alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
    16 thread(s)
    
    minimumweight = 0.000010
    autosubalignment = 0.000000
    nthread = 8
    randomseed = 0
    blosum 62 / kimura 200
    poffset = 0
    niter = 16
    sueff_global = 0.100000
    nadd = 16
    Loading 'hat3' ... done.
    generating a scoring matrix for nucleotide (dist=200) ... done
    
       30 / 40
    Segment   1/  1    1-1772
    001-0074-0 (thread    2) identical     
    Converged.
    done
    dvtditr (nuc) Version 7.407
    alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
    8 thread(s)
    
    
    Strategy:
     L-INS-i (Probably most accurate, very slow)
     Iterative refinement method (<16) with LOCAL pairwise alignment information
    
    If unsure which option to use, try 'mafft --auto input > output'.
    For more information, see 'mafft --help', 'mafft --man' and the mafft page.
    
    The default gap scoring scheme has been changed in version 7.110 (2013 Oct).
    It tends to insert more gaps into gap-rich regions than previous versions.
    To disable this change, add the --leavegappyregion option.
    
    mafft --thread -1 --maxiterate 1000 --localpair human-HA-sequences-by-flu-season/HA-H1N1-flu_season-2009-2010.fa > human-HA-sequences-by-flu-season/mafft-msa/HA-H1N1-flu_season-2009-2010.fa
    OS = linux
    The number of physical cores =  28
    outputhat23=16
    treein = 0
    compacttree = 0
    stacksize: 8192 kb
    generating a scoring matrix for nucleotide (dist=200) ... done
    All-to-all alignment.
    tbfast-pair (nuc) Version 7.407
    alg=L, model=DNA200 (2), 2.00 (6.00), -0.10 (-0.30), noshift, amax=0.0
    28 thread(s)
    
    outputhat23=16
    Loading 'hat3.seed' ... 
    done.
    Writing hat3 for iterative refinement
    generating a scoring matrix for nucleotide (dist=200) ... done
    Gap Penalty = -1.53, +0.00, +0.00
    tbutree = 1, compacttree = 0
    Constructing a UPGMA tree ... 
       40 / 51
    done.
    
    Progressive alignment ... 
    STEP    50 /50 (thread   13) 
    done.
    tbfast (nuc) Version 7.407
    alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
    16 thread(s)
    
    minimumweight = 0.000010
    autosubalignment = 0.000000
    nthread = 8
    randomseed = 0
    blosum 62 / kimura 200
    poffset = 0
    niter = 16
    sueff_global = 0.100000
    nadd = 16
    Loading 'hat3' ... done.
    generating a scoring matrix for nucleotide (dist=200) ... done
    
       40 / 51
    Segment   1/  1    1-1764
    001-0098-1 (thread    3) identical     
    Converged.
    done
    dvtditr (nuc) Version 7.407
    alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
    8 thread(s)
    
    
    Strategy:
     L-INS-i (Probably most accurate, very slow)
     Iterative refinement method (<16) with LOCAL pairwise alignment information
    
    If unsure which option to use, try 'mafft --auto input > output'.
    For more information, see 'mafft --help', 'mafft --man' and the mafft page.
    
    The default gap scoring scheme has been changed in version 7.110 (2013 Oct).
    It tends to insert more gaps into gap-rich regions than previous versions.
    To disable this change, add the --leavegappyregion option.
    
    mafft --thread -1 --maxiterate 1000 --localpair human-HA-sequences-by-flu-season/HA-H3N2-flu_season-2017-2018.fa > human-HA-sequences-by-flu-season/mafft-msa/HA-H3N2-flu_season-2017-2018.fa
    OS = linux
    The number of physical cores =  28
    outputhat23=16
    treein = 0
    compacttree = 0
    stacksize: 8192 kb
    generating a scoring matrix for nucleotide (dist=200) ... done
    All-to-all alignment.
    tbfast-pair (nuc) Version 7.407
    alg=L, model=DNA200 (2), 2.00 (6.00), -0.10 (-0.30), noshift, amax=0.0
    28 thread(s)
    
    outputhat23=16
    Loading 'hat3.seed' ... 
    done.
    Writing hat3 for iterative refinement
    generating a scoring matrix for nucleotide (dist=200) ... done
    Gap Penalty = -1.53, +0.00, +0.00
    tbutree = 1, compacttree = 0
    Constructing a UPGMA tree ... 
       40 / 50
    done.
    
    Progressive alignment ... 
    STEP    49 /49 (thread    1) 
    done.
    tbfast (nuc) Version 7.407
    alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
    16 thread(s)
    
    minimumweight = 0.000010
    autosubalignment = 0.000000
    nthread = 8
    randomseed = 0
    blosum 62 / kimura 200
    poffset = 0
    niter = 16
    sueff_global = 0.100000
    nadd = 16
    Loading 'hat3' ... done.
    generating a scoring matrix for nucleotide (dist=200) ... done
    
       40 / 50
    Segment   1/  1    1-1699
    001-0095-0 (thread    6) identical     
    Converged.
    done
    dvtditr (nuc) Version 7.407
    alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
    8 thread(s)
    
    
    Strategy:
     L-INS-i (Probably most accurate, very slow)
     Iterative refinement method (<16) with LOCAL pairwise alignment information
    
    If unsure which option to use, try 'mafft --auto input > output'.
    For more information, see 'mafft --help', 'mafft --man' and the mafft page.
    
    The default gap scoring scheme has been changed in version 7.110 (2013 Oct).
    It tends to insert more gaps into gap-rich regions than previous versions.
    To disable this change, add the --leavegappyregion option.
    
    mafft --thread -1 --maxiterate 1000 --localpair human-HA-sequences-by-flu-season/HA-H3N2-flu_season-2010-2011.fa > human-HA-sequences-by-flu-season/mafft-msa/HA-H3N2-flu_season-2010-2011.fa
    OS = linux
    The number of physical cores =  28
    outputhat23=16
    treein = 0
    compacttree = 0
    stacksize: 8192 kb
    generating a scoring matrix for nucleotide (dist=200) ... done
    All-to-all alignment.
    tbfast-pair (nuc) Version 7.407
    alg=L, model=DNA200 (2), 2.00 (6.00), -0.10 (-0.30), noshift, amax=0.0
    28 thread(s)
    
    outputhat23=16
    Loading 'hat3.seed' ... 
    done.
    Writing hat3 for iterative refinement
    generating a scoring matrix for nucleotide (dist=200) ... done
    Gap Penalty = -1.53, +0.00, +0.00
    tbutree = 1, compacttree = 0
    Constructing a UPGMA tree ... 
       60 / 62
    done.
    
    Progressive alignment ... 
    STEP    61 /61 (thread    7) 
    done.
    tbfast (nuc) Version 7.407
    alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
    16 thread(s)
    
    minimumweight = 0.000010
    autosubalignment = 0.000000
    nthread = 8
    randomseed = 0
    blosum 62 / kimura 200
    poffset = 0
    niter = 16
    sueff_global = 0.100000
    nadd = 16
    Loading 'hat3' ... done.
    generating a scoring matrix for nucleotide (dist=200) ... done
    
       60 / 62
    Segment   1/  1    1-1754
    001-0120-0 (thread    4) identical     
    Converged.
    done
    dvtditr (nuc) Version 7.407
    alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
    8 thread(s)
    
    
    Strategy:
     L-INS-i (Probably most accurate, very slow)
     Iterative refinement method (<16) with LOCAL pairwise alignment information
    
    If unsure which option to use, try 'mafft --auto input > output'.
    For more information, see 'mafft --help', 'mafft --man' and the mafft page.
    
    The default gap scoring scheme has been changed in version 7.110 (2013 Oct).
    It tends to insert more gaps into gap-rich regions than previous versions.
    To disable this change, add the --leavegappyregion option.
    
    mafft --thread -1 --maxiterate 1000 --localpair human-HA-sequences-by-flu-season/HA-H1N1-flu_season-2015-2016.fa > human-HA-sequences-by-flu-season/mafft-msa/HA-H1N1-flu_season-2015-2016.fa
    OS = linux
    The number of physical cores =  28
    outputhat23=16
    treein = 0
    compacttree = 0
    stacksize: 8192 kb
    generating a scoring matrix for nucleotide (dist=200) ... done
    All-to-all alignment.
    tbfast-pair (nuc) Version 7.407
    alg=L, model=DNA200 (2), 2.00 (6.00), -0.10 (-0.30), noshift, amax=0.0
    28 thread(s)
    
    outputhat23=16
    Loading 'hat3.seed' ... 
    done.
    Writing hat3 for iterative refinement
    generating a scoring matrix for nucleotide (dist=200) ... done
    Gap Penalty = -1.53, +0.00, +0.00
    tbutree = 1, compacttree = 0
    Constructing a UPGMA tree ... 
      190 / 198
    done.
    
    Progressive alignment ... 
    STEP   132 /197 (thread    7) 
    Reallocating (by thread 15) ..done. *alloclen = 3096
    STEP   197 /197 (thread    7) 
    done.
    tbfast (nuc) Version 7.407
    alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
    16 thread(s)
    
    minimumweight = 0.000010
    autosubalignment = 0.000000
    nthread = 8
    randomseed = 0
    blosum 62 / kimura 200
    poffset = 0
    niter = 16
    sueff_global = 0.100000
    nadd = 16
    Loading 'hat3' ... done.
    generating a scoring matrix for nucleotide (dist=200) ... done
    
      190 / 198
    Segment   1/  1    1-1065
    
    Converged.
    done
    dvtditr (nuc) Version 7.407
    alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
    8 thread(s)
    
    
    Strategy:
     L-INS-i (Probably most accurate, very slow)
     Iterative refinement method (<16) with LOCAL pairwise alignment information
    
    If unsure which option to use, try 'mafft --auto input > output'.
    For more information, see 'mafft --help', 'mafft --man' and the mafft page.
    
    The default gap scoring scheme has been changed in version 7.110 (2013 Oct).
    It tends to insert more gaps into gap-rich regions than previous versions.
    To disable this change, add the --leavegappyregion option.
    
    mafft --thread -1 --maxiterate 1000 --localpair human-HA-sequences-by-flu-season/HA-H3N2-flu_season-2014-2015.fa > human-HA-sequences-by-flu-season/mafft-msa/HA-H3N2-flu_season-2014-2015.fa
    OS = linux
    The number of physical cores =  28
    outputhat23=16
    treein = 0
    compacttree = 0
    stacksize: 8192 kb
    generating a scoring matrix for nucleotide (dist=200) ... done
    All-to-all alignment.
    tbfast-pair (nuc) Version 7.407
    alg=L, model=DNA200 (2), 2.00 (6.00), -0.10 (-0.30), noshift, amax=0.0
    28 thread(s)
    
    outputhat23=16
    Loading 'hat3.seed' ... 
    done.
    Writing hat3 for iterative refinement
    generating a scoring matrix for nucleotide (dist=200) ... done
    Gap Penalty = -1.53, +0.00, +0.00
    tbutree = 1, compacttree = 0
    Constructing a UPGMA tree ... 
      130 / 135
    done.
    
    Progressive alignment ... 
    STEP   134 /134 (thread    7) 
    done.
    tbfast (nuc) Version 7.407
    alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
    16 thread(s)
    
    minimumweight = 0.000010
    autosubalignment = 0.000000
    nthread = 8
    randomseed = 0
    blosum 62 / kimura 200
    poffset = 0
    niter = 16
    sueff_global = 0.100000
    nadd = 16
    Loading 'hat3' ... done.
    generating a scoring matrix for nucleotide (dist=200) ... done
    
      130 / 135
    Segment   1/  1    1-1713
    
    Converged.
    done
    dvtditr (nuc) Version 7.407
    alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
    8 thread(s)
    
    
    Strategy:
     L-INS-i (Probably most accurate, very slow)
     Iterative refinement method (<16) with LOCAL pairwise alignment information
    
    If unsure which option to use, try 'mafft --auto input > output'.
    For more information, see 'mafft --help', 'mafft --man' and the mafft page.
    
    The default gap scoring scheme has been changed in version 7.110 (2013 Oct).
    It tends to insert more gaps into gap-rich regions than previous versions.
    To disable this change, add the --leavegappyregion option.
    
    mafft --thread -1 --maxiterate 1000 --localpair human-HA-sequences-by-flu-season/HA-H3N2-flu_season-2016-2017.fa > human-HA-sequences-by-flu-season/mafft-msa/HA-H3N2-flu_season-2016-2017.fa
    OS = linux
    The number of physical cores =  28
    outputhat23=16
    treein = 0
    compacttree = 0
    stacksize: 8192 kb
    generating a scoring matrix for nucleotide (dist=200) ... done
    All-to-all alignment.
    tbfast-pair (nuc) Version 7.407
    alg=L, model=DNA200 (2), 2.00 (6.00), -0.10 (-0.30), noshift, amax=0.0
    28 thread(s)
    
    outputhat23=16
    Loading 'hat3.seed' ... 
    done.
    Writing hat3 for iterative refinement
    generating a scoring matrix for nucleotide (dist=200) ... done
    Gap Penalty = -1.53, +0.00, +0.00
    tbutree = 1, compacttree = 0
    Constructing a UPGMA tree ... 
      160 / 162
    done.
    
    Progressive alignment ... 
    STEP   161 /161 (thread   14) 
    done.
    tbfast (nuc) Version 7.407
    alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
    16 thread(s)
    
    minimumweight = 0.000010
    autosubalignment = 0.000000
    nthread = 8
    randomseed = 0
    blosum 62 / kimura 200
    poffset = 0
    niter = 16
    sueff_global = 0.100000
    nadd = 16
    Loading 'hat3' ... done.
    generating a scoring matrix for nucleotide (dist=200) ... done
    
      160 / 162
    Segment   1/  1    1-1702
    
    Converged.
    done
    dvtditr (nuc) Version 7.407
    alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
    8 thread(s)
    
    
    Strategy:
     L-INS-i (Probably most accurate, very slow)
     Iterative refinement method (<16) with LOCAL pairwise alignment information
    
    If unsure which option to use, try 'mafft --auto input > output'.
    For more information, see 'mafft --help', 'mafft --man' and the mafft page.
    
    The default gap scoring scheme has been changed in version 7.110 (2013 Oct).
    It tends to insert more gaps into gap-rich regions than previous versions.
    To disable this change, add the --leavegappyregion option.
    



```python
!mkdir human-HA-sequences-by-flu-season/concat
```


```python
!cat human-HA-sequences-by-flu-season/HA-H1N1*.fa > human-HA-sequences-by-flu-season/concat/HA-H1N1-all-seasons.fa
!cat human-HA-sequences-by-flu-season/HA-H3N2*.fa > human-HA-sequences-by-flu-season/concat/HA-H3N2-all-seasons.fa
```


```python
!mkdir human-HA-sequences-by-flu-season/concat/mafft-msa
```


```python
!parallel -v mafft --thread -1 --maxiterate 1000 --localpair {} ">" human-HA-sequences-by-flu-season/concat/mafft-msa/{/} ::: human-HA-sequences-by-flu-season/concat/*.fa
```

    mafft --thread -1 --maxiterate 1000 --localpair human-HA-sequences-by-flu-season/concat/HA-H1N1-all-seasons.fa > human-HA-sequences-by-flu-season/concat/mafft-msa/HA-H1N1-all-seasons.fa
    OS = linux
    The number of physical cores =  28
    outputhat23=16
    treein = 0
    compacttree = 0
    stacksize: 8192 kb
    generating a scoring matrix for nucleotide (dist=200) ... done
    All-to-all alignment.
    tbfast-pair (nuc) Version 7.407
    alg=L, model=DNA200 (2), 2.00 (6.00), -0.10 (-0.30), noshift, amax=0.0
    28 thread(s)
    
    outputhat23=16
    Loading 'hat3.seed' ... 
    done.
    Writing hat3 for iterative refinement
    generating a scoring matrix for nucleotide (dist=200) ... done
    Gap Penalty = -1.53, +0.00, +0.00
    tbutree = 1, compacttree = 0
    Constructing a UPGMA tree ... 
      340 / 350
    done.
    
    Progressive alignment ... 
    STEP   342 /349 (thread    8) 
    Reallocating (by thread 6) ..done. *alloclen = 4548
    STEP   349 /349 (thread   14) 
    done.
    tbfast (nuc) Version 7.407
    alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
    16 thread(s)
    
    minimumweight = 0.000010
    autosubalignment = 0.000000
    nthread = 8
    randomseed = 0
    blosum 62 / kimura 200
    poffset = 0
    niter = 16
    sueff_global = 0.100000
    nadd = 16
    Loading 'hat3' ... done.
    generating a scoring matrix for nucleotide (dist=200) ... done
    
      340 / 350
    Segment   1/  1    1-1786
    
    Converged.
    done
    dvtditr (nuc) Version 7.407
    alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
    8 thread(s)
    
    
    Strategy:
     L-INS-i (Probably most accurate, very slow)
     Iterative refinement method (<16) with LOCAL pairwise alignment information
    
    If unsure which option to use, try 'mafft --auto input > output'.
    For more information, see 'mafft --help', 'mafft --man' and the mafft page.
    
    The default gap scoring scheme has been changed in version 7.110 (2013 Oct).
    It tends to insert more gaps into gap-rich regions than previous versions.
    To disable this change, add the --leavegappyregion option.
    
    mafft --thread -1 --maxiterate 1000 --localpair human-HA-sequences-by-flu-season/concat/HA-H3N2-all-seasons.fa > human-HA-sequences-by-flu-season/concat/mafft-msa/HA-H3N2-all-seasons.fa
    OS = linux
    The number of physical cores =  28
    outputhat23=16
    treein = 0
    compacttree = 0
    stacksize: 8192 kb
    generating a scoring matrix for nucleotide (dist=200) ... done
    All-to-all alignment.
    tbfast-pair (nuc) Version 7.407
    alg=L, model=DNA200 (2), 2.00 (6.00), -0.10 (-0.30), noshift, amax=0.0
    28 thread(s)
    
    outputhat23=16
    Loading 'hat3.seed' ... 
    done.
    Writing hat3 for iterative refinement
    generating a scoring matrix for nucleotide (dist=200) ... done
    Gap Penalty = -1.53, +0.00, +0.00
    tbutree = 1, compacttree = 0
    Constructing a UPGMA tree ... 
      470 / 478
    done.
    
    Progressive alignment ... 
    STEP   477 /477 (thread    7) 
    done.
    tbfast (nuc) Version 7.407
    alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
    16 thread(s)
    
    minimumweight = 0.000010
    autosubalignment = 0.000000
    nthread = 8
    randomseed = 0
    blosum 62 / kimura 200
    poffset = 0
    niter = 16
    sueff_global = 0.100000
    nadd = 16
    Loading 'hat3' ... done.
    generating a scoring matrix for nucleotide (dist=200) ... done
    
      470 / 478
    Segment   1/  1    1-1754
    
    Converged.
    done
    dvtditr (nuc) Version 7.407
    alg=A, model=DNA200 (2), 1.53 (4.59), -0.00 (-0.00), noshift, amax=0.0
    8 thread(s)
    
    
    Strategy:
     L-INS-i (Probably most accurate, very slow)
     Iterative refinement method (<16) with LOCAL pairwise alignment information
    
    If unsure which option to use, try 'mafft --auto input > output'.
    For more information, see 'mafft --help', 'mafft --man' and the mafft page.
    
    The default gap scoring scheme has been changed in version 7.110 (2013 Oct).
    It tends to insert more gaps into gap-rich regions than previous versions.
    To disable this change, add the --leavegappyregion option.
    



```python

```
