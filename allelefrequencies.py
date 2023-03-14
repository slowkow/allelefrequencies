"""
allelefrequencies.py
====================

Kamil Slowikowski
2023-03-14


Overview
--------

This script downloads allele frequencies for HLA, KIR, MIC, and cytokine genes
from the Allele Frequency Net Database:

http://allelefrequencies.net

Then, the script writes .tsv files like this:

    ls allelefrequencies.net/hla/

        787K  A.tsv
        1.5M  B.tsv
        454K  C.tsv
        8.5K  DPA1.tsv
        201K  DPB1.tsv
         88K  DQA1.tsv
        291K  DQB1.tsv
        1.2M  DRB1.tsv

Each file looks like this:

    head allelefrequencies.net/A.tsv

     allele                              population indivs_over_n alleles_over_2n   n
    A*01:01                  Argentina Rosario Toba          15.1          0.0760  86
    A*01:01                Armenia combined Regions                        0.1250 100
    A*01:01 Australia Cape York Peninsula Aborigine                        0.0530 103
    A*01:01      Australia Groote Eylandt Aborigine                        0.0270  75
    A*01:01     Australia New South Wales Caucasian                        0.1870 134
    A*01:01            Australia Yuendumu Aborigine                        0.0080 191


Citation
--------

Please cite the latest manuscript about Allele Frequency Net Database:

https://pubmed.ncbi.nlm.nih.gov/31722398

Gonzalez-Galarza FF, McCabe A, Santos EJMD, Jones J, Takeshita L, Ortega-Rivera
ND, et al. Allele frequency net database (AFND) 2020 update: gold-standard data
classification, open access genotype data and new query tools. Nucleic Acids
Res. 2020;48: D783–D788. doi:10.1093/nar/gkz1029


Acknowledgments
---------------

Thanks to David A. Wells for sharing [scrapeAF][1], which inspired me to work
on this project.

[1]: https://github.com/DAWells/scrapeAF


License
-------

MIT License

Copyright (c) 2023 Kamil Slowikowski

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import os
import pandas as pd
import re
import requests
from bs4 import BeautifulSoup
from urllib.parse import urlencode
from tqdm import tqdm
import gzip

def main():
    download_hla()
    download_kir()
    download_cyt()
    download_mic()

def download_hla():
    out = 'hla'
    base = 'http://www.allelefrequencies.net/hla6006a.asp'
    loci = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1']
    for locus in loci:
        tsv_file = f'allelefrequencies.net/{out}/{locus}.tsv'
        print(tsv_file)
        if not os.path.exists(tsv_file):
            params = {
                'hla_locus': locus,
                'hla_locus_type': 'Classical',
                'hla_level': '2'
            }
            url = f'{base}?{urlencode(params)}'
            df = get_all_pages(url)
            df = df[[1,3,4,5,7]]
            df.columns = ['allele', 'population',
                    'indivs_over_n', 'alleles_over_2n', 'n']
            write_tsv(df, tsv_file)

def download_kir():
    out = 'kir'
    base = 'http://www.allelefrequencies.net/kir6002a.asp'
    loci = ['2DL1', '2DL2', '2DL3', '2DL4', '2DL5', '2DL5A', '2DL5B', '2DP1',
            '2DS1', '2DS2', '2DS3', '2DS4', '2DS5', '3DL1', '3DL2', '3DL3',
            '3DP1', '3DS1']
    for locus in loci:
        tsv_file = f'allelefrequencies.net/{out}/{locus}.tsv'
        print(tsv_file)
        if not os.path.exists(tsv_file):
            params = {
                'kir_locus': locus
            }
            url = f'{base}?{urlencode(params)}'
            df = get_all_pages(url)
            df = df[[1,3,5,6,8]]
            df.columns = ['allele', 'population',
                    'indivs_over_n', 'alleles_over_2n', 'n']
            write_tsv(df, tsv_file)

def download_cyt():
    out = 'cyt'
    base = 'http://www.allelefrequencies.net/cyt6001a.asp'
    loci = ['AIF-1/', 'bFGF/', 'EGF/', 'GM-CSF/', 'IFNgamma/', 'IGF-1/',
            'IL-10/', 'IL-12/', 'IL-12p40/', 'IL-13/', 'IL-15/', 'IL-18/',
            'IL-1alpha/', 'IL-1beta/', 'IL1RA/', 'IL1R/', 'IL-2/', 'IL-4/',
            'IL-4R alpha/', 'IL-4R/', 'IL-6/', 'NGF/', 'PDGF A/', 'PDGF B/',
            'RANTES/', 'TGFbeta1/', 'TNFalpha/', 'TNFbeta/', 'VEGF/']
    for locus in loci:
        tsv_file = f'allelefrequencies.net/{out}/{safe(locus)}.tsv'
        print(tsv_file)
        if not os.path.exists(tsv_file):
            params = {
                'cyt_gene': locus
            }
            url = f'{base}?{urlencode(params)}'
            df = get_all_pages(url)
            if not df.empty:
                df = df[[1,3,4,5,6]]
                df.columns = ['allele', 'population',
                        'indivs_over_n', 'alleles_over_2n', 'n']
                write_tsv(df, tsv_file)

def download_mic():
    out = 'mic'
    base = 'http://www.allelefrequencies.net/mic6001a.asp'
    loci = ['MICA', 'MICB']
    for locus in loci:
        tsv_file = f'allelefrequencies.net/{out}/{locus}.tsv'
        print(tsv_file)
        if not os.path.exists(tsv_file):
            params = {
                'mic_locus': locus,
                'mic_locus_type': 'Classical',
                'mic_order': 'order_1'
            }
            url = f'{base}?{urlencode(params)}'
            df = get_all_pages(url)
            df = df[[1,3,5,6,8]]
            df.columns = ['allele', 'population',
                    'indivs_over_n', 'alleles_over_2n', 'n']
            write_tsv(df, tsv_file)

def write_tsv(d, tsv_file):
    print(f'Writing {tsv_file}')
    os.makedirs(os.path.dirname(tsv_file), exist_ok=True)
    d.to_csv(tsv_file, sep='\t', index=False)

def safe(x):
    return re.sub(r"[/\\?%*:|\"<>\x7F\x00-\x1F ]", "-", x)

def get_url(url):
    if not os.path.exists('cache'):
        os.makedirs('cache', exist_ok=True)
    cache_file = f'cache/{safe(url)}.html.gz'
    if os.path.exists(cache_file):
        text = gzip.open(cache_file, 'rt').read()
    else:
        req = requests.get(url)
        text = req.text
        with gzip.open(cache_file, 'wt') as out:
            out.write(text)
    return text

def get_all_pages(url):
    text = get_url(url)
    bs = BeautifulSoup(text, 'html.parser')
    # Number of pages of results
    n_pages = bs(text = re.compile(r' of \d'))
    if n_pages:
        n_pages = int(n_pages[0].strip()[3:].replace(',', ''))
    else:
        n_pages = 1
    if bs(text = re.compile(r'we did not find any results')):
        return pd.DataFrame(None)
    print(f'{n_pages} pages of results')
    tables = []
    for i in tqdm(range(1, n_pages + 1)):
        text = get_url(f'{url}&page={i}')
        bs = BeautifulSoup(text, 'html.parser')
        table = get_df(bs)
        tables.append(table)
    df = pd.concat(tables)
    return df 

def get_df(bs):
    table = bs.find('table', {'class': 'tblNormal'})
    trs = table.find_all('tr')
    th = trs[0].find_all('th')
    columns = [x.get_text(strip=True) for x in th]
    rows = []
    for tr in trs[1:]:
        fields = [td.get_text(strip=True) for td in tr.find_all('td')]
        rows.append(fields)
    df = pd.DataFrame(rows)
    return df

if __name__ == '__main__':
    main()

