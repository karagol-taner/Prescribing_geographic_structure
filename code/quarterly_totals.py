import os as _os, sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path.insert(0, _HERE)
_ROOT = _os.path.dirname(_HERE)
_os.makedirs(_os.path.join(_ROOT, 'figures'), exist_ok=True)
_os.chdir(_os.path.join(_ROOT, 'data'))

"""Quarterly totals by area from the raw monthly Prescription Cost Analysis files (data/temporal_data.csv).
For each quarter and area: dispensed items, net ingredient cost (GBP) and the number of distinct BNF chemical
substances (the largest monthly count in the quarter). Reads the monthly files pca_YYYYMM.csv for January 2021 to
December 2025 one at a time; the area and item columns are named differently in the earlier (STP) and later (ICB)
files and are mapped to common names. The raw files (about 15 GB) are not included; see
code/01_aggregate_pca_data.ipynb for where to obtain them.
Usage: python code/quarterly_totals.py <folder with the monthly files>   (writes data/temporal_data.csv)"""
import sys
from pathlib import Path
import pandas as pd

SRC = Path(_ROOT, sys.argv[1]) if len(sys.argv) > 1 else None
if SRC is None or not SRC.is_dir():
    sys.exit('usage: python code/quarterly_totals.py <folder with pca_YYYYMM.csv files>')
files = [f for f in sorted(SRC.glob('pca_*.csv')) if '202101' <= f.stem.split('_')[1] <= '202512']
if len(files) != 60:
    sys.exit(f'expected the 60 monthly files for January 2021 to December 2025, found {len(files)}')

KEEP = {'STP_CODE', 'STP_NAME', 'ICB_CODE', 'ICB_NAME', 'BNF_CHEMICAL_SUBSTANCE', 'CHEMICAL_BNF_NAME',
        'GENERIC_BNF_EQUIVALENT_NAME', 'GENERIC_BNF_NAME', 'SUPPLIER_NAME', 'ITEMS', 'PRESCRIBED_ITEMS',
        'NIC', 'ACTUAL_COST', 'YEAR_MONTH'}
RENAME = {'ICB_CODE': 'STP_CODE', 'ICB_NAME': 'STP_NAME', 'BNF_CHEMICAL_SUBSTANCE': 'CHEMICAL_BNF_NAME',
          'GENERIC_BNF_EQUIVALENT_NAME': 'GENERIC_BNF_NAME', 'PRESCRIBED_ITEMS': 'ITEMS', 'ACTUAL_COST': 'NIC'}
parts = []
for f in files:
    df = pd.read_csv(f, usecols=lambda c: c in KEEP, dtype={'STP_CODE': str, 'ICB_CODE': str, 'YEAR_MONTH': str})
    df = df.rename(columns=RENAME)
    missing = [c for c in ['STP_CODE', 'STP_NAME', 'CHEMICAL_BNF_NAME', 'ITEMS', 'NIC', 'YEAR_MONTH'] if c not in df.columns]
    if missing:
        sys.exit(f'{f.name}: columns not found: {missing}')
    df['ITEMS'] = pd.to_numeric(df['ITEMS'], errors='coerce').fillna(0)
    df['NIC'] = pd.to_numeric(df['NIC'], errors='coerce').fillna(0)
    month = df['YEAR_MONTH'].str[4:].astype(int)
    df['Period'] = df['YEAR_MONTH'].str[:4].astype(int).astype(str) + '-Q' + ((month - 1) // 3 + 1).astype(str)
    parts.append(df.groupby(['Period', 'STP_CODE']).agg({'ITEMS': 'sum', 'NIC': 'sum', 'CHEMICAL_BNF_NAME': 'nunique'}).reset_index())
    print(f.name, len(df))
T = pd.concat(parts).groupby(['Period', 'STP_CODE']).agg({'ITEMS': 'sum', 'NIC': 'sum', 'CHEMICAL_BNF_NAME': 'max'}).reset_index()
T = T.rename(columns={'ITEMS': 'Total_Items', 'NIC': 'Total_Cost', 'CHEMICAL_BNF_NAME': 'Unique_Chemicals'})
T.to_csv('temporal_data.csv', index=False)
print('saved temporal_data.csv', T.shape)
