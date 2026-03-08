"""Get tumor variants
Only returns consensus calls from strelka and mutect in the following form
dict of patient to set
set contains tuples
tuples contain (chrom, pos, ref, alt, blood alt counts, tumor alt counts, gnomAD population AF, dbSNP population AF)
"""

from cyvcf2 import VCF
import gc
from glob import glob
import numpy as np
import os
import pandas as pd
import pysam
from tqdm import tqdm


main_path = '/path/to/data_analysis/'
patients = ['300105', '300103', '300101', '300102', '300104', '300109', '300110', '300108', '300106']

paths = ["/path/to/sequencing_run_1/Analysis/bwa/sarek/",
        "/path/to/sequencing_run_2/Analysis/bwa/sarek/",
        "/path/to/sequencing_run_3/Analysis/bwa/sarek/"]
patient_to_vcf_path = {p: glob(f"{m}sarek_{p}/outs/variant_calling/mutect2/*vs*/*.mutect2.filtered.vcf.gz")[0]
                    for m in paths for p in [x.split('sarek_')[-1] for x in glob(f"{m}*sarek*")]}
patient_to_strelka_path = {p: glob(f"{m}sarek_{p}/outs/variant_calling/strelka/*vs*/*.strelka.somatic_snvs.vcf.gz")[0]
                    for m in paths for p in [x.split('sarek_')[-1] for x in glob(f"{m}*sarek*")]}

patient_to_maf_path = {}
for m in paths:
    for patient_dir in glob(f"{m}*sarek*"):
        p = patient_dir.split('sarek_')[-1]
        pattern = f"{m}sarek_{p}/outs/variant_calling/mutect2/*vs*/*.mutect2.filtered.funcotated.maf"
        matches = glob(pattern)
        
        if matches:
            patient_to_maf_path[p] = matches[0]
        else:
            print(f"❌ No files found for patient {p}")
            print(f"   Looked in: {pattern}")
            # Debug: check what files actually exist
            print(f"   Available: {glob(f'{m}sarek_{p}/outs/variant_calling/mutect2/*vs*/')}")

patient_to_vcf_path = {p:v for p,v in patient_to_vcf_path.items() if p in patients}
patient_to_strelka_path = {p:v for p,v in patient_to_strelka_path.items() if p in patients}
patient_to_maf_path = {p:v for p,v in patient_to_maf_path.items() if p in patients}


def _read_and_tag_maf(f):
    dat = pd.read_csv(
        f, 
        sep="\t", 
        dtype=str, 
        comment="#", 
        quoting=3,  # csv.QUOTE_NONE = 3
        engine='python'
    )
    dat['file_path'] = f
    return dat

patient_to_maf = {}
for p, f in patient_to_maf_path.items():
    try:
        df = _read_and_tag_maf(f)
        if len(df) > 0:
            count_cols = ['n_ref_count', 'n_alt_count', 't_ref_count', 't_alt_count']
            df[count_cols] = df[count_cols].astype(int)
            df['POPAF'] = df['POPAF'].astype(float)
            df = df[df.t_alt_count > 0]
            df.loc[df.POPAF == 6.0, 'POPAF'] = np.nan
            df['dbSNP_AF'] = np.round(pd.to_numeric(df['dbSNP_TOPMED'].str.split(',').str[1], errors='coerce') * 100, 2)

            df['variant'] = list(zip(df['Chromosome'], df['Start_Position'].astype(int), df['Reference_Allele'], df['Tumor_Seq_Allele2']))
            df['full_variant_code'] = list(zip(df['Chromosome'], df['Start_Position'].astype(int),
                                               df['Reference_Allele'], df['Tumor_Seq_Allele2'], df['n_alt_count'],
                                               df['t_alt_count'], df['POPAF'], df['dbSNP_AF']))
            patient_to_maf[p] = df

    except Exception as e:
        print(f"Error reading {p}: {e}")

print(f"\nLoaded {len(patient_to_maf)} MAF files")

# The below code is redundant; already loading in variants from maf
'''
patient_to_mutect_variants = {}
print('Processing tumor variants...')
for patient_id, vcf_path in tqdm(patient_to_vcf_path.items()):
    tumor_vcf = VCF(vcf_path)
    maf_filtered_variants = set(patient_to_maf[patient_id]['variant'])

    patient_to_mutect_variants[patient_id] = {
        (rec.CHROM, rec.POS, rec.REF, alt)
        for rec in tumor_vcf
        for alt in rec.ALT
        if (rec.FILTER is None or str(rec.FILTER).lower() == "pass")
        and (rec.CHROM, rec.POS, rec.REF, alt) in maf_filtered_variants
    }
print('Done!')
print(f"Loaded mutect variants for {len(patients)} patients")
'''

patient_to_strelka_variants = {}
print('Processing tumor variants...')
for patient_id, vcf_path in tqdm(patient_to_strelka_path.items()):
    tumor_vcf = VCF(vcf_path)
    patient_to_strelka_variants[patient_id] = {
        (rec.CHROM, rec.POS, rec.REF, alt)
        for rec in tumor_vcf
        for alt in rec.ALT
        if rec.FILTER is None or str(rec.FILTER).lower() == "pass"
    }
print('Done!')
print(f"Loaded strelka variants for {len(patients)} patients")

for p, maf_df in patient_to_maf.items():
    patient_to_maf[p] = maf_df[maf_df.variant.isin(patient_to_strelka_variants[p])]
# patient_to_tumor_variants = {p:patient_to_strelka_variants[p].intersection(patient_to_mutect_variants[p]) for p in patients}
patient_to_tumor_variants = {p: set(patient_to_maf[p]['full_variant_code']) for p in patients}
del patient_to_strelka_variants, patient_to_maf
gc.collect()