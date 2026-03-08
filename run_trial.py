from match_mutations import calculate_f1_score
from get_cell_variants import master_dict_all
from get_tumor_variants import patients, patient_to_tumor_variants
from get_epithelial_cells import pt_cells, ep_cells

import pandas as pd


def apply_filters(variant_list, min_dp, min_alt_reads, min_af):
    """
    Apply filtering thresholds to a list of variants.
    Input: list of (chrom, pos, ref, alt, dp, alt_depth, af)
    Output: set of (chrom, pos, ref, alt)
    """
    filtered = set()
    for chrom, pos, ref, alt, dp, alt_depth, af in variant_list:
        if dp >= min_dp and alt_depth >= min_alt_reads and af >= min_af:
            filtered.add((chrom, pos, ref, alt))
    return filtered

def apply_filters_to_all(master_dict_all, min_dp, min_alt_reads, min_af):
    """
    Apply filters to all cells across all patients.
    Input: master_dict_all[patient][cell] -> list of (chrom, pos, ref, alt, dp, alt_depth, af)
    Output: master_dict_filtered[patient][cell] -> set of (chrom, pos, ref, alt)
    """
    master_dict_filtered = {}
    for patient, cells_dict in master_dict_all.items():
        master_dict_filtered[patient] = {}
        for cell_name, variant_list in cells_dict.items():
            master_dict_filtered[patient][cell_name] = apply_filters(variant_list, min_dp, min_alt_reads, min_af)
    return master_dict_filtered

def apply_filters_tumors(patient_to_tumor_variants, max_blood_alt_counts, min_tumor_alt_counts, max_gnomAD_freq, max_dbSNP_freq):
    patient_to_tumor_variants_filtered = {}
    for patient, variant_set in patient_to_tumor_variants.items():
        patient_to_tumor_variants_filtered[patient] = set()
        for chrom, pos, ref, alt, blood_alt_counts, tumor_alt_counts, gnomAD_freq, dbSNP_freq in variant_set:
            if (blood_alt_counts < max_blood_alt_counts and 
                tumor_alt_counts > min_tumor_alt_counts and 
                (pd.isna(gnomAD_freq) or gnomAD_freq < max_gnomAD_freq) and 
                (pd.isna(dbSNP_freq) or dbSNP_freq < max_dbSNP_freq)):
                patient_to_tumor_variants_filtered[patient].add((chrom, pos, ref, alt))
    return patient_to_tumor_variants_filtered

def objective(trial):
    # Suggest hyperparameters
    min_dp = trial.suggest_int('min_dp', 1, 10)
    min_alt_reads = trial.suggest_int('min_alt_reads', 1, 10)
    min_af = trial.suggest_float('min_af', 0.1, 0.75, step=0.05)
    max_blood_alt_counts = trial.suggest_int('max_blood_alt_counts', 0, 10)
    min_tumor_alt_counts = trial.suggest_int('min_tumor_alt_counts', 1, 10)
    max_gnomAD_freq = trial.suggest_float('max_gnomAD_freq', 0.0, 5, step=0.1)
    max_dbSNP_freq = trial.suggest_float('max_dbSNP_freq', 0.0, 5, step=0.1)
    
    # Apply filters
    master_dict_filtered = apply_filters_to_all(master_dict_all, min_dp, min_alt_reads, min_af)
    patient_to_tumor_variants_filtered = apply_filters_tumors(patient_to_tumor_variants, max_blood_alt_counts, min_tumor_alt_counts, max_gnomAD_freq, max_dbSNP_freq)
    
    # Calculate F1 score
    f1_score, precision, recall, tp, fp, tn, fn = calculate_f1_score(master_dict_filtered, patient_to_tumor_variants_filtered, ep_cells, pt_cells, patients)
     
    # Log additional metrics
    trial.set_user_attr('precision', precision)
    trial.set_user_attr('recall', recall)
    trial.set_user_attr('true_pos', tp)
    trial.set_user_attr('false_pos', fp)
    trial.set_user_attr('true_neg', tn)
    trial.set_user_attr('false_neg', fn)
    
    return f1_score