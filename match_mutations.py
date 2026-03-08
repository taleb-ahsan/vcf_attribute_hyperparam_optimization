from get_cell_variants import master_dict_all
from get_tumor_variants import patient_to_tumor_variants
from get_blood_variants import patient_to_blood_variants

def calculate_f1_score(master_dict_all, patient_to_tumor_variants, ep_cells, pt_cells, patients):
    """
    Calculate F1 score for a given set of filtered variants.
    master_dict_all: dict[patient][cell] -> set of (chrom, pos, ref, alt)
    """
    # Separate by tissue type
    master_dict_bm = {}
    master_dict_pt = {}
    
    for patient, cells_dict in master_dict_all.items():
        bm_dict = {}
        pt_dict = {}
        for cell_name, variants in cells_dict.items():
            if cell_name in pt_cells:
                pt_dict[cell_name] = variants
            else:
                bm_dict[cell_name] = variants
        master_dict_bm[patient] = bm_dict
        master_dict_pt[patient] = pt_dict
    
    # Calculate matches
    def match_mutations(x, y, master_dict_all):
        if x == y:
            target = patient_to_tumor_variants[y]
        else:
            target = patient_to_tumor_variants[y] - patient_to_tumor_variants[x]
        
        target = target - patient_to_blood_variants[x]

        num_matching_variants = []
        if x == y:
            cell_names = [e for e in list(master_dict_pt[x]) if e in ep_cells]
        else:
            cell_names = list(master_dict_all[x])
        
        for cell_name in cell_names:
            cell_set = master_dict_all[x][cell_name]
            if not cell_set:
                continue
            num_matching = len(cell_set & target)
            num_matching_variants.append(num_matching)
        
        return num_matching_variants
    
    same_pairs = [(p, p) for p in patients]
    diff_pairs = [(a, b) for a in patients for b in patients if a != b]
    
    num_same_matching = []
    num_diff_matching = []
    
    for x, y in same_pairs:
        num = match_mutations(x, y, master_dict_all)
        num_same_matching.extend(num)
    
    for x, y in diff_pairs:
        num = match_mutations(x, y, master_dict_all)
        num_diff_matching.extend(num)
    
    true_pos = len([e for e in num_same_matching if e > 0])
    false_pos = len([e for e in num_diff_matching if e > 0])
    true_neg = len([e for e in num_diff_matching if e == 0])
    false_neg = len([e for e in num_same_matching if e == 0])
    
    precision = true_pos / (true_pos + false_pos) if (true_pos + false_pos) > 0 else 0.0
    recall = true_pos / (true_pos + false_neg) if (true_pos + false_neg) > 0 else 0.0
    f1_score = (2 * precision * recall / (precision + recall)) if (precision + recall) > 0 else 0.0
    
    return f1_score, precision, recall, true_pos, false_pos, true_neg, false_neg