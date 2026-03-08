from get_tumor_variants import main_path, patients

from concurrent.futures import ThreadPoolExecutor, as_completed
from cyvcf2 import VCF
import gc
from glob import glob
import os
from tqdm import tqdm

max_workers = int(os.environ.get('VCF_WORKERS', 16))

def process_vcf_all(vcf_path):
    """
    Parse VCF and store ALL variants with metadata (no filtering yet).
    Returns: cell_name, list of (chrom, pos, ref, alt, dp, alt_depth, af) tuples
    """
    cell_name = vcf_path.split('/')[-1].split('.')[0]
    variant_list = []
    cell_vcf = VCF(vcf_path)

    for record in cell_vcf:
        if record.CHROM == "chrM":
            continue

        flt = record.FILTER
        if flt is not None and flt not in (".", "") and flt.lower() != "pass":
            continue

        gt_types = record.gt_types
        gt_type = gt_types[0]
        if gt_type == 0:  # skip homozygous ref
            continue

        ad = record.format("AD")
        dp = record.format("DP")

        if ad is None or dp is None:
            continue

        ad = ad[0]
        dp = dp[0]
        dp_val = int(dp[0]) if hasattr(dp, "__len__") else int(dp)

        if dp_val == 0:
            continue

        for alt_idx, alt in enumerate(record.ALT, start=1):
            if alt_idx >= len(ad):
                continue
            alt_depth = int(ad[alt_idx])
            af = alt_depth / dp_val if dp_val > 0 else 0.0
            variant_list.append((record.CHROM, record.POS, record.REF, alt, dp_val, alt_depth, af))

    return cell_name, variant_list

# Parse all VCFs in parallel — submit ALL patients at once for better utilization
master_dict_all = {p: {} for p in patients}
cell_vcf_path = '/path/to/cell_vcfs/'

# Collect all (vcf_path, patient) pairs upfront
all_jobs = []
for patient in patients:
    fpath = f'{cell_vcf_path}{patient}/'
    vcf_paths = glob(f"{fpath}*.vcf")
    for p in vcf_paths:
        all_jobs.append((p, patient))

print(f"Parsing {len(all_jobs)} VCFs across {len(patients)} patients with {max_workers} workers...")

with ThreadPoolExecutor(max_workers=max_workers) as ex:
    future_to_info = {
        ex.submit(process_vcf_all, vcf_path): (vcf_path, patient)
        for vcf_path, patient in all_jobs
    }

    for future in tqdm(as_completed(future_to_info), total=len(future_to_info), desc="Cells", smoothing=0):
        vcf_path, patient = future_to_info[future]
        try:
            cell_name, variant_list = future.result()
            master_dict_all[patient][cell_name] = variant_list
        except Exception as e:
            print(f'Failed {vcf_path}: {e}')

print("Done parsing VCFs!")

# Report statistics
total_cells = sum(len(cells) for cells in master_dict_all.values())
total_variants = sum(len(variants) for cells in master_dict_all.values() for variants in cells.values())
print(f"Parsed {total_cells} cells with {total_variants} total variant records")

gc.collect()
# Remove duplicate cell names
for patient, cell_dict in master_dict_all.items():
    to_delete = []
    for cell in cell_dict:
        base = cell.replace('_Plate_', '_')
        if '_Plate_' in cell and base in cell_dict:
            to_delete.append(cell)
        if (len(cell.split('_')) > 3) and ('_'.join(cell.split('_')[:3]) in cell_dict):
            to_delete.append(cell)

    for cell in to_delete:
        del cell_dict[cell]

print("Cell names cleaned")
