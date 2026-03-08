from get_tumor_variants import patients
from get_cell_variants import master_dict_all

import anndata as ad
from glob import glob
import pandas as pd

adata = ad.read_h5ad('/path/to/data_analysis/anndata_objects/anndata_file.h5ad')
print(f"Loaded AnnData: {adata.n_obs} cells × {adata.n_vars} genes")

wo_promoter = [e for e in adata.var_names if e.startswith('IGHD') or e.startswith('IGHJ') or e.startswith('IGKJ') or e.startswith('IGLJ')]
wo_promoter.remove('IGHD')
adata.obs['dj_positive'] = (adata[:, wo_promoter].layers['decontx_counts'].toarray() >= 1).any(axis=1)
bcells = adata[adata.obs.dj_positive]
print(f"Identified {len(bcells)} B cells")

batch_df_paths = glob('/path/to/batch_metadata/*combined*.xlsx')

cell_to_type = {}
for pathname in batch_df_paths:
    if 'batch' not in pathname.lower():
        continue
    batch = int(pathname.split('tch_')[-1].split('_')[0])
    if batch == 4:
        df1 = pd.read_excel(pathname, engine='openpyxl', sheet_name='Execution plan 4a')
        df2 = pd.read_excel(pathname, engine='openpyxl', sheet_name='Execution plan 4b')
        df = pd.concat([df1, df2])
        del df1, df2
    else:
        try:
            df = pd.read_excel(pathname, engine='openpyxl')
        except:
            continue

    df.columns = df.columns.str.lower()
    if not {'patient id', 'plate #', 'old well position'}.issubset(set(df.columns)):
        continue
    if 'cell details' in df.columns:
        cell_to_type.update(df.set_index(['patient id', 'plate #', 'old well position'])['cell details'].to_dict())
    else:
        idx = df.set_index(['patient id', 'plate #', 'old well position']).index
        cell_to_type.update(dict.fromkeys(idx, 'bone marrow'))

cell_to_type = {'_'.join(map(str, k)):v for k,v in cell_to_type.items()}
cell_to_type = {k: 'primary tumor' if 'tumor' in str(v).lower() and 'marrow' not in str(v).lower() else 'bone marrow' for k, v in cell_to_type.items()}

print(f"Loaded cell type metadata for {len(cell_to_type)} cells")

ep_cells = set(adata[(adata.obs.wu_label.isin({'Normal Epithelial', 'Cancer Epithelial'})) & (adata.obs.dj_positive == False)].obs.index)
bm_cells, pt_cells = set(), set()
for patient in master_dict_all.keys():
    for cell in master_dict_all[patient].keys():
        abrv_cell_name = ('_').join(cell.split('_')[:3])
        if abrv_cell_name in cell_to_type:
            if cell_to_type[abrv_cell_name] == 'primary tumor':
                if cell in ep_cells:
                    pt_cells.add(cell)
            else:
                bm_cells.add(cell)

# Add cells from plates 5 and 6 to PT if not already classified
for patient in patients:
    pt_cells.update([e for e in master_dict_all[patient] if e not in bm_cells.union(pt_cells) and ('_5_' in e or '_6_' in e)])

print(f"Classified {len(bm_cells)} BM cells and {len(pt_cells)} PT cells")