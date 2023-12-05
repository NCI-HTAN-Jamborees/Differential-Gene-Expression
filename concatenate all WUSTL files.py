import os
import scanpy as sp
import pandas as pd
sp.settings.dask_distribute = True



# Create an empty AnnData object to store the combined data
combined_adata = None
data_dir = '../project-files/HTAN_WUSTL/'

# Loop through each sample and load the data
for sample_name in os.listdir(data_dir):
    if sample_name.endswith('.mtx.gz'):
        # Construct file paths for each sample

        base_sample_name = '-'.join(sample_name.split(".")[0].split('-')[0:2])
        print(f"reading {base_sample_name}")
        mtx_file = os.path.join(data_dir, sample_name)

        if (os.path.isfile(os.path.join(data_dir, f'{base_sample_name}-peaks.bed'))):
            print(f"skipping {sample_name} as it is scATAC-Seq")
            continue
        features_file = os.path.join(data_dir, f'{base_sample_name}_features.tsv')
        barcodes_file = os.path.join(data_dir, f'{base_sample_name}_barcodes.tsv')

        # Read the Matrix Market (.mtx) file using scipy
        adata = sp.read_mtx(mtx_file).transpose()

        # Read features and barcodes using scanpy
        adata.var_names = pd.read_csv(features_file, header=None, index_col=0).index
        adata.obs_names = pd.read_csv(barcodes_file, header=None, index_col=0).index

        # Optionally, you can read additional metadata files if available
        # metadata_file = os.path.join(data_dir, f'{sample_name.split(".")[0]}_metadata.tsv')
        # adata.obs['metadata'] = pd.read_csv(metadata_file, sep='\t', index_col=0)

        # Concatenate the current sample to the combined AnnData object
        if combined_adata is None:
            combined_adata = adata
        else:
            combined_adata = combined_adata.concatenate(adata, index_unique=None, join='outer')

print('done reading all files')
print(combined_adata)
