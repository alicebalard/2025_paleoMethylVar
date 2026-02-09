#!/usr/bin/env python3
import h5py
import numpy as np
import pandas as pd
import gc

# Paths
file_path = "/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/10X/all_matrix_noscale.h5"
metadata_path = "/SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/05_hvCpGalgorithm/dataPrev/SupTab1_Loyfer2023.csv"
out_path = "/SAN/ghlab/epigen/Alice/hvCpG_project/data/WGBS_human/AtlasLoyfer/10X/mean_cpg_std_across_all_datasets.csv"

# Load metadata
print("Loading metadata...")
sample_metadata = pd.read_csv(metadata_path)
sample_metadata["dataset"] = (
    sample_metadata["Source Tissue"].astype(str) + " - " +
    sample_metadata["Cell type"].astype(str)
)

# Within each dataset, keep only one row per PatientID
meta_unique_pat = sample_metadata.drop_duplicates(
    subset=["dataset", "PatientID"], keep="first"
)
print(f"Found {len(meta_unique_pat)} unique dataset-patient combinations")

# Initialize accumulators (float32 & int16 to minimize memory)
with h5py.File(file_path, 'r') as h5file:
    h5_samples = h5file['samples'][:].astype(str)
    cpg_names = h5file['cpg_names'][:]
    matrix = h5file['matrix']
    n_rows = matrix.shape[0]
    
    print(f"Processing {n_rows:,} CpGs across {len(h5_samples)} samples")
    
    # Memory-efficient accumulators (20M * 4B + 20M * 2B = ~120MB total)
    means = np.zeros(n_rows, dtype='f4')  # Sum of std values
    counts = np.zeros(n_rows, dtype='i2') # Number of valid datasets per CpG
    
    batch_size = 20000  # Tune based on available RAM
    
    # Process each dataset
    for i, (dataset_label, df_ds) in enumerate(meta_unique_pat.groupby("dataset")):
        print(f"[{i+1}] Processing '{dataset_label}': {len(df_ds)} samples")
        
        # Extract Z... sample IDs (from your prior debugging)
        ds_samples_raw = df_ds["Sample name"].str.split('-Z', expand=True)[1]
        ds_samples = ['Z' + x if pd.notna(x) else '' for x in ds_samples_raw]
        ds_samples = [s.strip() for s in ds_samples if s]
        
        # Find matching columns in HDF5
        mask = np.isin(h5_samples, ds_samples)
        sample_indices = np.where(mask)[0]
        
        if sample_indices.size == 0:
            print(f"  Warning: no matching samples found, skipping")
            continue
        
        print(f"  Computing std for {sample_indices.size} matching samples...")
        
        # Compute std in batches (key memory optimization)
        std_values = np.full(n_rows, np.nan, dtype='f4')
        for start in range(0, n_rows, batch_size):
            end = min(start + batch_size, n_rows)
            chunk = matrix[start:end, sample_indices]  # h5py lazy loads
            std_values[start:end] = np.nanstd(chunk, axis=1, dtype='f4')
        
        # Accumulate properly without nan_to_num dtype issue
        valid_mask = ~np.isnan(std_values)
        means[valid_mask] += std_values[valid_mask]
        counts[valid_mask] += 1
        
        # Cleanup
        del std_values, chunk, mask
        gc.collect()
    
    print("Finalizing results...")
    # Compute final means only for CpGs with data
    valid_cpgs = counts > 0
    mean_std_values = means[valid_cpgs] / counts[valid_cpgs]
    
    # Handle cpg_names properly (avoid decode issues)
    valid_cpg_names = []
    for cpg in cpg_names[valid_cpgs]:
        if isinstance(cpg, bytes):
            valid_cpg_names.append(cpg.decode('utf-8'))
        else:
            valid_cpg_names.append(str(cpg))
    
    # Create final DataFrame (much smaller than full 20M)
    results = pd.DataFrame({
        'CpG_Name': valid_cpg_names,
        'Mean_Standard_Deviation': mean_std_values,
        'Datasets_with_data': counts[valid_cpgs]
    })
    
    print(f"Saving {len(results):,} CpGs with data to {out_path}")
    results.to_csv(out_path, index=False)
    print("Done!")
