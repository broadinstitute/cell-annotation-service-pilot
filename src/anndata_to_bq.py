# https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csr_matrix.html
# https://stackoverflow.com/questions/4319014/iterating-through-a-scipy-sparse-vector-or-matrix
# https://github.com/theislab/anndata2ri/blob/master/src/anndata2ri/scipy2ri/py2r.py

import argparse
import anndata as ad
import gzip
import numpy as np
import time


def current_milli_time():
    return round(time.time() * 1000)


def dump_core_matrix(x, row_lookup, col_lookup):
    counter = 0
    start = current_milli_time()
    with gzip.open('cas_raw_counts.tsv.gz', 'wt') as f:
        f.write("cas_cell_index\tcas_feature_index\tdata\n")

        for r in range(x.shape[0]):
            row = x[r]
            for c in row.indices:
                cas_cell_index = row_lookup[r]
                cas_feature_index = col_lookup[c]
                v = x[r, c]
            
                # Todo -- how can we ensure this is safe/right?
                v_int = int(v)
                f.write(f"{cas_cell_index}\t{cas_feature_index}\t{v_int}\n")
                counter = counter + 1
                if counter % 1000 == 0:
                    end = current_milli_time()
                    print(f"    Processed {counter} rows... in {end-start} ms")
                    start = end


def process(input_file, cas_cell_index_start, cas_feature_index_start):
    print("Loading data...")
    adata = ad.read(input_file, backed='r')

    # dump out cell info (obs) -- cell_id is index (26 columns)
    print("Processing cell/observation metadata...")
    adata.obs['original_cell_id'] = adata.obs.index    
    adata.obs['cas_cell_index'] = np.arange(cas_cell_index_start, cas_cell_index_start + len(adata.obs))

    row_index_to_cas_cell_index = [None] * len(adata.obs)
    for row in range(0, len(adata.obs)):
        row_index_to_cas_cell_index[row] = adata.obs['cas_cell_index'].iloc[[row]][0]
    
    adata.obs[['cas_cell_index', 'original_cell_id', 'cell_type']].to_csv(
        'cas_cell_info.tsv.gz', index=False, sep='\t', compression='gzip')
    
    # dump out feature info -- feature_id is index (12 columns)
    print("Processing feature/gene/variable metadata...")
    adata.var['original_feature_id'] = adata.var.index
    
    adata.var['cas_feature_index'] = np.arange(cas_feature_index_start, cas_feature_index_start + len(adata.var))
    col_index_to_cas_feature_index = [None] * len(adata.var)
    for col in range(0, len(adata.var)):
        col_index_to_cas_feature_index[col] = adata.var['cas_feature_index'].iloc[[col]][0]
    
    adata.var[['cas_feature_index', 'original_feature_id', 'feature_name']].to_csv(
        'cas_feature_info.tsv.gz', index=False, sep='\t', compression='gzip')
    
    print("Processing core data...")
    
    # recode the indexes to be the indexes of obs/var or should obs/var include these indices?
    dump_core_matrix(adata.raw.X, row_index_to_cas_cell_index, col_index_to_cas_feature_index)
        
# TODO: naive dump, compare
#    rows,cols = adata.X.nonzero()
#    for row,col in zip(rows,cols):
#        cas_cell_index = row_index_to_cas_cell_index[row]
#        cas_feature_index = col_index_to_cas_feature_index[col]
#        v = adata.X[row,col]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description='Convert AnnData Single Cell Expression Data into format for loading into BQ')

    parser.add_argument('--input', type=str, help='AnnData format input file', required=True)
    parser.add_argument('--cas_cell_index_start', type=int, help='starting number for cell index', required=False,
                        default=0)
    parser.add_argument('--cas_feature_index_start', type=int, help='starting number for feature index', required=False,
                        default=0)

    # Execute the parse_args() method
    args = parser.parse_args()

    process(args.input, args.cas_cell_index_start, args.cas_feature_index_start)
