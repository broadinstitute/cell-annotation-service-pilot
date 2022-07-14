# https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csr_matrix.html
# https://stackoverflow.com/questions/4319014/iterating-through-a-scipy-sparse-vector-or-matrix
# https://github.com/theislab/anndata2ri/blob/master/src/anndata2ri/scipy2ri/py2r.py

import anndata as ad
import argparse
from fastavro import writer, parse_schema
import numpy as np
import time


def current_milli_time():
    return round(time.time() * 1000)


def write_avro(generator, parsed_schema, filename, flush_batch_size, progress_batch_size):
    start = current_milli_time()
    with open(filename, 'a+b') as out:
        records = []
        counter = 0
        for record in generator():
            counter = counter + 1
            records.append(record)

            if counter % flush_batch_size == 0:
                writer(out, parsed_schema, records)
                records = []

            if counter % progress_batch_size == 0:
                end = current_milli_time()
                print(f"    Processed {counter} rows... in {end-start} ms")
                start = end

        if len(records) > 0:
            writer(out, parsed_schema, records)


def dump_core_matrix(x, row_lookup, col_lookup, basename, flush_batch_size, progress_batch_size):
    schema = {
        'doc': 'A raw datum indexed by cell and feature id in the CAS BigQuery schema',
        'namespace': 'cas',
        'name': 'CellFeature',
        'type': 'record',
        'fields': [
            {'name': 'cas_cell_index', 'type': 'int'},
            {'name': 'cas_feature_index', 'type': 'int'},
            {'name': 'raw_counts', 'type': 'int'}
        ]
    }
    parsed_schema = parse_schema(schema)

    def raw_counts_generator():
        cx = x.tocoo(copy=False)

        for i, j, v in zip(cx.row, cx.col, cx.data):
            cas_cell_index = row_lookup[i]
            cas_feature_index = col_lookup[j]
            
            # Todo -- how can we ensure this is safe/right?
            v_int = int(v)
            yield {
                u'cas_cell_index': cas_cell_index.item(),
                u'cas_feature_index': cas_feature_index.item(),
                u'raw_counts': v_int
            }

    write_avro(raw_counts_generator, parsed_schema, f'{basename}_raw_counts.avro',
               flush_batch_size, progress_batch_size)


def dump_cell_info(adata, basename, cas_cell_index_start, flush_batch_size, progress_batch_size):
    # Dump out cell info (obs) -- cell_id is index (26 columns).
    adata.obs['original_cell_id'] = adata.obs.index
    adata.obs['cas_cell_index'] = np.arange(cas_cell_index_start, cas_cell_index_start + len(adata.obs))

    row_index_to_cas_cell_index = [None] * len(adata.obs)
    for row in range(0, len(adata.obs)):
        row_index_to_cas_cell_index[row] = adata.obs['cas_cell_index'].iloc[[row]][0]

    schema = {
        'doc': 'CAS Cell Info',
        'namespace': 'cas',
        'name': 'CellInfo',
        'type': 'record',
        'fields': [
            {'name': 'cas_cell_index', 'type': 'int'},
            {'name': 'original_cell_id', 'type': 'string'},
            {'name': 'cell_type', 'type': 'string'}
        ]
    }

    parsed_schema = parse_schema(schema)

    def cell_generator():
        cells = adata.obs[['cas_cell_index', 'original_cell_id', 'cell_type']]
        for _, r in cells.iterrows():
            yield {
                u'cas_cell_index': r['cas_cell_index'],
                u'original_cell_id': r['original_cell_id'],
                u'cell_type': r['cell_type']
            }

    write_avro(cell_generator, parsed_schema, f'{basename}_cell_info.avro',
               flush_batch_size, progress_batch_size)

    return row_index_to_cas_cell_index


def dump_feature_info(adata, basename, cas_feature_index_start, flush_batch_size, progress_batch_size):
    # Dump out feature info -- feature_id is index (12 columns).
    adata.var['original_feature_id'] = adata.var.index

    adata.var['cas_feature_index'] = np.arange(cas_feature_index_start, cas_feature_index_start + len(adata.var))
    col_index_to_cas_feature_index = [None] * len(adata.var)
    for col in range(0, len(adata.var)):
        col_index_to_cas_feature_index[col] = adata.var['cas_feature_index'].iloc[[col]][0]

    schema = {
        'doc': 'CAS Feature Info',
        'namespace': 'cas',
        'name': 'FeatureInfo',
        'type': 'record',
        'fields': [
            {'name': 'cas_feature_index', 'type': 'int'},
            {'name': 'original_feature_id', 'type': 'string'},
            {'name': 'feature_name', 'type': 'string'}
        ]
    }

    parsed_schema = parse_schema(schema)

    def feature_generator():
        features = adata.var[['cas_feature_index', 'original_feature_id', 'feature_name']]
        for _, row in features.iterrows():
            yield {
                u'cas_feature_index': row['cas_feature_index'],
                u'original_feature_id': row['original_feature_id'],
                u'feature_name': row['feature_name']
            }

    write_avro(feature_generator, parsed_schema, f'{basename}_feature_info.avro',
               flush_batch_size, progress_batch_size)
    return col_index_to_cas_feature_index


def process(input_file, cas_cell_index_start, cas_feature_index_start, basename, flush_batch_size, progress_batch_size):
    basename = "cas" if not basename else basename
    flush_batch_size = 10000 if not flush_batch_size else flush_batch_size
    progress_batch_size = 1000000 if not progress_batch_size else progress_batch_size
    print("Loading data...")
    adata = ad.read(input_file)

    print("Processing cell/observation metadata...")
    row_index_to_cas_cell_index = dump_cell_info(adata, basename, cas_cell_index_start,
                                                 flush_batch_size, progress_batch_size)
    # dump out feature info -- feature_id is index (12 columns)
    print("Processing feature/gene/variable metadata...")
    col_index_to_cas_feature_index = dump_feature_info(adata, basename, cas_feature_index_start,
                                                       flush_batch_size, progress_batch_size)

    print("Processing core data...")
    # recode the indexes to be the indexes of obs/var or should obs/var include these indices?
    dump_core_matrix(adata.raw.X, row_index_to_cas_cell_index, col_index_to_cas_feature_index, basename,
                     flush_batch_size, progress_batch_size)
        
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
    parser.add_argument('--basename', type=str,
                        help='Base filename to use for output files, e.g. <base>_cell_info.avro etc', required=False)
    parser.add_argument('--cas_cell_index_start', type=int, help='starting number for cell index', required=False,
                        default=0)
    parser.add_argument('--cas_feature_index_start', type=int, help='starting number for feature index', required=False,
                        default=0)
    parser.add_argument('--flush_batch_size', type=int, help='max size of Avro batches to flush', required=False)
    parser.add_argument('--progress_batch_size', type=int, help='size of batches on which to report progress',
                        required=False)

    # Execute the parse_args() method
    args = parser.parse_args()

    process(args.input, args.cas_cell_index_start, args.cas_feature_index_start, args.basename,
            args.flush_batch_size, args.progress_batch_size)
