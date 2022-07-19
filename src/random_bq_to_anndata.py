from google.cloud import bigquery
import argparse
import random
import numpy as np
import anndata as ad
from scipy.sparse import csr_matrix

# assumes that cas_cell_info.cas_feature_index values are a contiguous list of ints
def get_random_ids(project, dataset, client, num_cells):
    query = client.query(f"SELECT MIN(cas_cell_index) AS min_cas_cell_index, MAX(cas_cell_index) AS max_cas_cell_index FROM `{project}.{dataset}.cas_cell_info`")
    row = list(query.result())[0]
    min_cas_cell_index, max_cas_cell_index = row.min_cas_cell_index, row.max_cas_cell_index
    print(f"Getting {num_cells} random IDs between {min_cas_cell_index} and {max_cas_cell_index}...")
    cell_ids = list(range(min_cas_cell_index, max_cas_cell_index + 1))
    random.shuffle(cell_ids)
    del cell_ids[num_cells:]
    print(f"Random IDs: {cell_ids}")
    return cell_ids


def get_cell_data(project, dataset, client, num_cells):
    random_ids = get_random_ids(project, dataset, client, num_cells)
    in_clause = f" matrix.cas_cell_index IN ({','.join(map(str, random_ids))})"

    # at some point, we will probably want create temp table of cell_ids and then JOIN on it
    # instead of an IN clause
    sql = f"SELECT matrix.cas_cell_index, original_cell_id, cell_type, matrix.cas_feature_index, original_feature_id, feature_name, raw_counts AS count FROM `{project}.{dataset}.cas_cell_info` AS cell, `{project}.{dataset}.cas_feature_info` AS feature, `{project}.{dataset}.cas_raw_count_matrix` AS matrix WHERE matrix.cas_cell_index = cell.cas_cell_index AND matrix.cas_feature_index = feature.cas_feature_index AND" + in_clause + " ORDER BY matrix.cas_cell_index, matrix.cas_feature_index"
    print(f"Getting {num_cells} random cells' data from {project}.{dataset}...")
    query = client.query(sql)
    return query.result()


def random_bq_to_anndata(project, dataset, num_cells, output_file_prefix):
    client = bigquery.Client(project=project)
    # Note that this method requires that the result set returned by get_cell_data be sorted by cas_cell_index (or, could also be sorted by original_cell_id)
    cell_data = get_cell_data(project, dataset, client, num_cells)

    (index_ptr, indices, data, cell_names, cell_types, feature_ids, feature_names) = generate_sparse_matrix(cell_data)

    # Create the matrix from the sparse data representation generated above.
    counts = csr_matrix((data, indices, index_ptr), dtype=np.float32)
    adata = ad.AnnData(counts)
    adata.obs_names = cell_names
    adata.obs["cell_type"] = cell_types
    adata.var_names = feature_ids
    adata.var["feature_name"] = feature_names

    # See https://anndata.readthedocs.io/en/latest/generated/anndata.AnnData.raw.html?highlight=raw#anndata.AnnData.raw
    #  for why we set 'raw' thusly: "The raw attribute is initialized with the current content of an object by setting:"
    adata.raw = adata
    adata.write(f'{output_file_prefix}.h5ad', compression="gzip")

def generate_sparse_matrix(cell_data):
    column_count = 0
    feature_to_column = {}
    last_cell_type = None

    # For representation of the sparse matrix
    index_ptr = [0]
    indices = []
    data = []

    # For storage of the metadata
    cell_names = []
    cell_types = []
    feature_ids = []
    feature_names = []

    # Builds the sparse matrix for storing the original_cell_ids (as obs) and features (as var).
    # From: https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csr_matrix.html
    # If we want to store the following sparse matrix:
    # OBS \ VAR
    #       ENS1 ENS2 ENS3 ENS4
    # AAA   1    7     2   n/a
    # AAC   0    n/a   4   n/a
    # AAG   3    n/a   0   1
    #
    # we represent it as:
    # index_ptr = [0, 3, 5, 8]
    # indices   = [0, 1, 2, 0, 2, 0, 2, 3]
    # data      = [1, 7, 2, 0, 4, 3, 0, 1]

    last_original_cell_id = None
    for row in cell_data:
        # print("cas_cell_index={}, original_cell_id={}, cell_type={}, cas_feature_index={}, original_feature_id={}, feature_name={}, count={}".format(row["cas_cell_index"], row["original_cell_id"], row["cell_type"], row["cas_feature_index"], row["original_feature_id"], row["feature_name"], row["count"]))
        # print(f"{row}")
        original_cell_id = row["original_cell_id"]
        cell_type = row["cell_type"]
        original_feature_id = row["original_feature_id"]
        count = row["count"]

        if original_cell_id != last_original_cell_id:
            if last_original_cell_id is not None:
                # We have just started reading data for a new 'original_cell_id'.
                index_ptr.append(len(indices))
                cell_names.append(last_original_cell_id)
                cell_types.append(last_cell_type)
                print(f"finishing record for original_cell_id: {last_original_cell_id}")
                print(f"index_ptr: {index_ptr}")
            last_original_cell_id = original_cell_id
            last_cell_type = cell_type

        if original_feature_id not in feature_to_column:
            # We have not seen this 'original_feature_id' before.
            feature_to_column[original_feature_id] = column_count
            feature_ids.append(original_feature_id)
            feature_names.append(row["feature_name"])
            column_count += 1
        col_num = feature_to_column[original_feature_id]

        indices.append(col_num)
        data.append(count)

    #TODO - catch the edge case if the last row is different from last-1
    # And add a test - especially for last-1
    # Deal with the last row.
    index_ptr.append(len(indices))
    cell_names.append(last_original_cell_id)
    cell_types.append(last_cell_type)
    print(f"finishing record for original_cell_id: {last_original_cell_id}")
    print(f"index_ptr: {index_ptr}")
    # print(f"indices:   {indices}")
    # print(f"data:      {data}")

    return index_ptr, indices, data, cell_names, cell_types, feature_ids, feature_names


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False, description='Query CASP tables for random cells')
    parser.add_argument('--project', type=str, help='BigQuery Project', required=True)
    parser.add_argument('--dataset', type=str, help='BigQuery Dataset', required=True)
    parser.add_argument('--num_cells', type=int, help='Number of cells to return', required=True)
    parser.add_argument('--output_file_prefix', type=str, help='The prefix of the anndata (.h5ad) file that will be created', required=True)

    args = parser.parse_args()
    random_bq_to_anndata(args.project, args.dataset, args.num_cells, args.output_file_prefix)
