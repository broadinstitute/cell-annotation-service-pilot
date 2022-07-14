from google.cloud import bigquery
import argparse
import random
import numpy as np
import anndata as ad
from scipy.sparse import csr_matrix

# assumes that cas_cell_info.cas_feature_index values are a contiguous list of ints starting at 0
def get_random_ids(project, dataset, client, num_cells):
    query = client.query(f"SELECT MIN(cas_cell_index) AS min_table_number, MAX(cas_cell_index) AS max_table_number FROM `{project}.{dataset}.cas_cell_info`")
    min_cell_id = int([row.min_table_number for row in list(query.result())][0])
    max_cell_id = int([row.max_table_number for row in list(query.result())][0])
    print(f"Getting {num_cells} random IDs between {min_cell_id} and {max_cell_id}...")
    cell_ids = list(range(min_cell_id, max_cell_id + 1))
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
    cell_data = get_cell_data(project, dataset, client, num_cells)

    col_count = 0
    feature_to_column = {}
    row_num = -1
    last_original_cell_id = None
    last_cell_type = None

    indptr = [0]
    indices = []
    data = []
    cell_names = []
    cell_types = []
    feature_ids = []
    feature_names = []

    for row in list(cell_data):
        # print("cas_cell_index={}, original_cell_id={}, cell_type={}, cas_feature_index={}, original_feature_id={}, feature_name={}, count={}".format(row["cas_cell_index"], row["original_cell_id"], row["cell_type"], row["cas_feature_index"], row["original_feature_id"], row["feature_name"], row["count"]))
        original_cell_id = row["original_cell_id"]
        cell_type = row["cell_type"]
        original_feature_id = row["original_feature_id"]
        feature_name = row["feature_name"]
        count = row["count"]

        if original_cell_id != last_original_cell_id:
            if last_original_cell_id is not None:
                # We have just started a New row (and it's not the very first), update indptr to point to the next row
                row_num += 1
                cell_names.append(last_original_cell_id)
                cell_types.append(last_cell_type)
                indptr.append(len(indices))
                print(f"Row Number: {row_num}, original_cell_id: {last_original_cell_id}")
                print(f"Indptr: {indptr}")
            last_original_cell_id = original_cell_id
            last_cell_type = cell_type

        if original_feature_id not in feature_to_column:
            feature_to_column[original_feature_id] = col_count
            feature_ids.append(original_feature_id)
            feature_names.append(feature_name)
            col_count += 1
        col_num = feature_to_column[original_feature_id]

        indices.append(col_num)
        data.append(count)

    # Deal with the last row.
    row_num += 1
    cell_names.append(last_original_cell_id)
    cell_types.append(last_cell_type)
    indptr.append(len(indices))
    print(f"Row Number: {row_num}, original_cell_id: {last_original_cell_id}")
    print(f"Indptr: {indptr}")

    counts = csr_matrix((data, indices, indptr), dtype=np.float32)
    adata = ad.AnnData(counts)
    adata.obs_names = cell_names
    adata.obs["cell_type"] = cell_types
    adata.var_names = feature_ids
    adata.var["feature_name"] = feature_names
    adata.write(f'{output_file_prefix}.h5ad', compression="gzip")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False, description='Query CASP tables for random cells')
    parser.add_argument('--project', type=str, help='BigQuery Project', required=True)
    parser.add_argument('--dataset', type=str, help='BigQuery Dataset', required=True)
    parser.add_argument('--num_cells', type=int, help='Number of cells to return', required=True)
    parser.add_argument('--output_file_prefix', type=str, help='The prefix of the anndata (.h5ad) file that will be created', required=True)

    args = parser.parse_args()
    random_bq_to_anndata(args.project, args.dataset, args.num_cells, args.output_file_prefix)
