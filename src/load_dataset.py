import argparse
from google.api_core.exceptions import Conflict
from google.cloud import bigquery
from google.cloud import storage
import os


def create_table(client, project, dataset, tablename, schema, clustering_fields):
    table_id = f"{project}.{dataset}.{tablename}"

    table = bigquery.Table(table_id, schema=schema)
    if clustering_fields:
        table.clustering_fields = clustering_fields

    try:
        _ = client.create_table(table)  # Make an API request.
        print(f"Created '{table_id}'.")
    except Conflict:
        print(f"Table '{table_id}' exists, continuing.")


def create_dataset(client, project, dataset, location):
    dataset_id = f"{project}.{dataset}"

    # Construct a full Dataset object to send to the API.
    dataset = bigquery.Dataset(dataset_id)
    dataset.location = location

    # Send the dataset to the API for creation, with an explicit timeout.
    # Raises google.api_core.exceptions.Conflict if the Dataset already
    # exists within the project.
    try:
        _ = client.create_dataset(dataset, timeout=30)  # Make an API request.
        print(f"Created dataset {dataset_id}.")
    except Conflict:
        print(f"Dataset {dataset_id} exists, continuing.")


def check_avro_files_exist(avro_prefix):
    file_types = ['cell_info', 'feature_info', 'raw_counts']
    filenames = [f'{avro_prefix}_{file_type}.avro' for file_type in file_types]
    missing = list(filter(lambda f: not os.path.exists(f), filenames))
    if len(missing) > 0:
        raise ValueError(
            f"Missing Avro files for loading to BigQuery: {', '.join(missing)}")
    return filenames


def bucket_and_prefix(project, gcs_prefix):
    import re
    client = storage.Client(project=project)
    (bucket_name, object_prefix) = re.search(r'gs://([^/]+)/(.*)$', gcs_prefix).groups()
    bucket = client.bucket(bucket_name)
    return bucket, object_prefix.rstrip('/')


def process(project, dataset, avro_prefix, gcs_prefix, force_bq_append):
    client = bigquery.Client(project=project)
    create_dataset(client, project, dataset, "US")

    create_table(client, project, dataset, "cas_cell_info",
                 [
                     bigquery.SchemaField("cas_cell_index", "INTEGER", mode="REQUIRED"),
                     bigquery.SchemaField("original_cell_id", "STRING", mode="REQUIRED"),
                     bigquery.SchemaField("cell_type", "STRING", mode="REQUIRED")
                 ],
                 []
                 )

    create_table(client, project, dataset, "cas_feature_info",
                 [
                     bigquery.SchemaField("cas_feature_index", "INTEGER", mode="REQUIRED"),
                     bigquery.SchemaField("original_feature_id", "STRING", mode="REQUIRED"),
                     bigquery.SchemaField("feature_name", "STRING", mode="REQUIRED")
                 ],
                 []
                 )

    create_table(client, project, dataset, "cas_raw_count_matrix",
                 [
                     bigquery.SchemaField("cas_cell_index", "INTEGER", mode="REQUIRED"),
                     bigquery.SchemaField("cas_feature_index", "INTEGER", mode="REQUIRED"),
                     bigquery.SchemaField("raw_counts", "INTEGER", mode="REQUIRED")
                 ],
                 ["cas_cell_index"]
                 )

    (cell_filename, feature_filename, raw_counts_filename) = check_avro_files_exist(avro_prefix)
    # noinspection SqlResolve
    query = f"""select table_name, total_rows from `{dataset}.INFORMATION_SCHEMA.PARTITIONS` where total_rows > 0"""
    job = client.query(query)
    tables_with_data = [r[0] for r in list(job.result())]

    (bucket, object_prefix) = bucket_and_prefix(project, gcs_prefix)

    def stage_file(file_to_stage):
        blob = bucket.blob(f'{object_prefix}/{file_to_stage}')
        blob.upload_from_filename(file_to_stage)
        print(f"Staged '{file}' to '{gcs_prefix}/{file}'.")

    def unstage_file(file_to_unstage):
        blob = bucket.blob(f'{object_prefix}/{file_to_unstage}')
        blob.delete()
        print(f"Removed staged file '{gcs_prefix}/{file}'.")

    staged_files = []
    gcs_prefix = gcs_prefix.rstrip('/')
    pairs = [
        ("cell_info", cell_filename),
        ("feature_info", feature_filename),
        ("raw_count_matrix", raw_counts_filename)
    ]
    job_config = bigquery.LoadJobConfig(source_format=bigquery.SourceFormat.AVRO)

    for table, file in pairs:
        table = f"cas_{table}"
        table_id = f'{project}.{dataset}.{table}'

        if table in tables_with_data and not force_bq_append:
            print(f"Table '{table_id}' contains data and `--force_bq_append` not specified, skipping load of '{file}'.")
            continue

        stage_file(file)
        staged_files.append(file)

        uri = f'{gcs_prefix}/{file}'
        load_job = client.load_table_from_uri(
            uri, table_id, job_config=job_config
        )  # Make an API request.
        load_job.result()  # Waits for the job to complete.
        destination_table = client.get_table(table_id)
        print(f"{destination_table.num_rows} rows loaded into BigQuery table '{table_id}'.")

    for f in staged_files:
        unstage_file(f)

    print("Done.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(allow_abbrev=False, description='Initialize CASP tables')

    parser.add_argument('--project', type=str, help='BigQuery Project', required=True)
    parser.add_argument('--dataset', type=str, help='BigQuery Dataset', required=True)
    parser.add_argument('--avro_prefix', type=str, help='Prefix with which Avro files are named', required=True)
    parser.add_argument('--gcs_prefix', type=str, help='GCS prefix to which Avro files should be staged', required=True)
    parser.add_argument('--force_bq_append', type=bool,
                        help='Append data to BigQuery tables even if data some data is already loaded', required=False)

    args = parser.parse_args()

    process(args.project, args.dataset, args.avro_prefix, args.gcs_prefix, args.force_bq_append)
