# Cell Annotation Service Pilot

This repository is intended to hold the functional prototype for a single cell annotation service

### Prerequisites / Installation

 - Python 3.7+
 
### Testing Data

For local testing, we are using a small data set from CellXGene.

[Horizontal Cells in Human Retina](https://cellxgene.cziscience.com/collections/af893e86-8e9f-41f1-a474-ef05359b1fb7) - 7,348 cells

And a larger data set

[Tabula Sapiens Endothelial](https://cellxgene.cziscience.com/collections/e5f58829-1a66-40b5-a624-9046778e74f5) -  31,691 cells


### Convert AnnData to Avro Load Format

First convert from AnnData to the Avro format that will be used to load BigQuery:

```
python src/anndata_to_avro.py --input data/horizontal-cells-in-human-retina.h5ad --avro_prefix retina --cas_cell_index_start 1000 --cas_feature_index_start 5000
```

NOTE: if multiple data sets are being loaded into the same dataset, offsets for the cell and feature indexes must be provided via the `cas_cell_index_start` and `cas_feature_index_start` parameters respectively.

### Initialize dataset and ingest Avro files into BigQuery

Continuing this example:

```
python src/load_dataset.py --project <google-project> --dataset <dataset-name> --avro_prefix retina --gcs_prefix <gcs-prefix-to-avro-files>
```

All of these parameters are required.

* `project` specifies the Google project id that will contain the BigQuery dataset.
* `dataset` specifies the dataset name that will contain the BigQuery tables.
* `avro_prefix` specifies the prefix used to name the output Avro from the AnnData conversion in the previous step.
* `gcs_prefix` specifies the GCS prefix to which the Avro files should be staged for BigQuery ingestion.

This step:

* Creates the dataset in this specified project
* Creates the CASP tables in the specified dataset
* Stages the specified Avro files to the specified GCS path
* Loads the specified tables from the staged Avro files
* Cleans up the staged Avro files from GCS


### Extract Random Subset

This will randomly get a specified number of cells' data from your BigQuery dataset:

```
python src/random_bq_to_anndata.py --project <google-project> --dataset <dataset-name> --num_cells <num>
```

TODO: use data returned from BigQuery to output an AnnData file.
