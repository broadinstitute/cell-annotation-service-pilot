import unittest
import numpy

from random_bq_to_anndata import generate_sparse_matrix

class MockCasRow(object):
    original_cell_id = ""
    cell_type = ""

    # The class "constructor" - It's actually an initializer
    def __init__(self, original_cell_id, cell_type):
        self.original_cell_id = original_cell_id
        self.cell_type = cell_type

def make_cas_row(original_cell_id, cell_type):
    return MockCasRow(original_cell_id, cell_type)


class TestRandomBqToAnndata(unittest.TestCase):
    def test_generate_sparse_matrix(self):

        #       ENS1 ENS2 ENS3 ENS4
        # AAA   1    7     2   n/a
        # AAC   0    n/a   4   n/a
        # AAG   3    n/a   0   1
        # AAT   n/a  n/a   n/a 9
        #
        # we represent it as:
        # index_ptr = [0, 3, 5, 8, 9]
        # indices   = [0, 1, 2, 0, 2, 0, 2, 3, 3]
        # data      = [1, 7, 2, 0, 4, 3, 0, 1, 9]

        cell_ids = ["AAA", "AAC", "AAG", "AAT"]
        feature_ids = ["ENS1", "ENS2", "ENS3", "ENS4"]

        matrix = numpy.full([len(cell_ids), len(feature_ids)], None)
        matrix[0][0] = 1
        matrix[0][1] = 7
        matrix[0][2] = 2

        matrix[1][0] = 0
        matrix[1][2] = 4

        matrix[2][0] = 3
        matrix[2][2] = 0
        matrix[2][3] = 1

        matrix[3][3] = 9
        print(f"{matrix[2][1]}")
        print(f"{matrix[2][2]}")

        cas_rows = []
        for i in range (0, len(cell_ids)):
            cell_id = cell_ids[i]
            cell_type = "CT_" + cell_id
            for j in range (0, len(feature_ids)):
                if (matrix[i][j] is not None):
                    cas_row = {}
                    cas_row["original_cell_id"] = cell_id
                    cas_row["cell_type"] = cell_type
                    cas_row["original_feature_id"] = feature_ids[j]
                    cas_row["feature_name"] = "FT_" + feature_ids[j]
                    cas_row["count"] = matrix[i][j]
                    cas_rows.append(cas_row)

        expected_index_ptr = [0, 3, 5, 8, 9]
        expected_indices   = [0, 1, 2, 0, 2, 0, 2, 3, 3]
        expected_data      = [1, 7, 2, 0, 4, 3, 0, 1, 9]
        expected_cell_ids  = ["AAA", "AAC", "AAG", "AAT"]
        expected_cell_types  = ["CT_AAA", "CT_AAC", "CT_AAG", "CT_AAT"]
        expected_feature_ids = ["ENS1", "ENS2", "ENS3", "ENS4"]
        expected_feature_names = ["FT_ENS1", "FT_ENS2", "FT_ENS3", "FT_ENS4"]

        (index_ptr, indices, data, cell_names, cell_types, feature_ids, feature_names) = generate_sparse_matrix(cas_rows)
        self.assertEqual(index_ptr, expected_index_ptr)
        self.assertEqual(indices, expected_indices)
        self.assertEqual(data, expected_data)
        self.assertEqual(cell_names, expected_cell_ids)
        self.assertEqual(cell_types, expected_cell_types)
        self.assertEqual(feature_ids, expected_feature_ids)
        self.assertEqual(feature_names, expected_feature_names)



        # self.assertEqual(value, 7)

def main():
    unittest.main()

if __name__ == "__main__":
    main()