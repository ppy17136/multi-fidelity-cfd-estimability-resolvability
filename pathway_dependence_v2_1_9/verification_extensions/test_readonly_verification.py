import sys
from pathlib import Path
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from compare_reconstruction import compare_values

class ComparisonTests(unittest.TestCase):
    def test_last_digit_accepted(self):
        compare_values({'diameter':0.002541864613862918},{'diameter':0.0025418646138629154})
    def test_large_difference_rejected(self):
        with self.assertRaises(AssertionError):
            compare_values({'diameter':0.0025},{'diameter':0.0026})
    def test_certificate_cost_not_rounded(self):
        with self.assertRaises(AssertionError):
            compare_values({'exact_diameter_cost':0.5},{'exact_diameter_cost':1.0})
    def test_witness_not_ignored(self):
        with self.assertRaises(AssertionError):
            compare_values({'partition':[[0,1],[2]]},{'partition':[[0],[1,2]]})
    def test_counts_and_flags_exact(self):
        with self.assertRaises(AssertionError):
            compare_values({'qualifies':True},{'qualifies':False})
    def test_schema_changes_rejected(self):
        with self.assertRaises(AssertionError):
            compare_values({'diameter':1.0},{'diameter':1.0,'extra':0})

if __name__=='__main__':
    unittest.main(verbosity=2)
