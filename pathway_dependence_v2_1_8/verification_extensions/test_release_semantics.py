"""Release routing checks; numerical correctness is tested by existing suites."""
import sys,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import run_all
class ReleaseRouting(unittest.TestCase):
    def test_default_is_manuscript_only(self):
        self.assertEqual(len(run_all.selected_outputs()),5)
        self.assertFalse(any('ROBUST_UNION' in n for n in run_all.selected_outputs()))
    def test_auxiliary_is_opt_in(self):
        self.assertEqual(len(run_all.selected_outputs(True)),13)
        self.assertTrue(all(n.startswith('auxiliary/robust_union/') for n in run_all.AUXILIARY_DERIVED))
    def test_declared_outputs_exist(self):
        self.assertTrue(all((ROOT/n).is_file() for n in run_all.selected_outputs(True)))
    def test_main_figure_sources(self):
        for name in ('make_expansion_figures.py','plot_juq_provenance_cost_profiles.py'):
            text=(ROOT/name).read_text(encoding='utf-8')
            self.assertIn('37_REAL_CFD_COMMON_PARTITION_PROFILES',text)
            self.assertNotIn('35_REAL_CFD_EDGE_COST_SENSITIVITY',text)
            self.assertNotIn('ROBUST_UNION',text)
if __name__=='__main__':unittest.main(verbosity=2)
