import itertools,json,unittest
from pathlib import Path
import numpy as np
from weight_stability import envelope_bounds

class WeightStabilityTests(unittest.TestCase):
    def setUp(self):
        self.d=np.array([[1,0,1,0],[0,0,1,1],[0,1,1,0]])
        self.lo=[.15]*4;self.hi=[.35]*4
    def test_closed_form_matches_enumeration(self):
        rng=np.random.default_rng(20260906)
        for w in rng.dirichlet(np.ones(4),size=500):
            self.assertAlmostEqual(min(self.d@w),w[2]+min(w[[0,1,3]]))
    def test_bounds_and_commitment_gap(self):
        r=envelope_bounds(self.d,self.lo,self.hi)
        self.assertAlmostEqual(r['minimum'],3/10)
        self.assertAlmostEqual(r['maximum_reoptimized'],17/30)
        self.assertAlmostEqual(r['minimum_fixed_witness_worst_cost'],7/10)
        self.assertAlmostEqual(r['commitment_gap'],2/15)
    def test_all_witnesses_can_be_optimal_none_uniformly(self):
        r=envelope_bounds(self.d,self.lo,self.hi)
        for a in r['witness_regions']:
            self.assertTrue(a['optimal_somewhere']);self.assertAlmostEqual(a['worst_regret'],.2)
    def test_concavity_monotonicity_homogeneity(self):
        rng=np.random.default_rng(63)
        for _ in range(300):
            a=rng.random(4);b=rng.random(4);t=rng.random()
            value=lambda w:float(min(self.d@w))
            self.assertGreaterEqual(value(t*a+(1-t)*b)+1e-12,t*value(a)+(1-t)*value(b))
            self.assertGreaterEqual(value(a+b)+1e-12,value(a))
            self.assertAlmostEqual(value(3*a),3*value(a))
            self.assertLessEqual(abs(value(a)-value(b)),sum(abs(a-b))+1e-12)
    def test_exact_witness_region(self):
        rng=np.random.default_rng(25)
        for w in rng.dirichlet(np.ones(4),size=300):
            for d in self.d:
                region=np.all((d-self.d)@w<=1e-12)
                self.assertEqual(region, bool(d@w<=min(self.d@w)+1e-12))
    def test_domain_rejects_zero_and_infeasible_bounds(self):
        for lo,hi in [([0]*4,[1]*4),([.3]*4,[.4]*4),([.1]*4,[.2]*4)]:
            with self.assertRaises(ValueError):envelope_bounds(self.d,lo,hi)
    def test_single_witness_has_no_commitment_gap(self):
        r=envelope_bounds([[1,0,0,1]],self.lo,self.hi)
        self.assertAlmostEqual(r['commitment_gap'],0)
    def test_cfd_complete_enumeration_structure(self):
        x=json.loads((Path(__file__).parent/'weight_stability_results.json').read_text())
        self.assertEqual(x['audit_count'],8)
        for a in x['results']:
            self.assertEqual(len(a['partitions']),12)
            d=np.asarray(a['pareto_cuts']);self.assertEqual(d.shape,(3,4))
            anchor=2 if a['support']=='fine_F100' else 3
            self.assertTrue(np.all(d[:,anchor]==1));self.assertTrue(np.all(d.sum(axis=1)==2))
            self.assertEqual(set(np.where(d.sum(axis=0)==1)[0]),set(range(4))-{anchor})
            self.assertAlmostEqual(a['bounds']['maximum_reoptimized'],17/30)

if __name__=='__main__':unittest.main()
