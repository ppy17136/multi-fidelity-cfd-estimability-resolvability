import itertools
import unittest
from fractions import Fraction
from types import SimpleNamespace
from unittest.mock import patch
import numpy as np
from graph_constrained_coupling import two_label_graph_relaxation, exact_decision_relaxation_frontier
from graph_coupling_milp import milp_decision_destroying_relaxation


class NumericalCertificateGuards(unittest.TestCase):
    def test_large_weight_forcing(self):
        for w in (2.**53, 2.**54, 1e100, 1e308):
            c=two_label_graph_relaxation([[(0,0)],[(1,1)]],[(0,1,w)])
            self.assertEqual(c.labels,((0,0),(1,1)))
            self.assertEqual(c.cost,w)

    def test_sum_overflow_but_finite_optimum(self):
        c=two_label_graph_relaxation([[(0,0)],[(0,0),(1,1)],[(1,1)]],[(0,1,1e308),(1,2,1e308)])
        self.assertEqual(c.cost,1e308)
        self.assertEqual(c.labels[0],(0,0));self.assertEqual(c.labels[2],(1,1))

    def test_unrepresentable_cost_raises(self):
        with self.assertRaises(OverflowError):
            two_label_graph_relaxation([[(0,0)],[(1,1)],[(0,0)]],[(0,1,1e308),(1,2,1e308)])

    def test_subnormal_and_large_capacity_combinations(self):
        tiny=np.nextafter(0.,1.)
        c=two_label_graph_relaxation([[(0,0)],[(0,0),(1,1)],[(1,1)]],[(0,1,tiny),(1,2,1e100)])
        self.assertEqual(c.cost,tiny)
        self.assertEqual(c.labels[1],(1,1))

    def test_two_label_against_exact_rational_enumeration(self):
        rng=np.random.default_rng(20260907)
        a=(0,0);b=(1,1);allowed=[(a,),(a,b),(a,b),(b,)]
        for _ in range(30):
            es=[(i,j,float(rng.choice([.1,.25,2.**53,1e-100]))) for i,j in itertools.combinations(range(4),2)]
            costs=[(sum((Fraction.from_float(w) for i,j,w in es if ls[i]!=ls[j]),Fraction()),ls) for ls in itertools.product(*allowed)]
            optimum=min(c for c,_ in costs)
            cert=two_label_graph_relaxation(allowed,es)
            actual=sum((Fraction.from_float(w) for i,j,w in es if cert.labels[i]!=cert.labels[j]),Fraction())
            self.assertEqual(actual,optimum)

    def test_near_boundary_infeasible_above(self):
        z=np.array([[1e-8,1.],[1.,1e-8]])
        with patch('graph_coupling_milp.milp',side_effect=AssertionError('must not call solver')):
            state,c=milp_decision_destroying_relaxation(z,[(0,1,1.)],1e-8)
        self.assertEqual(state,'resolved_above');self.assertEqual(c.status,'infeasible_exact_endpoint')

    def test_near_boundary_infeasible_below(self):
        z=-np.array([[1e-8,1.],[1.,1e-8]])
        state,c=milp_decision_destroying_relaxation(z,[(0,1,1.)],-1e-8)
        self.assertEqual(state,'resolved_below');self.assertEqual(c.status,'infeasible_exact_endpoint')

    def test_exact_crossing_boundary(self):
        state,c=milp_decision_destroying_relaxation(np.array([[1.,3.],[3.,1.]]),[(0,1,1.)],2.)
        self.assertEqual(c.status,'optimal');self.assertEqual(c.attained_value,2.)

    def test_initially_indeterminate_zero_cost(self):
        z=np.array([[0.,1.],[0.,1.]])
        state,c=milp_decision_destroying_relaxation(z,[(0,1,1.)],1.)
        direct=exact_decision_relaxation_frontier(z,[(0,1,1.)],1.)
        self.assertEqual(c.status,'initially_indeterminate');self.assertEqual(c.objective,0.)
        self.assertIsNone(c.labels);self.assertEqual(direct.minimum_relaxation_cost,0.)

    def test_cancellation_baseline_uses_exact_stored_inputs(self):
        # Rounded sequential sum can lose the unit contribution.
        z=np.array([[1e16,1e16],[1.,1.],[-1e16,-1e16]])
        state,c=milp_decision_destroying_relaxation(z,[(0,1,1.),(1,2,1.)],0.)
        self.assertEqual(state,'resolved_above');self.assertEqual(c.status,'infeasible_exact_endpoint')

    def test_nonfinite_threshold_rejected(self):
        for t in (np.nan,np.inf,-np.inf):
            with self.assertRaises(ValueError):
                milp_decision_destroying_relaxation(np.ones((2,2)),[(0,1,1.)],t)

    def test_solver_success_with_noncrossing_witness_rejected(self):
        z=np.array([[0.,1.],[1.,0.]])
        fake=SimpleNamespace(x=np.array([1.,0.,1.,0.,0.]),fun=0.,success=True,mip_dual_bound=0.,mip_gap=0.,mip_node_count=0)
        with patch('graph_coupling_milp.milp',return_value=fake):
            with self.assertRaisesRegex(RuntimeError,'crossing inequality'):
                milp_decision_destroying_relaxation(z,[(0,1,1.)],.5)

    def test_solver_infeasibility_cannot_override_endpoint(self):
        fake=SimpleNamespace(x=None,message='infeasible',mip_dual_bound=None,mip_gap=None,mip_node_count=0)
        with patch('graph_coupling_milp.milp',return_value=fake):
            with self.assertRaisesRegex(RuntimeError,'contradicts'):
                milp_decision_destroying_relaxation(np.array([[0.,1.],[1.,0.]]),[(0,1,1.)],.5)

if __name__=='__main__': unittest.main()
