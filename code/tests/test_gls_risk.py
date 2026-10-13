import math
import unittest
from gls_risk import solve_axis, tail_bound


class RiskSolverTests(unittest.TestCase):
    def test_single_gaussian(self):
        # Exact known two-sided Gaussian tail at 3 sigma.
        result=solve_axis([(1,0,.02)],math.erfc(3/math.sqrt(2)),0)
        self.assertAlmostEqual(result,.06,places=12)

    def test_fault_mass_reduction(self):
        a=solve_axis([(.99,0,.002),(.01,.1,.05)],1e-6,0)
        b=solve_axis([(.9999,0,.002),(.0001,.1,.05)],1e-6,0)
        self.assertLess(b,a)

    def test_unmodeled_floor(self):
        self.assertIsNone(solve_axis([(.99999,0,.01)],1e-6,1e-5))

    def test_reject_invalid_partition(self):
        with self.assertRaises(ValueError):solve_axis([(.9,0,.01)],1e-6,0)

    def test_separation_and_zero_sigma(self):
        self.assertEqual(tail_bound(.1,.2,.01),1)
        self.assertEqual(tail_bound(.3,.2,0),0)
        self.assertGreater(solve_axis([(1,.1,.02)],1e-6,0),.1)


if __name__=='__main__':unittest.main()
