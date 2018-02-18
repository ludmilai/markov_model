from fail_recover_graph import *

def nk_fail_graph(N, R, e, r):
	"""Build failure-recovery graph for N-K storage model given N, redundancy R=N-K, error rate e and recovery rate r"""
	g = StateGraph()
	g.add_state(0, (fault(1, N*e),), (), is_root = True)
	for i in range(1, R+1):
		g.add_state(i, (fault(i+1, (N-i)*e),), (recover(i-1, r),))
	g.add_term_state(R+1)
	return g

if __name__ == '__main__':
	from sympy import symbols
	R = 2
	N, e, r = symbols('N, e, r')
	g = nk_fail_graph(N, R, e, r)
	assert g.fault_resilence() == R
	print g
	print g.mttff_exact()
	print g.mttff_asymptotic()

