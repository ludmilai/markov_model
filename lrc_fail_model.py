from fail_recover_graph import *

def lrc2_fail_graph(M, R, e, rl, rg):
	"""
	Build failure-recovery graph for Azure LRC storage model with 2 local groups.
	The M is the number of information blocks in each local group. The R is the number
	of global parity blocks. There are one additional local parity block in each local
	group so there are R + 2 parity blocks in total in addition to the 2*M information
	blocks.	Other model parameters are error rate e, local recovery rate rl and global
	recovery rate rg.
	"""
	graph = StateGraph()
	# Label each state by (e1, e2, l1, l2, gp) tuple where:
	#  e1, e2 (0..M) - the number of failed information blocks in each group
	#  l1, l2 (0..1) - the number of failed local parity blocks in each group
	#  gp     (0..R) - the number of failed global parity blocks
	# We are going to generate all possible states in such a way that the
	# boundary of the generated id space will be covered by terminating states
	for e1 in range(R+3):
		for e2 in range(R+3):
			for l1 in range(2):
				for l2 in range(2):
					for gp in range(R+1):
						# The block id
						id = (e1, e2, l1, l2, gp)
						# Find the number of parity blocks useful for data reconstruction
						parities = R - gp # global ones
						if e1 and not l1: parities += 1
						if e2 and not l2: parities += 1
						# Check if errors are theoretically recoverable
						if e1 + e2 > parities:
							# Terminating state
							graph.add_term_state(id)
							continue
						# Information bocks failures
						e1fault = fault((e1 + 1, e2, l1, l2, gp), (M-e1)*e)
						e2fault = fault((e1, e2 + 1, l1, l2, gp), (M-e2)*e)						
						# Local parity blocks failures
						l1fault = fault((e1, e2, 1, l2, gp), e)
						l2fault = fault((e1, e2, l1, 1, gp), e)
						# Global parity blocks failures
						gpfault = fault((e1, e2, l1, l2, gp + 1), (R-gp)*e)
						# Local recovery
						l1recover = recover((0, e2, 0, l2, gp), rl)
						l2recover = recover((e1, 0, l1, 0, gp), rl)
						# Global recovery (right to the healthy state for simplicity)
						g_recover = recover((0, 0, 0, 0, 0), rg)
						# Build transition lists
						faults, recovers = [e1fault, e2fault], []
						if not l1: faults.append(l1fault)
						if not l2: faults.append(l2fault)
						if gp < R: faults.append(gpfault)
						# Consider local recovery first
						if e1 + l1 == 1: recovers.append(l1recover)
						if e2 + l2 == 1: recovers.append(l2recover)
						# If local recovery is impossible but we still have failed blocks
						# consider global recovery
						if not recovers and e1 + e2 + gp > 0:
							recovers.append(g_recover)
						# Add state to the graph
						graph.add_state(id, faults, recovers, is_root = id == (0, 0, 0, 0, 0))
	return graph

if __name__ == '__main__':
	from sympy import symbols
	R = 2
	N, e, r = symbols('N, e, r')
	M = (N - R)/2 - 1
	g = lrc2_fail_graph(M, R, e, r, r)
	assert g.fault_resilence() == R + 1
	print g.mttff_asymptotic()

