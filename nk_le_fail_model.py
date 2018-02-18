from fail_recover_graph import *

def nk_le_fail_graph(N, R, e, r, el, rl):
	"""
	Build failure-recovery graph for N-K storage model with disk and sector (latent) errors given N, redundancy R=N-K,
	error rate e, recovery rate r, sector error rate es and sector error recovery rate rs.
	"""
	g = StateGraph()
	# Label states by (err cnt, sector err cnt) tuple
	for d in range(0, R+2):
		for s in range(0, R+2):
			disk_fault = fault((d+1, s), (N-d-s)*e)
			sect_fault = fault((d, s+1), (N-d-s)*es)
			bad_sect_disk_fault = fault((d+1, s-1), s*e)
			disk_recover = recover((d-1, s), r)
			sect_recover = recover((d, s-1), rs)
			if d + s > R:
				# Terminating state
				g.add_term_state((d, s))
			elif not d and not s:
				# Root, no error, state
				g.add_state(
					(0, 0),
					(disk_fault, sect_fault),
					(), is_root = True
				)
			elif not s:
				# No sector errors
				g.add_state(
					(d, s),
					(disk_fault, sect_fault),
					(disk_recover,)
				)
			elif not d:
				# No disk errors
				g.add_state(
					(d, s),
					(disk_fault, sect_fault, bad_sect_disk_fault),
					(sect_recover,)
				)
			else:
				g.add_state(
					(d, s),
					(disk_fault, sect_fault, bad_sect_disk_fault),
					(disk_recover, sect_recover)
				)
	return g

if __name__ == '__main__':
	from sympy import symbols
	R = 2
	N, e, r, es, rs = symbols('N, e, r, es, rs')
	g = nk_le_fail_graph(N, R, e, r, es, rs)
	assert g.fault_resilence() == R
	print g
	# print g.mttff_exact()
	print g.mttff_asymptotic()

