"""
The Markov states model of fail-recovery.
Given the state graph it calculates the exact solution
for the mean time to first failure (MTTFF)
and asymptotic estimation in symbolic form.
"""

from sympy import Matrix, zeros, factor
from sympy.core.numbers import Zero, One
from collections import namedtuple, defaultdict

# State graph node
State = namedtuple('State', ('id', 'faults', 'recovers'))

# Graph edge corresponding to state transition
Transition = namedtuple('Transition', ('target', 'rate', 'label'))

def fault(target, rate, label=None):
	"""Create fault transition"""
	return Transition(target, rate, label if label is not None else 'fault')

def recover(target, rate, label=None):
	"""Create recovery transition"""
	return Transition(target, rate, label if label is not None else 'recover')

class StateGraph:
	"""Graph representing Markov process states and transitions"""
	def __init__(self):
		"""Create empty graph"""
		self.root = None
		self.nodes = dict()
		self.term_nodes = set()

	def add_state(self, id, faults, recovers, is_root = False):
		"""
		Add node to the state graph. The id may be an arbitrary hashable object.
		The faults and recovers are tuples whose elements are Transition tuples
		corresponding to outgoing transitions.
		"""
		n = State(id, faults, recovers)
		self.nodes[id] = n
		if is_root:
			assert self.root is None
			assert not recovers
			self.root = n
		return n

	def add_term_state(self, id):
		"""Add state node with given id as terminating one"""
		assert id not in self.nodes
		self.term_nodes.add(id)

	def get_states(self):
		"""
		Returns the list of nodes. The first element of the list is the root node.
		Remaining are sorted in order of increasing ids.
		"""
		assert self.root is not None
		nodes = sorted(self.nodes.values())
		return [self.root] + [n for n in nodes if n is not self.root]

	def _tx_fmt(self, tx):
		if tx.target in self.term_nodes:
			target = '!' + str(tx.target)
		else:
			target = '#' + str(tx.target)
		return '%s[%s]->%s' % (tx.label, tx.rate, target)

	def _state_fmt(self, n):
		s = '#' + str(n.id)
		for e in n.faults: 
			s += ' ' + self._tx_fmt(e)
		for e in n.recovers:
			s += ' ' + self._tx_fmt(e)
		return s

	def __str__(self):
		"""Print states annotated with outgoing transitions"""
		return '\n'.join(self._state_fmt(n) for n in self.get_states())

	def transition_intensity_matrix(self):
		"""Returns transition intensity matrix for non-terminating states"""
		nodes = self.get_states()
		nodes_map = {n.id: i for i, n in enumerate(nodes)}
		M = zeros(len(nodes))
		for i, n in enumerate(nodes):
			total = 0
			for x in n.faults:
				total += x.rate
				if x.target not in self.term_nodes:
					M[i, nodes_map[x.target]] = x.rate
			for x in n.recovers:
				total += x.rate
				M[i, nodes_map[x.target]] = x.rate
			M[i,i] = -total
		return M

	def mttff_exact(self):
		"""Returns exact equation for the mean time to first failure"""
		M = self.transition_intensity_matrix()
		D = M.det()
		for i in range(M.rows):
			M[i, 0] = 1
		D1 = M.det()
		return factor(-D1/D)

	def find_node(self, id):
		"""Lookup node by its id"""
		if id not in self.term_nodes:
			return self.nodes[id]
		else:
			return None

	def visit_nodes(self, cb):
		"""Visit graph fault edges in BFS order"""
		visited, discovered, next = set(), [(self.root, 0)], 0
		# Do BFS along fault edges
		while next < len(discovered):
			n, lvl = discovered[next]
			next += 1
			id = n.id if n is not None else None
			if id in visited: # node is already visited
				continue
			visited.add(id)
			cb(n, lvl)
			if n:
				for t in n.faults:
					discovered.append((self.find_node(t.target), lvl + 1))
		# The terminating states are represented by single None key
		assert visited == set(self.nodes.keys() + [None])

	def fault_level_map(self):
		"""
		Returns mapping between state id and state's fault level. The fault level is the length
		of the shortest failure path from the root to the given state. The terminating states
		are represented by single map entry with id equal to None.
		"""
		level_map = dict()
		def visitor(n, lvl):
			level_map[n.id if n is not None else None] = lvl
		self.visit_nodes(visitor)
		return level_map

	def fault_resilence(self):
		"""Returns the maximum number of faults the system can tolerate without going to terminating state"""
		return self.fault_level_map()[None] - 1

	def visit_major_edges(self, cb):
		"""
		Visit major edges of the graph in topological order so all outgoing edges
		will be processed after all incoming edges.
		"""
		# We use level map to find minor transitions
		level_map = self.fault_level_map()
		# The minor transition may be ignored if the corresponding flow is asymptotically negligible.
		# By removing them we make graph acyclic and simplify solution. 
		def is_fault_minor(lvl, t):
			tlvl = level_map[t.target if t.target not in self.term_nodes else None]
			if tlvl <= lvl:
				return True
			assert tlvl == lvl + 1
			return False

		def is_recovery_minor(lvl, t):
			tlvl = level_map[t.target]
			assert tlvl <= lvl
			return tlvl < lvl

		# Build incoming degree map which will be used for topological sorting
		inc_degree = defaultdict(int)
		def degree_counter(n, lvl):
			if n is None:
				return
			for t in n.faults:
				if is_fault_minor(lvl, t):
					continue
				inc_degree[t.target] += 1
			for t in n.recovers:
				if is_recovery_minor(lvl, t):
					continue
				inc_degree[t.target] += 1

		self.visit_nodes(degree_counter)
		# Visit major edges in topological order. Every outgoing edge will be visited
		# after all incoming major edges been processed.
		assert inc_degree[self.root.id] == 0
		queue, next = [self.root], 0
		def visit_transition(n, t):
			target = self.find_node(t.target)
			cb(target, n, t)
			if target is None:
				return # Terminating transition
			inc_degree[t.target] -= 1
			if not inc_degree[t.target]:
				queue.append(target)

		while next < len(queue):
			n, next = queue[next], next + 1
			lvl = level_map[n.id]
			for t in n.faults:
				if not is_fault_minor(lvl, t):
					visit_transition(n, t)
			for t in n.recovers:
				if not is_recovery_minor(lvl, t):
					visit_transition(n, t)

		# Check we have processed all nodes
		assert next == len(self.nodes.keys())

	def mttff_asymptotic(self):
		"""
		Returns asymptotic equation for the mean time to first failure
		in the approximation of fault rates << recovery rates.
		"""
		class visitor:
			def __init__(this, root):
				this.term_rate = Zero()
				this.population_map = defaultdict(Zero)
				this.population_map[root.id] = One()

			def visit(this, n, parent, tr):
				if n is None: # transition to terminating state
					this.term_rate += this.population_map[parent.id] * tr.rate
				else: # transition to internal state
					recovery_rate = Zero()
					for r in n.recovers:
						recovery_rate += r.rate
					this.population_map[n.id] += this.population_map[parent.id] * tr.rate / recovery_rate

		v = visitor(self.root)
		self.visit_major_edges(v.visit)
		return factor(1/v.term_rate)

if __name__ == '__main__':
	from sympy import symbols
	e, r = symbols('e, r')
	# Simple mirrored disk model
	g = StateGraph()
	g.add_state(0, (fault(1, 2*e),), (), is_root = True)
	g.add_state(1, (fault(2, e),), (recover(0, r),))
	g.add_term_state(2)
	assert g.fault_resilence() == 1
	print g
	print g.transition_intensity_matrix()
	print g.mttff_exact()
	print g.mttff_asymptotic()
	print g.fault_level_map()

