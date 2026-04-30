def czerosgame(g,F=[],B=[]):
    S=set(F) #suspicuous vertices
    Black_vertices=set(F) # current black vertices
    again=1 # iterate again or not
    V=set(g.vertices())
    while again==1:
        #print("S=",S)
        again=0
        for y in V.difference(Black_vertices):
            N=set(g.neighbors(y))
            D=N.difference(Black_vertices) # set of white neighbors
            #print(x,len(D))
            if len(D)==0:
                Black_vertices.add(y)
                again=1
                #print("vertex ",y," forced itself")
                break
        for x in S:
            N=set(g.neighbors(x))
            
            D=N.difference(Black_vertices) # set of white neighbors
            #print(x,len(D))
            if len(D)==0:
                S.remove(x)
                Black_vertices.add(x)
                again=1
                #print("vertex ",x," forced itself")
                break
            if len(D)==1:
                for v in D:
                    y=v # the only white neighbor
                    if (((x,y) in B)==False) and (((y,x) in B)==False):
                        again=1
                        S.remove(x)
                        S.add(y)
                        Black_vertices.add(y)
                        #print(x," forced ",y)
                        break
    #print(Black_vertices)                    
    return(Black_vertices)


def gZ_leq(graph, support=[], bannedset=[],i=None):
	"""
	For a given graph with support and banned set, if there is a zero forcing set of size i then return it; otherwise return False.

	Input:
		graph: a simple graph
		support: a list of vertices of g
		bannedset: a list of tuples representing banned edges of graph
		i: an integer, the function check gZ <= i or not

	Output:
		if F is a zero forcing set of size i and support is a subset of F, then return F
		False otherwise

	Examples:
		sage: gZ_leq(graphs.PathGraph(5),[],[],1)
		set([0])
		sage: gZ_leq(graphs.PathGraph(5),[],[(0,1)],1) 
		False
	"""
	if i < len(support):
#		print 'i cannot less than the cardinality of support'
		return False
	j=i-len(support) # additional number of black vertices
	VX=graph.vertices()
	order=graph.order()
	for y in support:
		VX.remove(y)
	# VX is the vertices outside support now
	for subset in Subsets(VX,j):
		test_set=set(support).union(subset) # the set is tested to be a zero forcing set
		outcome=czerosgame(graph, test_set, bannedset)
		if len(outcome)==order:
			return test_set
	return False



def is_coupled_zero_forcing_set(B,G,matching):
    n=G.order()
    loops =[]
    H=G.copy()
    for j in matching:
        H.add_edge(j)
        
    #H.show(layout='circular')
    if len(czerosgame(H,B))==n:
        return True
    return False       

def coupled_Z(G,matching,all_sets=False):
    n=G.order()
    V=G.vertices()
    czf_sets=[]
    complete=False
    for k in range(n):
        subsets=Subsets(V,k)
        for s in subsets:
            #print(s)
            if is_coupled_zero_forcing_set(s,G,matching):
                czf_sets.append(s)
                complete=True
        if complete==True and all_sets==True:
            return czf_sets
        elif complete==True:
            return k
    return False       

def _build_partner_map(G, matching):
    """
    Build and validate a partner map c from a prescribed matching.

    The matching must form a perfect involution on G.vertices():
    each vertex appears exactly once, no self-pairs are allowed,
    and every vertex of G must be covered.

    Input:
        G       : a simple graph
        matching: a list of 2-tuples [(u1, v1), (u2, v2), ...] covering
                  every vertex of G exactly once

    Output:
        A dict c such that c[u] = v and c[v] = u for every (u, v) in matching.

    Raises:
        ValueError if matching is not a perfect involution on G.vertices().

    Examples:
        sage: G = graphs.PathGraph(4)
        sage: _build_partner_map(G, [(0, 1), (2, 3)])
        {0: 1, 1: 0, 2: 3, 3: 2}
        sage: _build_partner_map(G, [(0, 1)])
        Traceback (most recent call last):
            ...
        ValueError: matching has 1 pair(s) but graph has 4 vertices; matching must cover every vertex exactly once
    """
    vertices = set(G.vertices())
    n = G.order()

    if 2 * len(matching) != n:
        raise ValueError(
            "matching has {} pair(s) but graph has {} vertices; "
            "matching must cover every vertex exactly once".format(len(matching), n)
        )

    c = {}
    for pair in matching:
        if len(pair) != 2:
            raise ValueError(
                "each matching element must be a 2-tuple; got {!r}".format(pair)
            )
        u, v = pair
        if u == v:
            raise ValueError(
                "self-pair ({}, {}) is not allowed in matching".format(u, v)
            )
        if u not in vertices:
            raise ValueError(
                "vertex {!r} in matching is not a vertex of G".format(u)
            )
        if v not in vertices:
            raise ValueError(
                "vertex {!r} in matching is not a vertex of G".format(v)
            )
        if u in c:
            raise ValueError(
                "vertex {!r} appears more than once in matching".format(u)
            )
        if v in c:
            raise ValueError(
                "vertex {!r} appears more than once in matching".format(v)
            )
        c[u] = v
        c[v] = u

    for v in vertices:
        if v not in c:
            raise ValueError(
                "vertex {!r} is not covered by the matching".format(v)
            )

    return c


def qzerosgame(G, initial_blue, c):
    """
    Run the quantum zero forcing game on graph G with the given initial blue
    set and partner map c, returning the final set of blue vertices.

    The three forcing rules are applied iteratively until no new vertex
    can be forced:

    Rule 1 (standard forcing): A blue vertex v forces the unique white
        vertex w to blue if w is the only white vertex in N(v) ∪ {c(v)}.

    Rule 2 (white-vertex forcing): A white vertex w forces itself to blue
        if every vertex in N(w) ∪ {c(w)} is already blue.

    Rule 3 (double forcing): A blue vertex v forces two white vertices
        w and c(w) simultaneously if:
          (a) w and c(w) are the only white vertices in N(v) ∪ {c(v)}, AND
          (b) either w has no white neighbor, or c(w) is the unique white
              neighbor of w  (checked for both assignments of the two whites
              as "w" and "c(w)").

    Input:
        G           : a simple graph
        initial_blue: an iterable of initially blue vertices
        c           : a dict partner map with c[v] = mate of v

    Output:
        A set of all blue vertices after the game terminates.

    Examples:
        sage: G = graphs.PathGraph(4)
        sage: c = {0: 1, 1: 0, 2: 3, 3: 2}
        sage: qzerosgame(G, [0, 2], c) == {0, 1, 2, 3}
        True
    """
    Blue = set(initial_blue)
    V = set(G.vertices())

    changed = True
    while changed:
        changed = False

        # Rule 1: standard forcing
        for v in Blue:
            Nv = set(G.neighbors(v)) | {c[v]}
            white_in_Nv = Nv - Blue
            if len(white_in_Nv) == 1:
                Blue.add(next(iter(white_in_Nv)))
                changed = True
                break
        if changed:
            continue

        # Rule 2: white-vertex forcing
        for w in V - Blue:
            Nw = set(G.neighbors(w)) | {c[w]}
            if Nw.issubset(Blue):
                Blue.add(w)
                changed = True
                break
        if changed:
            continue

        # Rule 3: double forcing
        for v in Blue:
            Nv = set(G.neighbors(v)) | {c[v]}
            white_in_Nv = Nv - Blue
            if len(white_in_Nv) == 2:
                p, q = list(white_in_Nv)
                if c[p] == q:
                    # Try p as "w" and q as "c(w)"
                    white_nbrs_p = set(G.neighbors(p)) - Blue
                    cond_p = (len(white_nbrs_p) == 0 or
                              (len(white_nbrs_p) == 1 and q in white_nbrs_p))
                    # Try q as "w" and p as "c(w)"
                    white_nbrs_q = set(G.neighbors(q)) - Blue
                    cond_q = (len(white_nbrs_q) == 0 or
                              (len(white_nbrs_q) == 1 and p in white_nbrs_q))
                    if cond_p or cond_q:
                        Blue.add(p)
                        Blue.add(q)
                        changed = True
                        break
        # if changed is still False the loop terminates

    return Blue


def is_quantum_zero_forcing_set(B, G, matching):
    """
    Return True if B is a quantum zero forcing set of G with respect to the
    prescribed perfect matching, and False otherwise.

    Input:
        B       : an iterable of initially blue vertices (a subset of G.vertices())
        G       : a simple graph
        matching: a list of 2-tuples defining a perfect matching on G.vertices()

    Output:
        True if B is a quantum zero forcing set; False otherwise.

    Raises:
        ValueError if matching is not a perfect involution on G.vertices().

    Examples:
        sage: G = graphs.PathGraph(4)
        sage: is_quantum_zero_forcing_set([0, 2], G, [(0, 1), (2, 3)])
        True
        sage: is_quantum_zero_forcing_set([0], G, [(0, 1), (2, 3)])
        False
    """
    c = _build_partner_map(G, matching)
    return len(qzerosgame(G, B, c)) == G.order()


def quantum_Z(G, matching, all_sets=False):
    """
    Compute the quantum zero forcing number of G with respect to the
    prescribed perfect matching.

    Input:
        G        : a simple graph
        matching : a list of 2-tuples defining a perfect matching on G.vertices()
        all_sets : if False (default), return the minimum size k of a quantum
                   zero forcing set; if True, return a list of all quantum zero
                   forcing sets of that minimum size.

    Output:
        The minimum size of a quantum zero forcing set (all_sets=False), or
        a list of all minimum-size quantum zero forcing sets (all_sets=True).
        Returns False if no zero forcing set exists (should not happen for a
        non-empty graph).

    Raises:
        ValueError if matching is not a perfect involution on G.vertices().

    Examples:
        sage: G = graphs.PathGraph(4)
        sage: quantum_Z(G, [(0, 1), (2, 3)])
        2
        sage: quantum_Z(G, [(0, 1), (2, 3)], all_sets=True)
        [{0, 2}, {0, 3}, {1, 2}, {1, 3}]
    """
    c = _build_partner_map(G, matching)
    n = G.order()
    V = G.vertices()
    for k in range(n + 1):
        if all_sets:
            qzf_sets = [s for s in Subsets(V, k)
                        if len(qzerosgame(G, s, c)) == n]
            if qzf_sets:
                return qzf_sets
        else:
            for s in Subsets(V, k):
                if len(qzerosgame(G, s, c)) == n:
                    return k
    return False


def Jmatrix(n):
    if n%2==1:
        return False
    else:
        p=n/2
        I=identity_matrix(ZZ,p)
        O= zero_matrix(ZZ,p,p)
        J=block_matrix([[O,I],[-I,O]])
        return J*identity_matrix(ZZ,n)

def symp_prod(x,y):
    n=len(x)
    if n!=len(y) or n%2==1:
        return False
    else:
        return x.dot_product(Jmatrix(n)*y)


def symp_evalues(A):
    if A.nrows()%2 ==1 or A.nrows()!=A.ncols():
        return False
    else:
        R=Jmatrix(A.nrows())*A
        eig = R.eigenvalues()
        seig=[]
        for i in range(len(eig)):
            seig.append(eig[i].imag().abs())
        return sorted(list(set(seig)))
