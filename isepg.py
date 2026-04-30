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


def qzerosgame(G, initial_blue, matching):
    """
    Plays the quantum zero forcing game on graph G starting from initial_blue,
    using the partner map induced by matching.  Applies the three quantum zero
    forcing rules until no further vertices can be coloured blue, then returns
    the final blue set.

    Rules applied (iterated to a fixed point):
      1. Standard forcing: blue vertex v forces white vertex w to blue if w is
         the only white vertex in N(v) union {c(v)}.
      2. White vertex forcing: white vertex w forces itself to blue if every
         vertex in N(w) union {c(w)} is blue.
      3. Double forcing: blue vertex v forces white vertices w and c(w) to blue
         if w and c(w) are the only white vertices in N(v) union {c(v)}, and
         either w has no white neighbor in G, or c(w) is the unique white
         neighbor of w in G.

    Input:
        G: a simple graph
        initial_blue: a collection of initially blue vertices
        matching: a list of edges (2-tuples) giving the partner map c;
                  each edge (u, v) means c(u)=v and c(v)=u

    Output:
        the set of blue vertices after no more forces can be applied

    Examples::

        sage: G = graphs.CompleteGraph(4)
        sage: matching = [(0,1),(2,3)]
        sage: qzerosgame(G, [0,1], matching) == set([0,1,2,3])
        True
        sage: qzerosgame(G, [0], matching) == set([0])
        True
    """
    c = {}
    for edge in matching:
        u, v = edge[0], edge[1]
        c[u] = v
        c[v] = u

    Blue = set(initial_blue)
    V = set(G.vertices())

    again = 1
    while again == 1:
        again = 0
        White = V.difference(Blue)

        # Rule 2: white vertex forcing
        for w in White:
            Nw = set(G.neighbors(w))
            Nw_ext = Nw | {c[w]} if w in c else Nw
            if Nw_ext.issubset(Blue):
                Blue.add(w)
                again = 1
                break

        if again == 1:
            continue

        # Rules 1 and 3: blue vertex forcing
        for v in Blue:
            Nv = set(G.neighbors(v))
            Nv_ext = Nv | {c[v]} if v in c else Nv
            white_in_Nv = Nv_ext.intersection(White)

            # Rule 1: standard forcing
            if len(white_in_Nv) == 1:
                w = next(iter(white_in_Nv))
                Blue.add(w)
                again = 1
                break

            # Rule 3: double forcing
            if len(white_in_Nv) == 2:
                w1, w2 = list(white_in_Nv)
                if not (w1 in c and c[w1] == w2):
                    continue
                # check the condition for either labeling of (w, c(w))
                wnb1 = set(G.neighbors(w1)).intersection(White)
                wnb2 = set(G.neighbors(w2)).intersection(White)
                cond = (
                    (len(wnb1) == 0 or (len(wnb1) == 1 and w2 in wnb1))
                    or
                    (len(wnb2) == 0 or (len(wnb2) == 1 and w1 in wnb2))
                )
                if cond:
                    Blue.add(w1)
                    Blue.add(w2)
                    again = 1
                    break

    return Blue


def is_quantum_zero_forcing_set(B, G, matching):
    """
    Determines whether B is a quantum zero forcing set for graph G with the
    given matching.

    Input:
        B: a set or list of vertices of G (the initially blue vertices)
        G: a simple graph
        matching: a list of edges (2-tuples) representing a perfect matching on G

    Output:
        True if B is a quantum zero forcing set, False otherwise

    Examples::

        sage: G = graphs.CompleteGraph(4)
        sage: matching = [(0,1),(2,3)]
        sage: is_quantum_zero_forcing_set([0,1], G, matching)
        True
        sage: is_quantum_zero_forcing_set([0], G, matching)
        False
    """
    n = G.order()
    return len(qzerosgame(G, B, matching)) == n


def quantum_Z(G, matching, all_sets=False):
    """
    Computes the quantum zero forcing number of G with respect to the given
    matching.

    Input:
        G: a simple graph
        matching: a list of edges (2-tuples) representing a perfect matching on G
        all_sets: if False (default), return the minimum size k of a quantum
                  zero forcing set; if True, return all minimum quantum zero
                  forcing sets as a list

    Output:
        If all_sets=False: an integer k (the quantum zero forcing number)
        If all_sets=True: a list of all minimum quantum zero forcing sets
        False if no quantum zero forcing set is found

    Examples::

        sage: G = graphs.CompleteGraph(4)
        sage: matching = [(0,1),(2,3)]
        sage: quantum_Z(G, matching)
        2
        sage: quantum_Z(G, matching, all_sets=True)
        [{0, 1}, {2, 3}]
    """
    n = G.order()
    V = G.vertices()
    qzf_sets = []
    complete = False
    for k in range(n):
        subsets = Subsets(V, k)
        for s in subsets:
            if is_quantum_zero_forcing_set(s, G, matching):
                qzf_sets.append(s)
                complete = True
        if complete == True and all_sets == True:
            return qzf_sets
        elif complete == True:
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
