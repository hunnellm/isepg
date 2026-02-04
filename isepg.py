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
