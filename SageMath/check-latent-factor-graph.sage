import numpy as np

# Main function (relies on utils.sage)
def check_graph(G, GammaAdjMatrix, basis_Omega, it=3):
    
    # Boolean adjacency matrix
    adjMatrix = G['adjMatrix']
    nrow = sqrt(len(G['adjMatrix']))
    adjMatrix = np.reshape(np.array(adjMatrix), (nrow,nrow)).astype(bool)

    # Check which variables we need in ring
    variables = ["delta"]
    for i in range(nrow):
        for j in range(nrow):
            if adjMatrix[i, j]:
                variables.append("l"+str(i+1)+str(j+1))

    # Define ring
    order = f"deglex(1),deglex({len(variables)-1})" # block-monomial order
    R = PolynomialRing(QQ, variables, order=order,implementation="singular")

    # Define matrix Lambda
    Lambda = matrix(R, nrow)
    count = 1
    for i in range(nrow):
        for j in range(nrow):
            if adjMatrix[i, j]:
                Lambda[i, j] = R.gens()[count]
                count = count + 1
    
    # Check identifiability
    res_id = [False] * it
    res_fin2one = [False] * it
    for k in range(it):
        # Random Sigma
        sample_space = [1,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97, 
                        101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199]
        sample_space = sample_space + [-x for x in sample_space]
        Id = identity_matrix(R, nrow)
        rLambda = Lambda.apply_map(lambda x : np.random.choice(sample_space) if x!=0 else x)
        rOmega_diag = diagonal_matrix(np.random.choice(sample_space, size=(int(nrow)), replace=True))
        rGamma = matrix(GammaAdjMatrix).apply_map(lambda x : np.random.choice(sample_space) if x!=0 else x)
        rOmega = rGamma * rGamma.T + rOmega_diag
        rSigma = (Id-rLambda).T.inverse() * rOmega * (Id-rLambda).inverse()
        
        # Generate ideal I and get Groebner basis
        generators = get_constraints((Id-Lambda).T * rSigma * (Id-Lambda), basis_Omega)
        generators.append(R.gens()[0] * det((Id-Lambda))-1)
        I = R.ideal(generators)
        
        # Check identifiability
        if (I.dimension() == 0):
            variety_CC = len(I.variety(ring=CC))
            variety_RR = len(I.variety(ring=RR))
            if (variety_CC == 1):
                res_id[k] = True
                res_fin2one[k] = True
            if (variety_CC > 1):
                if (variety_RR > 1):
                    res_fin2one[k] = True
                else: 
                    res_fin2one[k] = True
                    print("Generically finite-to-one over CC but generically identifiable over RR.")
            if (len(I.variety()) == 0):
                print("ERROR, zero points in variety of dimension zero.")
    
    if any(res_fin2one):
        G["finite-to-one"] = True
        nr_id = sum(res_id)
        nr_notid = sum(res_fin2one) - sum(res_id) 
        if (nr_id > nr_notid):
            G["identifiable"] = True
        elif (nr_id < nr_notid):
            G["identifiable"] = False
        else:
            G["identifiable"] = None
            print("ERROR, cannot decide identifiability!")
    else:
        G["identifiable"] = False
        G["finite-to-one"] = False
  
    return G