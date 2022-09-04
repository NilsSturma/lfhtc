import math
import numpy as np
from itertools import combinations


def binom(n, k):
    return math.factorial(n) // math.factorial(k) // math.factorial(n - k)


def get_tetrads(A):  
    tetrads = []
    for c in combinations(range(A.nrows()), 4):
        t1 = A[c[0],c[1]] * A[c[2],c[3]] - A[c[0],c[2]] * A[c[1],c[3]]
        tetrads.append(t1)
        t2 = A[c[0],c[3]] * A[c[1],c[2]] - A[c[0],c[2]] * A[c[1],c[3]]
        tetrads.append(t2)
    return tetrads
  
  
def get_upper_tri(A, k=0):
    m = A.nrows()
    utri_ind = np.triu_indices(int(m), k=k)
    upper_tri_list = []
    for i in range(len(utri_ind[0])):
        row_ind = utri_ind[0][i]
        col_ind = utri_ind[1][i]
        upper_tri_list.append(A[row_ind, col_ind])
    return(upper_tri_list)
    
   
def get_vanishing_ideal_Omega(GammaAdjMatrix):
    
    # GammaAdjMatrix is mxk matrix
    m = GammaAdjMatrix.shape[0]
    k = GammaAdjMatrix.shape[1]    
    GammaAdjMatrix = GammaAdjMatrix.astype(bool)

    # Check which variables we need in ring
    var_gamma = []
    for i in range(m):
        for j in range(k):
            if GammaAdjMatrix[i, j]:
                var_gamma.append("g"+str(i+1)+str(j+1))
    var_diag = []
    for i in range(m):
        var_diag.append("d"+str(i+1))
    var_omega = []
    for i in range(m):
        for j in range(i, m):
            var_omega.append("w"+str(i+1)+str(j+1))
    variables = var_gamma + var_diag + var_omega

    # Define ring
    order = f"deglex"
    #f"deglex({len(var_gamma)}),deglex({len(var_diag)}),deglex({len(var_omega)}))" # block-monomial order
    R = PolynomialRing(QQ, variables, order=order,implementation="singular")

    # Define matrices Gamma, Omega_diag and Omega
    count = 0
    Gamma = matrix(R,m,k)
    for i in range(m):
        for j in range(k):
            if GammaAdjMatrix[i, j]:
                Gamma[i, j] = R.gens()[count]
                count = count + 1
    diag = diagonal_matrix(R, variables[count:(count+m)])
    count = count+m
    Omega = matrix(R,m)
    for i in range(m):
        for j in range(i,m):
            Omega[i,j] = variables[count]
            Omega[j,i] = variables[count]
            count = count + 1

    # Compute vanishing ideal by implicitization
    J = R.ideal(get_upper_tri(Omega - (Gamma * Gamma.T + diag)))
    to_eliminate = R.gens()[:(len(var_gamma)+len(var_diag))]
    I = J.elimination_ideal(to_eliminate, algorithm="libsingular:eliminate")
    basis = I.groebner_basis(algorithm="libsingular:groebner")

    # Interpret the ideal in smaller ring (only omega as variables)
    littleR= PolynomialRing(QQ, var_omega, order="deglex")
    littleI = littleR.ideal(basis)
    littleBasis = littleI.groebner_basis(algorithm="libsingular:groebner")
    
    return littleI, littleBasis
    
    
def get_constraints(A, basis):
    upper_tri = get_upper_tri(A, k=0)
    constraints = []
    for i in range(len(basis)):
        constraints.append(basis[i](upper_tri))
    return(constraints)