from sage.all import *

def QuaternionOrderBasis(alpha, O):
    assert alpha in O
    Obasis = O.basis()
    M_O = Matrix(QQ, [b.coefficient_tuple() for b in Obasis]).transpose()
    vec_alpha = vector(alpha.coefficient_tuple())
    coeffs = M_O.inverse() * vec_alpha
    coeffs = [ZZ(c) for c in coeffs]
    #assert alpha == sum([c * beta for c, beta in zip(coeffs, Obasis)])
    return coeffs

def MakePrimitive(alpha, O):
    v_alpha = QuaternionOrderBasis(alpha, O)
    d = gcd(v_alpha)
    assert alpha/d in O
    return alpha/d

def gram_matrix(basis):
    M = []
    for a in basis:
        M.append([QQ(2) * a.pair(b) for b in basis])
    return Matrix(QQ, M)

def ReducedBasis(basis):
    def _matrix_to_gens(M, B):
        return [sum(c * g for c, g in zip(row, B)) for row in M]

    G = gram_matrix(basis)
    U = G.LLL_gram().transpose()
    reduced_basis_elements = _matrix_to_gens(U, basis)
    return reduced_basis_elements