from sage.all import *

from theta_structures.couple_point import CouplePoint
from theta_isogenies.product_isogeny_sqrt import EllipticProductIsogenySqrt

    
class IsogenyProduct:
    def __init__(self, P1, P2, Q1, Q2, deg, iota=None):

        assert all([2**(deg-1)*R for R in [P1, P2, Q1, Q2]])
        assert all([not 2**(deg)*R for R in [P1, P2, Q1, Q2]])

        assert Q1.curve() == P1.curve()
        assert Q2.curve() == P2.curve()

        self.domain = (P1.curve(), P2.curve())
        self.endomorphisms = []
        self.diagonal_isogenies = []
        self.iota = iota
        self.deg = 2**deg
        index = 0
        if self.iota:
            P1_, P2_, Q1_, Q2_ = [2**(deg-1)*R for R in [P1, P2, Q1, Q2]]
            type_end = self.CheckEndomorphism(P1_, P2_, Q1_, Q2_)

            if type_end:
                self.endomorphisms.append(type_end)

            while type_end:
                index += 1
                P1, P2 = self.EvalEndomorphismType(type_end, P1, P2)
                Q1, Q2 = self.EvalEndomorphismType(type_end, Q1, Q2)
                P1_, P2_, Q1_, Q2_ = [2**(deg-1-index)*R for R in [P1, P2, Q1, Q2]]
                if not all([P1_, P2_, Q1_, Q2_]):
                    break
                type_end = self.CheckEndomorphism(P1_, P2_, Q1_, Q2_)
                # We need to remember the endomorphisms for images
                if type_end:
                    self.endomorphisms.append(type_end)
            
            flag, R1, R2 = self.CheckDiagonalIsogeny(P1_, P2_, Q1_, Q2_)
            while flag:
                index += 1
                phi_1, phi_2 = self.DiagonalIsogeny(R1, R2)
                F1 = phi_1.codomain()
                F2 = phi_2.codomain()
                P1, P2 = phi_1(P1), phi_2(P2)
                Q1, Q2 = phi_1(Q1), phi_2(Q2)
                temp = (phi_1, phi_2)
                self.diagonal_isogenies.append(temp)
                P1_, P2_, Q1_, Q2_ = [2**(deg-1-index)*R for R in [P1, P2, Q1, Q2]]
                flag, R1, R2 = self.CheckDiagonalIsogeny(P1_, P2_, Q1_, Q2_)

        K1, K2 = CouplePoint(P1, P2), CouplePoint(Q1, Q2)
        self.Phi = EllipticProductIsogenySqrt((K1, K2), deg-index)


    def CheckEndomorphism(self, P1, P2, Q1, Q2):
        if (P1 == P2) and (Q1 == Q2):
            return 1
        if (self.iota(P1) == P2) and (self.iota(Q1) == Q2):
            return 2
        return 0
    
    def codomain(self):
        return self.Phi.codomain()
    
    def EvalEndomorphismType(self, type_end, T1, T2):
        # Endomorphisms of E0xE0
        if type_end == 1:
            return T1 + T2, T1 - T2
        else:
            iota_T1 = self.iota(T1)
            return iota_T1 + T2, iota_T1 - T2
        
    def CheckDiagonalIsogeny(self, P1, P2, Q1, Q2):

        if P1 == 0 or Q2 == 0:
            return True, Q1, P2
        elif P2 == 0 or Q1 == 0:
            return True, P1, Q2
        else:
            return False, 0, 0
        
    def DiagonalIsogeny(self, R1, R2):
        E1 = R1.curve()
        E2 = R2.curve()

        phi_1 = E1.isogeny(R1, degree=2, model="montgomery")
        phi_2 = E2.isogeny(R2, degree=2, model="montgomery")

        return phi_1, phi_2
    
    def TwoTorsionImage(self, P, Q, e, u):
        assert all([2**(e-1)*R for R in [P, Q]])
        assert all([not 2**e*R for R in [P, Q]])

        ePQ = P.weil_pairing(Q, 2**e)

        PmQ = P-Q
        imP = self.EvalEmbedded(P)
        imQ = self.EvalEmbedded(Q)
        imPQ = self.EvalEmbedded(PmQ)

        E1, E2 = self.codomain()

        imP1 = imP[0]
        imP2 = imP[1]

        imQ1 = imQ[0]
        imQ2 = imQ[1]

        imPQ1 = imPQ[0]
        imPQ2 = imPQ[1]

        imP1, imQ1 = fix_minus_sign(imP1, imQ1, imPQ1)
        imP2, imQ2 = fix_minus_sign(imP2, imQ2, imPQ2)
        if ePQ**u == imP1.weil_pairing(imQ1, 2**e):
            return E1, imP1, imQ1
        
        if ePQ**u == imP2.weil_pairing(imQ2, 2**e):
            return E2, imP2, imQ2
        
        raise ValueError

    def EvalEmbedded(self, P):

        PP = (P, self.domain[1](0))
        for type_end in self.endomorphisms:
            PP = self.EvalEndomorphismType(type_end, PP[0], PP[1])

        for phi_1, phi_2 in self.diagonal_isogenies:
            PP = (phi_1(PP[0]), phi_2(PP[1]))

        return self.Phi(CouplePoint(PP[0], PP[1]))


def fix_minus_sign(P, Q, PQ):
    if P-Q == PQ or P-Q == -PQ:
        return P, Q
    else:
        assert P+Q == PQ or P+Q == -PQ
        return P, -Q