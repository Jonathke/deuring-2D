from sage.all import *
from special_extremal import SpecialExtremalOrder
from quaternions import ReducedBasis, SuccessiveMinima, reduced_ideal
from product_isogeny import IsogenyProduct
from utilities.discrete_log import BiDLP
try:
    from qlapoti.qlapoti import solve as qlapoti
except ModuleNotFoundError:
    from qlapoti import solve as qlapoti

class Deuring2D:

    class Failure(Exception): pass

    def __init__(self, p=None, E0=None, O0=None, iota=None):

        if p is not None:
            if not all(t is None for t in (E0, O0, iota)):
                raise ValueError('need either (just) p or (all of) E0,O0,iota')

            if p % 4 != 3:
                raise NotImplementedError

            self.O0 = SpecialExtremalOrder(p, 1)

            self.F = GF((p, 2), name="i", modulus=[1,0,1])
            self.E0 = EllipticCurve(self.F, [1, 0])

            self.iota = self.E0.automorphisms()[-1]
            self.pi = self.E0.frobenius_isogeny()

        else:
            if any(t is None for t in (E0, O0, iota)):
                raise ValueError('need either (just) p or (all of) E0,O0,iota')

            self.E0 = E0
            self.F = self.E0.base_field()
            q, p = (-i for i in O0.quaternion_algebra().invariants())
            assert q == iota.degree()
            assert p == E0.base_field().characteristic()
            self.O0 = SpecialExtremalOrder(p, q, order=O0)
            self.iota = iota
            self.pi = self.E0.frobenius_isogeny()

        #generate 2-torsion basis
        self.e = Integer(p+1).valuation(2)  # can add -1 here to avoid field extensions
        self.P, self.Q = self.E0.torsion_basis(2**self.e)


    def EvalEndomorphism(self, alpha, P, ord):
        P.set_order(ord)
        d = lcm(c.denominator() for c in alpha)
        endo = sum(ZZ(c)*phi for c,phi in zip(alpha*d, (self.E0.identity_morphism(), self.iota, self.pi, self.iota * self.pi)))
        endo._degree = ZZ((alpha * d).reduced_norm())
        from sage.schemes.elliptic_curves.hom_fractional import EllipticCurveHom_fractional
        endo = EllipticCurveHom_fractional(endo, d, check=False)
        endo._degree = ZZ(alpha.reduced_norm())
        return endo._eval(P)


    def FixedDegreeIsogeny(self, u, *, theta=None):
        # Input: Odd integer u, basis P, Q of E_0[2^e]
        # Output: phi: E_0 -> E, a u-isogeny embedded in a (2^e, 2^e)-isogeny
        assert u%2 == 1
        assert u < 2**self.e

        if theta is None:
            try:
                theta = self.O0.RepresentInteger(u * (2**self.e - u))
            except self.Failure:
                pass
        if not theta:
            raise self.Failure(f'could not represent integer {u}*(2^{self.e}-{u})')

        assert theta.reduced_norm() == u*(2**self.e - u)
        thetaP = self.EvalEndomorphism(theta, self.P, 2**self.e)
        thetaQ = self.EvalEndomorphism(theta, self.Q, 2**self.e)

        Phi = IsogenyProduct(u*self.P, thetaP, u*self.Q, thetaQ, self.e, iota=self.iota)
        
        Ei, uP, uQ = Phi.TwoTorsionImage(self.P, self.Q, self.e, u)

        return Phi, uP, uQ


    def IdealToIsogeny(self, I, *, suitable=None):
        # Input: Ideal I
        # Output: E_I
        N_I = I.norm()
        J, beta_ij = reduced_ideal(I, return_elt=True)
        N_J = J.norm()
        if suitable is None:
            print("Running qlapoti...")
            suitable = qlapoti(J, 2**self.e, transporters=True)#, verbose=True)
        if suitable is None:
            raise self.Failure('could not find any suitable ideals')

        print('found a suitable ideal!', file=sys.stderr)

        beta_1, beta_2 = suitable

        d1 = beta_1.reduced_norm()/N_J
        theta = (beta_2*beta_1.conjugate())/N_J
        thetaP = self.EvalEndomorphism(theta, self.P, 2**self.e)
        thetaQ = self.EvalEndomorphism(theta, self.Q, 2**self.e)

        x_P, y_P = BiDLP(thetaP, self.P, self.Q, 2**self.e)
        x_Q, y_Q = BiDLP(thetaQ, self.P, self.Q, 2**self.e)

        assert thetaP == x_P*self.P + y_P*self.Q
        assert thetaQ == x_Q*self.P + y_Q*self.Q

        Phi = IsogenyProduct(d1*self.P, thetaP, d1*self.Q, thetaQ, self.e, self.iota)
        E_I, phi_1P, phi_1Q = Phi.TwoTorsionImage(self.P, self.Q, self.e, d1)

        correction_quat = (beta_1*beta_ij/N_J)
        beta_1P = self.EvalEndomorphism(correction_quat, self.P, 2**self.e)
        beta_1Q = self.EvalEndomorphism(correction_quat, self.Q, 2**self.e)

        xx_P, yy_P = BiDLP(beta_1P, self.P, self.Q, 2**self.e)
        xx_Q, yy_Q = BiDLP(beta_1Q, self.P, self.Q, 2**self.e)

        phi_IP = inverse_mod(d1, 2**self.e)*(xx_P*phi_1P + yy_P*phi_1Q)
        phi_IQ = inverse_mod(d1, 2**self.e)*(xx_Q*phi_1P + yy_Q*phi_1Q)

        assert self.P.weil_pairing(self.Q, 2**self.e)**N_I == phi_IP.weil_pairing(phi_IQ, 2**self.e)

        return E_I, phi_IP, phi_IQ
