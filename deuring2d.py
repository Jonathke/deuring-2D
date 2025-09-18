from sage.all import *
from special_extremal import SpecialExtremalOrder
from quaternions import ReducedBasis, SuccessiveMinima, reduced_ideal
from product_isogeny import IsogenyProduct
from utilities.discrete_log import BiDLP
from qlapoti.qlapoti import solve as qlapoti

class Deuring2D:

    class Failure(Exception): pass

    def __init__(self, p):
        if p % 4 != 3:
            raise NotImplementedError
        self.p = p

        self.O0 = SpecialExtremalOrder(self.p, 1)

        self.F = GF((p, 2), name = "i", modulus=[1,0,1])
        self.E0 = EllipticCurve(self.F, [1, 0])
        self.e = Integer(p+1).valuation(2) #can add -1 here to avoid field extensions
        
        #generate 2-torsion basis
        cofac = (p+1)/(2**self.e)
        while True:
            P = self.E0.random_point()*cofac
            P_small = P*2**(self.e-1)
            if P_small:
                break
        while True:
            Q = self.E0.random_point()*cofac
            Q_small = Q*2**(self.e-1)
            if Q_small and Q_small != P_small:
                break

        self.P = P
        self.Q = Q

        self.sqrtm1 = self.F(-1).sqrt()

        def endo_i(P):
            Ebig = P.curve()
            Fbig = Ebig.base_field()
            x, y = P.xy()
            return Ebig(-x, Fbig(self.sqrtm1)*y)
        
        def endo_j(P):
            Ebig = P.curve()
            pi = Ebig.base_field().frobenius_endomorphism()
            x,y = P.xy()
            return Ebig(pi(x), pi(y))
        
        self.iota = endo_i
        self.pi = endo_j


    def EvalEndomorphism(self, alpha, P, ord):

        assert not self.E0.defining_polynomial()(*P)
        assert P*ord == 0

        d = lcm(c.denominator() for c in alpha)

        if gcd(d, ord) == 2:
            #TODO should probably make sure this doesnt need to be called
            alpha = 2*alpha
            try:
                P = P.division_points(2)[0]
            except IndexError:
                Fbig, _ = self.F.extension(4,'A').objgen()
                Ebig = self.E0.base_extend(Fbig)
                P = Ebig(P).division_points(2)[0]
        elif gcd(d, ord) > 1:
            raise NotImplementedError('denominators > 2 are currently unsupported')

        iP = self.iota(P)
        jP = self.pi(P)
        kP = self.iota(jP)
        coeffs = [coeff % (ord*d) for coeff in alpha]
        return self.E0(sum(c*Q for c, Q in zip(coeffs, [P, iP, jP, kP])))


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
