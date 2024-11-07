from sage.all import *
from special_extremal import SpecialExtremalOrder
from quaternions import ReducedBasis, SuccessiveMinima
from product_isogeny import IsogenyProduct
from utilities.discrete_log import BiDLP

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
        return P.curve()(sum(c*Q for c, Q in zip(coeffs, [P, iP, jP, kP])))


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


    def SuitableIdeals(self, I):
        basis_I = ReducedBasis(I.basis())
        N_I = I.norm()
        d = ceil(300*sqrt(self.p))
        Bs = [isqrt(d/(alpha.reduced_norm()/N_I))//4 for alpha in basis_I]
        print([RR(log(lam/N_I, self.p)) for lam in SuccessiveMinima(I)])
        print(Bs)
        temp = basis_I[0]*basis_I[1].conjugate()
        for _ in range(100000):
            xs = [randint(-B, B) for B in Bs]
            ys = [randint(-B, B) for B in Bs]
            beta_1 = sum(x*alpha for x, alpha in zip(xs, basis_I))
            beta_2 = sum(y*alpha for y, alpha in zip(ys, basis_I))
            d1, d2 = beta_1.reduced_norm()/N_I, beta_2.reduced_norm()/N_I
            if d1 % 2 == 0 or d2 % 2 == 0 or gcd(d1, d2) > 1:
                continue
            u = (2**self.e*inverse_mod(d1, d2)) % d2
            v = (2**self.e - u*d1)/d2
            if v <= 0:
                continue
            assert u*d1 + v*d2 == 2**self.e
            e1 = min(u.valuation(2), v.valuation(2))
            u = Integer(u/(2**e1))
            v = Integer(v/(2**e1))
            f = self.e - e1

            d1 = beta_1.reduced_norm()/N_I
            theta = (beta_2*beta_1.conjugate())/N_I

            try:
                theta_u = self.O0.RepresentInteger(u * (2**self.e - u))
                theta_v = self.O0.RepresentInteger(v * (2**self.e - v))
            except self.Failure:
                continue

            yield beta_1, beta_2, u, v, f, theta_u, theta_v

    def IdealToIsogeny(self, I, *, suitable=None):
        # Input: Ideal I
        # Output: E_I

        N_I = I.norm()

        if suitable is None:
            try:
                suitable = next(self.SuitableIdeals(I))
            except StopIteration:
                pass
        if suitable is None:
            raise self.Failure('could not find any suitable ideals')

        beta_1, beta_2, u, v, f, theta_u, theta_v = suitable

        d1 = beta_1.reduced_norm()/N_I
        theta = (beta_2*beta_1.conjugate())/N_I

        Phi_u, uP, uQ = self.FixedDegreeIsogeny(u, theta=theta_u)
        Phi_v, vP, vQ = self.FixedDegreeIsogeny(v, theta=theta_v)

        e1 = 2**(self.e-f)

        thetaP = self.EvalEndomorphism(theta, self.P, 2**self.e)
        thetaQ = self.EvalEndomorphism(theta, self.Q, 2**self.e)

        x_P, y_P = BiDLP(thetaP, self.P, self.Q, 2**self.e)
        x_Q, y_Q = BiDLP(thetaQ, self.P, self.Q, 2**self.e)

        assert thetaP == x_P*self.P + y_P*self.Q
        assert thetaQ == x_Q*self.P + y_Q*self.Q

        Phi = IsogenyProduct(e1*d1*uP, e1*(x_P*vP + y_P*vQ), e1*d1*uQ, e1*(x_Q*vP + y_Q*vQ), f)
        E_I, phi_1P, phi_1Q = Phi.TwoTorsionImage(uP, uQ, self.e, d1*u)
        phi_1P = inverse_mod(u, 2**self.e)*phi_1P
        phi_1Q = inverse_mod(u, 2**self.e)*phi_1Q

        beta_1P = self.EvalEndomorphism(beta_1, self.P, 2**self.e)
        beta_1Q = self.EvalEndomorphism(beta_1, self.Q, 2**self.e)
        xx_P, yy_P = BiDLP(beta_1P, self.P, self.Q, 2**self.e)
        xx_Q, yy_Q = BiDLP(beta_1Q, self.P, self.Q, 2**self.e)

        phi_IP = inverse_mod(d1, 2**self.e)*(xx_P*phi_1P + yy_P*phi_1Q)
        phi_IQ = inverse_mod(d1, 2**self.e)*(xx_Q*phi_1P + yy_Q*phi_1Q)

        print(self.P.weil_pairing(self.Q, 2**self.e)**N_I)
        print(phi_IP.weil_pairing(phi_IQ, 2**self.e))

        return E_I, phi_IP, phi_IQ