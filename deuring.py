from sage.all import *
from special_extremal import SpecialExtremalOrder
from quaternions import ReducedBasis
from product_isogeny import IsogenyProduct
from utilities.discrete_log import BiDLP

class Deuring2D:
    def __init__(self, p):
        if p % 4 != 3:
            raise NotImplementedError
        self.p = p

        self.O0 = SpecialExtremalOrder(self.p, 1)

        self.F = GF((p, 2), name = "i", modulus=[1,0,1])
        self.E0 = EllipticCurve(self.F, [1, 0])
        self.e = Integer(p+1).valuation(2) #the -1 is there to avoid field extensions
        
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

        assert P.curve() == self.E0
        assert P*ord == 0

        d = lcm(c.denominator() for c in alpha)

        if gcd(d, ord) == 2:
            #TODO should probably make sure this doesnt need to be called
            alpha = d*alpha
            Fbig, _ = self.F.extension(4,'A').objgen()
            Ebig = self.E0.base_extend(Fbig)
            P = Ebig(P).division_points(2)[0]

        iP = self.iota(P)
        jP = self.pi(P)
        kP = self.iota(jP)
        coeffs = [coeff % (ord*d) for coeff in alpha]
        return self.E0(sum(c*Q for c, Q in zip(coeffs, [P, iP, jP, kP])))


    def FixedDegreeIsogeny(self, u):
        # Input: Odd integer u, basis P, Q of E_0[2^e]
        # Output: phi: E_0 -> E, a u-isogeny embedded in a (2^e, 2^e)-isogeny
        assert u%2 == 1
        assert u < 2**self.e

        theta = self.O0.FullRepresentInteger(u*(2**self.e - u))
        assert theta.reduced_norm() == u*(2**self.e - u)
        thetaP = self.EvalEndomorphism(theta, self.P, 2**self.e)
        thetaQ = self.EvalEndomorphism(theta, self.Q, 2**self.e)

        Phi = IsogenyProduct(u*self.P, thetaP, u*self.Q, thetaQ, self.e, iota=self.iota)
        
        return Phi, Phi.EvalTopLeft(self.P), Phi.EvalTopLeft(self.Q)


    def SuitableIdeals(self, I):
        basis_I = ReducedBasis(I.basis())
        N_I = I.norm()
        d = ceil(100*sqrt(self.p))
        Bs = [floor(sqrt(d/(alpha.reduced_norm()/N_I))/4) for alpha in basis_I]
        print(Bs)
        for _ in range(10000):
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
            e1 = min(u.valuation(2), v.valuation(2))
            u = Integer(u/(2**e1))
            v = Integer(v/(2**e1))
            f = self.e - e1
            return beta_1, beta_2, u, v, f
        
        raise ValueError

    def IdealToIsogeny(self, I):
        # Input: Ideal I
        # Output: E_I
        beta_1, beta_2, u, v, f = self.SuitableIdeals(I)
        N_I = I.norm()
        d1 = beta_1.reduced_norm()/N_I
        theta = (beta_2*beta_1.conjugate())/N_I

        Phi_u, uP, uQ = self.FixedDegreeIsogeny(u)
        Phi_v, vP, vQ = self.FixedDegreeIsogeny(v)

        print(theta)
        e1 = 2**(self.e-f)
        print(f"{f=}")
        print(f"{e1=}")

        thetaP = self.EvalEndomorphism(theta, self.P, 2**self.e)
        thetaQ = self.EvalEndomorphism(theta, self.Q, 2**self.e)

        beta1P = self.EvalEndomorphism(beta_1, self.P, 2**self.e)
        beta2P = self.EvalEndomorphism(beta_2, self.P, 2**self.e)
        beta1Q = self.EvalEndomorphism(beta_1, self.Q, 2**self.e)
        beta2Q = self.EvalEndomorphism(beta_2, self.Q, 2**self.e)

        x1_P, y1_P = BiDLP(beta1P, self.P, self.Q, 2**self.e)
        x2_P, y2_P = BiDLP(beta2P, self.P, self.Q, 2**self.e)

        x1_Q, y1_Q = BiDLP(beta1Q, self.P, self.Q, 2**self.e)
        x2_Q, y2_Q = BiDLP(beta2Q, self.P, self.Q, 2**self.e)
        #assert thetaP == x_P*self.P + y_P*self.Q
        #assert thetaQ == x_Q*self.P + y_Q*self.Q

        Phi = IsogenyProduct(e1*d1*uP, e1*Phi_v.EvalTopLeft(thetaP), e1*d1*uQ, e1*(Phi_v.EvalTopLeft(thetaQ)), f)
        #Phi = IsogenyProduct(e1*(Phi_u.EvalTopLeft(beta1P)), e1*(Phi_v.EvalTopLeft(beta2P)), e1*(Phi_u.EvalTopLeft(beta1Q)), e1*(Phi_v.EvalTopLeft(beta2Q)), f)
        #from theta_structures.couple_point import CouplePoint
        #from theta_isogenies.product_isogeny_sqrt import EllipticProductIsogenySqrt

        #Phi = EllipticProductIsogenySqrt((CouplePoint(e1*d1*uP, e1*(x_P*vP + y_P*vQ)), CouplePoint(e1*d1*uQ, e1*(x_Q*vP + y_Q*vQ))), f)

        E_I, E_aux = Phi.codomain()
        print(E_I)
        print(E_aux)