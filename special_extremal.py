from sage.all import *
from quaternions import MakePrimitive

# proof.all(False)
# Class for a special p-extremal order, with corresponding functions
class SpecialExtremalOrder:
    def __init__(self, p, q):
        if not (p % 4 == 3 and q == 1):
            raise NotImplementedError
        self.p = p
        self.q = q
        B = QuaternionAlgebra(-q, -p)
        i, j, k = B.gens()
        self.i = i
        self.j = j
        self.k = k
        self.order = B.quaternion_order([1, i, (i+j)/2, (1+k)/2])
        self.QF = BinaryQF([1,0,q])
        
    def Cornacchia(self, N):
        assert N > 0
        smallfacs = [l**e for l,e in factor(N, limit=1000) if l < 1000]
        N_prime = prod(smallfacs)
        if len(smallfacs) > 5 or not is_pseudoprime(N // N_prime):
            # self.QF.solve_integer() will be expensive -- bail out
            return
        algo = ('general', 'cornacchia')[is_pseudoprime(N)]
        return self.QF.solve_integer(N, algorithm=algo)
    
    def RepresentInteger(self, N):
        if N < self.p:
            if not (sol := self.Cornacchia(N)):
                return
            x, y = sol
            return x + y*self.i

        m1 = isqrt(N / self.p)
        assert N - self.p * m1**2 > 0

        for _ in range(1000):
            z = randint(-m1, m1)
            m2 = isqrt(N / self.p - z**2)
            assert 4*N - self.p * self.QF(z, m2) > 0
            w = randint(-m2, m2)
            Mm = N - self.p * self.QF(z, w)
            if not (sol := self.Cornacchia(Mm)):
                continue
            x, y = sol
            gamma = x + y*self.i + z*self.j + w*self.k
            gamma = MakePrimitive(gamma, self.order)
            if gamma.reduced_norm() == N:
                return gamma

    # Is there any reason to use this for us? Want to avoid using to avoid point division.
    def FullRepresentInteger(self, N):
        if N < self.p:
            if not (sol := self.Cornacchia(N)):
                raise ValueError
            x, y = sol
            return x + y*self.i
            
        m1 = isqrt(4*N / self.p)
        assert 4*N - self.p * m1**2 > 0

        for _ in range(1000):
            z = randint(-m1, m1)
            m2 = isqrt(4*N / self.p - z**2)
            assert 4*N - self.p * self.QF(z, m2) > 0
            w = randint(-m2, m2)
            Mm = 4*N - self.p*self.QF(z, w)
            if not (sol := self.Cornacchia(Mm)):
                continue
            x, y = sol
            if x % 2 == w % 2 and y % 2 == z % 2 and gcd([x,y,z,w]) == 1:
                gamma = x + y*self.i + z*self.j + w*self.k
                gamma = MakePrimitive(gamma,self.order)
                if gamma.reduced_norm() == N:
                    return gamma

