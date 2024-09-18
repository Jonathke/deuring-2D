from sage.all import *
from quaternions import MakePrimitive

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
        N_prime = prod([l**e for l, e in factor(N, limit=100) if l < 100])
        if N/N_prime > 2**100 and not is_pseudoprime(N/N_prime):
            return None, None, False
        sol = self.QF.solve_integer(N)
        if not sol:
            return None, None, False
        return sol[0], sol[1], True
    
    def RepresentInteger(self, N):
        if N < self.p:
            x, y, found = self.Cornacchia(N)
            if found:
                return x + y*self.i
            else:
                raise ValueError
        m1 = max(floor(sqrt(round((N)/self.p, 5))), 100)
        for _ in range(1000):
            z = randint(-m1, m1)
            m2 = floor(sqrt(round((N-z**2)/self.p, 5)))
            w = randint(-m2, m2)
            Mm = N - self.p*self.QF(z,w)
            x, y, found = self.Cornacchia(Mm)
            if not found:
                continue
            gamma = x + y*self.i + z*self.j + w*self.k
            gamma = MakePrimitive(gamma,self.order)
            if gamma.reduced_norm() == N:
                return gamma

    # Is there any reason to use this for us? Want to avoid using to avoid point division.
    def FullRepresentInteger(self, N):
        if N < self.p:
            x, y, found = self.Cornacchia(N)
            if found:
                return x + y*self.i
            else:
                raise ValueError
            
        m1 = max(floor(sqrt(round((4*N)/self.p, 5))), 100)
        for _ in range(1000):
            z = randint(-m1, m1)
            m2 = floor(sqrt(round((4*N-z**2)/self.p, 5)))
            w = randint(-m2, m2)
            Mm = 4*N - self.p*self.QF(z,w)
            x, y, found = self.Cornacchia(Mm)
            if not found:
                continue
            if x % 2 == w % 2 and y % 2 == z % 2 and gcd([x,y,z,w]) == 1:
                gamma = x + y*self.i + z*self.j + w*self.k
                gamma = MakePrimitive(gamma,self.order)
                if gamma.reduced_norm() == N:
                    return gamma
        