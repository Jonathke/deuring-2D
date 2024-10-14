from deuring2d import *
from special_extremal import *
from quaternions import *

def TestSpecialExtremal():
    p = 2**50*11 - 1
    O0 = SpecialExtremalOrder(p, 1)
    N = randint(floor(p**1.1), ceil(2*p**1.1))
    gamma = O0.FullRepresentInteger(N)
    assert gamma
    assert gamma.reduced_norm() == N
    print("TestSpecialExtremal passed!")

def TestFixedDegreeIsogeny():
    p = 2**50*11 - 1
    ctx = Deuring2D(p)

    u = randint(2**35, 2**36)
    while u%2 == 0:
        u = randint(2**35, 2**36)
    Phi, imP, imQ = ctx.FixedDegreeIsogeny(u)
    print(Phi.codomain()[0])
    print(Phi.codomain()[0].j_invariant())
    print(Phi.codomain()[1])
    print(Phi.codomain()[1].j_invariant())

def TestIdealToIsogeny():
    p = 2**50*11 - 1
    ctx = Deuring2D(p)
    l = next_prime(randint(2**50, 2**51))
    #ll = next_prime(l)
    #while not (ll % 4 == l % 4 == 1):
    #    l = next_prime(randint(2**20, 2**21))
    #    ll = next_prime(randint(2**20, 2**21))
    ll = next_prime(randint(2**50, 2**51))
    alpha = ctx.O0.FullRepresentInteger(l*ll)
    I = ctx.O0.order * l + ctx.O0.order * alpha
    E_I, phi_IP, phi_IQ = ctx.IdealToIsogeny(I)

    print(E_I)
    print(E_I.j_invariant())

def TestIdealToIsogenyBig():
    print("---- Testing big stuff ----")
    p = 2**500*27 - 1
    ctx = Deuring2D(p)
    l = next_prime(randint(p, 100*p))
    #ll = next_prime(l)
    #while not (ll % 4 == l % 4 == 1):
    #    l = next_prime(randint(2**20, 2**21))
    #    ll = next_prime(randint(2**20, 2**21))
    ll = next_prime(randint(p, 100*p))
    alpha = ctx.O0.FullRepresentInteger(l*ll)
    I = ctx.O0.order * l + ctx.O0.order * alpha
    E_I, phi_IP, phi_IQ = ctx.IdealToIsogeny(I)

    print(E_I)
    print(E_I.j_invariant())

if __name__=="__main__":
    proof.all(False)
    #TestSpecialExtremal()
    TestFixedDegreeIsogeny()
    TestIdealToIsogeny()
    TestIdealToIsogenyBig()