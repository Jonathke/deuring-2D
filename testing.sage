from deuring import *
from special_extremal import *

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


def TestIdealToIsogeny():
    p = 2**50*11 - 1
    ctx = Deuring2D(p)
    l = next_prime(randint(2**60, 2**61))
    alpha = ctx.O0.FullRepresentInteger(l*next_prime(l))
    I = ctx.O0.order * l + ctx.O0.order * alpha
    ctx.IdealToIsogeny(I)

if __name__=="__main__":
    #TestSpecialExtremal()
    TestFixedDegreeIsogeny()
    TestIdealToIsogeny()