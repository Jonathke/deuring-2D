#!/usr/bin/env sage

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
    print("> Testing ideal to isogeny")
    p = 2**50*11 - 1
    ctx = Deuring2D(p)
    l = next_prime(randint(2**50, 2**51))
    #ll = next_prime(l)
    #while not (ll % 4 == l % 4 == 1):
    #    l = next_prime(randint(2**20, 2**21))
    #    ll = next_prime(randint(2**20, 2**21))
    ll = next_prime(randint(2**50, 2**51))
    alpha = ctx.O0.FullRepresentInteger(l*ll)
    if alpha is None:
        raise Exception('test failed: FullRepresentInteger() did not find alpha')
    I = ctx.O0.order * l + ctx.O0.order * alpha
    E_I, phi_IP, phi_IQ = ctx.IdealToIsogeny(I)

    print(E_I)
    print(E_I.j_invariant())
    print(">>> Success!")

def TestIdealToIsogenyOwn():
    print("> Testing ideal to isogeny (with externally supplied order)")
    p = 2**127 - 1
    F.<i> = GF((p, 2), modulus=[1,0,1])
    E0 = EllipticCurve(F, [73300411938425210204353420537549415058, 96840771522044021527333883178334690669])
    q = 23
    iota = E0.isogeny([155415534247160515642272280914556724100, 148981912015670309004775159867160076038, 61025058819502448980768543275782209691, 151212779892207632179602417493632419553, 71209404954542770202672836886535034920, 147669429963117838129792671368121776970, 9398173484427566967687441310036216901, 134509300987373454650287121289856282090, 163854451224215204029367366925403956845, 12159765651977588628590571553508757087, 152735015202933157213038210414289363095, 1], codomain=E0)
    assert iota.degree() == q
    Quat.<i,j,k> = QuaternionAlgebra(-q, -p)
    O0 = Quat.quaternion_order([Quat.one(), 1/2 + 1/2*i, j, 6/23*i + 1/2*j + 1/46*k])
    ctx = Deuring2D(E0=E0, O0=O0, iota=iota)
    l = next_prime(randint(2**50, 2**51))
    I = O0.random_ideal('left', l)
    E_I, phi_IP, phi_IQ = ctx.IdealToIsogeny(I)
    print(">>> Success!")

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
    if alpha is None:
        raise Exception('test failed: FullRepresentInteger() did not find alpha')
    I = ctx.O0.order * l + ctx.O0.order * alpha
    E_I, phi_IP, phi_IQ = ctx.IdealToIsogeny(I)

    print(E_I)
    print(E_I.j_invariant())
    print(">>> Success!")

if __name__=="__main__":
    proof.all(False)
    #TestSpecialExtremal()
    #TestFixedDegreeIsogeny()
    TestIdealToIsogeny()
    TestIdealToIsogenyOwn()
    TestIdealToIsogenyBig()

