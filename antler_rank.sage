
from sage.crypto.boolean_function import BooleanFunction
from sage.rings.polynomial.pbori.pbori import *

B.<x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10> = BooleanPolynomialRing()

def random_S(n):
    """
    Sample a negacyclic lookup table uniformly at random.
    """
    L = []
    
    for i in range(2**(n-1)):
        L.append(randrange(2**n))
        
    for i in range(2**(n-1)):
        L.append((-L[i]) % 2**n)
    return L

def antler(x, y, z, L1, L2, L3, n):
    """
    Evaluate an antler function (with constants added to the input if needed).
    """
    return L3[(z + L1[x] + L2[y]) % 2**n]

def antler_one_var(x_compl, L1, L2, L3, n):
    """
    Compute the antler function using only one input variable.
    """
    x =  x_compl %  2**n
    y = (x_compl // 2**n) % 2**n
    z =  x_compl // 2**(2*n)

    return antler(x, y, z, L1, L2, L3, n) % 2

def anf_antler(L1, L2, L3, n):
    """
    Compute the algebraic normal form of an Antler function defined with the 3 NLUTs L1, L2, L3.
    """
    truth_table = []
    
    for x_compl in range(2**(3*n-1)):
        truth_table.append(antler_one_var(x_compl, L1, L2, L3, n))
    return BooleanFunction(truth_table).algebraic_normal_form()



def get_rank(n, N_iterations):
    """
    Get the rank of the full family of Antler functions.
    """
    l = []

    for k in range(N_iterations):
        L1, L2, L3 = random_S(n), random_S(n), random_S(n)
        l.append(anf_antler(L1, L2, L3, n))

    seq  = Sequence(l, B)
    A, _ = seq.coefficient_matrix()
    return A.rank() , A.ncols()


def test_3_family():
    """
    Run the experiments for n=3.
    """
    n = 3
    N_iterations = 500

    family_rank , n_monomials = get_rank(n, N_iterations)

    print("The rank of the space of Antler functions for n=3 is: {}".format(family_rank))
    print("The number of monomials for n=3 is: {}".format(n_monomials))


def test_4_family():
    """
    Run the experiments for n=4.
    """
    n = 4
    N_iterations = 1500

    family_rank , n_monomials = get_rank(n, N_iterations)

    print("The rank of the space of Antler functions for n=4 is: {}".format(family_rank))
    print("The number of monomials for n=4 is: {}".format(n_monomials))


if __name__=="__main__":
    test_3_family()
    test_4_family()

