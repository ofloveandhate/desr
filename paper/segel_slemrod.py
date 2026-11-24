"""
Michaelis-Menten kinetics on two timescales: the analysis of Segel and Slemrod.

Their "outer" equations (8a-b) and their reduction (21a-d), each reached as a
nondimensionalisation, with the alternatives obtained by column operations on the
Hermite multiplier.  It is worked out on the documentation page
`docs/source/examples/segel_slemrod.rst`, and this prints every computation on that
page.

L. A. Segel and M. Slemrod, The quasi-steady-state assumption: a case study in
perturbation, SIAM Review 31 (1989) 446-477, https://doi.org/10.2307/2031405.
"""
from _from_page import print_page

if __name__ == '__main__':
    print_page('segel_slemrod.rst')
