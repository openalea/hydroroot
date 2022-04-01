# Tests of the simulated architecture
import sys
sys.path.insert(0, '../src')

from hydroroot import radius, length, flux
from hydroroot.generator import markov # 21-12-14: FB __init__.py in src not doing job

length_data = [0., 0.03, 0.05, 0.16], [0., 0., 0.01, 0.13]

def test_linear():
    g = markov.linear(5)

    assert len(g) == 6
    assert g.nb_vertices(scale=g.max_scale()) == 5
    return g

def get_orders(g):
    orders = g.property('order')
    max_order = max(orders.values())
    order_dict = {order: sum(1 for v in orders if orders[v]==order) for order in range(max_order+1)}
    return order_dict

def check_mtg(g, _orders):

    orders = get_orders(g)

    assert orders == _orders, orders


def test_markov(n=600, lr=0.1,  length=1e-4):
    """ Generate an MTG based on a 1st order markov model.

    The test check if the invariance of the number of vertices in each order.
    """
    # topology
    g = markov.markov_binary_tree(nb_vertices=n, branching_variability=lr, seed=2)

    order = {0:n, 1:3550, 2:1755}
    check_mtg(g,order)

    return g

def test_law(segment_length=1e-4):
    """ Test the validity of the law.
    The length law is an interpolation of data.
    Plus a scaling property.
    Check if everything is correct.
    """
    segment_length = segment_length
    xl, yl = length_data
    law = length.fit_law(xl, yl, scale=segment_length)

    for i in range(len(xl)):
        x, y = xl[i], yl[i]
        assert int(law(x/segment_length)) == int(y/segment_length)
    return law

def test_markov_with_length_law(n=600, beta=0.298):
    law = test_law()

    g = markov.markov_binary_tree(nb_vertices=n, branching_variability=beta,
                                  length_law= law, seed=2)
    order = {0: 600, 1: 1188}
    check_mtg(g, order)
    return g


#def test_extract length_law(n=600, beta=0.298):
def length_law():
    from hydroroot import markov, radius, length, flux

    length_data = [0., 0.03, 0.05, 0.16], [0., 0., 0.01, 0.13]
    n=1600; beta=0.298
    segment_length = 1e-4
    # test_law
    segment_length = segment_length
    xl, yl = length_data
    law = length.fit_law(xl, yl, scale=segment_length)

    g = markov.markov_binary_tree(nb_vertices=n, branching_variability=beta,
                                  length_law=law)

    X, Y = flux.ramification_length_law(g,root=1, dl=segment_length)
    length_law = length.fit_law(X, Y, ext=2)

    return length_law, X, Y
