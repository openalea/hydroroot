from millet.architecture import millet_mtg


def test_mtg():
    g = millet_mtg(nb_vertices=100, seed=1)
    seminal = g.Trunk(g.root, Scale=1)
    assert len(seminal) == 102

    # test if every vertex has a label
    labels = g.property('label')
    assert len(labels) == len(g) - 1

    # test if we have both crown and seminal
    diff_labels = set(labels.itervalues())
    assert len(diff_labels) == 3
    print diff_labels
test_mtg()

