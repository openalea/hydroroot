""" Test architectural model built from data.

"""
import pandas

from openalea.hydroroot.generator.measured_root import mtg_from_aqua_data # Added F. Bauget 2019-12-16


def test_reconstruct_from_aqua_data():
    """ Added F. Bauget 2019-12-16
    test the reconstruction from a data given in the aquaporin format see hydroroot.generator.measured_root
    the file to test is very simple 1 primary, 2 1st order lateral and 2 2d order ones
    the real total length is 3.1 mm. It has 4 tips so the maximum length with e discretization of segment dl
    is 3.1 + 4*dl
    We test the total length of this against the total length from MTG

    """
    segment_length = 1.0e-4
    fn='data/test_reconstruct_from_aqua_data.txt'
    df = pandas.read_csv(fn, sep = '\t')
    df['db'] = df['distance_from_base_(mm)'] / 1.e3
    df['lr'] = df['lateral_root_length_(mm)'] / 1.e3
    g = mtg_from_aqua_data(df, segment_length=segment_length)

    total_length = g.nb_vertices(scale = 1) * segment_length

    n = len(df[df.lr > 0]['lr']) + 1 # nb of tips from data nb of lateral + the PR
    real_total_length = max(df.db[df.order == "1"]) + sum(df['lr'])
    max_total_length = real_total_length + n * segment_length

    assert real_total_length <= total_length <= max_total_length, "error on total length"
