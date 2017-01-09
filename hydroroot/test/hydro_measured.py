"""

"""
import pandas as pd
import numpy as np

from openalea.deploy.shared_data import shared_data

import random

from openalea.mtg import MTG, fat_mtg
from openalea.mtg.traversal import post_order2
from openalea.mtg import algo

import hydroroot
from hydroroot.length import fit_law
from hydroroot import radius, markov, flux, conductance, measured_root



# order 0
#order0 = df[df.order=='1']
#print order0

def read_data(fn):
    df = pd.read_csv(fn, sep='\t')
    df['db'] = df['distance_from_base_(mm)']/1.e3
    df['lr'] = df['lateral_root_length_(mm)']/1.e3

    df.path = df.order.str.split('-')

    n = max(df.path.map(len))

    for i in range(1,n):
        df['order_%d'%i] = df.path.map(lambda x: 0 if len(x) ==1 else int(x[i]))
    return df

############################

def reconstruct(df, segment_length=1e-4):

    g = MTG()
    rid = vid = g.add_component(g.root)

    rnid = g.node(rid)
    rnid.base_length = 0.
    rnid.length = segment_length
    rnid.label = 'S'
    rnid.order = 0

    # Primary root

    # g.node(vid).tip_length = total_length
    prev_len = 0.
    ramifs = {}

    df_order = df[df.order=='1']
    length_base = df_order.index

    path = [1]
    count = 0
    for i in length_base:
        len_base = df_order.iloc[i].db
        code='1'
        while len_base-prev_len > 0:
            prev_len += segment_length
            vid = g.add_child(vid, edge_type='<', label='S', base_length=prev_len, length=segment_length, order=0, code=code)
        vid = g.add_child(vid, edge_type='<', label='S', base_length=len_base, length=segment_length, order=0)
        prev_len = len_base
        len_lateral = df_order.iloc[i].lr
        if len_lateral > 0.:
            count += 1
            p = tuple([1, count])
            ramifs.setdefault(p, []).append((vid, len_lateral))

    new_ramifs = {}
    len_base = 0

    for path in ramifs:
        vid, lr = ramifs[path][0]
        order = '-'.join(map(str,path))
        df_order = df[df.order==order]
        length_base = df_order.index
        code = order
        len_base = lr

        prev_len = 0.
        _root_id = vid
        parent_base = g.node(vid).base_length

        if df_order.empty:
            while len_base-prev_len > 0:
                prev_len += segment_length
                edge_type = '+' if _root_id == vid else '<'
                vid = g.add_child(vid, edge_type=edge_type, label='S', base_length=parent_base+prev_len, length=segment_length, order=1, code=code)
        else:
            for i in length_base:
                len_base = df_order.db[i]
                edge_type = '+' if _root_id == vid else '<'
                while len_base-prev_len > 0:
                    prev_len += segment_length
                    vid = g.add_child(vid, edge_type=edge_type, label='S', base_length=parent_base+prev_len, length=segment_length, order=1, code=code)
                    edge_type = '<'
                vid = g.add_child(vid, edge_type='<', label='S', base_length=parent_base+len_base+segment_length, length=segment_length, order=1, code=code)
                prev_len = len_base
                len_lateral = df_order.lr[i]
                if len_lateral > 0.:
                    count += 1
                    p = tuple([2, count])
                    new_ramifs.setdefault(p, []).append((vid, len_lateral))

    ramifs, new_ramifs = new_ramifs, {}
    len_base = 0
    Order = 2
    # Copy & Paste
    for path in ramifs:
        vid, lr = ramifs[path][0]
        order = '-'.join(map(str,path))
        df_order = df[df.order==order]
        length_base = df_order.index
        code = order

        len_base = lr

        prev_len = 0.
        _root_id = vid
        parent_base = g.node(vid).base_length
        new_ramifs = {}
        if df_order.empty:
            while len_base-prev_len > 0:
                prev_len += segment_length
                edge_type = '+' if _root_id == vid else '<'
                vid = g.add_child(vid, edge_type=edge_type, label='S', base_length=parent_base+prev_len, length=segment_length, order=Order, code=code)
        else:
            for i in length_base:
                len_base = df_order.db[i]
                edge_type = '+' if _root_id == vid else '<'
                while len_base-prev_len > 0:
                    prev_len += segment_length
                    vid = g.add_child(vid, edge_type=edge_type, label='S', base_length=parent_base+prev_len, length=segment_length, order=Order, code=code)
                    edge_type = '<'
                vid = g.add_child(vid, edge_type='<', label='S', base_length=parent_base+len_base+segment_length, length=segment_length, order=Order, code=code)
                prev_len = len_base
                len_lateral = df_order.lr[i]
                if len_lateral > 0.:
                    count += 1
                    p = tuple([2, count])
                    new_ramifs.setdefault(p, []).append((vid, len_lateral))

    return g

def my_flux(g,
            segment_length=1e-4,
            ref_radius=1e-4,
            order_decrease_factor=0.7,
            k0=300,
            Jv=0.1,
            psi_e=0.4,
            psi_base=0.1,
            axial_conductivity_data=None,
            radial_conductivity_data=None,
            ):

    xa, ya = axial_conductivity_data
    ya = list(np.array(ya) * (segment_length / 1.e-4))
    axial_conductivity_law = fit_law(xa, ya)

    xr, yr = radial_conductivity_data
    radial_conductivity_law = fit_law(xr, yr)

    # compute radius property on MTG
    g = radius.ordered_radius(g, ref_radius=ref_radius, order_decrease_factor=order_decrease_factor)

    # compute length property and parametrisation
    g = radius.compute_length(g, segment_length)
    g = radius.compute_relative_position(g)

    # Compute K using axial conductance data
    g = conductance.fit_property_from_spline(g, axial_conductivity_law, 'position', 'K')

    g, surface = radius.compute_surface(g)
    g, volume = radius.compute_volume(g)

    # Compute the flux

    g = conductance.fit_property_from_spline(g, radial_conductivity_law, 'position', 'k0')
    g = conductance.compute_k(g, k0='k0')

    # TODO: return Keq base and Jv
    g = flux.flux(g, Jv, psi_e, psi_base, invert_model=True)

    Keqs = g.property('Keq')
    v_base = g.component_roots_at_scale_iter(g.root, scale=1).next()

    Keq = Keqs[v_base]
    Jv_global = Keq * (psi_e - psi_base)

    return g, surface, volume, Keq, Jv_global

def run(fn):
    ref_radius = 1e-4 # in m
    order_decrease_factor = 1.

    # parameters
    k0 = 400.
    Jv = 0.1
    psi_e = 0.4
    psi_base = 0.

    # laws
    acol = axial_conductivity_data = (
        [0., 0.03,  0.06, 0.09, 0.12, 0.15, 0.18],
        [2.9e-4, 34.8e-4, 147.4e-4, 200.3e-4,292.6e-4,262.5e-4,511.1e-4]
    )
    def radial(v=300, scale=1):
        xr = [0., 0.015, 0.03, 0.045, 0.06, 0.075, 0.09, 0.105, 0.135, 0.15, 0.16]
        yr = [v*scale]*len(xr)
        return xr, yr

    def axial(scale=1):
        x, y = acol
        y = [a*scale for a in y]
        return x, y

    rcol = radial_conductivity_data = radial(300)

    axfold=1
    radfold=1


    df = read_data(fn)
    g = reconstruct(df)
    g, surface, volume, Keq, Jv_global = my_flux(g,
            ref_radius=ref_radius, order_decrease_factor=order_decrease_factor,
            k0=k0, Jv=Jv, psi_e=psi_e, psi_base=psi_base,
            axial_conductivity_data=acol, radial_conductivity_data=rcol)

    return g, surface, volume, Keq, Jv_global

def plot(g):
    from IPython import get_ipython

    ipython = get_ipython()
    if ipython:
        ipython.magic('gui qt')

    from openalea.plantgl.all import Viewer
    from hydroroot.display import plot as _plot



    Viewer.display(_plot(g))


if __name__ == '__main__':
    files = shared_data(hydroroot, share_path='share', pattern='hydroroot_measured*.txt')
    fn = files[0]

    g, surface, volume, Keq, Jv_global= run(fn)
