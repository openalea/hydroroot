"""
This file is part of the OpenAlea project, a platform for
visual programming and simulation in plant science.

Revision 11/07/2025
- Depends only on mtg
- The function can be generalised to compute age of any root type
- The module can be added to hydroroot as is
"""

# -*- coding: utf-8 -*-
from openalea.mtg import algo
from openalea.mtg.algo import orders, axis
from openalea.mtg.traversal import pre_order2_with_filter, post_order2



def compute_age_with_constant_growth_speed(g,sr_age=15,delta_lr_appareance=5,delta_cr_appareance=6):
    """ compute age of different root types from growth laws"""

    max_scale = g.max_scale()
    _orders = algo.orders(g, scale=max_scale)
    ramifs = [vid for vid in g.vertices(scale=max_scale) if g.edge_type(vid)=='+']
    age = {}
    sr= g.Trunk(g.root, max_scale)
    len_sr = len(sr)
    delta_sr = sr_age / float(len_sr)


    age={}
    root = next(g.component_roots_at_scale_iter(g.root, scale=max_scale))
    for v in post_order2(g, root):
        if g.is_leaf(v):
            age[v] = sr_age

        if _orders[v]==0:
            ch = g.children(v)
            if ch:
                age[v] = min([age[cid] for cid in ch if g.edge_type(cid)=='<'])-delta_sr
            #age[v] = max(age.get(cid,sr_age) for cid in g.children(v))-delta_sr
    for v in ramifs:
        if g.label(v)=='Seminal':
            lr= algo.local_axis(g,v)
            lr=list(lr)
            len_lr = len(lr)
            p = g.parent(v)
            age_p = age[p]
            age[v] = age_p+ delta_lr_appareance
            age_first_vertex_lr = age[v]
            delta_lr = (sr_age-age_first_vertex_lr) /float(len_lr)
            for vid in lr:
                age[vid] = age_first_vertex_lr
                age_first_vertex_lr += delta_lr
        else:
            date = delta_cr_appareance
            cr = algo.axis(g, v)
            cr = list(cr)
            len_cr = len(cr)
            p= g.parent(v)
            age[v] = age[p]+ delta_cr_appareance
            age_first_vertex_cr = age[v]
            speed = (sr_age-age_first_vertex_cr) / float(len_cr)
            for cid in cr:
                age[cid] = age_first_vertex_cr
                age_first_vertex_cr += speed

    g.properties()['age'] = age
    return g


def compute_age_with_constant_growth_speed_cpl(g,
                                               sr_age=15,
                                               delta_lr_appearance=5,
                                               delta_cr_appearance=6):
    """Compute age of different root types from growth laws.

    Parameters:
    ===========
      - g: MTG
      - sr_age: seminal age
      - delta...
    """
    max_scale = g.max_scale()
    _orders = algo.orders(g, scale=max_scale)

    ramifs = [vid for vid in g.vertices(scale=max_scale) if g.edge_type(vid)!='<']

    order_axis = {}  # order : list of ramifications
    for ramif in ramifs:
        order_axis.setdefault(_orders[ramif], []).append(ramif)

    max_order = max(order_axis)

    # Compute age by traversing the axes order by order
    age = {}

    for order in range(max_order+1):
        for ramif in order_axis[order]:
            axis = g.Axis(ramif)

            len_axis = len(axis)
            pid = g.parent(ramif)
            age_base = 0.
            if pid is None:
                age_base = 0.
            elif g.label(pid)=='Crown':
                age_base = age[pid] + delta_cr_appearance
            else:
                age_base = age[pid] + delta_lr_appearance

            age[ramif] = age_base

            # hypothesis
            age_tip = sr_age

            delta_age = (age_tip - age_base) / len_axis

            if delta_age < 0.:
                print('WARNING: age is negative. Please increase seminal age or decrease delta lr appearance.')

            _age = delta_age
            for vid in axis[1:]:
                _age = _age + delta_age
                age[vid] = _age

    g.properties()['age'] = age

    return g
