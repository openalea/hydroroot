def mtg_out_flux(g):
    """
    Flux at the root base in uL/s

    :param g: MTG
    :return j
    """
    j = g.property('J_out')[1]; 

    # return outputs
    return j,
