
def axial_curve_generator(x, y_min, y_max):
    """

    :param x: 
    :param y_min: 
    :param y_max: 

    """
    y_interval = list(zip(y_min, y_max))

    def f(v):
        """

        :param v: 

        """
        y = [ym + v * (yM - ym) for ym, yM in y_interval] 
        return x, y
    return f
