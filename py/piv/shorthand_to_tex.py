__author__ = 'Jwely'


def _character_to_symbol(character):
    """ converts a single character to the mathematic symbol, does not apply all tex formatting """
    if character == "T":
        return r"v_{\theta}"
    elif character == "W":
        return r"v_{z}"
    elif character == "R":
        return r"v_{r}"
    elif character == "U":
        return r"v_{x}"
    elif character == "V":
        return r"v_{y}"
    elif character == "t":
        return r"v_{\theta}^{\prime}"
    elif character == "w":
        return r"v_{z}^{\prime}"
    elif character == "r":
        return r"v_{r}^{\prime}"
    elif character == "u":
        return r"v_{x}^{\prime}"
    elif character == "v":
        return r"v_{y}^{\prime}"

    elif character == "M":
        return r"v_{magnitude}"
    elif character == "P":
        return r"v_{in-plane}"
    else:
        return character


def _overline(symbol_list):
    """ applies the overline function """
    if not isinstance(symbol_list, list):
        symbol_list = [symbol_list]
    return r"\overline{{{0}}}".format(" ".join(symbol_list))


def _tex(string):
    """ applies the dollar signs needed to bound math expressions """
    return "${0}$".format(string)


def shorthand_to_tex(component):
        """
        Converts component symbols used extensively in the code syntax into TeX formatted math expressions.

        :param component: some shorthand component such as 'U', 'V', 'rr', 'rw', 'ctke', 'r_mesh', etc...
        :return: a tex formatted string that will parse to a mathematical expression for that components
        shorthand
        """


        if component == "turb_visc":        # turbulent viscosity
            return r"$\nu_{t}$"

        elif component == "ctke":           # turbulent energy
            return r"$k$"

        elif component == "num":            # number of measurements N
            return r"$N$"

        elif component == "dPdr":
            return r"$\frac{d\bar{P}}{dr}$"

        # handles single components and double correlation components.
        elif len(component) == 1:
            return _tex(_overline(_character_to_symbol(component)))
        elif len(component) == 2:
            return _tex(_overline([_character_to_symbol(component[0]), _character_to_symbol(component[1])]))

        # handles partial derivatives of format 'dxdy'
        elif len(component) == 4 and "d" in component:
            comp1 = component[1]
            comp2 = component[3]
            fmt = r"\frac{{\partial {0}}}{{\partial {1}}}"
            return _tex(_overline(fmt.format(_character_to_symbol(comp1), comp2)))

        # meshgrid attributes
        elif "mesh" in component:
            return "${0}$".format(component.upper().split("_")[0])

        # give the user their component string back if nothing matches
        else:
            Warning("Don't recognize input {0}".format(component))
            return component


# testing
if __name__ == "__main__":

    from py.piv import AxialVortex
    a = AxialVortex()
    keys = a.meshgrid.keys() + a.mean_set.keys() + a.derivative_set.keys() + a.equation_terms.keys()
    for key in keys:
        print(key, shorthand_to_tex(key))

