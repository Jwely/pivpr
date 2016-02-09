__author__ = 'Jwely'


def shorthand_to_tex(component):
        """
        Converts component symbols used extensively in the code syntax into TeX formatted math expressions.

        :param component: some shorthand component such as 'U', 'V', 'rr', 'rw', 'ctke', 'r_mesh', etc...
        :return: a tex formatted string that will parse to a mathematical expression for that components
        shorthand
        """

        # turbulent viscosity is a special case
        if component == "turb_visc":
            return "$\\nu_{t}$"

        # turbulent energy is a special case
        elif component == "ctke":
            return "$k$"

        elif component == "num":
            return "$N$"

        # handles single components
        elif len(component) == 1:

            # stable components
            if component.isupper():
                return "$\overline{{{comp}}}$".format(comp=component.lower())

            # fluctuating components
            else:
                return "$\overline{{{comp}^\prime}}$".format(comp=component)

        # double correlations, reynolds stresses.
        elif len(component) == 2:
            # note the use of double enclosing brackets, this is required to eval literal {}
            return "$\overline{{{comp0}^\prime {comp1}^\prime}}$".format(comp0=component[0], comp1=component[1])

        # handles partial derivatives of format 'dxdy'
        elif len(component) == 4 and "d" in component:
            return "$\\frac{{\\partial {0}}}{{\\partial {1}}}$".format(component[1], component[3])

        # meshgrid attributes
        elif "mesh" in component:
            return "${0}$".format(component.upper().split("_")[0])

        # give the user their component string back if nothing matches
        else:
            Warning("Don't recognize input {0}".format(component))
            return component