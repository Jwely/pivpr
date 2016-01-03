__author__ = 'Jwely'


def dynamic_import(module):
    """
    Simple dynamic importer that can load modules from strings
    :param module: any string that would eval as a module import
    """

    components = module.split(".")
    mod = __import__(components[0])
    for comp in components[1:]:
        mod = getattr(mod, comp)
    return mod


def create_all_final_figures():
    """
    Iterates through every figure generating class in the list below. Take a look at the
    FinalFigure class for better understanding of requirements for each entry in the
    figs_class_list.
    """

    figs_class_list = ["final_figures.ExamplePlots55",
                       ]

    for fig_class_name in figs_class_list:
        fig_class = dynamic_import(fig_class_name)
        fig = fig_class()
        fig.generate()


if __name__ == "__main__":
    create_all_final_figures()

