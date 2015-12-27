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

    figs_class_list = ["ExamplePlots55.ExamplePlots55"]

    for fig_class_name in figs_class_list:
        fig_class = dynamic_import(fig_class_name)
        fig = fig_class()
        fig.generate()









if __name__ == "__main__":
    create_all_final_figures()

