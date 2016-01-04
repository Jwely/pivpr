__author__ = 'Jwely'

from matplotlib import cm
from py.texmanager import TeXRunFigurePage
from py.utils import shorthand_to_tex as stt



def create_appendix_figs(run_id):
    """
    creates the example figures used in the PIV section
    """

    def contour_component_plotter(components):
        """ subfunction for iterating through component contour plots """
        for component in components:
            figdoc.add_contour_plot(
                component, "Plot of {0} for run ID {1}".format(stt(component), run_id), contour_width)

    contour_width = '5in'
    scatter_width = '7in'
    quiver_width = '5in'
    stream_width = '5in'

    figdoc = TeXRunFigurePage("../../texdocs/main.tex",
                              "appendix_run_id_{0}".format(run_id),     # note the color scheme tag
                              run_id,
                              force_recalc=True)

    # populate the figdoc with content
    figdoc.add_text("\subsection{{Run ID {0}}}".format(run_id))
    figdoc.add_quiver_plot("Quiver plot of run ID {0}.".format(run_id), quiver_width)
    figdoc.add_stream_plot("Stream plot of run ID {0}.".format(run_id), stream_width)
    contour_component_plotter(['num', 'ctke'])

    figdoc.add_text("\subsubsection{Cylindrical Coordinates}")
    contour_component_plotter(['R', 'T', 'W', 'rt', 'rw', 'tw', 'rr', 'tt', 'ww'])

    figdoc.add_text("\subsubsection{Cartesian Coordinates}")
    contour_component_plotter(['P', 'U', 'V', 'W', 'uv', 'vw', 'uw', 'uu', 'vv', 'ww'])

    figdoc.write()


def create_all_appendix_figs():
    run_ids = range(1, 71)
    for run_id in run_ids:
        create_appendix_figs(run_id)



if __name__ == "__main__":
    #create_appendix_figs(7)
    create_all_appendix_figs()
