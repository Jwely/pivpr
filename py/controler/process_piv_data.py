__author__ = 'Jwely'

from matplotlib import cm
from py.tex import TeXRunFigurePage
from py.utils import shorthand_to_tex as stt
from py.config import *



def process_piv_data(run_id):
    """
    Processes all vortex PIV data. At this point, everything has become a little convoluted, but running this
    function will complete all processing for input run_id and produce all the relevant figures which embody the
    dataset and and place them in an appendix `.tex` file ready for inclusion in the master thesis document main.

    The hierarchy at this point is:
        TeXRunFigurePage extends
        TexWriter invokes
        TexFigureGenerator wraps
        construct_experiments returns
        Experiment with attribute
        AxialVortex which has plotting methods to produce figures of interest.
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

    figdoc = TeXRunFigurePage(TEX_MAIN_PATH,
                              "appendix_run_id_{0}".format(run_id),     # note the color scheme tag
                              run_id,
                              force_recalc=False)

    # populate the figdoc with content
    figdoc.add_text("\subsection{{Run ID {0}}}".format(run_id))
    #figdoc.add_quiver_plot("Quiver plot of run ID {0}.".format(run_id), quiver_width)
    figdoc.add_stream_plot("Stream plot of run ID {0}.".format(run_id), stream_width)
    contour_component_plotter(['num', 'ctke'])

    contour_component_plotter(['R', 'T', 'W', 'rt', 'rw', 'tw', 'rr', 'tt', 'ww'])


    figdoc.write()


def process_all_piv_data():
    run_ids = range(1,71)
    for run_id in run_ids:
        process_piv_data(run_id)



if __name__ == "__main__":
    process_all_piv_data()
