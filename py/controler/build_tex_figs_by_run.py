__author__ = 'Jwely'

from matplotlib import cm
from py.tex import TeXRunFigurePage
from py.piv import shorthand_to_tex as stt
from py.config import *



def build_tex_figs_by_run(run_id):
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
                              "appendix_run_id_{0}".format(run_id),
                              run_id,
                              force_recalc=False)
    av = figdoc.axial_vortex

    # populate the figdoc with content
    figdoc.add_text("\subsection{{Run ID {0}}}".format(run_id))

    # add contour plots
    for component in ['R', 'T', 'W', 'rt', 'rw', 'tw', 'rr', 'tt', 'ww', 'ctke']:

        caption_fmt = "Contour plot of {0} for run {1} at $z/c$={2}, $V_{{free}}$={3}"
        caption = caption_fmt.format(stt(component), run_id, av.z_location, av.velocity_fs)
        figdoc.add_contour_plot(component, caption, contour_width)

    # add scatter plots
    scatter_kwargs = {"t_range": (10, 80),
                      "r_range": (0, 100),
                      "symmetric": True}
    #t_kwargs = scatter_kwargs.update()
    #figdoc.add_scatter_plot('r_mesh', 'T', t_kwargs)

    figdoc.write()


if __name__ == "__main__":
    run_ids = range(1,71)
    for run_id in run_ids:
        build_tex_figs_by_run(run_id)
