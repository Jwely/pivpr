__author__ = 'Jwely'

from py.tex import TeXRunFigurePage
from py.piv import shorthand_to_tex as stt
from py.config import *
from py.utils import merge_dicts



def build_tex_figs_by_run(run_id, force_recalc=False):
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

    contour_width = '5in'
    scatter_width = '7in'
    quiver_width = '5in'
    stream_width = '5in'

    figdoc = TeXRunFigurePage(TEX_MAIN_PATH,
                              "appendix_run_id_{0}".format(run_id),
                              run_id,
                              force_recalc=force_recalc)
    av = figdoc.axial_vortex

    # populate the figdoc with content
    figdoc.add_text("\subsection{{Run ID {0}}}".format(run_id))

    station_id = ((run_id + 9) % 10)
    z_location = av.z_location / 101.6

    # add stream plot
    caption = "Stream plot at $z/c$={0}, $V_{{free}}$={1}, station {2}.".format(z_location, av.velocity_fs, station_id)
    figdoc.add_stream_plot(caption, stream_width)

    # full contour plots with no radial or angular subseting
    for component in ['R', 'T', 'W', 'rt', 'rw', 'tw', 'rr', 'tt', 'ww', 'ctke', 'num']:
        contour_kwargs = {"r_range": ('0r', '6r')}
        caption_fmt = "Contour plot of {0} at $z/c$={2}, $V_{{free}}$={3}, station {1}."
        caption = caption_fmt.format(stt(component), station_id, z_location, av.velocity_fs)
        figdoc.add_contour_plot(component, caption, contour_width, create_kwargs=contour_kwargs)


    # radius scatter plots with kwargs
    scatter_kwargs = {"x_range": (0, 4)}

    t_kwargs = merge_dicts(scatter_kwargs, {"title": "Azimuthal Velocity vs Radius"})
    caption = "Scatter plot of azimuthal velocity vs radius at $z/c$={0}, $V_{{free}}$={1}, station{2}".format(
        z_location, av.velocity_fs, station_id)
    figdoc.add_scatter_plot('r_mesh', 'T', caption, scatter_width, create_kwargs=t_kwargs)

    k_kwargs = merge_dicts(scatter_kwargs, {"title": "Turbulent Kinetic Energy"})
    caption = "Scatter plot of turbulent kinetic energy vs radius at $z/c$={0}, $V_{{free}}$={1}, station{2}".format(
        z_location, av.velocity_fs, station_id)
    figdoc.add_scatter_plot('r_mesh', 'ctke', caption, scatter_width, create_kwargs=k_kwargs)


    # logarithmic plots wtih kwqargs
    log_kwargs = merge_dicts(scatter_kwargs, {"log_y": True,
                                              "y_range": (1e3, 1e9),
                                              "y_label": " ",
                                              "t_range": (10, 80),
                                              "symmetric": True})

    kwargs = merge_dicts(log_kwargs, {"title": r"$\frac{1}{r^2}\frac{d}{dr}[r^2 \overline{t^\prime r^\prime}]$"})
    caption = "Scatter plot of reynolds stress term vs radius at $z/c$={0}, $V_{{free}}$={1}, station{2}".format(
        z_location, av.velocity_fs, station_id)
    figdoc.add_scatter_plot('r_mesh', 'turb_visc_ettap_top', caption, scatter_width, create_kwargs=kwargs)

    kwargs = merge_dicts(log_kwargs, {"title": r"$\frac{d^2\bar{t}}{dr^2} + \frac{d}{dr}\Big(\frac{\bar{t}}{r}\Big)$"})
    caption = "Scatter plot of velocity gradient term vs radius at $z/c$={0}, $V_{{free}}$={1}, station{2}".format(
        z_location, av.velocity_fs, station_id)
    figdoc.add_scatter_plot('r_mesh', 'turb_visc_ettap_bot', caption, scatter_width, create_kwargs=kwargs)



    figdoc.write()


if __name__ == "__main__":
    run_ids = range(1, 71)
    for run_id in run_ids:
        build_tex_figs_by_run(run_id, force_recalc=True)
