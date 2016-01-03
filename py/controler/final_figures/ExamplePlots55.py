__author__ = 'Jwely'

from matplotlib import cm
import os
from FinalFigure import FinalFigure


class ExamplePlots55(FinalFigure):
    """
    Builds plots for use in example section of final thesis document
    """

    def __init__(self):
        # invoke parent class init
        FinalFigure.__init__(self, experiment_ids=[55], sub_dir="example_vortex_figs")


    def generate(self):
        av = self.experiments[0].axial_vortex

        av.quiver_plot(title="Quiver Plot",
                       outpath=os.path.join(self.fig_dir, "example_quiver.png"))

        av.stream_plot(title="Stream Plot Colored by In-plane Velocities",
                       outpath=os.path.join(self.fig_dir, "example_stream.png"))
    
        av.contour_plot('num', title="Number of samples ($N$)",
                        outpath=os.path.join(self.fig_dir, "example_num_contour.png"))
    
        av.contour_plot('U', title=r"$\overline{u}$",
                        outpath=os.path.join(self.fig_dir, "example_U_contour.png"))
    
        av.contour_plot('V', title=r"$\overline{v}$",
                        outpath=os.path.join(self.fig_dir, "example_V_contour.png"))
    
        av.contour_plot('T', title=r"$\overline{t}$",
                        outpath=os.path.join(self.fig_dir, "example_T_contour.png"))
    
        av.contour_plot('R', title=r"$\overline{r}$",
                        outpath=os.path.join(self.fig_dir, "example_R_contour.png"))
    
        av.contour_plot('W', title=r"$\overline{w}$",
                        outpath=os.path.join(self.fig_dir, "example_W_contour.png"))
    
        av.contour_plot('rr', title=r"$\overline{r^\prime r^\prime}$",
                        outpath=os.path.join(self.fig_dir, "example_rr_contour.png"))
    
        av.contour_plot('tt', title=r"$\overline{t^\prime t^\prime}$",
                        outpath=os.path.join(self.fig_dir, "example_tt_contour.png"))
    
        av.contour_plot('ww', title=r"$\overline{w^\prime w^\prime}$",
                        outpath=os.path.join(self.fig_dir, "example_ww_contour.png"))
    
        av.contour_plot('rt', title=r"$\overline{r^\prime t^\prime}$",
                        outpath=os.path.join(self.fig_dir, "example_rt_contour.png"))
    
        av.contour_plot('rw', title=r"$\overline{r^\prime w^\prime}$",
                        outpath=os.path.join(self.fig_dir, "example_rw_contour.png"))
    
        av.contour_plot('tw', title=r"$\overline{t^\prime w^\prime}$",
                        outpath=os.path.join(self.fig_dir, "example_tw_contour.png"))
    
        av.contour_plot('ctke', title=r"Turbulent Kinetic Energy ($k$)",
                        outpath=os.path.join(self.fig_dir, "example_ctke_contour.png"))
    
        av.scatter_plot('r_mesh', 'T', 'ctke', cmap=cm.jet,
                        title="Tangential Velocity and Turbulent Kinetic Energy",
                        x_label="Distance to vortex core (mm)",
                        y_label="Tangential velocity (m/s)",
                        c_label="Turbulent kinetic energy (TKE)",
                        x_range=(0, 60),
                        y_range=(0, 6.5),
                        tight=True,
                        outpath=os.path.join(self.fig_dir, "example_TscatterTKE.png"))
