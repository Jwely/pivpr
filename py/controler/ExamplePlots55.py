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
    
        av.contour_plot('num', titles="Number of samples",
                        outpath=os.path.join(self.fig_dir, "example_num_contour.png"))
    
        av.contour_plot('U', titles="Stable Horizontal Velocity (U)",
                        outpath=os.path.join(self.fig_dir, "example_U_contour.png"))
    
        av.contour_plot('V', titles="Stable Vertical Velocity (V)",
                        outpath=os.path.join(self.fig_dir, "example_V_contour.png"))
    
        av.contour_plot('T', titles="Stable Tangential Velocity (T)",
                        outpath=os.path.join(self.fig_dir, "example_T_contour.png"))
    
        av.contour_plot('R', titles="Stable Radial Velocity (R)",
                        outpath=os.path.join(self.fig_dir, "example_R_contour.png"))
    
        av.contour_plot('W', titles="Stable Axial Velocity (W)",
                        outpath=os.path.join(self.fig_dir, "example_W_contour.png"))
    
        av.contour_plot('rr', titles="Radial Fluctuation Squared (rr)",
                        outpath=os.path.join(self.fig_dir, "example_rr_contour.png"))
    
        av.contour_plot('tt', titles="Tangential Fluctuation Squared (tt)",
                        outpath=os.path.join(self.fig_dir, "example_tt_contour.png"))
    
        av.contour_plot('ww', titles="Axial Fluctuation Squared (ww)",
                        outpath=os.path.join(self.fig_dir, "example_ww_contour.png"))
    
        av.contour_plot('rt', titles="Radial/Tangential Reynolds Stress (rt)",
                        outpath=os.path.join(self.fig_dir, "example_rt_contour.png"))
    
        av.contour_plot('rw', titles="Radial/Axial Reynolds Stress (rw)",
                        outpath=os.path.join(self.fig_dir, "example_rw_contour.png"))
    
        av.contour_plot('tw', titles="Tangential/Axial Reynolds Stress (tw)",
                        outpath=os.path.join(self.fig_dir, "example_tw_contour.png"))
    
        av.contour_plot('ctke', titles="Turbulent Kinetic Energy (k)",
                        outpath=os.path.join(self.fig_dir, "example_ctke_contour.png"))
    
        av.scatter_plot('r_mesh', 'T', 'ctke', cmap=cm.jet,
                        title="Tangential Velocity and Turbulent Kinetic Energy",
                        x_label="Distance to vortex core (mm)",
                        y_label="Tangential velocity (m/s)",
                        c_label="Turbulent kinetic energy (TKE)",
                        xrange=(0, 60),
                        yrange=(0, 6.5),
                        tight=True,
                        outpath=os.path.join(self.fig_dir, "example_TscatterTKE.png"))
