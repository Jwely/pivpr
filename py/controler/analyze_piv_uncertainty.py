__author__ = 'Jwely'

import os
from py.uncertainty import ArtificialVecField
from py.tex import TeXWriter
from py.config import *


def analyze_piv_uncertainty():

    for fold_name in os.listdir(SYNTHESIZED_PIV_DIR):
        v3d_path = os.path.join(SYNTHESIZED_PIV_DIR, fold_name, "{0}.v3d".format(fold_name))
        param_path = os.path.join(SYNTHESIZED_PIV_DIR, fold_name, "params.json")
        tex_path = os.path.join(TEX_FIGURE_DIR, "tex", "{0}.tex".format(fold_name))

        avf = ArtificialVecField(v3d_path, param_path)

        texwriter = TeXWriter(TEX_MAIN_PATH, tex_path)
        # figure output directory
        fig_dir = os.path.join(TEX_FIGURE_DIR, fold_name)
        if not os.path.exists(fig_dir):
            os.mkdir(fig_dir)

        # for each component
        for component in ['U', 'V', 'W']:
            fig_name = "uncertainty_{0}_{1}.png".format(fold_name, component)
            fig_path = os.path.join(fig_dir, fig_name)

            results = avf.plot_histogram(component, title="", outpath=fig_path)
            caption_fmt = "Histogram of error in ${c}$ measurements at station {s}. " \
                          "Simulated conditions $u={u} m/s$, $v={v} m/s$, $w={w} m/s$, $dt={dt} \mu s$. " \
                          "Mean bias=${mb}%$ precision=$\pm {p}%$"
            caption = caption_fmt.format(c=component.lower(), s=fold_name.split("-")[1],
                                         u=avf.piv_params['u'],
                                         v=avf.piv_params['v'],
                                         w=avf.piv_params['w'],
                                         dt=avf.piv_params['dt'],
                                         mb=results['mean_er'],
                                         p=results['max_er'])
            print caption
            #texwriter.add_figure(fig_path, caption, '5in')


# test area
if __name__ == "__main__":
    analyze_piv_uncertainty()

