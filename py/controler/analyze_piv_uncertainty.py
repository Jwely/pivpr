__author__ = 'Jwely'

import os
import pandas as pd
from py.utils import csv_to_tex
from py.uncertainty import ArtificialVecField
from py.tex import TeXWriter
from py.config import *


def calculate_uncertainty(name, n_measurements=200):
    """
    Finds the files in the artificial_images folder of the uncertainty package and calculates
    uncertainty values for each image set. The controller function `synthesize_piv_uncertainty_images` must
    have already been run, and the resulting images must have been processed using PIV software with
    resulting v3d files now present.

    :param name: the names of the v3d files to process such as 'Ely_May28th01001'
    :param n_measurements:  the number of measurements should be 1, or 200ish
    :return:
    """

    v3d_path = os.path.join(SYNTHESIZED_PIV_DIR, "{0}.v3d".format(name))
    param_path = os.path.join(SYNTHESIZED_PIV_DIR, "{0}.json".format(name))
    tex_path = os.path.join(TEX_FIGURE_DIR, "tex", "{0}.tex".format(name))

    # throw error if required files are not present
    if not os.path.exists(v3d_path):
        raise Exception("Missing {0}".format(v3d_path))

    if not os.path.exists(param_path):
        raise Exception("Missing {0}".format(param_path))

    avf = ArtificialVecField(v3d_path, param_path)

    # initalize the tex writing requirements
    texwriter = TeXWriter(TEX_MAIN_PATH, tex_path)
    fig_dir = os.path.join(TEX_FIGURE_DIR, name)
    if not os.path.exists(fig_dir):
        os.mkdir(fig_dir)

    # for each component
    for component in ['U', 'V', 'W']:
        fig_name = "uncertainty_{0}_{1}.jpg".format(name, component)
        fig_path = os.path.join(fig_dir, fig_name)

        results = avf.plot_histogram(component, n_measurements, title="", outpath=fig_path)
                      # r"$\beta_{c}=\pm {mb:1.4f}$, $P_1{c}=\pm {p1:1.4f}$"\
        caption_fmt = r"Histogram of ${c}$ measurements at station {s}. " \
                      r"Simulated conditions $(u,v,w)=({u}, {v}, {w})$, $dt={dt} \mu s$, " \
                      r"$U_1_{c}=\pm {u1:1.3f}$, $U_200_{c}=\pm {u200:1.3f}$"
        caption = caption_fmt.format(c=component.upper(), s=name[12],
                                     u=avf.piv_params['u'],
                                     v=avf.piv_params['v'],
                                     w=avf.piv_params['w'],
                                     dt=avf.piv_params['dt'],
                                     mb=results['bias'],
                                     p1=results['precision_1'],
                                     u1=results['uncertainty_1'],
                                     p200=results['precision_n'],
                                     u200=results['uncertainty_n']
                                     )
        print caption
        texwriter.add_figure(fig_path, caption, '6in')

    texwriter.write()
    return avf


def make_csv_uncertainty_tables(stations, conditions, csv_name, verbose=False):

    if not isinstance(stations, list):
        stations = [stations]

    if not isinstance(conditions, list):
        conditions = [conditions]

    def ff(myfloat):
        if isinstance(myfloat, float):
            return "{0:2.3f}".format(myfloat)
        else:
            return myfloat

    # empty lists
    u_list = []
    w_list = []
    v_list = []

    # the order of columns to output to dict with pandas dataframe.to_csv
    order = ["Station",
             "$dt$",
             "$u_{{sim}}$",
             "$v_{{sim}}$",
             "$w_{{sim}}$",
             ]

    for station in stations:

        names = ["Ely_May28th{0}{1}".format(str(station).zfill(2), str(condition).zfill(3))
                 for condition in conditions]

        for name in names:
            station = name[12]
            condition = name[15]
            avf = calculate_uncertainty(name)

            u_entry = {"Station": station,
                       "$dt$": avf.piv_params['dt'],
                       "$u_{sim}$": avf.piv_params['u'],
                       "$v_{sim}$": avf.piv_params['v'],
                       "$w_{sim}$": avf.piv_params['w'],
                       "$\\bar{u}$": ff(avf.error_data['U']['mean']),
                       "$\\beta_u$": ff(avf.error_data['U']['bias']),
                       "$P_{u^{\prime}}$": ff(avf.error_data['U']['precision_1']),
                       "$U_{u^{\prime}}$": ff(avf.error_data['U']['uncertainty_1']),
                       "$P_{\\bar{u}}$": ff(avf.error_data['U']['precision_n']),
                       "$U_{\\bar{u}}$": ff(avf.error_data['U']['uncertainty_n']),
                       }

            v_entry = {"Station": station,
                       "$dt$": avf.piv_params['dt'],
                       "$u_{sim}$": avf.piv_params['u'],
                       "$v_{sim}$": avf.piv_params['v'],
                       "$w_{sim}$": avf.piv_params['w'],
                       "$\\bar{v}$": ff(avf.error_data['V']['mean']),
                       "$\\beta_v$": ff(avf.error_data['V']['bias']),
                       "$P_{v^{\prime}}$": ff(avf.error_data['V']['precision_1']),
                       "$U_{v^{\prime}}$": ff(avf.error_data['V']['uncertainty_1']),
                       "$P_{\\bar{v}}$": ff(avf.error_data['V']['precision_n']),
                       "$U_{\\bar{v}}$": ff(avf.error_data['V']['uncertainty_n']),
                       }
            w_entry = {"Station": station,
                       "$dt$": avf.piv_params['dt'],
                       "$u_{sim}$": avf.piv_params['u'],
                       "$v_{sim}$": avf.piv_params['v'],
                       "$w_{sim}$": avf.piv_params['w'],
                       "$\\bar{w}$": ff(avf.error_data['W']['mean']),
                       "$\\beta_w$": ff(avf.error_data['W']['bias']),
                       "$P_{w^{\prime}}$": ff(avf.error_data['W']['precision_1']),
                       "$U_{w^{\prime}}$": ff(avf.error_data['W']['uncertainty_1']),
                       "$P_{\\bar{w}}$": ff(avf.error_data['W']['precision_n']),
                       "$U_{\\bar{w}}$": ff(avf.error_data['W']['uncertainty_n']),
                       }

            u_list.append(u_entry)
            v_list.append(v_entry)
            w_list.append(w_entry)

    # now write these lists of dictionaries to a pandas dataframe and csv file
    returnlist = []
    for comp in ["u", "v", "w"]:
        fpath = os.path.join(TEX_TABLE_DIR, "{0}_{1}.csv".format(csv_name, comp))
        df = pd.DataFrame(eval("{0}_list".format(comp)))
        # the bracketing here has gotten really tricky
        headers = order + ["$\\bar{{{c}}}$", "$\\beta_{c}$", "$P_{{{c}^{{\prime}}}}$",
                           "$P_{{\\bar{{{c}}}}}$", "$U_{{{c}^{{\prime}}}}$", "$U_{{\\bar{{{c}}}}}$"]
        headers = [head.format(c=comp) for head in headers]
        df = df[headers]
        df.to_csv(fpath, index=False, float_format='%.2f')
        returnlist.append(fpath)

        if verbose:
            print df

    return tuple(returnlist)


def perform_uncertainty_analysis():
    """
    Final setp in producing uncertainty tables and figures.

    :return:
    """

    # create comprehensive uncertainty table
    u_path, v_path, w_path = make_csv_uncertainty_tables([1, 2, 3, 4, 5, 6, 7], [1, 2], "uncertainties")
    csv_to_tex(u_path,
               caption="Uncertainty in $X$ direction velocity measurements. Unlabelled units are $m/s$.",
               horizontal_line_rows=[1])
    csv_to_tex(v_path,
               caption="Uncertainty in $Y$ direction velocity measurements. Unlabelled units are $m/s$.",
               horizontal_line_rows=[1])
    csv_to_tex(w_path,
               caption="Uncertainty in $Z$ direction velocity measurements. Unlabelled units are $m/s$.",
               horizontal_line_rows=[1])

# test area
if __name__ == "__main__":
    perform_uncertainty_analysis()
