__author__ = 'Jwely'

import os
import pandas as pd
from py.utils import csv_to_tex
from py.uncertainty import ArtificialVecField
from py.tex import TeXWriter
from py.config import *


def calculate_uncertainty(name, n_measurements=1):
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
        fig_name = "uncertainty_{0}_{1}.png".format(name, component)
        fig_path = os.path.join(fig_dir, fig_name)

        results = avf.plot_histogram(component, n_measurements, title="", outpath=fig_path)
        caption_fmt = r"Histogram of ${c}$ measurements at station {s}. " \
                      r"Simulated conditions $(u,v,w)=({u}, {v}, {w})$, $dt={dt} \mu s$, " \
                      r"$\beta_{c}=\pm {mb:1.4f}$, $P_{c}=\pm {p:1.4f}$, $U_{c}=\pm {uncert:1.4f}$"
        caption = caption_fmt.format(c=component.upper(), s=name[12],
                                     u=avf.piv_params['u'],
                                     v=avf.piv_params['v'],
                                     w=avf.piv_params['w'],
                                     dt=avf.piv_params['dt'],
                                     mb=results['bias'],
                                     p=results['stdev'],
                                     uncert=results['uncertainty'])
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
            return "{0:2.4f}".format(myfloat)
        else:
            return myfloat

    # empty lists
    u_list = []
    w_list = []
    v_list = []

    # the order of columns to output to dict with pandas dataframe.to_csv
    order = ["Station",
             "$dt$",
             "$U_{sim}$",
             "$V_{sim}$",
             "$W_{sim}$"]

    for station in stations:

        names = ["Ely_May28th{0}{1}".format(str(station).zfill(2), str(condition).zfill(3)) for condition in conditions]

        for name in names:
            station = name[12]
            condition = name[15]
            avf = calculate_uncertainty(name, n_measurements=1)

            u_entry = {"Station": station,
                       "$dt$": avf.piv_params['dt'],
                       "$U_{sim}$": avf.piv_params['u'],
                       "$V_{sim}$": avf.piv_params['v'],
                       "$W_{sim}$": avf.piv_params['w'],
                       "$\\bar{U}$": ff(avf.error_data['U']['mean']),
                       "$\\beta_U$": ff(avf.error_data['U']['bias']),
                       "$P_U$": ff(avf.error_data['U']['precision']),
                       "$U_U$": ff(avf.error_data['U']['uncertainty'])}

            v_entry = {"Station": station,
                       "$dt$": avf.piv_params['dt'],
                       "$U_{sim}$": avf.piv_params['u'],
                       "$V_{sim}$": avf.piv_params['v'],
                       "$W_{sim}$": avf.piv_params['w'],
                       "$\\bar{V}$": ff(avf.error_data['V']['mean']),
                       "$\\beta_V$": ff(avf.error_data['V']['bias']),
                       "$P_V$": ff(avf.error_data['V']['precision']),
                       "$U_V$": ff(avf.error_data['V']['uncertainty'])}

            w_entry = {"Station": station,
                       "$dt$": avf.piv_params['dt'],
                       "$U_{sim}$": avf.piv_params['u'],
                       "$V_{sim}$": avf.piv_params['v'],
                       "$W_{sim}$": avf.piv_params['w'],
                       "$\\bar{W}$": ff(avf.error_data['W']['mean']),
                       "$\\beta_W$": ff(avf.error_data['W']['bias']),
                       "$P_W$": ff(avf.error_data['W']['precision']),
                       "$U_W$": ff(avf.error_data['W']['uncertainty'])}

            u_list.append(u_entry)
            v_list.append(v_entry)
            w_list.append(w_entry)

    # now write these lists of dictionaries to a pandas dataframe and csv file
    u_file = os.path.join(TEX_TABLE_DIR, "{0}_u.csv".format(csv_name))
    udf = pd.DataFrame(u_list)
    udf = udf[order + ["$\\bar{U}$", "$\\beta_U$", "$P_U$", "$U_U$"]]
    udf.to_csv(u_file, index_label=" ")

    v_file = os.path.join(TEX_TABLE_DIR, "{0}_v.csv".format(csv_name))
    vdf = pd.DataFrame(v_list)
    vdf = vdf[order + ["$\\bar{V}$", "$\\beta_V$", "$P_V$", "$U_V$"]]
    vdf.to_csv(v_file, index_label=" ")

    w_file = os.path.join(TEX_TABLE_DIR, "{0}_w.csv".format(csv_name))
    wdf = pd.DataFrame(w_list)
    wdf = wdf[order + ["$\\bar{W}$", "$\\beta_W$", "$P_W$", "$U_W$"]]
    wdf.to_csv(w_file, index_label=" ")

    if verbose:
        print udf
        print vdf
        print wdf

    return u_file, v_file, w_file



def perform_uncertainty_analysis():
    """
    Final setp in producing uncertainty tables and figures.

    :return:
    """

    # create comprehensive uncertainty table
    u_path, v_path, w_path = make_csv_uncertainty_tables([1, 2, 3, 4, 5, 6, 7], [1, 2, 3, 4], "uncertainties")
    csv_to_tex(u_path,
               caption="Uncertainty in $X$ direction velocity measurements. Unlabelled units are $m/s$.",
               justification="|cccccccccc|",
               horizontal_line_rows=[1])
    csv_to_tex(v_path,
               caption="Uncertainty in $Y$ direction velocity measurements. Unlabelled units are $m/s$.",
               justification="|cccccccccc|",
               horizontal_line_rows=[1])
    csv_to_tex(w_path,
               caption="Uncertainty in $Z$ direction velocity measurements. Unlabelled units are $m/s$.",
               justification="|cccccccccc|",
               horizontal_line_rows=[1])

# test area
if __name__ == "__main__":
    perform_uncertainty_analysis()
