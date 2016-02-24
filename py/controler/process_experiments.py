__author__ = "Jwely"

import os
import json
import pandas as pd
from py.tex import csv_to_tex
from py.piv import construct_experiments
from py.config import *


def process_experiments(outname, table_caption, ids=None, min_points=DEFAULT_MIN_POINTS,
                        include_dynamic=False, force_recalc=False):

    if isinstance(ids, int):
        ids = [ids]
    elif ids is None:
        ids = range(1, 71)

    # build all the json files.
    for i in ids:
        construct_experiments(EXPERIMENT_TABLE_PATH, DATA_FULL_DIR,
                              i, min_points, include_dynamic, force_recalc)

    # now list all the json files
    table = {}
    json_paths = [os.path.join(PIV_PICKLE_DIR, fname) for fname in os.listdir(PIV_PICKLE_DIR) if "json" in fname]
    for j in json_paths:
        if any(["ID-{0}_Z".format(i) in j for i in ids]):
            atts = json.loads(open(j).read())
            table[atts["experiment_id"]] = atts

    # build a dataframe, rename the columns to tex formatting
    df = pd.DataFrame(table.values())
    df['pres_atm'] /= 1000
    df['eta_p'] = df['eta_p'].map(lambda x: "%.2f" % x)     # increase precision on eta_p args

    kw_names = ["experiment_id", "z_location_mm", "v_nominal", "dt", "velocity_free_stream", "pres_atm",
                "temp_tunnel", "rel_humid", "eta_p", "r_mesh_core", "T_max", "W_mean"]

    units = ["", "$mm$", "$m/s$", "$\\mu s$", "$m/s$", "$KPa$", "K", "$\\%$", "$\\mu s$", "$mm$", "$m/s$", "$m/s$"]

    tex_names = ["Run", "$I_Z$", "$V_{nom}$", "$dt$", "$V_{free}$", "$P_{atm}$", "$T_{tunnel}$",
                 "$\\phi$", "$\\eta_P$", "$R_{core}$", "$\\overline{v_{\\theta}}_{max}$", "$\\overline{v_{\\bar{z}}}$"]

    df = df[kw_names]
    df.columns = pd.MultiIndex.from_tuples(zip(tex_names, units))
    result_path = os.path.join(TEX_TABLE_DIR, outname)
    if os.path.exists(result_path):
        os.remove(result_path)

    df.sort(columns="Run", inplace=True)  # ensures they are in ascending order, which they aren't somehow
    df.to_csv(result_path, index=False, float_format='%.1f')

    # push experiment results summary table to TEX
    csv_to_tex(result_path,
               caption=table_caption,
               horizontal_line_rows=[2])


if __name__ == "__main__":
    # lets make tables split into groups by station

    for i in [1, 2, 3, 4, 5, 6, 7]:
        process_experiments(ids=range(10 * (i - 1) + 1, 10 * i + 1),
                            outname="experiment_results_{i}.csv".format(i=i),
                            table_caption="Summary of experimental measurements at station {i}".format(i=i),
                            include_dynamic=False, force_recalc=False)
