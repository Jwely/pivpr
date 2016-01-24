__author__ = "Jwely"

import os
import json
import pandas as pd
from py.piv import construct_experiments
from py.config import *


def process_experiments(ids=None, min_points=DEFAULT_MIN_POINTS,
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
        atts = json.loads(open(j).read())
        table[atts["experiment_id"]] = atts

    # build a dataframe, rename the columns to tex formatting
    df = pd.DataFrame(table.values())

    kw_names = ["experiment_id", "z_location", "v_nominal", "dt", "velocity_free_stream", "q", "pres_atm",
                "temp_tunnel", "rel_humid", "r_mesh_core", "T_max", "W_core"]
    tex_names = ["Run", "$I_Z$", "$V_nom$", "$dt$", "$V_{fs}$", "$Q$", "$P_{atm}$",
                 "$T_{tunnel}$", "$\phi$", "$R_{core}$", "$\\overline{t}_{max}$", "$\\overline{w}_{core}$"]
    df = df[kw_names]
    df.columns = tex_names
    df.to_csv("results_table.csv", index=False)




if __name__ == "__main__":
    process_experiments(include_dynamic=False, force_recalc=False)
