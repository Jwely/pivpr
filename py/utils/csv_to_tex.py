__author__ = 'Jwely'

import pandas as pd
import os


def csv_to_tex(csv_path, caption, justification, horizontal_line_rows=None):
    """
    Quick custom function to write Tex format tables from csv's, with just the most
    critical customization inputs.

    :param csv_path: filepath to csv to convert
    :param caption: caption to place in table
    :param justification: Tex style justification argument to show column spacing and lines
            example: "|ccc|" produces a three column table with all columns centered
    :param horizontal_line_rows: integer list of rows below which to draw horizontal line
            example: [1,2] draws a horizontal line below first row (headers) and second (units)
    """

    if isinstance(horizontal_line_rows, int):
        horizontal_line_rows = [horizontal_line_rows]

    print("loading from {0}".format(csv_path))
    df = pd.read_csv(csv_path, header=None)
    df = df.dropna(axis=1, how='all')           # drops last column of empty caused by trailing commas

    # write code at the top of the file
    texdata = ["\\renewcommand\\baselinestretch{1.3}\\selectfont",
               "\\begin{table}[H]",
               "\\begin{center}",
               "\\begin{tabular}{%s}" % justification,
               "\t\\hline"]

    if horizontal_line_rows is None:
        horizontal_line_rows = []

    for i, row in enumerate(df.as_matrix()):
        row = [element if element != "nan" else " " for element in map(str, row)]
        if any([i == ix for ix in horizontal_line_rows]):
            texdata.append("\t\\hline")
        texdata.append("\t" + " & ".join(row) + "\\\\")


    texdata += ["\t\\hline",
                "\\end{tabular}",
                "\\caption{%s}" % caption,
                "\\label{table:%s}" % os.path.basename(csv_path).replace(".csv", ""),
                "\\end{center}",
                "\\end{table}",
                "\\renewcommand\\baselinestretch{2}\\selectfont"]


    outpath = csv_path.replace(".csv", ".tex")
    with open(outpath, 'w+') as f:
        for line in texdata:
            print(line)
            f.write(line + "\n")

if __name__ == "__main__":

    csv_to_tex("../tables/test_matrix_table.csv",
               "|ccc||ccc||ccc|",
               "Table of experimental conditions for all 70 experiments",
               [2])



