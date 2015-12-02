__author__ = 'Jwely'

import pip
import os


def main():
    # just change this directory to the local filepath of the resource folder
    wheel_dir = r"C:\Users\Jeff\Downloads\python27_64bit_resources"

    # install microsoft visual studio. You may hit cancel without error if its already installed.
    os.system(os.path.join(wheel_dir, r"VCForPython27.msi"))

    # pip installation thingies
    pip.main(["install", "--upgrade", "pip"])
    pip.main(["install", os.path.join(wheel_dir, "numpy-1.9.3+mkl-cp27-none-win_amd64.whl")])
    pip.main(["install", "pandas"])
    pip.main(["install", "openpiv"])
    pip.main(["install", "plotly"])



