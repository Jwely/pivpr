__author__ = 'Jwely'

import pip
import os


def install_dependencies(wheel_dir):
    """
    Installs all the dependencies required to use pivpr/py. Fetch them from the assets release
    by downloading and extracting `python27_64bit_resources.7z` and pointing this function to
    the extracted directory. You can go get all these yourself, but this is probably easier.

    :param wheel_dir:   directory with python resource files in it
    """

    # install microsoft visual studio and graphviz msi files. If they are already installed
    # you may hit cancel without any negative effects.
    os.system(os.path.join(wheel_dir, r"VCForPython27.msi"))
    os.system(os.path.join(wheel_dir, r"graphviz-2.38.msi"))

    # pip installation thingies
    pip.main(["install", "--upgrade", "pip"])
    pip.main(["install", os.path.join(wheel_dir, "numpy-1.9.3+mkl-cp27-none-win_amd64.whl")])
    pip.main(["install", "pandas"])
    pip.main(["install", "matplotlib"])
    pip.main(["install", os.path.join(wheel_dir, "pygraphviz-1.3.1-cp27-none-win_amd64.whl")])


if __name__ == "__main__":
    install_dependencies(r"C:\Users\Jeff\Downloads\python27_64bit_resources")



