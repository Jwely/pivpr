__author__ = 'Jwely'

import pip
import os

# just change this directory to the local filepath of the resource folder
wheel_dir = r"C:\Users\Jeff\Downloads\python27_64bit_resources"

os.system(os.path.join(wheel_dir, r"VCForPython27.msi"))
pip.main(["install", "--upgrade", "pip"])
pip.main(["install", "pandas"])
pip.main(["install", "openpiv"])



