import sys
import subprocess

def install_package(package):
    subprocess.check_call([sys.executable, "-m","pip","install",package])

dependencies = ["numpy","pandas","statsmodels","scipy","matplotlib","DynaTMT-py"]

for i in dependencies:
    install_package(i)
print('Finished')