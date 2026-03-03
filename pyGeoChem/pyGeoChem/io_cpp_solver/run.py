# %%
import subprocess
import pathlib as pt


def run_simulator(input_file, windows=False, dir=''):
    '''
    input_file: GeoChemX input file
    windows: True/False use Windows or linux version
    dir: relative path to Geochem repo folder
    '''
    if windows:
        bin = 'build_win/bin/Debug/GeoChemX.exe'
    else:
        bin = 'build/bin/GeoChemX'
    p = pt.Path(dir + bin)
    if not p.exists():
        raise ValueError('Could not find executable: ' + str(p))
    cmd = [str(p) + ' EQSOLVER ' + input_file]
    return subprocess.run(cmd, shell=True)
