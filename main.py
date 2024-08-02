#!/usr/bin/env python3
import os
import shutil
import subprocess
import numpy as np
from time import time
from pathlib import Path
import matplotlib.pyplot as plt
from quasilinear_gs2 import quasilinear_estimate
from simsopt.mhd import Vmec
from simsopt.mhd.vmec_diagnostics import vmec_fieldlines
import matplotlib
import warnings
matplotlib.use('Agg') 
warnings.filterwarnings("ignore",category=matplotlib.MatplotlibDeprecationWarning)
import matplotlib.cbook
import argparse
from configurations import CONFIG
from gs2_functions import getgamma, eigenPlot, replace, gammabyky, output_to_csv, remove_files, to_gs2
from vmecPlot2 import main as vmecPlot2
parser = argparse.ArgumentParser()
parser.add_argument("--type", type=int, default=2)
args = parser.parse_args()
this_path = Path(__file__).parent.resolve()
#########
## To run this file use the following command:
## python3 run_gs2.py --type 2
## where type 1 is HSX, type 2 is W7-X, type 3 is QI, type 4 is QH, type 5 is QA
######## INPUT PARAMETERS ########
run_vmecPlot2 = True
run_nonlinear = False
home_directory = os.path.expanduser("~")
gs2_executable = f'{home_directory}/local/gs2/bin/gs2'
results_folder = 'output'
vmec_equilibria_folder = 'equilibria'
config = CONFIG[args.type]
PARAMS = config['params']
phi_GS2 = np.linspace(-PARAMS['nperiod']*np.pi, PARAMS['nperiod']*np.pi, PARAMS['nphi'])
########################################
# Go into the output directory
OUT_DIR_APPENDIX=config['output_dir']
OUT_DIR = os.path.join(this_path,results_folder,config['output_dir'])
os.makedirs(OUT_DIR, exist_ok=True)
os.chdir(OUT_DIR)
vmec = Vmec(os.path.join(os.path.join(this_path,vmec_equilibria_folder), config['wout']),verbose=False)
output_path_parameters=f"{OUT_DIR_APPENDIX}.csv"
output_csv = os.path.join(OUT_DIR,f'scan_ln_lt_{OUT_DIR_APPENDIX}.csv')

print(f'#######################')
print(f'Configuration: {OUT_DIR_APPENDIX}')
print(f'Output directory: {OUT_DIR}')
print(f'#######################')
print()

if run_vmecPlot2:
    print(f'#######################')
    print('Running vmecPlot2')
    print(f'#######################')
    vmecPlot2(os.path.join(os.path.join(this_path,vmec_equilibria_folder), config['wout']))
    print(f'#######################')
    print()

# Run GS2
gridout_file = os.path.join(OUT_DIR,f'grid_gs2.out')
to_gs2(gridout_file, vmec, PARAMS["s_radius"], PARAMS["alpha_fieldline"], phi1d=phi_GS2, nlambda=PARAMS["nlambda"])
fl1 = vmec_fieldlines(vmec, PARAMS["s_radius"], PARAMS["alpha_fieldline"], phi1d=phi_GS2, plot=True, show=False)
plt.savefig(f'geometry_profiles_s{PARAMS["s_radius"]}_alpha{PARAMS["alpha_fieldline"]}.png');plt.close()
def run_gs2(ln, lt, show_fig=True, save_fig=True):
    start_time_local = time()
    gs2_input_name = f"gs2Input-LN{ln:.1f}-LT{lt:.1f}"
    gs2_input_file = os.path.join(OUT_DIR,f'{gs2_input_name}.in')
    shutil.copy(os.path.join(this_path,'gs2Input-linear.in'),gs2_input_file)
    replace(gs2_input_file,' gridout_file = "grid.out"',f' gridout_file = "grid_gs2.out"')
    replace(gs2_input_file,' nstep = 150',f' nstep = {PARAMS["nstep"]}')
    replace(gs2_input_file,' delt = 0.4 ! Time step',f' delt = {PARAMS["dt"]} ! Time step')
    replace(gs2_input_file,' fprim = 1.0 ! -1/n (dn/drho)',f' fprim = {ln} ! -1/n (dn/drho)')
    replace(gs2_input_file,' tprim = 3.0 ! -1/T (dT/drho)',f' tprim = {lt} ! -1/T (dT/drho)')
    replace(gs2_input_file,' aky_min = 0.4',f' aky_min = {PARAMS["aky_min"]}')
    replace(gs2_input_file,' aky_max = 5.0',f' aky_max = {PARAMS["aky_max"]}')
    replace(gs2_input_file,' naky = 4',f' naky = {PARAMS["naky"]}')
    replace(gs2_input_file,' vnewk = 0.01 ! collisionality parameter',f' vnewk = {PARAMS["vnewk"]} ! collisionality parameter')
    replace(gs2_input_file,' ngauss = 3 ! Number of untrapped pitch-angles moving in one direction along field line.',
         f' ngauss = {PARAMS["ngauss"]} ! Number of untrapped pitch-angles moving in one direction along field line.')
    replace(gs2_input_file,' negrid = 10 ! Total number of energy grid points',
          f' negrid = {PARAMS["negrid"]} ! Total number of energy grid points')
    if run_nonlinear:
        replace(gs2_input_file,' nonlinear_mode = "off" ! Include nonlinear terms? ("on","off")',
                              f' nonlinear_mode = "on"  ! Include nonlinear terms? ("on","off")')
        replace(gs2_input_file, 'grid_option = "range" ! The general layout of the perpendicular grid.',
                                'grid_option = "box" ! The general layout of the perpendicular grid.')
        replace(gs2_input_file,' ny = 48',f' ny = {3*PARAMS["naky"]}')
        # replace(gs2_input_file,' nx = 48',f' nx = {PARAMS["naky"]}')
    bashCommand = f"{gs2_executable} {gs2_input_file}"
    p = subprocess.Popen(bashCommand.split(),stderr=subprocess.STDOUT)#,stdout=subprocess.DEVNULL)#stdout=fp)
    p.wait()
    file2read = os.path.join(OUT_DIR,f"{gs2_input_name}.out.nc")
    eigenPlot(file2read)
    growth_rate, omega, ky = getgamma(file2read,savefig=save_fig)
    kyX, growthRateX, realFrequencyX = gammabyky(file2read,savefig=save_fig)
    weighted_growth_rate = np.sum(quasilinear_estimate(file2read,show=show_fig,savefig=save_fig))/PARAMS["naky"]
    # output_to_csv(output_csv, growth_rate, omega, ky, weighted_growth_rate, ln, lt)
    print(f'  LN={ln:1f}, LT={lt:1f}, growth rate={growth_rate:1f}, omega={omega:1f}, ky={ky:1f}, weighted gamma={weighted_growth_rate:1f} took {(time()-start_time_local):1f}s')
    return growth_rate, omega, ky, weighted_growth_rate
print(f'#######################')
print('Starting GS2 run')
print(f'#######################')
start_time = time()
growth_rate, omega, ky, weighted_growth_rate = run_gs2(PARAMS['LN'],PARAMS['LT'])
print(f'#######################')
print(f'Running GS2 took {time()-start_time}s')

remove_files()