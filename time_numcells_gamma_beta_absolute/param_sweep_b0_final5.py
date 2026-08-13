# This script provides simple parameter exploration functionality. The script creates
# a new folder (subdirectory) for each set of parameters, makes changes to a default 
# configuration (.xml) file using specified parameter values (in an accompanying .txt file),
# copies the new config file into the new folder, then
# runs the simulation (in the background) which writes results into the new folder.
# 

import xml.etree.ElementTree as ET
from shutil import copyfile
import os
import sys
# import subprocess
import numpy as np

# print(len(sys.argv))
if (len(sys.argv) < 2):
  usage_str = "Usage: %s <exec_pgm>" % (sys.argv[0])
  print(usage_str)
  print("e.g.:  python param_sweep.py project")
  exit(1)
else:
   exec_pgm = sys.argv[1]

background_str = " &"  # works on Unix
background_str = " "  # works on Unix
if sys.platform == 'win32':
    background_str = ""

xml_file_in = 'config_numcells_gamma_beta.xml'
xml_file_out = 'temp.xml'
copyfile(xml_file_in, xml_file_out)

tree = ET.parse(xml_file_out)
xml_root = tree.getroot()
output_dirs = []

# Ideally want T=[15, 27, 68, 136, 271], but settle for:
#  18 ===> gamma: 0.0026
#  26 ===> gamma: 0.317
#  68 ===> gamma: 0.7266  # 0.727 ?
# 135 ===> gamma: 0.887
# 266 ===> gamma: 0.933

# gamma_final4 = [0.0026, 0.3170, 0.7266, 0.9505]
gamma_final5 = [0.0026, 0.317, 0.7266, 0.887, 0.933]
print("gamma_final5=",gamma_final5)

# maps to calibrated times, T:
calib_time = [15, 27, 68, 136, 271]
# where, in plot_cell_scalars_4states.py:
# self.cycle_duration = 443.5
# calib_time = mins / self.cycle_duration
cycle_duration = 443.5
num_samples = 200
# --> led to:
save_intervals = [30, 60, 150, 300, 600]  

for beta in [0.0]:
    # for gamma in gamma_final5:
    # for gamma, tval in zip(gamma_final5, calib_time):
    # for gamma, tval in zip(gamma_final5, calib_time):
    for gamma, save_interval in zip(gamma_final5, save_intervals):
        # save_interval = tval * cycle_duration / num_samples
        print(f"{gamma}: save_interval={save_interval}")
        folder_name = "final_cells_b" + str(beta) + "_g" + str(gamma)
        output_dirs.append(folder_name)
        if (not os.path.exists(folder_name)):
            print("--- mkdir ", folder_name)
            os.makedirs(folder_name)

        xml_file_out = os.path.join(folder_name, 'config.xml')  # copy config file into the output dir

        print('---write config file (and start sim): ', xml_file_out)
        log_file = folder_name + ".log"  
        cmd =  exec_pgm + " " + xml_file_out + " > " + log_file + " " + background_str

        try:
            xml_root.find('.//' + 'folder').text = str(folder_name)   # beware of rules folder!
        except:
            print("--- Error setting output folder")
            exit(-1)

        try:
            xml_root.find('.//' + 'beta_threshold').text = str(beta)
            xml_root.find('.//' + 'gamma_threshold').text = str(gamma)

            xml_root.find(".//full_data//interval").text = str(save_interval)
        except:
            print("--- Error setting beta or gamma")
            exit(-1)

        # if gamma > 0.9:  # make output interval larger
            # xml_root.find('.//' + 'interval').text = '14400'

        xml_file_out = os.path.join(folder_name, 'config.xml')  # copy config file into the output dir
        tree.write(xml_file_out)   # will create folder_name/config.xml

        print("----- cmd = ",cmd)
        os.system(cmd)   # <------ Execute the simulation


print("\n ------\n Your output results will appear in these directories:\n   ",output_dirs)
print("and check for a .log file of each name for your terminal output from each simulation.\n")
