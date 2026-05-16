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
# first_time = True
output_dirs = []

# NEW: 5/5/26 -------------------
# beta_vals=  [0.0, 0.988, 0.995, 0.999, 1.0]
# beta_vals=  [0.0]
# gamma_vals=  [0.0, 0.791, 0.893, 0.937, 0.975]
# gamma_vals=  [0.791, 0.893, 0.937, 0.975]
# gamma_vals=  [0.5]

# for beta in [0.0]:
for gamma in [0.0]:
    # for gamma in gamma_vals:
    # for gamma in np.arange(0.0, 1.005, 0.01):   # 100 samples
    # for gamma in [0.36]:   # 100 samples
        # if gamma == 0.0 or gamma == 0.5 or gamma == 1.0:
        #     continue
    # for gamma in np.arange(0.0, 1.005, 0.5):   # test
    for beta in np.arange(0.0, 1.005, 0.01):   # 100 samples
    # for beta in [0.52]:   # 100 samples
        folder_name = "out_num_cells_b" + str(beta) + "_g" + str(gamma)
        output_dirs.append(folder_name)
        if (not os.path.exists(folder_name)):
            print("--- mkdir ", folder_name)
            os.makedirs(folder_name)

        xml_file_out = os.path.join(folder_name, 'config.xml')  # copy config file into the output dir

        print('---write config file (and start sim): ', xml_file_out)
        # tree.write(xml_file_out)   # will create folder_name/config.xml
        log_file = folder_name + ".log"  
        cmd =  exec_pgm + " " + xml_file_out + " > " + log_file + " " + background_str
        # print("----- cmd = ",cmd)
        # os.system(cmd)   # <------ Execute the simulation
        # subprocess.Popen([exec_pgm, xml_file_out])
        # with open(log_file,"w") as outf:
        #     subprocess.Popen([exec_pgm, xml_file_out],stdout=outf)

        try:
            xml_root.find('.//' + 'folder').text = str(folder_name)   # beware of rules folder!
        except:
            print("--- Error setting output folder")
            exit(-1)

        try:
            xml_root.find('.//' + 'beta_threshold').text = str(beta)
            xml_root.find('.//' + 'gamma_threshold').text = str(gamma)
        except:
            print("--- Error setting beta or gamma")
            exit(-1)

        xml_file_out = os.path.join(folder_name, 'config.xml')  # copy config file into the output dir
        tree.write(xml_file_out)   # will create folder_name/config.xml

        print("----- cmd = ",cmd)
        os.system(cmd)   # <------ Execute the simulation

#---------------
if False:
    # beta_vals=  [0.0, 0.988, 0.995, 0.999, 1.0]
    beta_vals=  [0.988, 0.995, 0.999, 1.0]
    beta_vals=  [0.5]
    # gamma_vals=  [0.0, 0.791, 0.893, 0.937, 0.975]
    gamma_vals=  [0.0]

    for beta in beta_vals:
        for gamma in gamma_vals:
            folder_name = "out_num_cells_b" + str(beta) + "_g" + str(gamma)
            output_dirs.append(folder_name)
            if (not os.path.exists(folder_name)):
                print("--- mkdir ", folder_name)
                os.makedirs(folder_name)

            xml_file_out = os.path.join(folder_name, 'config.xml')  # copy config file into the output dir

            print('---write config file (and start sim): ', xml_file_out)
            # tree.write(xml_file_out)   # will create folder_name/config.xml
            log_file = folder_name + ".log"  
            cmd =  exec_pgm + " " + xml_file_out + " > " + log_file + " " + background_str
            # print("----- cmd = ",cmd)
            # os.system(cmd)   # <------ Execute the simulation
            # subprocess.Popen([exec_pgm, xml_file_out])
            # with open(log_file,"w") as outf:
            #     subprocess.Popen([exec_pgm, xml_file_out],stdout=outf)

            try:
                xml_root.find('.//' + 'folder').text = str(folder_name)   # beware of rules folder!
            except:
                print("--- Error setting output folder")
                exit(-1)

            try:
                xml_root.find('.//' + 'beta_threshold').text = str(beta)
                xml_root.find('.//' + 'gamma_threshold').text = str(gamma)
            except:
                print("--- Error setting beta or gamma")
                exit(-1)

            xml_file_out = os.path.join(folder_name, 'config.xml')  # copy config file into the output dir
            tree.write(xml_file_out)   # will create folder_name/config.xml

            print("----- cmd = ",cmd)
            os.system(cmd)   # <------ Execute the simulation

print("\n ------\n Your output results will appear in these directories:\n   ",output_dirs)
print("and check for a .log file of each name for your terminal output from each simulation.\n")
