
import sys
import os
import numpy as np
import glob
import csv
from pyMCDS_cells import pyMCDS_cells
import matplotlib.pyplot as plt

output_dirs = []

cell_cycle_duration = 443.5    # 5*T

# fig = plt.figure(figsize=(4,3), dpi=200, tight_layout=True)
fig = plt.figure(figsize=(6,4))
fig.subplots(1)
ax0 = fig.gca()

gamma_vals = np.arange(0.0, 1.0, 0.01)   # 100 samples
# gamma_vals = np.arange(0.0, 1.005, 0.5)  # test
# gamma_vals = np.arange(0.0, 0.10, 0.043)   # 100 samples
# print("type(gamma_vals)=", type(gamma_vals))
# print(gamma_vals)

tvals = []
gvals = []
for beta in [0.0]:
    # idx = -1
    # for gamma in np.arange(0.0, 1.05, 0.01):   # 100 samples
    # for gamma in np.arange(0.0, 0.10, 0.043):   # 100 samples
    for gamma in gamma_vals:
        folder_name = "out_num_cells_b" + str(beta) + "_g" + str(gamma)
        if (not os.path.exists(folder_name)):
            print("--- ERROR missing dir: ", folder_name)
            sys.exit(-1)

        label = "b:"+str(beta) + ",g:"+str(gamma)
        # idx += 1

        data_dir = "out_num_cells_b" + str(beta) + "_g" + str(gamma)
        # print('data_dir = ',data_dir)

        os.chdir(data_dir)
        xml_files = glob.glob('output*.xml')
        os.chdir('..')
        xml_files.sort()
        # print('xml_files = ',xml_files)

        ds_count = len(xml_files)
        # print("ds_count = ",ds_count)
        # ds_count = 192
        # print("----- ds_count = ",ds_count)
        mcds = [pyMCDS_cells(xml_files[i], data_dir) for i in range(ds_count)]

        tval = np.linspace(0, mcds[-1].get_time(), ds_count)

        tval /= cell_cycle_duration
        # print("tval= ",tval)
        final_time = tval[-1]
        # print(f'{data_dir} final time= {final_time}')

        if final_time > 0:
            tvals.append(final_time)
            gvals.append(gamma)
        else:
            print(f"  --- bogus time {final_time} in {data_dir} ")

file_out = f'gamma_time_10K.csv'
print("--> ",file_out)
with open(file_out, "w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(['gamma','time'])
    for jdx in range(len(gvals)):
        writer.writerow([gvals[jdx],tvals[jdx]])

# print("gvals=",gvals)
# print("tvals=",tvals)
ax0.plot(gvals,tvals)

ax0.set_title("beta=0", fontsize=12)
ax0.set_xlabel('gamma')
ax0.set_ylabel('Time (calibrated for 5T)')
# ax0.savefig(data_dir + '.png')
plt.show()

