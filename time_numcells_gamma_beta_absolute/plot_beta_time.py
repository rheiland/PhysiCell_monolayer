
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

# beta_vals = np.arange(0.0, 1.0, 0.01)   # 100 samples
# print("type(beta_vals)=", type(beta_vals))
# print(beta_vals)

# beta_vals = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95]
tvals = []
bvals = []
eps = 0.0001
for gamma in [0.0]:
    # idx = -1
    # for gamma in np.arange(0.0, 1.05, 0.01):   # 100 samples
    # for gamma in np.arange(0.0, 0.10, 0.043):   # 100 samples
    # for beta in beta_vals[:-1]:
    # for beta in beta_vals:
    # if primary_pass:
    for beta in np.arange(0.1, 0.95+eps, 0.01):
        b2 = int((beta+eps) * 100) / 100
        folder_name = "out_num_cells_b" + str(b2) + "_g" + str(gamma)
        if (not os.path.exists(folder_name)):
            print("--- ERROR missing dir: ", folder_name)
            sys.exit(-1)

        label = "b:"+str(b2) + ",g:"+str(gamma)
        # idx += 1

        data_dir = "out_num_cells_b" + str(b2) + "_g" + str(gamma)
        print('data_dir = ',data_dir)

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
        print(f'{data_dir} final time= {final_time}')

        if final_time > 0:
            tvals.append(final_time)
            bvals.append(b2)
        else:
            print(f"  --- bogus time {final_time} in {data_dir} ")


    # extra runs
    for beta in [0.996, 0.997, 0.998, 0.999]:
        b2 = beta
        folder_name = "out_num_cells_b" + str(b2) + "_g" + str(gamma)
        if (not os.path.exists(folder_name)):
            print("--- ERROR missing dir: ", folder_name)
            sys.exit(-1)

        label = "b:"+str(b2) + ",g:"+str(gamma)
        # idx += 1

        data_dir = "out_num_cells_b" + str(b2) + "_g" + str(gamma)
        print('data_dir = ',data_dir)

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
        print(f'{data_dir} final time= {final_time}')

        if final_time > 0:
            tvals.append(final_time)
            bvals.append(b2)
        else:
            print(f"  --- bogus time {final_time} in {data_dir} ")

file_out = f'beta_time_10K.csv'
print("--> ",file_out)
with open(file_out, "w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(['beta','time'])
    for jdx in range(len(bvals)):
        writer.writerow([bvals[jdx],tvals[jdx]])

# print("bvals=",bvals)
# print("tvals=",tvals)
# ax0.plot(bvals,tvals)
ax0.plot(bvals,tvals,'.-', color='k')
# ax0.plot(bvals[:-2],tvals[:-2])

ax0.set_title("gamma=0", fontsize=12)
ax0.set_xlabel('beta')
ax0.set_xlim(0., 1.1)
ax0.set_ylabel('Time (calibrated for 5T)')
# ax0.savefig(data_dir + '.png')
plt.show()
