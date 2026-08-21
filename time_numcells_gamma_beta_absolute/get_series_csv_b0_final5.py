#  Get CSV info for final (10K) cells of gamma sweep (beta=0):
#   x,y,r,inhibited
# 
import os
# import sys
import numpy as np
import csv
import glob
from pyMCDS_cells import pyMCDS_cells
from pathlib import Path


cell_cycle_duration = 443.5    # 5*T

gamma_final5 = [0.0026, 0.317, 0.7266, 0.887, 0.933]
# gamma_final5 = [0.0026]
print("gamma_final5=",gamma_final5)

for beta in [0.0]:
    # for gamma in gamma_final5:
    # for gamma, tval in zip(gamma_final5, calib_time):
    # for gamma, tval in zip(gamma_final5, calib_time):
    for gamma in gamma_final5:
        # save_interval = tval * cycle_duration / num_samples
        data_dir = "final_cells_b" + str(beta) + "_g" + str(gamma)


        csv_dir = "series_b" + str(beta) + "_g" + str(gamma)
        # print(f'      {csv_dir} --- calib time= {tval}')
        print(f'--- csv_dir= {csv_dir}')
        try:
            os.makedirs(csv_dir)
        except:
            print(f"... {csv_dir} already exists")
        csv_path = Path(csv_dir).resolve()

        # xml_file_out = os.path.join(data_dir, 'config.xml')  # copy config file into the output dir

        os.chdir(data_dir)
        xml_files = glob.glob('output*.xml')
        # xml_files = glob.glob(f'{data_dir}/output*.xml')
        # os.chdir('..')
        xml_files.sort()
        # print('xml_files = ',xml_files)

        ds_count = len(xml_files)
        print("----- ds_count = ",ds_count)

        for idx, xml_file in enumerate(xml_files):

            # mcds = [pyMCDS_cells(xml_files[i], data_dir) for i in range(ds_count)]
            # mcds = pyMCDS_cells(xml_file, data_dir)
            mcds = pyMCDS_cells(xml_file)

            # tval = np.linspace(0, mcds[-1].get_time(), ds_count)
            # tval = np.linspace(0, mcds.get_time(), ds_count)
            tval = mcds.get_time()

            # extract: x,y,r,a,f  (Dom will compute inhibition boolean based on these and render)
            # df_cells = mcds[-1].get_cell_df()
            df_cells = mcds.get_cell_df()
            xvals = df_cells['position_x']
            yvals = df_cells['position_y']
            radii = df_cells['cell_radius']
            a_i = df_cells['a_i']
            f_i = df_cells['f_i']
            inhibition_value = df_cells['beta_or_gamma']

            tval /= cell_cycle_duration
            # print("tval= ",tval)

            # file_out = f'physicell_10K_b0_g{gamma}.csv'

            # os.chdir(csv_dir)
            # directory = Path(csv_dir)
            # file_out = directory / f'pc_{idx}.csv'
            file_out = csv_path / f'pc_{idx:03d}.csv'
            if idx > 270:
                print("--> ",file_out)

            with open(file_out, "w", newline="") as file:
                writer = csv.writer(file)
                # writer.writerow(['x','y','inhibition'])
                writer.writerow(['x','y','r','a','f','inhibition_value'])
                for jdx in range(len(xvals)):
                    writer.writerow([xvals[jdx],yvals[jdx],radii[jdx],a_i[jdx],f_i[jdx],int(inhibition_value[jdx])])

        os.chdir('..')