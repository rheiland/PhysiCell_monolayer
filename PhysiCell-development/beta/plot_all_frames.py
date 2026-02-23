# cf. montage_monolayer.py

import os
#import subprocess

# output_dirs = ['out_b0.5_g0.5', 'out_b0.5_g0.6', 'out_b0.5_g0.7', 'out_b0.5_g0.8', 'out_b0.5_g0.9', 'out_b0.6_g0.5', 'out_b0.6_g0.6', 'out_b0.6_g0.7', 'out_b0.6_g0.8', 'out_b0.6_g0.9', 'out_b0.7_g0.5', 'out_b0.7_g0.6', 'out_b0.7_g0.7', 'out_b0.7_g0.8', 'out_b0.7_g0.9', 'out_b0.8_g0.5', 'out_b0.8_g0.6', 'out_b0.8_g0.7', 'out_b0.8_g0.8', 'out_b0.8_g0.9', 'out_b0.9_g0.5', 'out_b0.9_g0.6', 'out_b0.9_g0.7', 'out_b0.9_g0.8', 'out_b0.9_g0.9']
# output_dirs = ['out_b0.5_g0.5', 'out_b0.5_g0.6']
output_dirs = []

montage_cmd = "montage -geometry +0+0 -tile 9x9 "
montage_cmd = "montage -geometry +0+0 -tile 5x5 "

# Roman had beta in X (columns), gamma in Y (rows)
gamma_vals = [0.5, 0.6, 0.7, 0.8, 0.9]
gamma_vals = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95]
gamma_vals = [0.3, 0.5, 0.7, 0.8, 0.9]
print("type(gamma_vals)=", type(gamma_vals))
print(gamma_vals)
gamma_r = list(reversed(gamma_vals))
print("type(gamma_r)=", type(gamma_r))
print(gamma_r)
for gamma in gamma_r:
    # for beta in [0.5, 0.6, 0.7, 0.8, 0.9]:
    # for beta in [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95]:
    for beta in [ 0.6, 0.7, 0.8, 0.9, 0.95]:
        # output_dir = "out_new_b" + str(beta) + "_g" + str(gamma)
        output_dir = "out_cell_area_b" + str(beta) + "_g" + str(gamma)
        cmd = "python beta/plot_cell_scalars.py -o " + output_dir + " -s beta_or_gamma -f -1 -a -x0 -800 -x1 800 -y0 -800 -y1 800"
        os.system(cmd)

        montage_cmd += output_dir + "/keep.png "

# montage -geometry +0+0 -tile 5x5 run01/keep.png run02/keep.png run03/keep.png run04/keep.png run05/keep.png row1.png
montage_cmd += " montage.png"
print("montage_cmd= ",montage_cmd)