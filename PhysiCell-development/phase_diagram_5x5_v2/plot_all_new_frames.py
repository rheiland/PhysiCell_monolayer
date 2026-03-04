# cf. montage_monolayer.py

import os
#import subprocess

output_dirs = []

montage_cmd = "montage -geometry +0+0 -tile 5x5 "

# Roman had beta in X (columns), gamma in Y (rows)

# copy from param_sweep.py
gamma_vals = [0.0, 0.683, 0.816, 0.908, 0.937]
beta_vals  = [0.0, 0.970, 0.987, 0.996, 0.998]
print("g=",gamma_vals)
print("b=",gamma_vals)

gamma_r = list(reversed(gamma_vals))
# print("type(gamma_r)=", type(gamma_r))
print("g reversed=",gamma_r)
for gamma in gamma_r:
# for gamma in [0.937]:
    for beta in beta_vals:
    # for beta in [0.998]:
        output_dir = "out_cell_area_b" + str(beta) + "_g" + str(gamma)
        cmd = "python ../beta/plot_cell_scalars-2.py -o " + output_dir + " -s beta_or_gamma -f -1 -a -x0 -780 -x1 780 -y0 -780 -y1 780"   # --show_colorbar"
        os.system(cmd)

        montage_cmd += output_dir + "/keep.png "

# montage -geometry +0+0 -tile 5x5 run01/keep.png run02/keep.png run03/keep.png run04/keep.png run05/keep.png row1.png
montage_cmd += " montage.png"
print("montage_cmd= ",montage_cmd)