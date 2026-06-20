
python plot_cell_scalars_4states.py -s beta_or_gamma -o final_cells_b0.0_g0.0026 -f -1 -a -x0 -820 -x1 820 -y0 -820 -y1 820
python plot_cell_scalars_4states.py -s beta_or_gamma -o final_cells_b0.0_g0.317 -f -1 -a -x0 -820 -x1 820 -y0 -820 -y1 820
python plot_cell_scalars_4states.py -s beta_or_gamma -o final_cells_b0.0_g0.7266 -f -1 -a -x0 -820 -x1 820 -y0 -820 -y1 820
python plot_cell_scalars_4states.py -s beta_or_gamma -o final_cells_b0.0_g0.9505 -f -1 -a -x0 -820 -x1 820 -y0 -820 -y1 820

#open final_cells_b0.0_g0.0026/keep.png
#open final_cells_b0.0_g0.317/keep.png
#open final_cells_b0.0_g0.7266/keep.png
#open final_cells_b0.0_g0.9505/keep.png
montage -geometry +0+0 -tile 4x1 final_cells_b0.0_g0.0026/keep.png final_cells_b0.0_g0.317/keep.png final_cells_b0.0_g0.7266/keep.png final_cells_b0.0_g0.9505/keep.png  physicell_beta0_final4.png

