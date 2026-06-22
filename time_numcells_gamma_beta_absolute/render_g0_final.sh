
python plot_cell_scalars_4states.py -s beta_or_gamma -o final_cells_b0.663_g0.0 -f -1 -a -x0 -820 -x1 820 -y0 -820 -y1 820
python plot_cell_scalars_4states.py -s beta_or_gamma -o final_cells_b0.925_g0.0 -f -1 -a -x0 -820 -x1 820 -y0 -820 -y1 820

# python plot_cell_scalars_4states.py -s beta_or_gamma -o final_cells_b0.963_g0.0 -f -1 -a -x0 -820 -x1 820 -y0 -820 -y1 820
# python plot_cell_scalars_4states.py -s beta_or_gamma -o out_num_cells_b0.996_g0.0 -f -1 -a -x0 -820 -x1 820 -y0 -820 -y1 820
# python plot_cell_scalars_4states.py -s beta_or_gamma -o out_num_cells_b0.97_g0.0 -f -1 -a -x0 -820 -x1 820 -y0 -820 -y1 820
# python plot_cell_scalars_4states.py -s beta_or_gamma -o out_num_cells_b0.98_g0.0 -f -1 -a -x0 -820 -x1 820 -y0 -820 -y1 820
python plot_cell_scalars_4states.py -s beta_or_gamma -o out_num_cells_b0.99_g0.0 -f -1 -a -x0 -820 -x1 820 -y0 -820 -y1 820


python plot_cell_scalars_4states.py -s beta_or_gamma -o out_num_cells_b0.9991_g0.0 -f -1 -a -x0 -820 -x1 820 -y0 -820 -y1 820

python plot_cell_scalars_4states.py -s beta_or_gamma -o out_num_cells_b0.9976_g0.0 -f -1 -a -x0 -820 -x1 820 -y0 -820 -y1 820

# python plot_cell_scalars_4states.py -s beta_or_gamma -o out_num_cells_b0.9999_g0.0 -f -1 -a -x0 -820 -x1 820 -y0 -820 -y1 820

#open final_cells_b0.0_g0.0026/keep.png
#open final_cells_b0.0_g0.317/keep.png
#open final_cells_b0.0_g0.7266/keep.png
#open final_cells_b0.0_g0.9505/keep.png
# out_num_cells_b0.98_g0.0/keep.png \
# out_num_cells_b0.996_g0.0/keep.png
montage -geometry +0+0 -tile 5x1 final_cells_b0.663_g0.0/keep.png \
final_cells_b0.925_g0.0/keep.png \
out_num_cells_b0.99_g0.0/keep.png \
out_num_cells_b0.9976_g0.0/keep.png \
out_num_cells_b0.9991_g0.0/keep.png  physicell_gamma0_final.png

# montage -geometry +0+0 -tile 5x1 final_cells_b0.663_g0.0/keep.png \
# final_cells_b0.925_g0.0/keep.png \
# out_num_cells_b0.996_g0.0/keep.png \
# out_num_cells_b0.9976_g0.0/keep.png \
# out_num_cells_b0.9991_g0.0/keep.png  physicell_gamma0_final.png

echo "\n ---> physicell_gamma0_final.png"
