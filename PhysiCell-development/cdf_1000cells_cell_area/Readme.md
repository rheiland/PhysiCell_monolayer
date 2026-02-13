'''
(base) M1P~/git/monolayer_PhysiCell/cdf_1000cells_cell_area$ python ../beta/param_00_1000cells_cell_area.py ../project 


'''

To plot results in the Studio (but won't have the 4-state colormap):

'''
python ~/git/studio_dev/bin/studio.py -c bg00_cells1000_0/config.xml 
'''
and then in Plot tab, Select the output dir, then select the .mat for cells.


To plot the last frame (.mat), showing f_i or a_i::
'''
python ../beta/plot_cell_scalars-2.py -o bg00_cells1000_0 -f -1 -s a_i
or,
python ../beta/plot_cell_scalars-2.py -o bg00_cells1000_0 -f -1 -s f_i
args= Namespace(output_dir='bg00_cells1000_0', frame=-1, axes_fixed=False, colorbar_name=None, show_colorbar=False, scalar_name='f_i', xmin=None, xmax=None, ymin=None, ymax=None)
unknown= []
output_dir=bg00_cells1000_0, current_frame=-1, axes_fixed=False, colorbar=RdBu, show_colorbar=False, xmax=100.0
'''

To plot results using the 4-state colormap ("-s beta_or_gamma"; however, it's not a discrete colormap!):
'''
python ../beta/plot_cell_scalars-2.py -s beta_or_gamma --show_colorbar -o bg00_cells1000_0 -f -1 

python ../beta/plot_cell_scalars_4states.py -s beta_or_gamma --show_colorbar -o bg00_cells1000_0 -f -1

'''

To plot histograms and cumulative dist fns (CDFs), do so from the root dir:
'''
(base) M1P~/git/monolayer_PhysiCell$ 
python beta/all_CDF.py a_i
python beta/all_CDF_percentiles.py

python beta/all_CDF.py 1 cdf_1000cells_cell_area bg00_cells1000_ a_i
python beta/all_CDF.py 100 cdf_1000cells_cell_area bg00_cells1000_ a_i


- prepare proper subdir and .csv files for Dom's script:
(base) M1P~/git/monolayer_PhysiCell/cdf_1000cells_cell_area$ python ../analysis/gen_cdf_csv.py

(base) M1P~/git/monolayer_PhysiCell/cdf_1000cells_cell_area$ python ../analysis/dom.py

'''

