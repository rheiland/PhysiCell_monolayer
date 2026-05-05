python ../analysis/gen_cdf_time_series_data.py 100 5   #   <max_runs> <cell_radius>
python ../analysis/gen_cdf_idx_time_csv.py 100 443.5   #   <max_runs>  <5*T cycle duration>
rm -rf PhysiCell_MonolayerGrowth_1000_Data
mkdir PhysiCell_MonolayerGrowth_1000_Data
mv run* PhysiCell_MonolayerGrowth_1000_Data
zip -r PhysiCell_MonolayerGrowth_1000_Data.zip PhysiCell_MonolayerGrowth_1000_Data   # confirm < 100MB for github!
#mv PhysiCell_MonolayerGrowth_1000_Data.zip ../results/
