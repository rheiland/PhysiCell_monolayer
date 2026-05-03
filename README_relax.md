# PhysiCell model for a growing monolayer (an OpenVT reference model)

```
(base) M1P~/git/mechanics_relaxation$ make load PROJ=cells11_final
(base) M1P~/git/mechanics_relaxation$ make 
(base) M1P~/git/mechanics_relaxation$ project  
or:
(base) M1P~/git/mechanics_relaxation$ pcstudio  

(base) M1P~/git/mechanics_relaxation$ python analysis/plot_11cells_crop.py 88
(generates pc_plot_11cells.csv)
```
<img src="images/physicell_relax11.png" width=400/>

```
- assuming we have Chaste's results
(base) M1P~/git/mechanics_relaxation$ python analysis/plot_11cells_csv.py pc_plot_11cells.csv
```
<img src="images/physicell_vs_chaste.png" width=500/>

```
(base) M1P~/git/mechanics_relaxation$ python analysis/plot_21cells_crop.py 1440
(generates pc_plot_21cells_width.csv)
```
<img src="images/physicell_relax21.png" width=400/>

<img src="images/mechanics_cells21_relax_2026-04-22.gif" width=400/>


## Funding

National Science Foundation 2303695 and National Cancer Institute 1U24CA284156-01A1

