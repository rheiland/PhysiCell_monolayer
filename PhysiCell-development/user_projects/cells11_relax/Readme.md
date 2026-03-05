Reminder:
```
(base) M1P~/git/PhysiCell_monolayer/PhysiCell-development$ make load PROJ=cells11_relax   # or "save" if updating
(base) M1P~/git/PhysiCell_monolayer/PhysiCell-development$ make 
(base) M1P~/git/PhysiCell_monolayer/PhysiCell-development$ cp project project_11cells 

(base) M1P~/git/PhysiCell_monolayer/PhysiCell-development$ project_11cells config/test_11cells.xml 

- look for time_90pct.txt in the output folder, e.g.:
(base) M1P~/git/PhysiCell_monolayer/PhysiCell-development$ ty output_test_11cells/time_90pct.txt
88.3

or via Studio:
(base) M1P~/git/PhysiCell_monolayer/PhysiCell-development$ pcstudio -e project_11cells -c config/test_11cells.xml

```
