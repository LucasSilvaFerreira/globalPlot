# globalPlot
Extract coordinates and global plot  for pybigwig


```
!pip uninstall globalPlot -y
!pip install git+https://github.com/LucasSilvaFerreira/globalPlot.git
```


```python from globalPlot import extract_coordinates_nostranded 
global_extract = extract_coordinates_nostranded('bed_merged_all.bed',    H1_atac_bw_list)
global_values = global_extract.extract_positions()
```
