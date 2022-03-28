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
```python
from globalPlot import global_plot

for (g,c), v in df_merged_bed_centered.groupby(['group', 'cell']):
    print (g,c)
    
    v.to_csv('temp.bed', sep='\t', header=None, index=None)
    global_plot('temp.bed', modify_start_global=0, modify_end_global=0, vector_coordinates_global_unstranded= [H1_atac_bw_list[2]], color='blue', )
    global_plot('temp.bed', modify_start_global=0, modify_end_global=0, vector_coordinates_global_unstranded= [H2_atac_bw_list[2]], color='red')
    plt.title(f'ATAC-Seq {g} {c} \n  n: {v.shape[0]}', fontsize=20)
    plt.xticks([0,2000,4000], ['-2k', 'ATAC-center', '+2k'], fontsize=15)
    plt.ylim(0,2.6)
    plt.show()
```
