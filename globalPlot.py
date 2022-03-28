import pyBigWig
import pandas as pd
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from astropy.convolution import  Gaussian1DKernel, convolve

class extract_coordinates_nostranded:
    def __init__(self, bed_file, vector_coordinates, modify_start=0, modify_end=0, bed_header=False, force_strand=True,force_return_same_direction=False, set_tss=False ,set_tes=False, create_100_bins  = False):
        
        if bed_header == False:
            bed_header = None
        else:
            bed_header = 1
            
        self.bed = pd.read_csv( bed_file, sep='\t', header=bed_header)
        #print self.bed
        self.vector_coodinates = vector_coordinates
    
        self.modify_start = modify_start
        self.modify_end =   modify_end
        self.force_strand = True
        self.force_return_same_direction = force_return_same_direction
        self.set_tss=set_tss
        self.set_tes=set_tes
        self.create_100_bins = create_100_bins
        
    def extract_positions(self):
        extracted_coords = []
        for p in tqdm(self.bed.values.tolist(), desc='extracting...' ):
            p[1], p[2] = int(p[1]), int(p[2])
            
            if self.set_tss:
                if p[5] == '+':
                    p[1] = p[1] 
                    p[2] = p[1]+1
                else:
                    p[1] = p[2] 
                    p[2] = p[2]+1
            
            if self.set_tes:
                if p[5] == '+':
                    p[1] = p[2]
                    p[2] = p[2]+1
                else:
                    p[1] = p[1] 
                    p[2] = p[1]+1
            
            
            
            if self.force_strand:
                if p[5] == '+':
                    p[1] = p[1] - self.modify_start
                    p[2] = p[2] + self.modify_end
                else:
                    p[1] = p[1] -  self.modify_end
                    p[2] = p[2] +  self.modify_start 
                    
            
            return_vector = []
            for m in self.vector_coodinates :
                bigwig_file = pyBigWig.open(m)
                try:
                    if self.force_return_same_direction and p[5] == '-':
                        if self.create_100_bins:
                            remove_none  = np.nan_to_num( bigwig_file.stats(p[0], int(p[1]), int(p[2]), type="mean", nBins=100 ))[::-1]
                            remove_zero = [i if i is not None else 0 for i in remove_zero]
                            return_vector.append (np.array(remove_none))
                        else:
                            return_vector.append (np.nan_to_num( np.array (bigwig_file.values(p[0], int(p[1]), int(p[2]) ) ))[::-1] )
                    else:
                        if self.create_100_bins:
                            remove_zero  = np.nan_to_num( bigwig_file.stats(p[0], int(p[1]), int(p[2]), type="mean", nBins=100 ))
                            remove_zero = [i if i is not None else 0 for i in remove_zero]
                            return_vector.append (np.array(remove_zero))
                        else:
                            return_vector.append (np.nan_to_num( np.array (bigwig_file.values(p[0], int(p[1]), int(p[2]) ) )) )
                except:
                    #print(p[0], int(p[1]), int(p[2]), 'ISSUE REGION')
                    
                    if self.create_100_bins:
                        return_vector.append( np.array([0] *  100 ))

                    else:
                        return_vector.append( np.array([0] *  (int(p[2]) - int(p[1])) ))
                    
                    
                    
                bigwig_file.close()
            extracted_coords.append([return_vector, p])
            #print (extracted_coords)
        return extracted_coords  
    
    
    
    
class extract_coordinates_stranded:
    def __init__(self, bed_file, vector_positive, vector_negative,  modify_start=0, modify_end=0, bed_header=None, force_strand=True, force_return_same_direction=False, set_tss=False ,set_tes=False, create_100_bins  = False):
        
        if bed_header == False:
            bed_header = None
        else:
            bed_header = 1
            
            
        self.bed = pd.read_csv( bed_file, sep='\t', header=bed_header)
        #print self.bed
        self.vector_positive = vector_positive
        self.vector_negative = vector_negative
        
        self.modify_start = modify_start
        self.modify_end =   modify_end
        self.force_strand = True
        self.force_return_same_direction = force_return_same_direction
        self.set_tss=set_tss
        self.set_tes=set_tes
        self.create_100_bins = create_100_bins
        
        
    def extract_positions(self):
        extracted_coords = []
        for p in tqdm(self.bed.values.tolist(),desc='extracting..'):
            p[1], p[2] = int(p[1]), int(p[2])
            
            if self.set_tss:
                if p[5] == '+':
                    p[1] = p[1] 
                    p[2] = p[1]+1
                else:
                    p[1] = p[2] 
                    p[2] = p[2]+1
            
            if self.set_tes:
                if p[5] == '+':
                    p[1] = p[2]
                    p[2] = p[2]+1
                else:
                    p[1] = p[1] 
                    p[2] = p[1] + 1
            
            #print (p[:4], 'uma base?')
            if self.force_strand:
                if p[5] == '+':
                    p[1] = p[1] - self.modify_start
                    p[2] = p[2] + self.modify_end
                else:
                    p[1] = p[1] -  self.modify_end
                    p[2] = p[2] +  self.modify_start 
                    
            
            return_vector = []
            for plus, minus  in zip(self.vector_positive, self.vector_negative):
               
                if p[5] == '+':
                    modify_strand = 1
                    m = plus
                else:
                    modify_strand = -1
                    m = minus
                        
                
                #print (p[:4])
                bigwig_file = pyBigWig.open(m)
                try:
                    if self.force_return_same_direction and p[5] == '-':
                        if self.create_100_bins: 
                            remove_zero = np.nan_to_num( np.array (bigwig_file.stats(p[0], int(p[1]), int(p[2]),  type="mean", nBins=100 ) ))[::-1] * modify_strand
                            remove_zero = [i if i is not None else 0 for i in remove_zero]
                            return_vector.append (np.array(remove_zero))
                        else:
                            return_vector.append (np.nan_to_num( np.array (bigwig_file.values(p[0], int(p[1]), int(p[2]) ) ))[::-1] * modify_strand )
                    else:
                        if self.create_100_bins:
                            remove_zero = np.nan_to_num( np.array (bigwig_file.stats(p[0], int(p[1]), int(p[2]),  type="mean", nBins=100 ) )) * modify_strand 
                            remove_zero = [i if i is not None else 0 for i in remove_zero]
                            return_vector.append (np.array(remove_zero))
                            
                        else:
                            return_vector.append (np.nan_to_num( np.array (bigwig_file.values(p[0], int(p[1]), int(p[2]) ) )) * modify_strand )
                        
                except:
                    #print(p[0], int(p[1]), int(p[2]), 'ISSUE REGION')
                    if self.create_100_bins:
                        return_vector.append( np.array([0] *  100 ))

                    else:
                        return_vector.append( np.array([0] *  (int(p[2]) - int(p[1])) ))
                bigwig_file.close()
            extracted_coords.append([return_vector, p])
        return extracted_coords
    

    
    
    
# def split_same_len(a, n):
#     k, m = divmod(len(a), n)
#     return (a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in np.arange(n))




    
def create_percentage(numpy_values):
    
    if numpy_values.shape[0] < 100:
        print ('ERROR: coordinate range smallest than 100 bp')
    else:
        #print numpy_values.shape[0]
        return numpy_values



    

    
def global_plot(bed_file_global, modify_start_global=0 ,modify_end_global=0, vector_coordinates_global_unstranded=None, vector_coordinates_global_plus_samples=None ,     vector_coordinates_global_minus_samples = None, color='red', merge_operation = 'mean', view_values_distribuition=False, remove_outlier=False, bed_header=False, use_len_percentage=False, set_tss=False, set_tes=False, smooth_plot=0,heatmap_output=False, percent_add_ext_start_and_end=0, __supress_plot__=False):
    
    '''
    Usage example:
    
        global_plot(bed, vector_coordinates_global_plus_samples= [arquivos[0], arquivos[2]] , vector_coordinates_global_minus_samples=[arquivos[1], arquivos[3]], modify_end_global=150, modify_start_global=150, color='red')
        global_plot(bed, vector_coordinates_global_plus_samples= [arquivos_2[0], arquivos_2[2]] , vector_coordinates_global_minus_samples=[arquivos_2[1], arquivos_2[3]],  modify_end_global=150, modify_start_global=150, color='blue')
        
    
    
    
    '''
    boolean_create_100_bins=False
    if use_len_percentage:
        boolean_create_100_bins = True
    

    
    
    if vector_coordinates_global_unstranded == None and vector_coordinates_global_plus_samples == None and vector_coordinates_global_minus_samples == None:
        print ('Provide  a vector contaning the unstranded samples (vector_coordinates_global_unstranded) \n or \n Provide  positive (vector_coordinates_global_plus_samples) and negative vectors (vector_coordinates_global_minus_samples)')
        
        
        
    
    if vector_coordinates_global_unstranded != None:
    
        object_extracted = extract_coordinates_nostranded(bed_file=bed_file_global,
                                                           modify_end=modify_end_global,
                                                           modify_start=modify_start_global,
                                                          vector_coordinates=vector_coordinates_global_unstranded,
                                                          force_return_same_direction=True,
                                                          bed_header=bed_header,
                                                          set_tss=set_tss,
                                                          set_tes=set_tes,
                                                          create_100_bins=boolean_create_100_bins)


    
        regions = object_extracted.extract_positions()

    if vector_coordinates_global_plus_samples != None and  vector_coordinates_global_minus_samples != None:
        
        object_extracted = extract_coordinates_stranded(bed_file=bed_file_global,
                               modify_end=modify_end_global,
                               modify_start=modify_start_global,
                               vector_positive=vector_coordinates_global_plus_samples,
                               vector_negative  =vector_coordinates_global_minus_samples, 
                                                        force_return_same_direction=True,
                                                        bed_header=bed_header,
                                                           set_tss=set_tss,
                                                          set_tes=set_tes,
                                                        create_100_bins=boolean_create_100_bins)





    
        regions = object_extracted.extract_positions()

    plot_data = []
    

    for r in regions:
        bw_regions = r[0]
        if merge_operation == 'mean':
            
            plot_data.append(np.array([x.tolist() for x in  bw_regions]).mean(axis=0))
        elif  merge_operation == 'sum':
            plot_data.append(np.array([x.tolist() for x in  bw_regions]).sum(axis=0))
        else:
            print ("use a valid operation: ('mean' or 'sum') ")

    
    
    if use_len_percentage:
        
        plot_data = [create_percentage(p) for p in plot_data]

    
    
#     df_outlier =  pd.DataFrame()
#     df_outlier['values'] = np.array(plot_data).mean(axis=1)
#     df_outlier['order'] = df_outlier.index.values.tolist()
#     df_outlier = df_outlier.sort_values('values')
#     #print (df_outlier)
#     two_and_a_half_percent = int(df_outlier.shape[0] * 0.005)
#     #print (two_and_a_half_percent)
#     possible_outliers_set = set(df_outlier['order'].values.tolist()[:two_and_a_half_percent] + df_outlier['order'].values.tolist()[-two_and_a_half_percent:])
#     #print (possible_outliers_set)
    

    
   


    
    if remove_outlier:

        
        
        df_to_zscore = pd.DataFrame(plot_data)
        df_zscored = (df_to_zscore - df_to_zscore.mean())/df_to_zscore.std()
        c_sub = df_to_zscore.columns.values
        reconstruct_df = []
        for c in  tqdm(c_sub):
            df_ope = pd.DataFrame( {'zscored': df_zscored[c].values,
                           'original_values': df_to_zscore[c].values}
                                 )
            
            
            
            get_limit_zscore_cut = df_ope[5- df_ope.zscored < 5 ].sort_values('original_values').iloc[0]['original_values']
            df_to_zscore[c] = df_to_zscore[c].apply(lambda x : get_limit_zscore_cut if x >get_limit_zscore_cut else x)
    
    
        plot_data = df_to_zscore.values # smoothing outliers
    
    
    
    

    
    
    
    
    if view_values_distribuition: # Checking value distribuition Plot
        print ('alert! DISABLE view_values_distribuition in order to plot the global profile')
        plt.scatter(  range(0, len(vector_stats.tolist())), vector_stats.tolist() )
        plt.show()
        sns.heatmap(plot_data)
        plt.show()
        
        df_to_zscore = pd.DataFrame(plot_data)
        df_zscored = (df_to_zscore - df_to_zscore.mean())/df_to_zscore.std()
        sns.boxplot(x='variable', y='value',data=df_zscored.melt())
        plt.title('zscore distribuition before filter')
        plt.show()
        plt.plot(df_to_zscore.mean())
        plt.title('beforer_remove')
        plt.show()
        
        c_sub = df_to_zscore.columns.values
        reconstruct_df = []
        for c in  tqdm(c_sub):
            df_ope = pd.DataFrame( {'zscored': df_zscored[c].values,'original_values': df_to_zscore[c].values} )
            get_limit_zscore_cut = df_ope[5- df_ope.zscored < 5 ].sort_values('original_values').iloc[0]['original_values']
            df_to_zscore[c] = df_to_zscore[c].apply(lambda x : get_limit_zscore_cut if x >get_limit_zscore_cut else x)
    
    

        plt.plot(df_to_zscore.mean())
        plt.title('after_remove')
        plt.show()
        sns.heatmap(df_to_zscore)
        plt.title('after_remove')
        plt.show()


       
        for c in  tqdm(c_sub):
            df_ope = pd.DataFrame( {'zscored': df_zscored[c].values,'original_values': df_to_zscore[c].values} )
            get_limit_zscore_cut = df_ope[5- df_ope.zscored < 5 ].sort_values('original_values').iloc[0]['original_values']
            df_to_zscore[c] = df_to_zscore[c].apply(lambda x : get_limit_zscore_cut if x >get_limit_zscore_cut else x)
    
    
        plot_data = df_to_zscore.values # smoothing outliers
    
        
    

    
    else:
       
            
        #print ([len(x_data)  for x_data in plot_data])
        data_to_plot  = np.array(plot_data).mean(axis=0)
        if smooth_plot>0 and percent_add_ext_start_and_end==0:
            gauss_kernel = Gaussian1DKernel(smooth_plot)  
            data_to_plot = convolve(data_to_plot, gauss_kernel, boundary='extend')
        


        if percent_add_ext_start_and_end >0:
            
            p_ext = percent_add_ext_start_and_end
            to_add_start = global_plot(bed_file_global,
                        modify_start_global=p_ext,
                        modify_end_global=0,
                        vector_coordinates_global_unstranded=vector_coordinates_global_unstranded,
                        vector_coordinates_global_plus_samples=vector_coordinates_global_plus_samples,
                        vector_coordinates_global_minus_samples = vector_coordinates_global_minus_samples,
                        color=color,
                        merge_operation = merge_operation,
                        view_values_distribuition=False,
                        remove_outlier=remove_outlier,
                        bed_header=bed_header,
                        use_len_percentage=False,
                        set_tss=True, set_tes=False,
                        smooth_plot=smooth_plot,
                        heatmap_output=True,
                        percent_add_ext_start_and_end=0,
                        __supress_plot__=True)

            to_add_stop = global_plot(bed_file_global,
                        modify_start_global=0,
                        modify_end_global=p_ext,
                        vector_coordinates_global_unstranded=vector_coordinates_global_unstranded,
                        vector_coordinates_global_plus_samples=vector_coordinates_global_plus_samples,
                        vector_coordinates_global_minus_samples = vector_coordinates_global_minus_samples,
                        color=color,
                        merge_operation = merge_operation,
                        view_values_distribuition=False,
                        remove_outlier=remove_outlier,
                        bed_header=bed_header,
                        use_len_percentage=False,
                        set_tss=False, set_tes=True,
                        smooth_plot=smooth_plot,
                        heatmap_output=True,
                        percent_add_ext_start_and_end=0,
                        __supress_plot__=True)
            
        
            start_data_to_plot = to_add_start[0].mean(axis=0)
            stop_data_to_plot = to_add_stop[0].mean(axis=0)
            
            
            
            #print ('__EXT__ ultimo')
            #print (start_data_to_plot.shape,data_to_plot.shape, stop_data_to_plot.shape )
            
            start_data_percent       = start_data_to_plot.tolist()
            final_point_data_percent = data_to_plot.tolist()
            end_data_percent         =  stop_data_to_plot.tolist()
            
            y_plot= start_data_percent + final_point_data_percent + end_data_percent
            #print (y_plot)
            
            x_plot = [(float(100)/len(start_data_percent))*r for r in range(len(start_data_percent))] + [100+((float(800)/len(final_point_data_percent))) * r for r in range(len(final_point_data_percent))] + [ 900+ ((float(100)/len(end_data_percent)) *r) for r in range(len(end_data_percent))]
            plt.plot(x_plot, y_plot, color=color)
            _ = plt.xticks([0,100,366,366+266,900,1000], ('-{}'.format(len(start_data_percent)+1), '%0', '33%', '%66', '%100', str(len(end_data_percent)-1) ))


        elif __supress_plot__== False:
            plt.plot(data_to_plot, color=color)
        
    if heatmap_output:
        return plot_data, [i[-1] for i in regions]
    


    
    
    return regions
