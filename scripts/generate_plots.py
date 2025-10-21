#%% libraries
import numpy as np
from numpy import arange, add, array
from scipy import stats
from scipy.stats import sem
from numpy.random import uniform
from pandas import DataFrame, concat, Categorical
import plotnine
from plotnine import ggplot, aes, geom_point, geom_jitter, ylim, theme, geom_rect, xlab, ylab, labs, element_text, geom_vline, geom_segment, coord_cartesian, scale_fill_brewer, scale_fill_hue, qplot, ggtitle


from collections import Counter
import json
import os
from os.path import join
import glob


#%% get all partitions

#%% data

centers1 = [0.1, 0.9, 0.1, 0.9]
partition1 = [[1,2], [3,4]]
centers2 = [0.1, 0.9, 0.9, 0.9]
partition2 = [[1,2], [3], [4]]
centers3 = centers2
partition3 = [[1], [2], [3], [4]]
centers4 = [0.9] * 4
partition4 = partition3


p1s = arange(0.5, 1 + 1/40, 1/40)

centers = [centers1] * 12 + [centers2] * 1 + [centers2] * 1 + [centers2] * 2 + [centers2] * 1 + [centers2] * 3 + [centers4] * 1

partitions = [partition1] * 12 + [partition2] * 1 + [partition2] * 1 + [partition2] * 2 + [partition2] * 1 + [partition2] * 3 + [partition4] * 1

def unlist(list_of_lists):
    
    return [item for sublist in list_of_lists for item in sublist]

def make_plotting_mini_dataframe(p1, centers, partition, pch_width):
    fecundity = len(centers)
    p1s = [p1] * fecundity
    partitions = [[str(i)] * len(el) for i, el in enumerate(partition)]
    partitions = unlist(partitions)
    
    
    

    if len(partition) == 1:
        shift_ = [p1] * fecundity
    elif len(partition) == 2:
        shift_ = []
        for i, el in enumerate(partition):
            if i == 0:
                shift_ += [p1 - 0.5 * pch_width] * len(el)
            elif i == 1:
                shift_ += [p1 + 0.5 * pch_width] * len(el)
    elif len(partition) == 3:
        shift_ = []
        for i, el in enumerate(partition):
            if i == 0:
                shift_ += [p1 - pch_width] * len(el)
            elif i == 1:
                shift_ += [p1] * len(el)
            elif i == 2:
                shift_ += [p1 + pch_width] * len(el)
    elif len(partition) == 4:
        shift_ = []
        for i, el in enumerate(partition):
            if i == 0:
                shift_ += [p1 - 1.5 * pch_width] * len(el)
            elif i == 1:
                shift_ += [p1 - 0.5 * pch_width] * len(el)
            elif i == 2:
                shift_ += [p1 + 0.5 * pch_width] * len(el)
            elif i == 3:
                shift_ += [p1 + 1.5 * pch_width] * len(el)
    
    
    return DataFrame({"p1": p1s, "center": centers, "patch": partitions, "shift_": shift_})
    
def make_plotting_dataframe(p1s, centers, partitions, pch_width):
    mini_frames = [make_plotting_mini_dataframe(el_p1, el_cent, el_part, pch_width) for el_p1, el_cent, el_part in zip(p1s, centers, partitions)]
    
    data_frame = concat(mini_frames, ignore_index = True)
    
    return data_frame

data = make_plotting_dataframe(p1s, centers, partitions, pch_width = 0.006)
#%% rectangles
possible_partitions = [[1, 2, 3, 4]], [[1, 2, 3], [4]], [[1, 2], [3, 4]], [[1, 2], [3], [4]], [[1], [2], [3], [4]]
partition_strings = ["1/2/3/4", "1/2/3, 4", "1/2, 3/4", "1/2, 3, 4", "1, 2, 3, 4"]

unique_p1 = sorted(data['p1'].unique())
rectangles = DataFrame({
    'xmin': [x - 1/(2*40) for x in unique_p1],  
    'xmax': [x + 1/(2*40) for x in unique_p1],  
    'ymin': [0] * len(unique_p1),  
    'ymax': [1] * len(unique_p1),  
    'fill': [partition_strings[possible_partitions.index(el)] for el in partitions]  # Alternating colors
})





 


#%% read in data
os.chdir('/Users/avivbrokman/Documents/Kentucky/Grad School/ms_project/branching1')

def read_json(file_path):
    with open(file_path, 'r') as file:
        data = json.load(file)
    return data

p1s = arange(0.5, 1 + 1/40, 1/40)


#%% utilities
def unlist(list_of_lists):
    
    return [item for sublist in list_of_lists for item in sublist]
#%% 
def calculate_patch_shifts(num_patches, pch_width):
    
    return [(i - 0.5 * (num_patches - 1)) * pch_width for i in range(num_patches)]

def calculate_partition_shift(partition, pch_width):
    
    patch_shifts = calculate_patch_shifts(len(partition), pch_width)
    
    offspring_shifts = [[patch_shifts[i]] * len(el) for i, el in enumerate(partition)]
    
    return unlist(offspring_shifts)

def make_plotting_mini_dataframe(p1, centers, partition, pch_width):
    fecundity = len(centers)
    p1s = [p1] * fecundity
    partitions = [[str(i)] * len(el) for i, el in enumerate(partition)]
    partitions = unlist(partitions)
    offsets = calculate_partition_shift(partition, pch_width)
    
    return DataFrame({"p1": p1s, "center": centers, "patch": partitions, "shifted_p1": add(p1s, offsets)})

#%% rectangles
def partition2string(partition):
    string_list = ['/'.join([str(el) for el in patch]) for patch in partition]
    return ', '.join(string_list)
    
def make_rectangles(partitions, p1s):

    rectangles = DataFrame({
        'xmin': [x - 1/(2*40) for x in p1s],  
        'xmax': [x + 1/(2*40) for x in p1s],  
        'ymin': [0] * len(p1s),  
        'ymax': [1] * len(p1s),  
        'fill': [partition2string(el) for el in partitions]
    })
    
    partition_strings = ["1, 2, 3, 4", "1/2/3, 4", "1/2, 3/4", "1/2, 3, 4", "1/2/3/4"]
    rectangles['fill'] = Categorical(rectangles['fill'], categories = partition_strings, ordered = True)
     
    return rectangles

#%%
def make_datasets(pattern):
    files = glob.glob(pattern, recursive = True)
    mini_datasets = []
    p1_partitions = []
    extinction_probabilities = []
    for file in files:
        results = read_json(file)
        p1 = results['p1']
        centers = results['centers']
        partition = results['partition']
        mini_dataset = make_plotting_mini_dataframe(p1, centers, partition, 0.006)
        
        mini_datasets.append(mini_dataset)
        p1_partitions.append((p1, partition))
        extinction_probabilities.append(results['extinction_probability'])
        
    data = concat(mini_datasets, ignore_index = True)  
    
    p1s, partitions = zip(*p1_partitions)
    rectangles = make_rectangles(partitions, p1s)
    extinction_data = DataFrame({"p1": p1s, "extinction_probability": extinction_probabilities})
    
    
    return data, rectangles, extinction_data

def make_extinction_probability_dataset(pattern):
    files = glob.glob(pattern, recursive = True)
    extinction_probabilities = []
    partitions = []
    p1s = []
    for file in files:
        results = read_json(file)
# =============================================================================
#         print(results)
# =============================================================================
        for el in results:
            print(el)
            p1s.append(el['p1'])
            partitions.append(partition2string(el['partition']))
            extinction_probabilities.append(el['extinction_probability'])
        
    
    extinction_data = DataFrame({"p1": p1s, "extinction_probability": extinction_probabilities, "partition": partitions})
    
    
    return extinction_data

#%% make average_dataset
def make_plotting_mini_avg_dataframe(p1, centers, partition, pch_width):
    fecundity = len(centers)
    p1s = [p1] * fecundity
    partitions = [[str(i)] * len(el) for i, el in enumerate(partition)]
    partitions = unlist(partitions)
    offsets = calculate_partition_shift(partition, pch_width)
    
    return DataFrame({"p1": p1s, "center": centers, "patch": partitions, "shifted_p1": add(p1s, offsets)})

def make_avg_datasets(pattern):
    files = glob.glob(pattern, recursive = True)
    mini_datasets = []
    p1_partitions = []
    extinction_probabilities = []
    for file in files:
        results = read_json(file)
        p1 = results['p1']
        centers = results['centers']
        partition = results['partition']
        mini_dataset = make_plotting_mini_dataframe(p1, centers, partition, 0.006)
        
        mini_datasets.append(mini_dataset)
        p1_partitions.append((p1, partition))
        extinction_probabilities.append(results['extinction_probability'])
        
    data = concat(mini_datasets, ignore_index = True)  
    
    p1s, partitions = zip(*p1_partitions)
    rectangles = make_rectangles(partitions, p1s)
    extinction_data = DataFrame({"p1": p1s, "extinction_probability": extinction_probabilities})
    
    
    return data, rectangles, extinction_data
#%% plotting

def plot(data, rectangles):
    g = ggplot()
    g = g + ylim((0,1))
    g = g + theme(figure_size=(10, 4))
    
    g = g + geom_rect(data = rectangles, mapping = aes(xmin = 'xmin', xmax = 'xmax', ymin = 'ymin', ymax = 'ymax', fill = 'fill'), alpha = 1) #alpha = 0.5
    
    g = g + geom_point(data, aes(x = 'shifted_p1', y = 'center', shape = 'patch'), show_legend = False)
    
    g = g + xlab('$p_1$') + ylab('Center')
    g = g + labs(fill = "Distribution")
# =============================================================================
#     g = g + scale_fill_brewer(palette='Set1')
# =============================================================================
    
    g = g + scale_fill_hue()

    g = g + theme(
        axis_title = element_text(size = 14),
        axis_text = element_text(size = 12),
        legend_text = element_text(size = 12),
        legend_title = element_text(size = 14),
        plot_title = element_text(size = 16)
    )
    
    verts = rectangles.xmin.sort_values()
    verts = verts[1:]
    for el in verts:
        g = g + geom_segment(aes(x = el, xend = el, y = 0, yend = 1), color = "#B8B8B8")
    
# =============================================================================
#     my_colors = ["#e6a8a4", "#d5e6a4", "#a4e6c2", "#a4bae6", "#dda4e6"]
# =============================================================================
    my_colors = ["#FFA8A8", "#d5e6a4", "#a4e6c2", "#96E7FF", "#dda4e6"] #e6a8a4 #A1CDE6 #99FFFF
    g = g + scale_fill_manual(values = my_colors)
    
    print(g)
    return

def extinction_plot(data):
    g = ggplot(data, aes(x = "p1", y = "extinction_probability"))
    g = g + geom_point(aes(color = "partition"))
    print(g)
    return 

#%%
# pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/bimodal/alpha1_0.6__p1_*/output.json'
pattern = 'output/study/hierarchical/bimodal/seed_0/alpha1_0.6__p1_*/partition_output.json'

data, rectangles, extinction_data = make_datasets(pattern)
plot(data, rectangles)

pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/bimodal/alpha1_0.6__p1_*/partition_output.json'

extinction_data = make_extinction_probability_dataset(pattern)
extinction_plot(extinction_data)

#%%
pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/bimodal/alpha1_0.9__p1_*/output.json'

data, rectangles, extinction_data = make_datasets(pattern)
plot(data, rectangles)
plot_extinction_probability(extinction_data)

pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/bimodal/alpha1_0.9__p1_*/partition_output.json'

extinction_data = make_extinction_probability_dataset(pattern)
extinction_plot(extinction_data)

extinction_data = make_extinction_probability_dataset(pattern)
extinction_plot(extinction_data)

#%%
pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/central_mode__varying_width/alpha1_4__beta1_2__p1_*/output.json'

data, rectangles, extinction_data = make_datasets(pattern)
plot(data, rectangles)
plot_extinction_probability(extinction_data)

pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/central_mode__varying_width/alpha1_4__beta1_2__p1_*/partition_output.json'

extinction_data = make_extinction_probability_dataset(pattern)
extinction_plot(extinction_data)
#%%
pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/central_mode__varying_width/alpha1_7__beta1_3__p1_*/output.json'

data, rectangles, extinction_data = make_datasets(pattern)
plot(data, rectangles)
plot_extinction_probability(extinction_data)

pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/central_mode__varying_width/alpha1_7__beta1_3__p1_*/partition_output.json'

extinction_data = make_extinction_probability_dataset(pattern)
extinction_plot(extinction_data)

#%%
pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/central_mode__varying_width/alpha1_20__beta1_10__p1_*/output.json'

data, rectangles, extinction_data = make_datasets(pattern)
plot(data, rectangles)

pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/central_mode__varying_width/alpha1_20__beta1_10__p1_*/partition_output.json'

extinction_data = make_extinction_probability_dataset(pattern)
extinction_plot(extinction_data)

#%%
pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/central_mode__varying_width/alpha1_23__beta1_7__p1_*/output.json'

data, rectangles, extinction_data = make_datasets(pattern)
plot(data, rectangles)

pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/central_mode__varying_width/alpha1_23__beta1_7__p1_*/partition_output.json'

extinction_data = make_extinction_probability_dataset(pattern)
extinction_plot(extinction_data)

#%%
pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/extreme_mode__varying_mass_everywhere_else/beta1_0.05__p1_*/output.json'

data, rectangles, extinction_data = make_datasets(pattern)
plot(data, rectangles)

pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/extreme_mode__varying_mass_everywhere_else/beta1_0.05__p1_*/partition_output.json'

extinction_data = make_extinction_probability_dataset(pattern)
extinction_plot(extinction_data)

#%%
pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/extreme_mode__varying_mass_everywhere_else/beta1_0.95__p1_*/output.json'

data, rectangles, extinction_data = make_datasets(pattern)
plot(data, rectangles)

pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/extreme_mode__varying_mass_everywhere_else/beta1_0.95__p1_*/partition_output.json'

extinction_data = make_extinction_probability_dataset(pattern)
extinction_plot(extinction_data)

#%%
pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/extreme_mode__varying_slope/alpha1_1.8__p1_*/output.json'

data, rectangles, extinction_data = make_datasets(pattern)
plot(data, rectangles)

pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/extreme_mode__varying_slope/alpha1_1.8__p1_*/partition_output.json'

extinction_data = make_extinction_probability_dataset(pattern)
extinction_plot(extinction_data)

#%%
pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/extreme_mode__varying_slope/alpha1_2__p1_*/output.json'

data, rectangles, extinction_data = make_datasets(pattern)
plot(data, rectangles)

pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/extreme_mode__varying_slope/alpha1_2__p1_*/partition_output.json'

extinction_data = make_extinction_probability_dataset(pattern)
extinction_plot(extinction_data)

#%%
pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/extreme_mode__varying_slope/alpha1_3__p1_*/output.json'

data, rectangles, extinction_data = make_datasets(pattern)
plot(data, rectangles)

pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/extreme_mode__varying_slope/alpha1_3__p1_*/partition_output.json'

extinction_data = make_extinction_probability_dataset(pattern)
extinction_plot(extinction_data)

#%%
pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/extreme_mode__varying_slope/alpha1_4__p1_*/output.json'

data, rectangles, extinction_data = make_datasets(pattern)
plot(data, rectangles)

pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/extreme_mode__varying_slope/alpha1_4__p1_*/partition_output.json'

extinction_data = make_extinction_probability_dataset(pattern)
extinction_plot(extinction_data)


#%%
pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/extreme_mode__varying_slope/alpha1_5__p1_*/output.json'

data, rectangles, extinction_data = make_datasets(pattern)
plot(data, rectangles)

pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/extreme_mode__varying_slope/alpha1_5__p1_*/partition_output.json'

extinction_data = make_extinction_probability_dataset(pattern)
extinction_plot(extinction_data)

#%%
pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/extreme_mode__varying_slope/alpha1_7.5__p1_*/output.json'

data, rectangles, extinction_data = make_datasets(pattern)
plot(data, rectangles)

pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/extreme_mode__varying_slope/alpha1_7.5__p1_*/partition_output.json'

extinction_data = make_extinction_probability_dataset(pattern)
extinction_plot(extinction_data)


#%%
pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/extreme_mode__varying_slope/alpha1_10__p1_*/output.json'

data, rectangles, extinction_data = make_datasets(pattern)
plot(data, rectangles)

pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/extreme_mode__varying_slope/alpha1_10__p1_*/partition_output.json'

extinction_data = make_extinction_probability_dataset(pattern)
extinction_plot(extinction_data)


#%%


def add_fine(partition_results, raw):
    mean_dict = {key: np.mean(value) for key, value in raw.items()}
    
    for key, value in mean_dict.items():
        partition_results[key]["fine_extinction_probability"] = value
    

def get_best_partition(partition_results):
    best_extinction_probability = 2
    for partition, results in partition_results.items():
        if results["fine_extinction_probability"] < best_extinction_probability:
            best_extinction_probability = results["fine_extinction_probability"]
            best_partition_results = results
    return best_partition_results



#%%

def make_datasets(partition_pattern, raw_pattern):
    partition_files = glob.glob(partition_pattern, recursive = True)
    raw_files = glob.glob(raw_pattern, recursive = True)
    mini_datasets = []
    p1_partitions = []
    
    p1s = []
    partitions = []
    extinction_probabilities = []
    
    for partition_file, raw_file in zip(partition_files, raw_files):
        partition_results = read_json(partition_file)
        raw_results = read_json(raw_file)
        
        add_fine(partition_results, raw_results)
        
        best_results = get_best_partition(partition_results)
        
        p1 = best_results['p1']
        centers = best_results['centers']
        partition = best_results['partition']
        mini_dataset = make_plotting_mini_dataframe(p1, centers, partition, 0.006)
        
        for partition, results in partition_results.items():
            p1s.append(results['p1'])
            partitions.append(partition2string(results['partition']))
            extinction_probabilities.append(results['fine_extinction_probability'])
        
        
        mini_datasets.append(mini_dataset)
        p1_partitions.append((p1, best_results['partition']))
        #extinction_probabilities.append(results['fine_extinction_probability'])
        
    data = concat(mini_datasets, ignore_index = True)  
    
    rectangle_p1s, rectangle_partitions = zip(*p1_partitions)
    rectangles = make_rectangles(rectangle_partitions, rectangle_p1s)
    
    extinction_data = DataFrame({"p1": p1s, "extinction_probability": extinction_probabilities, "partition": partitions})
    #extinction_data = DataFrame({"p1": p1s, "extinction_probability": extinction_probabilities})
    
    return data, rectangles, extinction_data





#%%
base = 'output/repeated1/hierarchical/fecundity4/delta0.1/central_mode__varying_width/seed0/alpha1_20__beta1_10__p1_*/'

raw_pattern = base + 'robust_output.json'

partition_pattern = base + 'partition_output.json'

data, rectangles, extinction_data = make_datasets(partition_pattern, raw_pattern)
plot(data, rectangles)
extinction_plot(extinction_data)



#%%
base = 'output/repeated1/hierarchical/fecundity4/delta0.1/central_mode__varying_width/seed1/alpha1_20__beta1_10__p1_*/'

raw_pattern = base + 'robust_output.json'

partition_pattern = base + 'partition_output.json'

data, rectangles, extinction_data = make_datasets(partition_pattern, raw_pattern)
plot(data, rectangles)
extinction_plot(extinction_data)




#%%
base = 'output/repeated1/hierarchical/fecundity4/delta0.1/central_mode__varying_width/seed2/alpha1_20__beta1_10__p1_*/'

raw_pattern = base + 'robust_output.json'

partition_pattern = base + 'partition_output.json'

data, rectangles, extinction_data = make_datasets(partition_pattern, raw_pattern)
plot(data, rectangles)
extinction_plot(extinction_data)

#%%
base = 'output/repeated1/hierarchical/fecundity4/delta0.1/central_mode__varying_width/seed3/alpha1_20__beta1_10__p1_*/'

raw_pattern = base + 'robust_output.json'

partition_pattern = base + 'partition_output.json'

data, rectangles, extinction_data = make_datasets(partition_pattern, raw_pattern)
plot(data, rectangles)
extinction_plot(extinction_data)

#%%
base = 'output/repeated1/hierarchical/fecundity4/delta0.1/central_mode__varying_width/seed4/alpha1_20__beta1_10__p1_*/'

raw_pattern = base + 'robust_output.json'

partition_pattern = base + 'partition_output.json'

data, rectangles, extinction_data = make_datasets(partition_pattern, raw_pattern)
plot(data, rectangles)
extinction_plot(extinction_data)

#%%




#%%


#%%

#%%
base = 'output/repeatedLBFGS/hierarchical/fecundity4/delta0.1/central_mode__varying_width/seed0/alpha1_20__beta1_10__p1_*/'

raw_pattern = base + 'robust_output.json'

partition_pattern = base + 'partition_output.json'

data, rectangles, extinction_data = make_datasets(partition_pattern, raw_pattern)
plot(data, rectangles)
extinction_plot(extinction_data)


#%%
base = 'output/repeatedLBFGS/hierarchical/fecundity4/delta0.1/central_mode__varying_width/seed1/alpha1_20__beta1_10__p1_*/'

raw_pattern = base + 'robust_output.json'

partition_pattern = base + 'partition_output.json'

data, rectangles, extinction_data = make_datasets(partition_pattern, raw_pattern)
plot(data, rectangles)
extinction_plot(extinction_data)

#%%
base = 'output/repeatedLBFGS/hierarchical/fecundity4/delta0.1/central_mode__varying_width/seed2/alpha1_20__beta1_10__p1_*/'

raw_pattern = base + 'robust_output.json'

partition_pattern = base + 'partition_output.json'

data, rectangles, extinction_data = make_datasets(partition_pattern, raw_pattern)
plot(data, rectangles)
extinction_plot(extinction_data)

#%%
base = 'output/repeatedLBFGS/hierarchical/fecundity4/delta0.1/central_mode__varying_width/seed3/alpha1_20__beta1_10__p1_*/'

raw_pattern = base + 'robust_output.json'

partition_pattern = base + 'partition_output.json'

data, rectangles, extinction_data = make_datasets(partition_pattern, raw_pattern)
plot(data, rectangles)
extinction_plot(extinction_data)

#%%
base = 'output/repeatedLBFGS/hierarchical/fecundity4/delta0.1/central_mode__varying_width/seed4/alpha1_20__beta1_10__p1_*/'

raw_pattern = base + 'robust_output.json'

partition_pattern = base + 'partition_output.json'

data, rectangles, extinction_data = make_datasets(partition_pattern, raw_pattern)
plot(data, rectangles)
extinction_plot(extinction_data)




#%%%




#%%
def individual_merge(master, rectangles_for_a_seed, seed):
    rectangles_add = rectangles_for_a_seed[["xmin", "fill"]]
    rectangles_add.rename(columns={"fill": f"fill_{seed}"}, inplace = True)
    master = master.merge(rectangles_add, on = "xmin", how = "left")
    
    return master
    
def merge(rectangles_by_seed):
    
    rectangles = rectangles_by_seed[0]
    rectangles.rename(columns = {"fill": "fill_0"}, inplace = True)
    
    for i, el in enumerate(rectangles_by_seed):
        if i == 0:
            pass
        else:
            rectangles = individual_merge(rectangles, el, i)
    
    return rectangles    

# =============================================================================
# def merge(rectangles_by_seed):
#     
#     rectangles = rectangles_by_seed[0]
#     rectangles.rename(columns = {"fill": "fill_0"}, inplace = True)
# 
#     rectangles_add = rectangles_by_seed[1][["xmin", "fill"]]
#     rectangles_add.rename(columns={"fill": "fill_1"}, inplace = True)
#     rectangles = rectangles.merge(rectangles_add, on = "xmin", how = "left")
# 
#     rectangles_add = rectangles_by_seed[2][["xmin", "fill"]]
#     rectangles.rename(columns={"fill": "fill_2"}, inplace = True)
#     rectangles = rectangles.merge(rectangles_add, on = "xmin", how = "left")
#     
#     #rectangles = rectangles.merge(rectangles_add, on = "xmin", how = "left", suffixes = ('', '_2'))
#     
#     return rectangles
# =============================================================================

def has_majority(row):
    counter = Counter(row)
    return counter.most_common(1)[0][1] > row.shape[0]/2

def majority_vote(row):
    counter = Counter(row)
    return counter.most_common(1)[0][0]

def first_argmin(row):
    print(row)
    if row["fill_0"] == row["mode"]:
        return "fill_0"
    elif row["fill_1"] == row["mode"]:
        return "fill_1"
    else: 
        return "fill_2"
    
def first_argmin_idx(row):
    if row["fill_0"] == row["mode"]:
        return 0
    elif row["fill_1"] == row["mode"]:
        return 1
    else: 
        return 2

def get_multi_rectangles(rectangles_by_seed):
    
    num_seeds = len(rectangles_by_seed)
    
    multi_rectangles = merge(rectangles_by_seed)
    
    last_three_columns = multi_rectangles.iloc[:, -num_seeds:]
    
    multi_rectangles["has_mode"] = last_three_columns.apply(has_majority, axis = 1)
    multi_rectangles["mode"] = last_three_columns.apply(majority_vote, axis = 1)
    multi_rectangles["argmin"] = multi_rectangles.apply(first_argmin, axis = 1)
    multi_rectangles["argmin_idx"] = multi_rectangles.apply(first_argmin_idx, axis = 1)
    return multi_rectangles

def process(coarse_pattern, fine_pattern, num_seeds):
    rectangles_by_seed = list()
    data_by_seed = list()
    extinction_by_seed = list()
    
    for i in range(num_seeds):
        
# =============================================================================
#         base = f'output/LBFGS/hierarchical/fecundity4/delta0.1/{coarse_pattern}/seed_{i}/{fine_pattern}/'
# =============================================================================
        base = f'output/hierarchical1/fecundity4/delta0.1/{coarse_pattern}/seed_{i}/{fine_pattern}/'

        raw_pattern = base + 'robust_output.json'
        
        partition_pattern = base + 'partition_output.json'
        
        data, rectangles, extinction_data = make_datasets(partition_pattern, raw_pattern)
        
        rectangles_by_seed.append(rectangles)
        data_by_seed.append(data)
        extinction_by_seed.append(extinction_data)
        
    return rectangles_by_seed, data_by_seed, extinction_by_seed


def retrieve_best_dataset(row, data_by_seed):
    xmin = row.xmin
    xmax = row.xmax
    p1 = (xmin + xmax)/2
    best_idx = row.argmin_idx
    data = data_by_seed[best_idx]
    return data[data.p1 == p1]
    
    
def get_overall_dataset(data_by_seed, rectangles):
    
    mini_data_list = [retrieve_best_dataset(rectangles.loc[i], data_by_seed) for i in range(rectangles.shape[0])]
    
    data = concat(mini_data_list)
    
    return data
   
#%%
def plot(data, rectangles):
    g = ggplot()
    g = g + ylim((0,1))
    g = g + theme(figure_size=(10, 4))
    
    g = g + geom_rect(data = rectangles, mapping = aes(xmin = 'xmin', xmax = 'xmax', ymin = 'ymin', ymax = 'ymax', fill = 'mode'), alpha = 1) #alpha = 0.5
    
    g = g + geom_point(data, aes(x = 'shifted_p1', y = 'center', shape = 'patch'), show_legend = False)
    
    g = g + xlab('$p$') + ylab('Center')
    g = g + labs(fill = "Partition")
# =============================================================================
#     g = g + scale_fill_brewer(palette='Set1')
# =============================================================================
    
    g = g + scale_fill_hue()
    
    scaler = 1.5
    g = g + theme(
        axis_title = element_text(size = 14 * scaler),
        axis_text = element_text(size = 12 * scaler),
        legend_text = element_text(size = 12 * scaler),
        legend_title = element_text(size = 14 * scaler),
        plot_title = element_text(size = 16 * scaler, hjust = 0.5)
    )
    
    verts = rectangles.xmin.sort_values()
    verts = verts[1:]
    for el in verts:
        g = g + geom_segment(aes(x = el, xend = el, y = 0, yend = 1), color = "#B8B8B8")
    
# =============================================================================
#     my_colors = ["#e6a8a4", "#d5e6a4", "#a4e6c2", "#a4bae6", "#dda4e6"]
# =============================================================================
# =============================================================================
#     my_colors = ["#FFA8A8", "#d5e6a4", "#a4e6c2", "#96E7FF", "#dda4e6"] #e6a8a4 #A1CDE6 #99FFFF
# =============================================================================
   
    partition2color = {"1, 2, 3, 4": "#FFA8A8",
                       "1/2, 3, 4": "#d5e6a4",
                       "1/2, 3/4": "#a4e6c2",
                       "1/2/3, 4": "#96E7FF",
                       "1/2/3/4": "#dda4e6"}
    
    
    
    g = g + scale_fill_manual(values = partition2color)
    
    
    print(g)
    return g

def extinction_plot(data):
    g = ggplot(data, aes(x = "p1", y = "extinction_probability"))
    g = g + geom_point(aes(color = "partition"))
    
    partition2color = {"1, 2, 3, 4": "#FFA8A8",
                       "1/2, 3, 4": "#d5e6a4",
                       "1/2, 3/4": "#a4e6c2",
                       "1/2/3, 4": "#96E7FF",
                       "1/2/3/4": "#dda4e6"}
    
    
    g = g + scale_fill_manual(values = partition2color)
    
    base_size = 16
    g = g + theme(text = element_text(size = base_size),
                  axis_title=element_text(size=base_size * 1.2),
                  axis_text=element_text(size=base_size),
                  legend_title=element_text(size=base_size * 1.1),
                  legend_text=element_text(size=base_size * 0.9),
                  plot_title=element_text(size=base_size * 1.5),
                  strip_text=element_text(size=base_size * 1.2)
                  )
    
    print(g)
    return      
#%%


    

#%%
coarse_pattern = "bimodal"
fine_pattern = 'alpha1_0.6__p1_*'

rectangles_by_seed, data_by_seed, extinction_by_seed = process(coarse_pattern, fine_pattern, 5)

rectangles = get_multi_rectangles(rectangles_by_seed)

data = get_overall_dataset(data_by_seed, rectangles)

g = plot(data, rectangles)
extinction_plot(extinction_by_seed[0])

g = g + labs(title = "$(α = 0.6, β = 0.4)$")
print(g)


print(all(rectangles.has_mode))


# Analysis: When both distribtution common, two patches with both extreme specialists. As the right distribution increases in probability it transitions to one double-specialist patch, and two patches with right specialists. As p1 continues increasing, we transition to four right-specialists.

#%%
coarse_pattern = "bimodal"
fine_pattern = 'alpha1_0.9__p1_*'

rectangles_by_seed, data_by_seed, extinction_by_seed = process(coarse_pattern, fine_pattern, 5)

rectangles = get_multi_rectangles(rectangles_by_seed)

data = get_overall_dataset(data_by_seed, rectangles)

g = plot(data, rectangles)
extinction_plot(extinction_by_seed[0])

g = g + labs(title = "$(α = 0.9, β = 0.1)$")
print(g)

print(all(rectangles.has_mode))

# Analysis: Same pattern as above
#%%





#%%
coarse_pattern = "central_mode__varying_width"
fine_pattern = 'alpha1_4__beta1_2__p1_*'

rectangles_by_seed, data_by_seed, extinction_by_seed = process(coarse_pattern, fine_pattern, 7)

rectangles = get_multi_rectangles(rectangles_by_seed)

data = get_overall_dataset(data_by_seed, rectangles)

g = plot(data, rectangles)
extinction_plot(extinction_by_seed[0])

g = g + labs(title = "$(α = 4, β = 2)$")
print(g)

print(all(rectangles.has_mode))

# Analysis: Two patches of double-specialists. Transitions to a patch of triple specialists and a single moderate specialist

#%%
coarse_pattern = "central_mode__varying_width"
fine_pattern = 'alpha1_7__beta1_3__p1_*'

rectangles_by_seed, data_by_seed, extinction_by_seed = process(coarse_pattern, fine_pattern, 11)

rectangles = get_multi_rectangles(rectangles_by_seed)

data = get_overall_dataset(data_by_seed, rectangles)

g = plot(data, rectangles)
extinction_plot(extinction_by_seed[0])

g = g + labs(title = "$(α = 7, β = 3)$")
print(g)

print(all(rectangles.has_mode)) #!!!!!

# Analysis: Two patches of double-specialists. Transitions to one double-specialist patch and two single patch specialists. transitions to one 3-specialist patch and one single patch specialist.

#%%
coarse_pattern = "central_mode__varying_width"
fine_pattern = 'alpha1_20__beta1_10__p1_*'

rectangles_by_seed, data_by_seed, extinction_by_seed = process(coarse_pattern, fine_pattern, 7)

rectangles = get_multi_rectangles(rectangles_by_seed)

data = get_overall_dataset(data_by_seed, rectangles)

g = plot(data, rectangles)
extinction_plot(extinction_by_seed[0])

g = g + labs(title = "$(α = 20, β = 10)$")
print(g)

print(all(rectangles.has_mode))

# Analysis: One patch with four specialists.

#%%
coarse_pattern = "central_mode__varying_width"
fine_pattern = 'alpha1_23__beta1_7__p1_*'

rectangles_by_seed, data_by_seed, extinction_by_seed = process(coarse_pattern, fine_pattern, 11)

rectangles = get_multi_rectangles(rectangles_by_seed)

data = get_overall_dataset(data_by_seed, rectangles)

g = plot(data, rectangles)
extinction_plot(extinction_by_seed[0])

g = g + labs(title = "$(α = 23, β = 7)$")
print(g)


print(all(rectangles.has_mode)) #!!!!!!

# Analysis: One three-specialist patch, but with two of them focusing on the infrequent distribution, and one single-specialist alone. transitions to 

#%%


#%%
coarse_pattern = "extreme_mode__varying_mass_everywhere_else"
fine_pattern = 'beta1_0.95__p1_*'

rectangles_by_seed, data_by_seed, extinction_by_seed = process(coarse_pattern, fine_pattern, 5)

rectangles = get_multi_rectangles(rectangles_by_seed)

data = get_overall_dataset(data_by_seed, rectangles)

g = plot(data, rectangles)
extinction_plot(extinction_by_seed[0])

g = g + labs(title = "$(α = 1, β = 0.95)$")
print(g)

print(all(rectangles.has_mode))

# Analysis: similar to bimodal

#%%
coarse_pattern = "extreme_mode__varying_mass_everywhere_else"
fine_pattern = 'beta1_0.05__p1_*'

rectangles_by_seed, data_by_seed, extinction_by_seed = process(coarse_pattern, fine_pattern, 5)

rectangles = get_multi_rectangles(rectangles_by_seed)

data = get_overall_dataset(data_by_seed, rectangles)

g = plot(data, rectangles)
extinction_plot(extinction_by_seed[0])

g = g + labs(title = "$(α = 1, β = 0.05)$")
print(g)

print(all(rectangles.has_mode))

# Analysis: similar to bimodal

#%%



#%%
coarse_pattern = "extreme_mode__varying_slope"
fine_pattern = 'alpha1_1.8__p1_*'

rectangles_by_seed, data_by_seed, extinction_by_seed = process(coarse_pattern, fine_pattern, 7)

rectangles = get_multi_rectangles(rectangles_by_seed)

data = get_overall_dataset(data_by_seed, rectangles)

g = plot(data, rectangles)
extinction_plot(extinction_by_seed[0])

g = g + labs(title = "$(α = 1.8, β = 1)$")
print(g)

print(all(rectangles.has_mode))

# Analysis: two patches where each patch has two right-specialists, but staggered to cover more support. Transitions to splitting up one of those patches into two far-right specialists.

#%%
coarse_pattern = "extreme_mode__varying_slope"
fine_pattern = 'alpha1_2__p1_*'

rectangles_by_seed, data_by_seed, extinction_by_seed = process(coarse_pattern, fine_pattern, 7)

rectangles = get_multi_rectangles(rectangles_by_seed)

data = get_overall_dataset(data_by_seed, rectangles)

g = plot(data, rectangles)
extinction_plot(extinction_by_seed[0])

g = g + labs(title = "$(α = 2, β = 1)$")
print(g)

print(all(rectangles.has_mode))

# Analysis: All two-specialist patches as above

#%%
coarse_pattern = "extreme_mode__varying_slope"
fine_pattern = 'alpha1_3__p1_*'

rectangles_by_seed, data_by_seed, extinction_by_seed = process(coarse_pattern, fine_pattern, 7)

rectangles = get_multi_rectangles(rectangles_by_seed)

data = get_overall_dataset(data_by_seed, rectangles)

g = plot(data, rectangles)
extinction_plot(extinction_by_seed[0])

g = g + labs(title = "$(α = 3, β = 1)$")
print(g)


print(all(rectangles.has_mode))

# Analysis: three patches, each with far-right specialists and one patch also has a left specialist. Transitions to two patches, three spaced right-specialists and one right specialist.

#%%
coarse_pattern = "extreme_mode__varying_slope"
fine_pattern = 'alpha1_4__p1_*'

rectangles_by_seed, data_by_seed, extinction_by_seed = process(coarse_pattern, fine_pattern, 7)

rectangles = get_multi_rectangles(rectangles_by_seed)

data = get_overall_dataset(data_by_seed, rectangles)

g = plot(data, rectangles)
extinction_plot(extinction_by_seed[0])

g = g + labs(title = "$(α = 4, β = 1)$")
print(g)


print(all(rectangles.has_mode))


# Analysis: Two patches of both-extreme specialists. Transitions to splitting one patch to two patches of one right-specialist. Finally transitions to one patch of three staggered right-specialists, and one right-specialist.

#%%
coarse_pattern = "extreme_mode__varying_slope"
fine_pattern = 'alpha1_5__p1_*'

rectangles_by_seed, data_by_seed, extinction_by_seed = process(coarse_pattern, fine_pattern, 7)

rectangles = get_multi_rectangles(rectangles_by_seed)

data = get_overall_dataset(data_by_seed, rectangles)

g = plot(data, rectangles)
extinction_plot(extinction_by_seed[0])

g = g + labs(title = "$(α = 5, β = 1)$")
print(g)



print(all(rectangles.has_mode))

# Analysis: Same as above.

#%%
coarse_pattern = "extreme_mode__varying_slope"
fine_pattern = 'alpha1_7.5__p1_*'

rectangles_by_seed, data_by_seed, extinction_by_seed = process(coarse_pattern, fine_pattern, 7)

rectangles = get_multi_rectangles(rectangles_by_seed)

data = get_overall_dataset(data_by_seed, rectangles)

g = plot(data, rectangles)
extinction_plot(extinction_by_seed[0])

g = g + labs(title = "$(α = 7.5, β = 1)$")
print(g)



print(all(rectangles.has_mode))

# Analysis: same as above

#%%
coarse_pattern = "extreme_mode__varying_slope"
fine_pattern = 'alpha1_10__p1_*'

rectangles_by_seed, data_by_seed, extinction_by_seed = process(coarse_pattern, fine_pattern, 9)

rectangles = get_multi_rectangles(rectangles_by_seed)

data = get_overall_dataset(data_by_seed, rectangles)

g = plot(data, rectangles)
extinction_plot(extinction_by_seed[0])

g = g + labs(title = "$(α = 10, β = 1)$")
print(g)


print(all(rectangles.has_mode))

# Analysis: same as above

#%%
a = Series([0.8623, 0.8875], index = ["xmin", "xmax"])





#%%

for el in df:
    print(el)


#%%

def f(row, y):
    return row.probability + y




#%%

def plot_extinction_probability(extinction_data):
    
    q = qplot("p1", "extinction_probability", extinction_data)
    print(q)
    return

#%% raw
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multicomp import tukeyhsd


def do_tukey(data, alpha):
    
    partitions = []
    extinction_probabilities = []
    for i, el in enumerate(data):
        partitions += [partition2string(possible_partitions[i])] * len(el)
        extinction_probabilities += el
        
    data = {'partition': partitions, 'probability': extinction_probabilities}
    
    df = DataFrame(data)
    
    tukey_result = pairwise_tukeyhsd(endog = df['probability'], groups = df['partition'], alpha = alpha)
    
    

            
            
    
    
    data_df = DataFrame(data_dict)
    
    
    pairwise_tukeyhsd()
    return 
def does_overlap_workhorse(data):
    data = array(data)
    mean = np.mean(data, axis = 1)
    sterr = sem(data, axis = 1)


    (mean[2] - sterr[2], mean[2] + sterr[2])
    (mean[3] - sterr[3], mean[3] + sterr[3])
    return


#%% running raw
pattern = 'output/finer3/hierarchical/fecundity4/delta0.1/bimodal/alpha1_0.6__p1_*/robust_output.json'

files = glob.glob(pattern, recursive = True)
raw_results = [read_json(el) for el in files]



