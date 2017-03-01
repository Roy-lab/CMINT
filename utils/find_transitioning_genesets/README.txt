This directory has one C++ executable and two scripts to make figures.
1. findTransitionGenesets is used to generate input files for making heatmap figures.
2. makeHeatmaps_array.sh and makeHeatmaps_seq.sh are scripts that operate on the output of findTransitionGenesets to make .svg figures for each set of genes and a bird's eye figure with the average pattern figure for all
gene sets. 

findTransitionGenesets
----------------------
This program will take the merged cluster assignments and cluster profiles and produce multiple files, one for each set of genes with a specific pattern. The size of the gene set
is user defined to put a lower limit on the number of genes. It will produce two additional files, one is called all_assign.txt which has all the genes grouped by the specific patterns
they exhibit (including gene sets including 1 gene). The second file is called ordered_clusterset_means.txt that will have the mean cluster assignment profile and the mark profile for all gene sets
with more than the "minimum" number of genes. All these files are inputs to the makeHeatmap_array.sh or makeHeatmap_seq.sh files.

./findTransitionGenesets cmint_result_dir celltypeorder OGID_file srccelltype threshold outputdir mingenesetsize
cmint_result_dir: This is the location of the CMINT output directory
celltypeorder: This is the same input as for CMINT and specifies the order of the cell types in the OGID file
OGID_file: This file is the same input as for CMINT and has the names of regions, one line per region.
srccelltype: This is the root of the tree or the species with respect to which the cluster assignment files are created. Same as in CMINT
threshold: This is used in the hierarchical clustering to define clusters of regions with a specific pattern.
outputdir: This directory stores the result of running findTRansitionGenesets
mingenesetsize: This is the minimum number of genes a pattern should have to display its pattern in a separate clusterset file. The average pattern of this set is also outputted in the ordered_clusterset_means.txt file

Example usage with the reprogramming data is:
./findTransitionGenesets ../../outputs/cmint_modules_reprogramming/ ../../data/reprogramming/celltypeorder_reprogramming.txt ../../data/reprogramming/ogids_notfilterbyexp.txt ips 0.05 example_out2 10
Example usage with the seq data is:
./findTransitionGenesets  ../../outputs/cmint_modules_hemato_20k/ ../../data/hematopoesis/hemato_lineage.txt ../../data/hematopoesis/OGIDs_HematoLineage.txt LT 0.05 example_out3 10

Note 1:  when generating the average cluster profile for a cell type, there might be ties (we display the majority of the cluster assignment, but it is possible that for a cell type for those regions, there is no clear majority.). Then you will see a message: "Oops found a tie." This is a harmless message. The tie is resolved by using the max of the cluster ids.

Note 2: The directories in ../../outputs are tarred up right now. Please cd into ../../outputs to first untar to execute the above commands properly. 

makeHeatmaps_seq.sh and makeHeatmaps_array.sh
------------------------------------------------
These two shell scripts will use the files in the output directory of findTransitionGenesets and make heatmaps. The makeHeatmaps_seq.sh is for sequencing data and uses a sequential colormap. The min and max limits are set
for the CMINT results on the 20k all-observed regions for the hematopoetis analysis. To adjust, please change the COLORMAP variable in the script. The makeHeatmaps_seq.sh also does not output a figure for all_assign.svg
as this tends to become very big. This line is commented out. Please uncomment if you wish to make this figure.

The makeHeatmaps_array.sh is used for array data and uses a diverging scale for making the heatmaps.

Both scripts take  three arguments: outputdir, numberofclusters, heatmap_script_path
outputdir: this is the output location directory of findTransistionGenesets. All figures will be stored in a subdirectory called figs in outputdir
numberofclusters: this is the number of clusters CMINT outputs
heatmap_script_path: This is the location of the heatmap.awk script. We provide a version as part of the CMINT pacakge. This Heatmap.awk script is written by Pouya Kheradpour and was downloaded from http://compbio.mit.edu/pouyak/software/Heatmap.awk

Example usages
-------------
Example usage with makeHeatmaps_array.sh with the reprogramming data assuming we ran the findTransitionGenesets as described above:
./makeHeatmaps_array.sh example_out2 15 ../Heatmap.awk

Example usage with makeHeatmaps_seq.sh with the hematopoeisis data
./makeHeatmaps_seq.sh example_out3 16 ../Heatmap.awk


Acknowledgements
-----------------
We thank Pouya Kheradpour for writing the Heatmap.awk script. 


