OVERVIEW
--------
cmint_array is used for clustering chromatin mark profiles from array data where we do not have a lot of zeros and we can estimate the variance.
cmint_seq is used for clustering chromatin mark profiles from sequencing profiles. The input is assumed to be log transformed and normalized.
We do not estimate the variance within cmint_seq. cmint_seq also is more memory efficient compared to cmint_array. Both cmint_array and cmint_seq
have similar inputs and outputs.

USING CMINT
-----------
The cmint_array/cmint_seq executables take 10 different inputs can be run in the following manner.

./cmint_array celltype_order genegroup maxk celllineage clusterassignments rand[rseed|none] outputDir mode[learn|generate] srcnode inittype[uniform|branchlength] p_diagonal_nonleaf
The input arguments are as follows:

1. celltype_order: A file describing the order of the cell types and is needed to parse the genegroup file. Example: ../data/reprogramming/specorder.txt

2. genegroup: A file has the groups of genes/regions. The format of this file is the name of the group followed by a comma-separated unique identifier of a region in the cell type. The comma separates different celltypes and order of the cell types is specified by celltype_order. 
Example: ../data/ogids_notfilterexp.txt

3. maxk: The maximum number of clusters

4. celllineage: A file describing the tree of cell types. The format of this file is 2 column, tab-separated. First column is the name of the child and the second column is the name of the parent of each branch. Example: ../data/reprogramming/celltype_tree3_ancestor.txt

5. config : A file describing the locations of the initial cluster assignments and the mark data. 
The cluster assignment files can be used to specify which regions should be used for the CMINT algorithm. 
Example: ../data/reprogramming/config_k15.txt

6. outputDir: Location of results

7. model: run cmint in learn or sampling mode

8. srcnode: The celltype specifying a reference point with which the rows in some output files are specified.

9. inittype: An input argument specifying how the transition matrices should be initialized

10. p_diagonal_nonleaf: If inittype is uniform, p_diagonal_nonleaf is a number between 0-1 specifying the probability with which a region/gene maintains its module assignment. If inittype is branchlength, then p_diagonal_nonleaf is a file specifying the transition probabilities. WARNING, the option tested with CMINT is when inittype is uniform.

Example of using CMINT with the reprogramming data
--------------------------------------------------
The command line below assumes we are in the cmint_array directory.

./cmint_array $DATADIR/specorder.txt $DATADIR/ogids_notfilterbyexp.txt 15 $DATADIR/celltype_tree3_ancestor.txt $DATADIR/config_k15.txt none $RESULTS  learn ips uniform 0.8

DATADIR is set to ../../data/reprogramming/
RESULTS is set to ../../../results/reprogramming
This usage will run cmint on the reprogramming data with k=15 clusters. Results will be stored in the ../../results/reprogramming directory.


Example of using CMINT with the hematopoesis data
--------------------------------------------------
The command line below assumes we are in the cmint_seq directory. 
We have two sets of configuation files: one that is for the 1 million regions
and one that is for the 20K regions.

For the 20k regions, run cmint as follows:

./cmint_seq $DATADIR/hemato_lineage.txt $DATADIR/OGIDs_HematoLineage.txt 16 $DATADIR/hemato_lineage_tree.txt $DATADIR/config_k16_20k.txt none $RESULTS  learn  LT uniform 0.8
DATADIR is set to ../../data/hematopoesis/
RESULTS is set to ../../results/20k
NOTE: Running this command with the above inputs will need ~9GB of memory. The above command will produces results for the 20K regions for which we have non-zero measurements of all marks in all cell types. Results will be stored in ../../results/20k

For the 1 million regions for all 15 cell types, run cmint as follows:
./cmint_seq $DATADIR/hemato_lineage.txt $DATADIR/OGIDs_HematoLineage.txt 16 $DATADIR/hemato_lineage_tree.txt $DATADIR/config_k16_1million.txt none $RESULTS  learn  LT uniform 0.8
DATADIR is set to ../../data/hematopoesis/
RESULTS is set to ../../results/1million
NOTE: Running this command with the above inputs will need 90GB of memory. 
The high memory requirement is because of the large (15) number of cell types. The above command will produces results for all the 1million regions for which we have non-zero measurements of a mark in at least one cell type. Results will be stored in ../../results/1million


