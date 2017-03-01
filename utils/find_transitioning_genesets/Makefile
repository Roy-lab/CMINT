LFLAG = -lgsl -lgslcblas 
SRC = Distance.C  Framework.C  GeneExpManager.C  Heap.C  HierarchicalCluster.C  HierarchicalClusterNode.C  OptimalLeafOrder.C MappedOrthogroupReader.C MappedOrthogroup.C GeneMap.C

CC=g++
CFLAGS = -g

findTransitionGenesets: $(SRC)
	$(CC) $(SRC)  $(CFLAGS) -o findTransitionGenesets

