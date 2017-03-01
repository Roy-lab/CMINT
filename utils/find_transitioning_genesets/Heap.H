#ifndef _HEAP_
#define _HEAP_
#include <string>
using namespace std;
class Heap
{
	public:
	Heap();
	~Heap();
	struct Pair{
		int node1;
		int node2;
		double dist;
		int id;
	};
	Heap* insertToHeapNoHeapify(int n1, int n, double d);
	Heap* insertToHeapNoHeapify(int n1, int n, double d,int id);
	Heap* findNewPosition(Heap* currPos,Heap* findMe,int&);
	Heap* findNewPosition_iterative(Heap* currPos);
	int deleteFromHeap_getLeaf(Heap*);
	Heap* getLeaf(Heap* currPos);
	int heapifyDown(Heap* currPos);
	bool checkHeap();
	bool checkHeap(Heap*);
	bool checkPointers();
	bool checkPointers(Heap*);
	Heap::Pair* getMin();
	Heap* getRoot();
	int showHeap();
	int showHeap(Heap*);
	bool empty();
	Heap* left;
	Heap* right;
	Heap* root;
	Heap* parent;
	Pair p;
		
};
#endif
