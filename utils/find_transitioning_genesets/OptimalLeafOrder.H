#ifndef _OPTIMALLEAFORDER_
#define _OPTIMALLEAFORDER_
#include <map>
#include <string>
#include <vector>
using namespace std;
class HierarchicalClusterNode;
class OptimalLeafOrder
{
	public:
	OptimalLeafOrder();
	~OptimalLeafOrder();
	int setHierarchicalClusterNode(HierarchicalClusterNode*);
	int reorder(vector<string>&);
	int setDist(double**);
	private:
	typedef struct
	{
		int leftextreme;
		int rightextreme;
	}Pair;
	struct Node;
	struct Node
	{
		int id;
		string nodeName;
		map<int,map<int,double>*> childOrder;
		map<string,Pair*> extremes;
		Node* left;
		Node* right;
		Node* parent;
	};
	Node* root;
	map<string,Node*> minDist;
	int getordering(Node*,int,int,vector<string>&);
	int populateDist(HierarchicalClusterNode* node);
	double getSim(string&,string&);
	double getSim(int,int);
	int makeTree(HierarchicalClusterNode*,Node*);
	int reorder(Node* n);
	map<string,double> distance;
	HierarchicalClusterNode* hcRoot;
	double** dist;
	map<int,string> nodeIDNameMap;
};
#endif
