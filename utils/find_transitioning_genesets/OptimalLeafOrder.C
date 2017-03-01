#include <string.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include "HierarchicalClusterNode.H"
#include "OptimalLeafOrder.H"


OptimalLeafOrder::OptimalLeafOrder()
{
}

OptimalLeafOrder::~OptimalLeafOrder()
{
}

int 
OptimalLeafOrder::setHierarchicalClusterNode(HierarchicalClusterNode* aPtr)
{
	hcRoot=aPtr;
	return 0;
}

int
OptimalLeafOrder::reorder(vector<string>& ordering)
{
	root=new OptimalLeafOrder::Node;
	makeTree(hcRoot,root);
	if(root->left==NULL && root->right==NULL)
	{
		ordering.push_back(root->nodeName);
		return 0;
	}
	populateDist(hcRoot);
	reorder(root); 
	//Now get the highest scoring left and right child
	double min=10E7;
	string leftmost;
	string rightmost;
	int leftmostid=-1;
	int rightmostid=-1;
	for(map<int,map<int,double>*>::iterator cIter=root->childOrder.begin();cIter!=root->childOrder.end();cIter++)
	{
		map<int,double>* vals=cIter->second;
		for(map<int,double>::iterator vIter=vals->begin();vIter!=vals->end();vIter++)
		{
			//cout << cIter->first <<"\t" << vIter->first <<"\t" <<vIter->second << endl;
			if(vIter->second<min)
			{		
				min=vIter->second;
				leftmostid=cIter->first;
				rightmostid=vIter->first;
			}
		}
	}
	//cout<< "Mindist " <<min << endl;
	string key;
	char keychar[256];
	if(leftmostid<=rightmostid)
	{
		sprintf(keychar,"%d-%d",leftmostid,rightmostid);
	}
	else
	{
		sprintf(keychar,"%d-%d",rightmostid,leftmostid);
	}
	key.append(keychar);
	OptimalLeafOrder::Pair* p=root->extremes[key];
	if(root->left!=NULL && root->right!=NULL)
	{
		if(root->left->childOrder.find(leftmostid)!=root->left->childOrder.end())	
		{
			getordering(root->left,leftmostid,p->leftextreme,ordering);
			getordering(root->right,p->rightextreme,rightmostid,ordering);
		}
		//switch case
		else if(root->right->childOrder.find(leftmostid)!=root->right->childOrder.end())	
		{
			getordering(root->right,leftmostid,p->rightextreme,ordering);
			getordering(root->left,p->leftextreme,rightmostid,ordering);
		}
		double d=0;
		for(int i=0;i<ordering.size()-1;i++)
		{
			//d=d+getSim(ordering[i],ordering[i+1]);
		}
		//cout << "distance: " << d << endl;
	}
	return 0;
}
int
OptimalLeafOrder::setDist(double** d)
{
	dist=d;
	return 0;
}


int
OptimalLeafOrder::populateDist(HierarchicalClusterNode* node)
{
	if(node->left==NULL && node->right==NULL)
	{
		for(map<int,double>::iterator nIter=node->distToNeighbors.begin();nIter!=node->distToNeighbors.end();nIter++)
		{
			char key[1024];
			if(node->id<nIter->first)
			{
				sprintf(key,"%d-%d",node->id,nIter->first);
			}
			else
			{
				sprintf(key,"%d-%d",nIter->first,node->id);
			}
			if(nIter->second>10E7)
			{
				cout <<"Distance overflow error for "<<  key << " "<< nIter->second <<  endl;
				exit(0);
			}
			string keystr(key);
			distance[keystr]=nIter->second;
		}
	}
	if(node->left!=NULL)
	{
		populateDist(node->left);
	}
	if(node->right!=NULL)
	{
		populateDist(node->right);
	}
	
	return 0;
}

int
OptimalLeafOrder::makeTree(HierarchicalClusterNode* node,OptimalLeafOrder::Node* ooNode)
{
	ooNode->left=NULL;
	ooNode->right=NULL;
	ooNode->parent=NULL;
	ooNode->nodeName.append(node->nodeName);
	ooNode->id=node->id;
	nodeIDNameMap[ooNode->id]=node->nodeName;
	if(node->left!=NULL)
	{
		OptimalLeafOrder::Node* leftchild=new OptimalLeafOrder::Node;
		ooNode->left=leftchild;
		makeTree(node->left,leftchild);
	}
	if(node->right!=NULL)
	{
		OptimalLeafOrder::Node* rightchild=new OptimalLeafOrder::Node;
		ooNode->right=rightchild;
		makeTree(node->right,rightchild);
	}
	return 0;
}

int
OptimalLeafOrder::reorder(OptimalLeafOrder::Node* n)
{
	//If this is left node then this is the only node possible here
	if(n->left==NULL && n->right==NULL)
	{
		map<int,double>* v=new map<int,double>;
		n->childOrder[n->id]=v;
		(*v)[n->id]=0;
		return 0;
	}
	reorder(n->left);
	reorder(n->right);
	map<int,map<int,double>*>& leftleaves=n->left->childOrder;
	map<int,map<int,double>*>& rightleaves=n->right->childOrder;
	for(map<int,map<int,double>*>::iterator vlIter=leftleaves.begin();vlIter!=leftleaves.end();vlIter++)
	{
		//Possible ms in vl
		map<int,double>* mset=vlIter->second; 
		for(map<int,map<int,double>*>::iterator ulIter=rightleaves.begin();ulIter!=rightleaves.end();ulIter++)
		{
			//Possible ks in vr
			map<int,double>* nset=ulIter->second;
			double mindist=10E7;
			string xleft;
			string xright;
			int xleftid=vlIter->first;
			int xrightid=ulIter->first;
			for(map<int,double>::iterator mIter=mset->begin();mIter!=mset->end();mIter++)
			{
				for(map<int,double>::iterator nIter=nset->begin();nIter!=nset->end();nIter++)
				{
					//double score=getSim((string&)mIter->first,(string&)nIter->first)+mIter->second+nIter->second;
					double score=getSim(mIter->first,nIter->first)+mIter->second+nIter->second;
					if(score<mindist)
					{
						mindist=score;
						if(mIter->first!=vlIter->first)
						{
							xleftid=mIter->first;
						}
						if(nIter->first!=ulIter->first)
						{
							xrightid=nIter->first;
						}
					}
				}
			}
			if(mindist==1000)
			{
				cout <<"odd.. very odd. mindist is still small" <<endl;
			}
			map<int,double>* vals=NULL;
			if(n->childOrder.find(vlIter->first)==n->childOrder.end())
			{
				vals=new map<int,double>;
				n->childOrder[vlIter->first]=vals;
			}
			else
			{
				vals=n->childOrder[vlIter->first];
			}
			(*vals)[ulIter->first]=mindist;
			vals=NULL;
			if(n->childOrder.find(ulIter->first)==n->childOrder.end())
			{
				vals=new map<int,double>;
				n->childOrder[ulIter->first]=vals;
			}
			else
			{
				vals=n->childOrder[ulIter->first];
			}
			(*vals)[vlIter->first]=mindist;
			string key;
			char keychar[256];
			if(ulIter->first<=vlIter->first)
			{
				sprintf(keychar,"%d-%d",ulIter->first,vlIter->first);
			}
			else
			{
				sprintf(keychar,"%d-%d",vlIter->first,ulIter->first);
			}
			key.append(keychar);
			if(xleftid<0 && xrightid<0)
			{
				continue;
			}
			OptimalLeafOrder::Pair* p=new OptimalLeafOrder::Pair;
			if(xleftid>=0)
			{
				p->leftextreme=xleftid;
			}
			if(xrightid>=0)
			{
				p->rightextreme=xrightid;
			}
			n->extremes[key]=p;
		}
	}
	return 0;
}

int
OptimalLeafOrder::getordering(OptimalLeafOrder::Node* n, int leftextreme, int rightextreme, vector<string>& leafordering)
{
	if(n->left==NULL && n->right==NULL)
	{
		leafordering.push_back(n->nodeName);
		return 0;
	}
	else if(n->extremes.size()==0)
	{
		if(nodeIDNameMap.find(leftextreme)==nodeIDNameMap.end() || nodeIDNameMap.find(rightextreme)==nodeIDNameMap.end())
		{
			cout << "No node id " << leftextreme << " or " << rightextreme << endl;
			exit(0);
		}
		string& lfname=nodeIDNameMap[leftextreme];
		leafordering.push_back(lfname);
		string& rfname=nodeIDNameMap[rightextreme];
		leafordering.push_back(rfname);
	}
	else
	{
		//leafordering.push_back(leftextreme);
		//leafordering.push_back(leftextreme);
		string key;
		char keychar[256];
		if(leftextreme<=rightextreme)
		{
			sprintf(keychar,"%d-%d",leftextreme,rightextreme);
		}
		else
		{
			sprintf(keychar,"%d-%d",rightextreme,leftextreme);
		}
		key.append(keychar);
		if(n->extremes.find(key)==n->extremes.end())
		{
			cout <<"No key " << key << " at " << n->nodeName <<  endl;
			exit(0);
		}

		Pair* ex=n->extremes[key];
		//if(ex->leftextreme.size()>0)
		//{
			//because the left and rights can switch, it is important to check that we are
			//setting the boundaries of the right trees
			if(n->left->childOrder.find(leftextreme)!=n->left->childOrder.end())	
			{
				getordering(n->left,leftextreme,ex->leftextreme,leafordering);
				getordering(n->right,ex->rightextreme,rightextreme,leafordering);
			}
			//switch case
			else if(n->right->childOrder.find(leftextreme)!=n->right->childOrder.end())	
			{
				getordering(n->right,leftextreme,ex->rightextreme,leafordering);
				getordering(n->left,ex->leftextreme,rightextreme,leafordering);
			}
		//}
		//if(ex->rightextreme.size()>0)
		//{
		//}
	}

	return 0;
}

double
OptimalLeafOrder::getSim(string& s1,string& s2)
{
	string key;
	if(strcmp(s1.c_str(),s2.c_str())<=0)
	{
		key.append(s1.c_str());
		key.append("-");
		key.append(s2.c_str());
	}
	else
	{
		key.append(s2.c_str());
		key.append("-");
		key.append(s1.c_str());
	}
	if(distance.find(key)==distance.end())
	{
		cout <<"No distance for " << key <<endl;
	}
	double d=distance[key];
	return d;
}

double
OptimalLeafOrder::getSim(int s1, int s2)
{
	double d=dist[s1][s2];
	return d;
}
