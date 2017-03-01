#include <iostream>
#include <fstream>
#include <vector>
#include <string.h>
#include <stdlib.h>
#include "HierarchicalClusterNode.H"
#include "HierarchicalCluster.H"
#include "OptimalLeafOrder.H"
#include "GeneExpManager.H"

#include "GeneMap.H"
#include "MappedOrthogroup.H"
#include "MappedOrthogroupReader.H"

#include "Framework.H"

Framework::Framework()
{
}

Framework::~Framework()
{
}

/*This function will read in the name of the clusterassignment file and also the mark values that we will need. We need to store the mark values column wise as per
 * the columns of the data matrix */
int 
Framework::readDataMatrix(const char* aDirName)
{
	char aFName[1024];
	sprintf(aFName,"%s/allcelltypes_clusterassign_brk.txt",aDirName);
	ifstream inFile(aFName);
	char* buffer=NULL;
	int bufflen=0;
	string strbuff;
	int gid=1;
	int totalLoci=0;
	while(inFile.good())
	{
		getline(inFile,strbuff);
		if(strbuff.length()<=0)
		{
			continue;
		}
		if(bufflen<=strbuff.length())
		{
			if(buffer!=NULL)
			{
				delete[] buffer;
			}
			bufflen=strbuff.length()+1;
			buffer=new char[bufflen];
		}
		strcpy(buffer,strbuff.c_str());
		if(strstr(buffer,"Loci")!=NULL)
		{
			readColumns(buffer);
			continue;
		}
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;	
		string geneName;
		vector<double> valueSet;
		map<int,int> frequency;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				geneName.append(tok);
			}
			else
			{
				double v=atof(tok);
				/*if(v>0)
				{
					node->attrib[tokCnt-1]=1;
				}
				else
				{
					node->attrib[tokCnt-1]=0;
				}
				node->attrib[tokCnt-1]=v+1;*/
				valueSet.push_back(v+1);
				int intVal=(int) (v+1);	
				if(frequency.find(intVal)==frequency.end())
				{
					frequency[intVal]=1;
				}
				else
				{
					frequency[intVal]=frequency[intVal]+1;
				}
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		totalLoci++;
		if(frequency.size()==1)
		{
			frequency.clear();
			valueSet.clear();
			continue;
		}	
		HierarchicalClusterNode* node=new HierarchicalClusterNode;
		for(int i=0;i<valueSet.size();i++)
		{
			node->expr.push_back(valueSet[i]);
			//node->attrib[i]=valueSet[i];
		}
		node->nodeName.append(geneName.c_str());
		nodeSet[geneName]=node;
		backup[geneName]=node;
		nameIDMap[geneName]=gid;
		node->size=1;
		gid++;
		frequency.clear();
		valueSet.clear();
	}
	cout <<"Found" << nodeSet.size() << " loci that change cluster assignment of total " << totalLoci << endl;
	if(attribNameIDMap.size()==0)
	{
		cout <<"Did not find any cell types! "<<endl;
		exit(0);
	}
	inFile.close();
	//Now read the mark profiles. This should be read in the same order as the order of attribIDNameMap
	for(map<int,string>::iterator cIter=attribIDNameMap.begin();cIter!=attribIDNameMap.end();cIter++)
	{
		char geneFName[1024];
		sprintf(geneFName,"%s/%s_exprtab.txt",aDirName,cIter->second.c_str());
		GeneExpManager* geExp=new GeneExpManager;
		geExp->readExpression_Withheader(geneFName);
		markProfileSet[cIter->second]=geExp;
	}
	cout <<"Read " << markProfileSet.size() << " different datasets " << endl;
	createMarkHeaders();
	return 0;
}

int
Framework::readOGIDs(const char* celltypeOrder,const char* aFName)
{
	//mor.setSpeciesMapping(attribIDNameMap);
	mor.readSpeciesMapping(celltypeOrder);
	mor.readFile(aFName);
	return 0;
}

int
Framework::setSrcCellType(const char* src)
{
	srcCelltype.append(src);
	return 0;
}

//This function will generate one file which is the assignment of genes in all cell types and this will include everything. 
//The second set of files will be the clustersets that will include all the individual transitions. I think we should also generate the cluster mark profiles for the clusterset.
//It's better to output in heatmap.awk compatible format right here.
int 
Framework::generateTransitioningGeneSets(double threshold,const char* outdir,int mingenesetSize)
{
	map<int,map<string,int>*> modules;
	map<int,HierarchicalClusterNode*> attribs;
	cluster.setDistanceType(HierarchicalCluster::CITYBLOCK);
	cluster.cluster(modules,nodeSet,threshold,attribs);
	double** dist=cluster.getDist();
	char aFName[1024];
	sprintf(aFName,"%s/all_assign.txt",outdir);
	ofstream oFile(aFName);
	int c=0;
	int atcnt=0;
	olo.setDist(dist);
	for(map<int,HierarchicalClusterNode*>::iterator aIter=attribs.begin();aIter!=attribs.end();aIter++)
	{
		olo.setHierarchicalClusterNode(aIter->second);
		if(aIter->second->left!=NULL || aIter->second->right!=NULL)
		{
		//	cout <<"Stop here " << endl;
		}
		vector<string>* ordering=new vector<string>;
		olo.reorder(*ordering);
		/*if(ordering.size()<3) //|| ordering.size()>200)
		{
			continue;
		}*/
		for(int i=0;i<ordering->size();i++)
		{
			HierarchicalClusterNode* hc=backup[(*ordering)[i]];
			
			//for(map<int,double>::iterator aIter=hc->attrib.begin();aIter!=hc->attrib.end();aIter++)
			for(int j=0;j<hc->expr.size();j++)
			{
				string& colName=attribIDNameMap[j];
				oFile<<(*ordering)[i]<< "||6\t" <<colName <<"||8\t" << hc->expr[j]<<"|1|"<< hc->expr[j]<< endl;
			}
			//Now put a vertical spacer
			//we need this only once
			if(c==0)
			{
				oFile <<"|- MyVspacer1|Spacer"<< endl;
			}
			//Not print out out the mark profiles too
			//for(map<int,double>::iterator aIter=hc->attrib.begin();aIter!=hc->attrib.end();aIter++)
			for(int j=0;j<hc->expr.size();j++)
			{
				//string& colName=attribIDNameMap[aIter->first];
				string& colName=attribIDNameMap[j];
				GeneExpManager* geMgr=markProfileSet[colName];
				vector<string>& markNames=geMgr->getColNames();
				//here we have to be careful because we don't know what the locus name in the other cell type is
				const char* lociCelltype=getNameinCelltype((*ordering)[i].c_str(),colName.c_str());
				//vector<double>* expvals=geMgr->getExp((*ordering)[i]);
				vector<double>* expvals=geMgr->getExp(lociCelltype);
				for(int m=0;m<markNames.size();m++)
				{
					oFile<< (*ordering)[i]<<"||6\t"<<colName << "_" << markNames[m] << "||6\t" << (*expvals)[m]<<"|2" << endl;
				}
				if(c==0 && i==0 && (j<hc->expr.size()-1))
				{
					oFile <<"|- cellmarkerVspacer"<<(j+1)<<"|Spacer"<< endl;
				}
			}
			atcnt=hc->expr.size();	
		}
		//Now put a horizontal spacer		
		oFile <<"|Spacer||"<<c << " |-"<< endl;
		if((*ordering).size()>=mingenesetSize)
		{
			clusterset[c]=ordering;
			char cFName[1024];
			sprintf(cFName,"%s/clusterset%d.txt",outdir,c);
			ofstream coFile(cFName);
			//Show header
			for(int i=0;i<ordering->size();i++)
			{
				HierarchicalClusterNode* hc=backup[(*ordering)[i]];
			
				//for(map<int,double>::iterator aIter=hc->attrib.begin();aIter!=hc->attrib.end();aIter++)
				for(int j=0;j<hc->expr.size();j++)
				{
					coFile<<(*ordering)[i] <<"||6\t"<< attribIDNameMap[j] << "||8\t" << hc->expr[j]<<"|1|"<< hc->expr[j]<<endl;
				}
				coFile <<"|- MyVspacer1|Spacer"<< endl;
				//Not print out out the mark profiles too
				for(int j=0;j<hc->expr.size();j++)
				{
					string& colName=attribIDNameMap[j];
					GeneExpManager* geMgr=markProfileSet[colName];
					vector<string>& markNames=geMgr->getColNames();
					//vector<double>* expvals=geMgr->getExp((*ordering)[i]);
					const char* lociCelltype=getNameinCelltype((*ordering)[i].c_str(),colName.c_str());
					vector<double>* expvals=geMgr->getExp(lociCelltype);
					for(int m=0;m<markNames.size();m++)
					{
						coFile<< (*ordering)[i]<<"||6\t"<<colName << "_" << markNames[m] << "||6\t" << (*expvals)[m]<<"|2" << endl;
					}
					if(i==0 && j <hc->expr.size())
					{
						coFile <<"|- cellmarkerVspacer"<<(j+1)<<"|Spacer"<< endl;
					}
				}
			}
			coFile.close();
		}
		
		c++;
	}
	oFile.close();
	return 0;
}

int
Framework::generateOrderedClusterMeans(const char* outDirName)
{
	map<string,vector<double>*> meanC_Set;
	map<string,vector<int>*> meanD_Set;
	map<string,HierarchicalClusterNode*> meanNodeSet;
	int gid=0;
	//We need to compute the means, and reorder them for the clusters that are of reasonable size
	for(map<int,vector<string>*>::iterator aIter=clusterset.begin();aIter!=clusterset.end();aIter++)
	{
		vector<string>* clusterMembers=aIter->second;
		vector<double>* meanC=getMeanContinuous(clusterMembers);
		vector<int>*  meanD=getMeanDiscrete(clusterMembers);
		HierarchicalClusterNode* node=new HierarchicalClusterNode;
		for(int i=0;i<meanC->size();i++)
		{
			node->expr.push_back((*meanC)[i]);
		}
		char cname[1024];
		sprintf(cname,"clusterset%d",aIter->first);
		node->nodeName.append(cname);
		meanNodeSet[cname]=node;
		meanC_Set[cname]=meanC;
		meanD_Set[cname]=meanD;
		nameIDMap[cname]=gid;
		node->size=1;
		gid++;
	}
	map<int,map<string,int>*> modules;
	HierarchicalCluster clusterMeans;
	clusterMeans.setDistanceType(HierarchicalCluster::PEARSON);
	clusterMeans.cluster(modules,meanNodeSet,1);
	olo.setHierarchicalClusterNode(clusterMeans.getRoot());
	vector<string> ordering;
	olo.reorder(ordering);
	char aFName[1024];
	sprintf(aFName,"%s/ordered_clusterset_means.txt",outDirName);
	ofstream oFile(aFName);
	for(int i=0;i<ordering.size();i++)
	{
		HierarchicalClusterNode* hc=meanNodeSet[ordering[i]];
		//for(map<int,double>::iterator aIter=hc->attrib.begin();aIter!=hc->attrib.end();aIter++)
		/*for(map<int,string>::iterator aIter=attribIDNameMap_Profile.begin();aIter!=attribIDNameMap_Profile.end();aIter++)
		{
			//oFile<<ordering[i]<<"\t" << aIter->second<<"\t" << aIter->second << endl;
			oFile<<ordering[i]<<"||8\t" << aIter->second <<"||8\t" << hc->expr[aIter->first] <<"|2" << endl;
		}*/
		
		for(map<int,string>::iterator aIter=attribIDNameMap.begin();aIter!=attribIDNameMap.end();aIter++)
		{
			//Get the cellName 
			string& cellName=(string&)aIter->second; 
			//Get the marks
			GeneExpManager* geMgr=markProfileSet[cellName];
			vector<string>& colNames=geMgr->getColNames();
			char cellMark[1024];
			for(int j=0;j<colNames.size();j++)
			{
				sprintf(cellMark,"%s_%s",cellName.c_str(),colNames[j].c_str());
				string attribName(cellMark);
				int markID=attribNameIDMap_Profile[attribName];
				oFile<<ordering[i]<<"||8\t" << attribName <<"||8\t" << hc->expr[markID] <<"|2" << endl;
			}
			if(i==0 && aIter->first<attribIDNameMap.size())
			{
				oFile <<"|- cellmarkerVspacer"<<(aIter->first+1)<<"|Spacer"<< endl;
			}
		}
	//	oFile <<"|- MyVspacer1|Spacer"<< endl;
		vector<int>* meanD=meanD_Set[ordering[i]];
		for(int j=0;j<meanD->size();j++)
		{
			oFile << ordering[i] << "||8\t" << attribIDNameMap[j] <<"||8\t" << (*meanD)[j]<<"|1|"<<(*meanD)[j]<< endl;
		}
	}
	oFile.close();
	return 0;
}

//This will create a pseudo header using the ordering of the cell types in the cmint assignment set
int
Framework::createMarkHeaders()
{
	int markID=0;
	for(map<int,string>::iterator aIter=attribIDNameMap.begin();aIter!=attribIDNameMap.end();aIter++)
	{
		//Get the cellName 
		string& cellName=(string&)aIter->second; 
		//Get the marks
		GeneExpManager* geMgr=markProfileSet[cellName];
		vector<string>& colNames=geMgr->getColNames();
		char cellMark[1024];
		for(int i=0;i<colNames.size();i++)
		{
			sprintf(cellMark,"%s_%s",cellName.c_str(),colNames[i].c_str());
			string attribName(cellMark);
			attribIDNameMap_Profile[markID]=attribName;
			attribNameIDMap_Profile[attribName]=markID;
			markID++;	
		}
	}
	return 0;
}

vector<double>*
Framework::getMeanContinuous(vector<string>* geneSet)
{
	vector<double>* theMean=new vector<double>;	
	//get the cell types in the order of the cellNameIDMap
	for(map<int,string>::iterator aIter=attribIDNameMap.begin();aIter!=attribIDNameMap.end();aIter++)
	{
		string& colName=aIter->second;
		GeneExpManager* geMgr=markProfileSet[colName];
		vector<string>& colNames=geMgr->getColNames();
		for(int j=0;j<colNames.size();j++)
		{
			double s=0;
			for(int i=0;i<geneSet->size();i++)
			{
				const char* lociCelltype=getNameinCelltype((*geneSet)[i].c_str(),colName.c_str());
				vector<double>* vals=geMgr->getExp(lociCelltype);
				s=s+(*vals)[j];		
			}
			s=s/geneSet->size();
			theMean->push_back(s);
		}
	}
	return theMean;
}

vector<int>*
Framework::getMeanDiscrete(vector<string>* geneSet)
{
	vector<int>* theMean=new vector<int>;
	for(map<int,string>::iterator aIter=attribIDNameMap.begin();aIter!=attribIDNameMap.end();aIter++)
	{
		map<int,int> frequency;
		for(int i=0;i<geneSet->size();i++)
		{
			string& gname=(*geneSet)[i];
			HierarchicalClusterNode* hc=backup[gname];
			int val=(int)hc->expr[aIter->first];
			if(frequency.find(val)==frequency.end())
			{
				frequency[val]=1;
			}
			else
			{
				frequency[val]=frequency[val]+1;
			}
		}
		//Now get the majority. In case of a tie we will use the max of the assignments (hopefully there won't be ties!)
		int maxCnt=-1;
		int maxVal=-1;
		for(map<int,int>::iterator cIter=frequency.begin();cIter!=frequency.end();cIter++)
		{
			if(cIter->second>maxCnt)
			{
				maxCnt=cIter->second;
				maxVal=cIter->first;
			}
		}
		//Now check for ties
		for(map<int,int>::iterator cIter=frequency.begin();cIter!=frequency.end();cIter++)
		{
			if(cIter->second==maxCnt && maxVal!=cIter->first)
			{
				cout <<"Oops found a tie of " << maxVal<< " with " << cIter->first << endl;
				if(cIter->first>maxVal)
				{
					maxVal=cIter->first;
				}
			}
		}
		frequency.clear();
		theMean->push_back(maxVal);
	}
	return theMean;
}



int
Framework::readColumns(char* buffer)
{
	int tokCnt=0;
	char* tok=strtok(buffer,"\t");
	while(tok!=NULL)
	{
		if(tokCnt>0)
		{
			string colName(tok);
			attribNameIDMap[colName]=tokCnt-1;
			attribIDNameMap[tokCnt-1]=colName;
		}
		tok=strtok(NULL,"\t");
		tokCnt++;
	}
	return 0;
}

const char*
Framework::getNameinCelltype(const char* locus, const char* celltype)
{
	MappedOrthogroup* mgrp=mor.getMappedOrthogroup(locus,srcCelltype.c_str());
	if(mgrp==NULL)
	{
		cout <<"No region with name " << locus << " in " << srcCelltype << "!"<<endl;
		exit(0);
	}
	const char* name=mgrp->getSpeciesHit(celltype);
	if(name==NULL)
	{
		cout <<"No region representation of " << locus << "in " << celltype << "!"<< endl;
		exit(0);
	}
	return name;
}

int
main(int argc, const char** argv)
{
	if(argc!=8)
	{
		//This is configured for cmint and expects specific files. So make sure to have those files in the datadir
		cout <<"Usage: ./findTransitionGenesets cmint_result_dir celltypeorder OGID_file srccelltype threshold outputdir maxgenesetsize" << endl;
		return 0;
	}
	Framework fw;
	fw.readDataMatrix(argv[1]);
	fw.readOGIDs(argv[2],argv[3]);
	fw.setSrcCellType(argv[4]);
	fw.generateTransitioningGeneSets(atof(argv[5]),argv[6],atoi(argv[7]));
	fw.generateOrderedClusterMeans(argv[6]);
	return 0;
}

