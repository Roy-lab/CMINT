#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "GeneExpManager.H"


GeneExpManager::GeneExpManager()
{

}

GeneExpManager::~GeneExpManager()
{

}

int 
GeneExpManager::readExpression(const char* aFName)
{
	ifstream inFile(aFName);
	char* buffer=NULL;
	int bufflen=0;
	string buffstr;
	int linecnt=0;
	while(inFile.good())
	{
		getline(inFile,buffstr);
		if(buffstr.length()<=0)
		{
			continue;
		}
		if(linecnt>=500)
		{
			linecnt++;
			//continue;
		}

		if(bufflen <=buffstr.length())
		{
			bufflen=buffstr.length()+1;
			if(buffer!=NULL)
			{
				delete[] buffer;
			}
			buffer=new char[bufflen];
		}
		strcpy(buffer,buffstr.c_str());
		int tokCnt=0;
		char* tok=strtok(buffer,"\t");
		string geneName;
		vector<double>* expr=NULL;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				geneName.append(tok);
				expr=new vector<double>;
				exprSet[geneName]=expr;
			}
			else
			{
				expr->push_back(atof(tok));
			}
			tokCnt++;
			tok=strtok(NULL,"\t");
		}
		linecnt++;
	}	
	inFile.close();
	return 0;
}


int 
GeneExpManager::readExpression_Withheader(const char* aFName)
{
	ifstream inFile(aFName);
	char* buffer=NULL;
	int bufflen=0;
	string buffstr;
	int linecnt=0;
	while(inFile.good())
	{
		getline(inFile,buffstr);
		if(buffstr.length()<=0)
		{
			continue;
		}

		if(bufflen <=buffstr.length())
		{
			bufflen=buffstr.length()+1;
			if(buffer!=NULL)
			{
				delete[] buffer;
			}
			buffer=new char[bufflen];
		}
		strcpy(buffer,buffstr.c_str());
		//First line is header
		if(linecnt==0)
		{
			readColNames(buffer);
		}
		else
		{
			int tokCnt=0;
			char* tok=strtok(buffer,"\t");
			string geneName;
			vector<double>* expr=NULL;
			while(tok!=NULL)
			{
				if(tokCnt==0)
				{
					geneName.append(tok);
					expr=new vector<double>;
					exprSet[geneName]=expr;
				}
				else
				{
					expr->push_back(atof(tok));
				}
				tokCnt++;
				tok=strtok(NULL,"\t");
			}
		}		
		linecnt++;
	}
	inFile.close();
	return 0;
}

int
GeneExpManager::readColNames(char* buffer)
{
	char* tok=strtok(buffer,"\t");
	int tokCnt=0;
	while(tok!=NULL)
	{
		if(tokCnt>0)
		{
			string cName(tok);
			colNames.push_back(cName);
		}
		tok=strtok(NULL,"\t");
		tokCnt++;
	}
	return 0;
}

vector<string>&
GeneExpManager::getColNames()
{
	return colNames;
}

vector<double>* 
GeneExpManager::getExp(const string& geneKey)
{
	if(exprSet.find(geneKey)==exprSet.end())
	{
		return NULL;
	}
	vector<double>* expr=exprSet[geneKey];
	return expr;
}

map<string,vector<double>*>&
GeneExpManager::getGeneSet()
{
	return exprSet;
}

