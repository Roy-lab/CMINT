#ifndef _ORTHOGROUP_
#define _ORTHOGROUP_
#include <string>
#include <map>
using namespace std;

class GeneMap;
class MappedOrthogroup
{
	public:
	MappedOrthogroup();
	~MappedOrthogroup();
	int setMembers(map<string,string>&);
	int setID(int);
	int getID();
	int incrCnt();
	int getCnt();

	map<string,GeneMap*>& getOrthoMembers();
	map<string,string>& getOrthoMembers_names();
	GeneMap* getSpeciesHits(const char* specName);
	const char* getSpeciesHit(const char* specName);
	STRINTMAP* getSpeciesHitsForGene(const char* srcSpecName, const char* targetSpecName, const char* geneName);
	map<int,map<string,string>*>& getGeneSets();
	private:
	map<string,GeneMap*> orthoMembers;
	map<string,string> orthoMembers_names;
	map<string,int> scerHits;
	map<int,map<string,string>*> geneSets;
	map<string,int> specList;
	int oid;
	int cnts;
};
#endif
