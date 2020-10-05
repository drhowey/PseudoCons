/************************************************************************
 * PseudoCons, version 1.12
 * Copyright 2013-2014,
 * Richard Howey
 * Institute of Genetic Medicine, Newcastle University
 *
 * richard.howey@ncl.ac.uk
 * http://www.staff.ncl.ac.uk/richard.howey/
 *
 * This file is part of PseudoCons.
 *
 * PseudoCons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PseudoCons is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PseudoCons.  If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/


/*! \file Data.h
    \brief This file contains basic classes for storing genotype infomation.
    
*/

#ifndef __DATA
#define __DATA

#include <string>
#include <map>
#include <set>
#include <list>
#include <iostream>
#include <fstream>
#include <ostream>


using namespace std; // initiates the "std" or "standard" namespace


//! Maps the names of subjects and pedigrees given in file to an ordinal number.
class MapIds
{
private:
	map<string, unsigned int> idMap; //!< original name in file, counted id

public:

	MapIds() : idMap() {};

	//! Create a new ID for an item.
	unsigned int addItem(const string & name)
	{
		unsigned int newId = (unsigned)idMap.size() + 1;
		idMap[name] = newId;
		return newId; 
	};

	//! Look up ID from a name given from file.
	unsigned int getId(const string & name);
	//! Look up name from a ID
	string getName(const unsigned int & id);
};

//! Maps the allele names given in file to either a 0 or 1.

//! In the .ped pedigree file the SNP names may be given by a letter, e.g. A and G,
//! these letters are mapped to values 0 and 1 in the code which in turn map to
//! allele "1" and allele "2" respectively.
class MapSnpAlleleIds
{
private:
	map<unsigned int, map<bool, string> > mapSnpAlleles; //SNP ID no., allele id - 0 or 1, allele name

public:

	MapSnpAlleleIds() {};

	~MapSnpAlleleIds() {};

	void addSNPAlleleNames(unsigned int & snpID, map<bool, string> & names) {mapSnpAlleles[snpID] = names;};
	map<bool, string> getSnpAlleleNames(const unsigned int & snpId) const; 
};

void writeLastByteBeforeNextSNP(ofstream & psBedFile, unsigned int & psBitCount, int & aBit);

//! A subject with its attributes including parents.
class Subject
{
private:
	unsigned int id;
	unsigned int fatherId;
	unsigned int motherId;
	unsigned int sex; //1 = male, 2 = female, 0 = unknown
	bool affected;
	bool removed;
	
	string userIDName;

public:

	Subject(unsigned int & i, unsigned int & fi, unsigned int & mi, unsigned int & s, bool & a, string & uin): id(i), fatherId(fi), motherId(mi), sex(s), affected(a), removed(false), userIDName(uin) {};

	virtual ~Subject() {};
	
	virtual void addGenotype(const bool & all1, const bool & all2, const unsigned int & geno = 1) {};
	virtual pair<bool, bool> getGenotype(const unsigned int & geno = 1) const {return make_pair(true, false);};
	virtual bool genotypeDataPresent() const {return false;};
	
	bool getAffected() const {return affected;};
	unsigned int getFatherId() const {return fatherId;};
	unsigned int getMotherId() const {return motherId;};
	unsigned int getId() const {return id;};
	unsigned int getSexId() const {return sex;};
	string getIDName() {return userIDName;};
	void removeFromPedigree() {removed = true;};
	void restoreToPedigree() {removed = false;};
	bool isNotRemoved() {return !removed;};
	void outputDetails();
};

//! A subject with genotype data for only one SNP.
class SubjectOneSnp : public Subject
{
private:
	bool allele1, allele2; //let 0 = allele 1 in BIM file, 1 = allele 2 in BIM file

public:

	SubjectOneSnp(unsigned int & i, unsigned int & fi, unsigned int & mi, unsigned int & s, bool & a, string & n): Subject(i, fi, mi, s, a, n) {};
	
	~SubjectOneSnp() {};

	void addGenotype(const bool & all1, const bool & all2, const unsigned int & geno = 1) {allele1 = all1; allele2 = all2;};

	pair<bool, bool> getGenotype(const unsigned int & geno = 1) const
	{		
		return make_pair(allele1, allele2);
	};

	bool genotypeDataPresent() const {return (!allele1 || allele2);}; //genotype 1 / 0 denotes missing genotype	
};

//! A subject with genotype data for two SNPs.
class SubjectTwoSnps : public Subject
{
private:
	bool geno1Allele1, geno1Allele2; //let 0 = allele 1 in BIM file, 1 = allele 2 in BIM file
	bool geno2Allele1, geno2Allele2;

public:

	SubjectTwoSnps(unsigned int & i, unsigned int & fi, unsigned int & mi, unsigned int & s, bool & a, string & n): Subject(i, fi, mi, s, a, n) {};
	
	~SubjectTwoSnps() {};

	void addGenotype(const bool & all1, const bool & all2, const unsigned int & geno = 1)
	{
		if(geno == 1) {geno1Allele1 = all1; geno1Allele2 = all2;}
		else {geno2Allele1 = all1; geno2Allele2 = all2;};
	};

	pair<bool, bool> getGenotype(const unsigned int & geno = 1) const
	{			
		if(geno == 1) {return make_pair(geno1Allele1, geno1Allele2);}
		else {return make_pair(geno2Allele1, geno2Allele2);};
	};

	bool genotypeDataPresent() const {return ((!geno1Allele1 || geno1Allele2) && (!geno2Allele1 || geno2Allele2));}; //genotype 1 / 0 denotes missing genotype	
};

//! A simple grouping of a father, mother and child and pseudocontrols.
class Trio
{
private:
	string pedIdName;
	Subject * father;
	Subject * mother;
	Subject * child;
	list<Subject * > pseudoControls;

	//info for info file
	unsigned int pseudoConSetNo;
	list<bool> fatherGenoTypes; //0 = allele 1 in BIM file, 1 = allele 2 in BIM file, added in pairs
	list<bool> motherGenoTypes;
	list<bool> childGenoTypes;
	list<list<bool> > pseudoControlGenotypes;

public:
	Trio() : father(0), mother(0), child(0), pseudoControls() {};
	Trio(string & pId, Subject * f, Subject * m, Subject * c) : pedIdName(pId), father(f), mother(m), child(c), pseudoControls() {};

	~Trio()
	{
		for(list<Subject * >::iterator p = pseudoControls.begin(); p != pseudoControls.end(); ++p)
		{
			delete *p;
		};
	};

	void addFather(Subject * f) {father = f;};
	void addMother(Subject * m) {mother = m;};
	void addChild(Subject * c) {child = c;};
	Subject * getFather() {return father;};
	Subject * getMother() {return mother;};
	Subject * getChild() {return child;};
	void createPseudoControls(unsigned int & pseudoControlSetNo, unsigned int & noPseudoCons, bool & cepg, MapIds & subjectIds, ofstream & psFamFile, bool & interaction);
	void updatePsuedoControlSNPData(unsigned int & noPseudoCons, bool & cepg, ofstream & psBedFile, unsigned int & psBitCount, int & aBit);
	void updatePsuedoControlJointSNPData(unsigned int & noPseudoCons, bool & cepg, ofstream & psBedFile, unsigned int & psBitCount, int & aBit);
	void writePseudoControlSecondSNPData(ofstream & psBedFile, unsigned int & psBitCount, int & aBit);
	list<Subject * > getPseudoControls() const {return pseudoControls;};
	void storeSNPDataForInfoFile(const unsigned int & geno = 1);
	void outputInfoFile(ofstream & infoFile, MapSnpAlleleIds & mapSnpAlleleIds, bool & infoFileFatherGeno, bool & infoFileMotherGeno, bool & interaction, bool & cepg, unsigned int & snpNo1, unsigned int & snpNo2);
};

//! Contains a list of all subjects belonging to a pedigree.
class Pedigree
{
private:
	unsigned int pedId;
	string pedIdName;
	map<unsigned int, Subject *> subjects;//subject Id, pointer to subject object

public:

	Pedigree(unsigned int & p, string & n) : pedId(p), pedIdName(n), subjects() {};

	~Pedigree()
	{
		for(map<unsigned int, Subject *>::iterator s = subjects.begin(); s != subjects.end(); ++s)
		{
			delete s->second;
		};
	};

	void addSubject(unsigned int & subjectId, Subject * subject)
	{
		subjects[subjectId] = subject;
	};

	unsigned int getPedId() {return pedId;};
	Subject * getSubject(const unsigned int & id) const;
	Subject * getSubject(const unsigned int & id, bool & exists) const;
	Trio * findCaseParentTrio(bool & found);	
	Trio * findCaseParentTrio(bool & found, Subject * subject);
	
	unsigned int getNumberOfSubjects() const {return (unsigned)subjects.size();};
	void restoreAllSubjectsToPedigree();
}; 


//! Contains a list of all pedigrees.
class AllPedigreeData
{
private:
	

public:

	map<unsigned int, Pedigree *> pedigrees; // pedigree ID, pedigree

	AllPedigreeData() : pedigrees() {};

	//! Delete all the pedigrees, these belong to this class
	~AllPedigreeData()
	{
		for(map<unsigned int, Pedigree *>::iterator p = pedigrees.begin(); p != pedigrees.end(); ++p)
		{
			delete p->second;
		};
	};

	void addSubjectToPedigree(unsigned int & subjectId, Subject * subject, unsigned int & pedigreeId, string & pedIdName);
	void outputSummary();
	void restoreAllSubjectsToAllPedigrees();
	void subjectExistsAndAffected(const unsigned int & pedId, const unsigned int & subjectId, MapIds & subjectIds);  
};

//! Contains a list of trios
class AllTrioData
{
private:
	unsigned int noPseudoCons;	
	bool cepg;

public:

	list<Trio *> trios; //trios

	AllTrioData() : trios() {};

	//! Delete all the pedigrees, these belong to this class
	~AllTrioData()
	{
		for(list<Trio *>::iterator t = trios.begin(); t != trios.end(); ++t)
		{
			delete *t;
		};
	};

	void setNoPseudoConsAndCEPG(unsigned int & noPsCons, bool & c) {noPseudoCons = noPsCons; cepg = c;};
	void addTrio(Trio * aTrio);
	void createPseudoControls(MapIds & subjectIds, ofstream & psFamFile, bool & interaction);
	void updatePsuedoControlSNPData(ofstream & psBedFile, unsigned int & psBitCount, int & aBit);
	void updatePsuedoControlJointSNPData(ofstream & psBedFile, unsigned int & psBitCount, int & aBit);
	void outputPseudoControlInfo();
	void storeSNPDataForInfoFile(const unsigned int & geno = 1);
	void outputInfoFile(string & infoFilename, MapSnpAlleleIds & mapSnpAlleleIds, bool & infoFileFatherGeno, bool & infoFileMotherGeno, bool & interaction, bool & cepg, unsigned int & snpNo1, unsigned int & snpNo2);
};

#endif
