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


/*! \file Data.cpp
    \brief This file contains functions for processing basic genotype infomation.
       
*/

#include "Data.h"
#include "main.h"
#include "Analysis.h"

#include <string>
#include <map>
#include <set>
#include <math.h>
#include <iostream>
#include <fstream>
#include <ostream>
#include <sstream>
#include <cstdlib>

using namespace std; // initiates the "std" or "standard" namespace


//! Get the ID of a subject, family ID + subject ID maps to an ordinal number over all subjects. 
unsigned int MapIds::getId(const string & name)
{
  map<string, unsigned int>::const_iterator i = idMap.find(name);
  if(i != idMap.end()) return i->second;

  //if id does not exist create an id
  unsigned int id = (unsigned)idMap.size() + 1;
  idMap[name] = id;
  
  return id;
};

//! Get the name of a subject, given the id.
string MapIds::getName(const unsigned int & id)
{
	//loop thro' names to find the id, this is not too quick but is only used for reporting errors/warnings
	for(map<string, unsigned int>::const_iterator i = idMap.begin(); i != idMap.end(); ++i)
	{
		if(id == i->second) return i->first;
	};

	return "ID not found";
};


//! Get the allele names for a SNP
map<bool, string> MapSnpAlleleIds::getSnpAlleleNames(const unsigned int & snpId) const
{
	map<unsigned int, map<bool, string> >::const_iterator msa = mapSnpAlleles.find(snpId);

	if(msa != mapSnpAlleles.end()) return msa->second;
	else
	{
		map<bool, string> dummy;
		dummy[false] = "NoData ";
		dummy[true] = "NoData";
		return dummy;
	};
};

//! Adds a subject to the pedigree.
void AllPedigreeData::addSubjectToPedigree(unsigned int & subjectId, Subject * subject, unsigned int & pedigreeId, string & pedIdName)
{
	map<unsigned int, Pedigree *>::iterator p = pedigrees.find(pedigreeId);

	//if pedigree exists add the subject to it otherwise add subject to a new pedigree
	if(p != pedigrees.end())
	{
		p->second->addSubject(subjectId, subject);
	}
	else
	{
		Pedigree * pedigree = new Pedigree(pedigreeId, pedIdName);
		pedigree->addSubject(subjectId, subject);
		pedigrees[pedigreeId] = pedigree;
	};

};

//! Outputs basic stats of the pedigree to screen.
void AllPedigreeData::outputSummary()
{
	double noOfPedigrees = (double)pedigrees.size();
	double meanPedigreeSize = 0;
	double standardDeviationPedigreeSize = 0;

	for(map<unsigned int, Pedigree *>::const_iterator ped = pedigrees.begin(); ped != pedigrees.end(); ++ped)
	{
		meanPedigreeSize += ped->second->getNumberOfSubjects();
	};
	meanPedigreeSize /= noOfPedigrees;

	for(map<unsigned int, Pedigree *>::const_iterator ped = pedigrees.begin(); ped != pedigrees.end(); ++ped)
	{
		standardDeviationPedigreeSize += (ped->second->getNumberOfSubjects()- meanPedigreeSize)*(ped->second->getNumberOfSubjects()- meanPedigreeSize);
	};
	standardDeviationPedigreeSize /= (noOfPedigrees - 1);
	standardDeviationPedigreeSize = sqrt(standardDeviationPedigreeSize);

	out("Number of pedigrees: "); out(noOfPedigrees); out("\n");
	out("Mean pedigree size: "); out(meanPedigreeSize); out("\n");
	out("Standard deviation of pedigree size: "); out(standardDeviationPedigreeSize); out("\n\n");
};

//! Unmarks subjects in pedigree, used when allowing extras trios.
void AllPedigreeData::restoreAllSubjectsToAllPedigrees()
{
	for(map<unsigned int, Pedigree *>::iterator p = pedigrees.begin(); p != pedigrees.end(); ++p)
	{
		p->second->restoreAllSubjectsToPedigree();
	};

};

//! checks if a subject exists and if it does whether it is affected or not
void AllPedigreeData::subjectExistsAndAffected(const unsigned int & pedId, const unsigned int & subjectId, MapIds & subjectIds) 
{
	bool exists = false;
	bool affected = false;
	bool fatherExists = false;
	bool motherExists = false;

	Subject * subject;
	
	for(map<unsigned int, Pedigree *>::iterator p = pedigrees.begin(); p != pedigrees.end(); ++p)
	{
		subject = p->second->getSubject(subjectId, exists);
		if(exists)
		{
			affected = subject->getAffected();
			fatherExists = subject->getFatherId() != 0;
			motherExists = subject->getMotherId() != 0;
			break;
		};
	};

	if(!exists || !affected || !fatherExists || !motherExists)
	{
		string subjectName = subjectIds.getName(subjectId);
		
		if(!exists) {outErr("Warning: proband subject "); outErr(subjectName); outErr(" does not exist!\n");}
		else if(!affected) {outErr("Warning: proband subject "); outErr(subjectName); outErr(" is not affected!\n");}
		else if(!motherExists) {outErr("Warning: proband subject "); outErr(subjectName); outErr(" has no mother!\n");}
		else if(!fatherExists) {outErr("Warning: proband subject "); outErr(subjectName); outErr(" has no father!\n");};
	};
	
};

//! Adds trio to list and updates pseudocontrol basic info.
void AllTrioData::addTrio(Trio * aTrio)
{
	trios.push_back(aTrio);	
};	

//! Creates pseudocontrols for each trio.
void AllTrioData::createPseudoControls(MapIds & subjectIds, ofstream & psFamFile, bool & interaction)
{
	unsigned int pseudoControlSetNo = 1;
	for(list<Trio *>::iterator t = trios.begin(); t != trios.end(); ++t)
	{
		(*t)->createPseudoControls(pseudoControlSetNo, noPseudoCons, cepg, subjectIds, psFamFile, interaction);
		pseudoControlSetNo++;
	};
};

//! Updates SNP info for each pseudocontrol.
void AllTrioData::updatePsuedoControlSNPData(ofstream & psBedFile, unsigned int & psBitCount, int & aBit)
{
	for(list<Trio *>::iterator t = trios.begin(); t != trios.end(); ++t)
	{
		(*t)->updatePsuedoControlSNPData(noPseudoCons, cepg, psBedFile, psBitCount, aBit);
	};
};

//! Adds SNP data to trios to output for info file.
void AllTrioData::storeSNPDataForInfoFile(const unsigned int & geno)
{
	for(list<Trio *>::iterator t = trios.begin(); t != trios.end(); ++t)
	{
		(*t)->storeSNPDataForInfoFile(geno);
	};
};

//!Outputs info file.
void AllTrioData::outputInfoFile(string & infoFilename, MapSnpAlleleIds & mapSnpAlleleIds, bool & infoFileFatherGeno, bool & infoFileMotherGeno, bool & interaction, bool & cepg, unsigned int & snpNo1, unsigned int & snpNo2)
{
	ofstream infoFile(infoFilename.c_str());

	for(list<Trio *>::iterator t = trios.begin(); t != trios.end(); ++t)
	{
		(*t)->outputInfoFile(infoFile, mapSnpAlleleIds, infoFileFatherGeno, infoFileMotherGeno, interaction, cepg, snpNo1, snpNo2);
	};

	infoFile.close();
};

//! Write last byte and resets the Bit and bit counter
void writeLastByteBeforeNextSNP(ofstream & psBedFile, unsigned int & psBitCount, int & aBit)
{
	if(psBitCount == 0) return; //last byte may already be written if the number of pseudocontrols is a multiple of 4

	//shift right bits over
	while(psBitCount < 8)
	{
		aBit = aBit >> 1; //shift bits to the right
		psBitCount++;
	};
	
	//write to file
	char buffer[1];
	buffer[0] = aBit;
	psBedFile.write(buffer, 1);
	psBitCount = 0;
	aBit = 0;
};

//! Updates joint SNP info for each pseudocontrol.
void AllTrioData::updatePsuedoControlJointSNPData(ofstream & psBedFile, unsigned int & psBitCount, int & aBit)
{

	//update joint pseudocontrol data including writing data for the first SNP
	for(list<Trio *>::iterator t = trios.begin(); t != trios.end(); ++t)
	{
		(*t)->updatePsuedoControlJointSNPData(noPseudoCons, cepg, psBedFile, psBitCount, aBit);
	};

	//start new byte as starting a new SNP
	writeLastByteBeforeNextSNP(psBedFile, psBitCount, aBit);

	//now write the data for the second SNP
	for(list<Trio *>::iterator f = trios.begin(); f != trios.end(); ++f)
	{
		(*f)->writePseudoControlSecondSNPData(psBedFile, psBitCount, aBit);
	};
};

//! Outputs info about pseudocontrols to screen and log file.
void AllTrioData::outputPseudoControlInfo()
{
	list<Subject * > pseudoControls;
	unsigned int maleCountPseudo = 0;
	unsigned int femaleCountPseudo = 0;
	unsigned int unknownCountPseudo = 0;
	unsigned int maleCountCase = 0;
	unsigned int femaleCountCase = 0;
	unsigned int unknownCountCase = 0;
	unsigned int sex;

	for(list<Trio *>::iterator t = trios.begin(); t != trios.end(); ++t)
	{
		sex = (*t)->getChild()->getSexId();
		if(sex == 1) maleCountCase++;
		else if(sex == 2) femaleCountCase++;
		else unknownCountCase++;

		pseudoControls = (*t)->getPseudoControls();
		for(list<Subject * >::iterator p = pseudoControls.begin(); p != pseudoControls.end(); ++p)
		{
			sex = (*p)->getSexId();
			if(sex == 1) maleCountPseudo++;
			else if(sex == 2) femaleCountPseudo++;
			else unknownCountPseudo++;
		};
	};

	unsigned int totalNoOfCases = maleCountCase+femaleCountCase+unknownCountCase;

	double femalePercent, malePercent, unknownSexPercent;
	malePercent = (((double)(maleCountCase))/((double)(totalNoOfCases)))*100;
	femalePercent = (((double)(femaleCountCase))/((double)(totalNoOfCases)))*100;
	unknownSexPercent = (((double)(unknownCountCase))/((double)(totalNoOfCases)))*100;
	
	out("Number of cases: "); out(totalNoOfCases);
	out("\n          Males: "); out(maleCountCase); out(" ("); out(malePercent); out("%)");
	out("\n          Females: "); out(femaleCountCase); out(" ("); out(femalePercent); out("%)");
	out("\n          Unknown sex: "); out(unknownCountCase); out(" ("); out(unknownSexPercent); out("%)\n\n");
	
	unsigned int totalNoOfPsCons = maleCountPseudo+femaleCountPseudo+unknownCountPseudo;

	malePercent = (((double)(maleCountPseudo))/((double)(totalNoOfPsCons)))*100;
	femalePercent = (((double)(femaleCountPseudo))/((double)(totalNoOfPsCons)))*100;
	unknownSexPercent = (((double)(unknownCountPseudo))/((double)(totalNoOfPsCons)))*100;
	
	out("Number of pseudocontrols: "); out(totalNoOfPsCons);
	out("\n          Males: "); out(maleCountPseudo); out(" ("); out(malePercent); out("%)");
	out("\n          Females: "); out(femaleCountPseudo); out(" ("); out(femalePercent); out("%)");
	out("\n          Unknown sex: "); out(unknownCountPseudo); out(" ("); out(unknownSexPercent); out("%)\n\n");
	
};

//! Outputs details of a subject to screen.
void Subject::outputDetails()
{
	out("Subject id: "); out(id); out("\n");
	out("Father id: "); out(fatherId); out("\n");
	out("Mother id: "); out(motherId); out("\n");
	out("Sex id: "); out(sex); out("\n");
	out("Affected id: "); out(affected); out("\n\n");
};

//! Returns a subject in a pedigree
Subject * Pedigree::getSubject(const unsigned int & id) const
{
	map<unsigned int, Subject *>::const_iterator s = subjects.find(id);
	if(s != subjects.end()) return s->second;

	outErr("Error: subject "); outErr(" not found in pedigree!\n");
	exit(1);	       
};

//! Returns a subject in a pedigree and returns if it exists. 
Subject * Pedigree::getSubject(const unsigned int & id, bool & exists) const
{
	map<unsigned int, Subject *>::const_iterator s = subjects.find(id);
	if(s != subjects.end())
	{
		exists = true;
		return s->second;
	};

	exists = false;
	return 0;       
};

//! Checks if a subject is a child of a case parent trio in a pedigree.
Trio * Pedigree::findCaseParentTrio(bool & found, Subject * subject)
{
	bool fatherExists = false;
	bool motherExists = false;
	
	//check if subject is affected, has genotype data and has not been removed if allowing extra trios 
	if(subject->getAffected() /*&& subject->genotypeDataPresent(snpId)*/ && subject->isNotRemoved())
	{
		//check if parents exist and have genotype data
		if(subject->getFatherId() != 0 && subject->getMotherId() != 0)
		{
			Subject * father = getSubject(subject->getFatherId(), fatherExists);
			Subject * mother = getSubject(subject->getMotherId(), motherExists);
			if(fatherExists && motherExists /*&& father->isNotRemoved() && mother->isNotRemoved()*/)
			{
				//if(father->genotypeDataPresent(snpId) && mother->genotypeDataPresent(snpId))
				{ 
					Trio * aTrio = new Trio(pedIdName, father, mother, subject);					
					found = true;
					return aTrio;
				};
			};
		};
	};

	return 0;
};

//! Loops thro' subjects in pedigree to find a case parent trio.
Trio * Pedigree::findCaseParentTrio(bool & found)
{
	found = false;
	Trio * aTrio;	

	//loop through subjects looking for a case which has parents
	for(map<unsigned int, Subject *>::const_iterator s = subjects.begin(); (s != subjects.end() && !found); ++s)
	{
		 aTrio = findCaseParentTrio(found, s->second);	
	};

	return aTrio;
};

//! Create the pseudocontrols and add basic info but no genotype data yet.
void Trio::createPseudoControls(unsigned int & pseudoControlSetNo, unsigned int & noPseudoCons, bool & cepg, MapIds & subjectIds, ofstream & psFamFile, bool & interaction)
{
	Subject * aPseudoControl;
	string pseudoConIDName;
	unsigned int pseudoConID;
	unsigned int patID = 0, matID = 0;
	unsigned int sexID = child->getSexId();
	string childIDName = child->getIDName();
	bool affected = false; //not affected for control
	unsigned noPseudoConsUpdate = noPseudoCons; //number of pseudocontrols may need to be increased if assuming CEPG
	if(cepg)
	{
		if(noPseudoCons == 1) noPseudoConsUpdate = 3;
		else if(noPseudoCons == 3) noPseudoConsUpdate = 7;
		else if(noPseudoCons == 15) noPseudoConsUpdate = 31;
	};

	//output case to fam file
	psFamFile << pedIdName << " " << childIDName << " " << patID << " " << matID << " " << sexID << " 2\n";

	//for info file
	pseudoConSetNo = pseudoControlSetNo;
	
	for(unsigned int psNo = 1; psNo <= noPseudoConsUpdate; ++psNo)
	{
		pseudoConIDName = childIDName + "-pseudo-" + toString(psNo);
		pseudoConID = subjectIds.getId(pseudoConIDName);		

		if(interaction) aPseudoControl = new SubjectTwoSnps(pseudoConID, patID, matID, sexID, affected, pseudoConIDName);
		else aPseudoControl = new SubjectOneSnp(pseudoConID, patID, matID, sexID, affected, pseudoConIDName);

		pseudoControls.push_back(aPseudoControl);

		//output to fam file
		psFamFile << pedIdName << " " << pseudoConIDName << " " << patID << " " << matID << " " << sexID << " 1\n";
	};
	
};

//! Write binary data.
void writeBedGenotype(ofstream & psBedFile, unsigned int & psBitCount, int & aBit, const bool & allele1, const bool & allele2)
{

	aBit = aBit >> 1; //shift bits to the right
	if(allele1) aBit = aBit | 128; 
	aBit = aBit >> 1;  //shift bits to the right
	if(allele2) aBit = aBit | 128;
	
	psBitCount += 2;

	//write to file if byte is finished
	if(psBitCount == 8)
	{
		//write to file
		char buffer[1];
		buffer[0] = aBit;
		psBedFile.write(buffer, 1);
		psBitCount = 0;
		aBit = 0;
	};

};

//! Checks for a Medelian error
bool checkForMendelianError(pair<bool, bool> & fatherGeno, pair<bool, bool> & motherGeno, pair<bool, bool> & childGeno)
{
	//suppose that first child allele is from father and second from mother, then vice versa
	bool isOK = ((childGeno.first == fatherGeno.first || childGeno.first == fatherGeno.second) && (childGeno.second == motherGeno.first || childGeno.second == motherGeno.second))
		|| ((childGeno.second == fatherGeno.first || childGeno.second == fatherGeno.second) && (childGeno.first == motherGeno.first || childGeno.first == motherGeno.second));

	return !isOK;
};

//! Given trio genotypes returns 1 psuedo control based on the alleles that were not transmitted.
pair<bool, bool> getPseudoControlAllele(pair<bool, bool> & fatherGeno, pair<bool, bool> & motherGeno, pair<bool, bool> & childGeno)
{
	if((fatherGeno.first == 1 && fatherGeno.second == 0) || (motherGeno.first == 1 && motherGeno.second == 0)
		|| (childGeno.first == 1 && childGeno.second == 0) || checkForMendelianError(fatherGeno, motherGeno, childGeno)) //if mother or father data is missing then so is the pseudocontrol
	{
			return make_pair(true, false);
	};
	
	bool psAllele1, psAllele2;

	//add parent genotypes to set and then remove the ones that were transmitted
	multiset<bool> alleles;
	multiset<bool>::const_iterator i;
	alleles.insert(fatherGeno.first);
	alleles.insert(fatherGeno.second);
	alleles.insert(motherGeno.first);
	alleles.insert(motherGeno.second);

	i = alleles.find(childGeno.first);
	if(i != alleles.end()) alleles.erase(i); 
	i = alleles.find(childGeno.second);
	if(i != alleles.end()) alleles.erase(i);

	multiset<bool>::const_iterator a = alleles.begin();
	psAllele1 = *a; a++;
	psAllele2 = *a;

	return make_pair(psAllele1, psAllele2);
};


//! Given trio genotypes returns 3 psuedo controls assuming CEPG for -pc1.
list<pair<bool, bool> > getPseudoControlGenotypeCEPG3(pair<bool, bool> & fatherGeno, pair<bool, bool> & motherGeno, pair<bool, bool> & childGeno)
{
	list<pair<bool, bool> > ans;
	
	//add genotypes with alleles that were not passed on twice plus case genotypes again
	if((fatherGeno.first == 1 && fatherGeno.second == 0) || (motherGeno.first == 1 && motherGeno.second == 0)
		|| (childGeno.first == 1 && childGeno.second == 0) || checkForMendelianError(fatherGeno, motherGeno, childGeno)) //if mother or father data is missing then so is the pseudocontrol
	{
		for(unsigned int k = 1; k <= 3; ++k) ans.push_back(make_pair(true, false));
	}
	else
	{
		pair<bool, bool> aPair;

		//get untransmitted genotype firstly 
		aPair = getPseudoControlAllele(fatherGeno, motherGeno, childGeno);

		ans.push_back(aPair);

		//add case again as the pseudocontrol with reversed alleles - exept the order will not be known so just add it the same
		ans.push_back(childGeno);

		//add the untransmitted alleles with reversed alleles - exept the order will not be known so just add it the same 
		ans.push_back(aPair);
	};

	return ans;
};

//! Given trio genotypes returns 3 psuedo controls based on the 3 genotypes that were not transmitted.
list<pair<bool, bool> > getPseudoControlGenotype(pair<bool, bool> & fatherGeno, pair<bool, bool> & motherGeno, pair<bool, bool> & childGeno)
{
	list<pair<bool, bool> > ans;
	
	//add all possible genotypes that could be passed on and except the one that was
	if((fatherGeno.first == 1 && fatherGeno.second == 0) || (motherGeno.first == 1 && motherGeno.second == 0)
		|| (childGeno.first == 1 && childGeno.second == 0) || checkForMendelianError(fatherGeno, motherGeno, childGeno)) //if mother or father data is missing then so is the pseudocontrol
	{
		for(unsigned int k = 1; k <= 3; ++k) ans.push_back(make_pair(true, false));
	}
	else
	{
		//add genotypes in order and do not add the child genotype		
		bool skipped = false;
		pair<bool, bool> aPair;

		unsigned int allele2CountChild = childGeno.first + childGeno.second;
		unsigned int allele2Count = fatherGeno.first + motherGeno.first;

		if(allele2Count != allele2CountChild)
		{
				aPair = make_pair(fatherGeno.first, motherGeno.first);
				if(aPair.first && !aPair.second) aPair = make_pair(false, true);
				ans.push_back(aPair);				
		} else skipped = true;

		allele2Count = fatherGeno.first + motherGeno.second;

		if(skipped || (allele2Count != allele2CountChild))
		{
			    aPair = make_pair(fatherGeno.first, motherGeno.second);
				if(aPair.first && !aPair.second) aPair = make_pair(false, true);
				ans.push_back(aPair);						
		} else skipped = true;

		allele2Count = fatherGeno.second + motherGeno.first;

		if(skipped || (allele2Count != allele2CountChild))
		{
				aPair = make_pair(fatherGeno.second, motherGeno.first);
				if(aPair.first && !aPair.second) aPair = make_pair(false, true);
				ans.push_back(aPair);						
		} else skipped = true;

		allele2Count = fatherGeno.second + motherGeno.second;

		if(skipped || (allele2Count != allele2CountChild))
		{
				aPair = make_pair(fatherGeno.second, motherGeno.second);
				if(aPair.first && !aPair.second) aPair = make_pair(false, true);
				ans.push_back(aPair);				
		};

	};

	return ans;
};

//! Given trio genotypes returns 7 psuedo controls assuming CEPG for -pc1.
list<pair<bool, bool> > getPseudoControlGenotypeCEPG7(pair<bool, bool> & fatherGeno, pair<bool, bool> & motherGeno, pair<bool, bool> & childGeno)
{
	list<pair<bool, bool> > ans, somePairs;
	
	//add all possible genotypes that could be passed on and except the one that was
	if((fatherGeno.first == 1 && fatherGeno.second == 0) || (motherGeno.first == 1 && motherGeno.second == 0)
		|| (childGeno.first == 1 && childGeno.second == 0) || checkForMendelianError(fatherGeno, motherGeno, childGeno)) //if mother or father data is missing then so is the pseudocontrol
	{
		for(unsigned int k = 1; k <= 7; ++k) ans.push_back(make_pair(true, false));
	}
	else
	{
		//add the 3 pseudocontrols as for CPG firstly
		somePairs = getPseudoControlGenotype(fatherGeno, motherGeno, childGeno);
		ans = somePairs;

		//add the case again as the pseudocontrol with reversed alleles - exept the order will not be known so just add it the same
		ans.push_back(childGeno);

		//add the untransmitted genotypes with reversed alleles - exept the order will not be known so just add it the same		
		ans.splice(ans.end(), somePairs);
	};

	return ans;
};

//! Adds SNP data to the psuedo controls based on the trio.
void Trio::updatePsuedoControlSNPData(unsigned int & noPseudoCons, bool & cepg, ofstream & psBedFile, unsigned int & psBitCount, int & aBit)
{
	pair<bool, bool> fatherGeno = father->getGenotype();
	pair<bool, bool> motherGeno = mother->getGenotype();
	pair<bool, bool> childGeno = child->getGenotype();	
	list<Subject * >::iterator pc = pseudoControls.begin();

	list<pair<bool, bool> > pseudoControlsGeno;
	unsigned int noPseudoConsUpdate = noPseudoCons;

	//write case genotypes firstly
	writeBedGenotype(psBedFile, psBitCount, aBit, childGeno.first, childGeno.second);

	//get the pseudocontrol data	
	if(noPseudoCons == 1)
	{
		if(cepg)
		{
			noPseudoConsUpdate = 3;
			pseudoControlsGeno = getPseudoControlGenotypeCEPG3(fatherGeno, motherGeno, childGeno); 
		}
		else
		{
			pair<bool, bool> pseudoControlGeno = getPseudoControlAllele(fatherGeno, motherGeno, childGeno);
			pseudoControlsGeno.push_back(pseudoControlGeno);
		};
	}
	else // (noPseudoCons == 3)
	{
		if(cepg)
		{
			noPseudoConsUpdate = 7;
			pseudoControlsGeno = getPseudoControlGenotypeCEPG7(fatherGeno, motherGeno, childGeno); 
		}
		else
		{
			pseudoControlsGeno = getPseudoControlGenotype(fatherGeno, motherGeno, childGeno); 
		};
	};
		
	//add the data
	list<pair<bool, bool> >::const_iterator g = pseudoControlsGeno.begin();

	for(unsigned int psNo = 1; psNo <= noPseudoConsUpdate; psNo++)
	{
		if(pc == pseudoControls.end())
		{
			outErr("Failed to create pseudocontrol for subjects ");
			outErr(father->getIDName()); outErr(", ");
			outErr(mother->getIDName()); outErr(" and ");
			outErr(child->getIDName()); outErr("!\n");
			exit(1);
		};

		(*pc)->addGenotype(g->first, g->second);
		writeBedGenotype(psBedFile, psBitCount, aBit, g->first, g->second);			

		++pc; ++g;
	};

};

//! Removes SNP data into lists to output in info file.
void Trio::storeSNPDataForInfoFile(const unsigned int & geno)
{
	pair<bool, bool> fatherGeno = father->getGenotype(geno);
	pair<bool, bool> motherGeno = mother->getGenotype(geno);
	pair<bool, bool> childGeno = child->getGenotype(geno);

	fatherGenoTypes.push_back(fatherGeno.first); fatherGenoTypes.push_back(fatherGeno.second);
	motherGenoTypes.push_back(motherGeno.first); motherGenoTypes.push_back(motherGeno.second);
	childGenoTypes.push_back(childGeno.first); childGenoTypes.push_back(childGeno.second);

	pair<bool, bool> psGeno;

	//create list of ps control SNP data if it does not already exist
	if(pseudoControlGenotypes.empty())
	{
		list<bool> emptyList;
		for(unsigned int i = 1; i <= pseudoControls.size(); ++i)
		{
			pseudoControlGenotypes.push_back(emptyList);
		};
	};

	list<list<bool> >::iterator psg = pseudoControlGenotypes.begin();
	for(list<Subject * >::const_iterator ps = pseudoControls.begin(); ps != pseudoControls.end(); ++ps, ++psg)
	{
		psGeno = (*ps)->getGenotype(geno);
		(*psg).push_back(psGeno.first); (*psg).push_back(psGeno.second);
	};
};

//! Output genotype info.
void outputGenotypeInfo(ofstream & infoFile, MapSnpAlleleIds & mapSnpAlleleIds, const list<bool> & genotypes)
{
	unsigned int snpID = 1;
	map<bool, string> snpNames;
	string alleleName1, alleleName2;
	bool al1, al2;
	for(list<bool>::const_iterator cgt = genotypes.begin(); cgt != genotypes.end(); )
	{
		al1 = *cgt; ++cgt;
		al2 = *cgt; ++cgt;

		if(al1 && !al2) infoFile << " 0 0"; //missing data
		else
		{
			snpNames = mapSnpAlleleIds.getSnpAlleleNames(snpID);
			alleleName1 = snpNames.find(al1)->second;		
			alleleName2 = snpNames.find(al2)->second;
			infoFile << " " << alleleName1 << " " << alleleName2;
		};
	};
};

//! Outputs extra info for info file.
void Trio::outputInfoFile(ofstream & infoFile, MapSnpAlleleIds & mapSnpAlleleIds, bool & infoFileFatherGeno, bool & infoFileMotherGeno, bool & interaction, bool & cepg, unsigned int & snpNo1, unsigned int & snpNo2)
{
	//output info case/child firstly
	infoFile << pedIdName << " " << child->getIDName() << " " << child->getIDName() << " " << pseudoConSetNo;

	//output Father info
	if(infoFileFatherGeno) outputGenotypeInfo(infoFile, mapSnpAlleleIds, fatherGenoTypes);	

	//output Mother info
	if(infoFileMotherGeno) outputGenotypeInfo(infoFile, mapSnpAlleleIds, motherGenoTypes);	

	//output SNP data for child
	if(infoFileFatherGeno || infoFileMotherGeno) outputGenotypeInfo(infoFile, mapSnpAlleleIds, childGenoTypes);	
	
	infoFile << "\n";

	//now output info for the corresponding pseudocontrols
	list<list<bool> >::const_iterator pcg = pseudoControlGenotypes.begin();

	//if CEPG then the genotypes of the parents are switched, so switch them in the info file also
	unsigned int noPs = pseudoControls.size();
	unsigned int pseudoNo = 1;

	for(list<Subject * >::const_iterator ps = pseudoControls.begin(); ps != pseudoControls.end(); ++ps, ++pcg)
	{
		//output data for pseudocon
		infoFile << pedIdName << " " << (*ps)->getIDName() << " " << child->getIDName() << " " << pseudoConSetNo;

		//output Father info
		if(infoFileFatherGeno)
		{
			if(!cepg || (noPs==3 && pseudoNo == 1) || (noPs==7 && pseudoNo <= 3) || (noPs==31 && pseudoNo <= 15))
				outputGenotypeInfo(infoFile, mapSnpAlleleIds, fatherGenoTypes);	
			else outputGenotypeInfo(infoFile, mapSnpAlleleIds, motherGenoTypes); 
		};

		//output Mother info
		if(infoFileMotherGeno)
		{
			if(!cepg || (noPs==3 && pseudoNo == 1) || (noPs==7 && pseudoNo <= 3) || (noPs==31 && pseudoNo <= 15))
				outputGenotypeInfo(infoFile, mapSnpAlleleIds, motherGenoTypes);	
			else outputGenotypeInfo(infoFile, mapSnpAlleleIds, fatherGenoTypes); 
		};

		//output SNP data for pseudocon
		if(infoFileFatherGeno || infoFileMotherGeno) outputGenotypeInfo(infoFile, mapSnpAlleleIds, (*pcg));	

		infoFile << "\n";

		pseudoNo++;
	};
	
};

//! Adds SNP data to the psuedo controls based on the trio.
void Trio::updatePsuedoControlJointSNPData(unsigned int & noPseudoCons, bool & cepg, ofstream & psBedFile, unsigned int & psBitCount, int & aBit)
{
	pair<bool, bool> fatherGeno = father->getGenotype(1);
	pair<bool, bool> motherGeno = mother->getGenotype(1);
	pair<bool, bool> childGeno = child->getGenotype(1);
	pair<bool, bool> fatherGeno2 = father->getGenotype(2);
	pair<bool, bool> motherGeno2 = mother->getGenotype(2);
	pair<bool, bool> childGeno2 = child->getGenotype(2);

	list<Subject * >::iterator pc = pseudoControls.begin();

	//write first case SNP now and write second SNP later after all the data for the first SNP has been written
	writeBedGenotype(psBedFile, psBitCount, aBit, childGeno.first, childGeno.second);	

	list<unsigned int > genotypes1; //allele 2 counts
	list<unsigned int > genotypes2;

	//if mother or father data is missing then so is the pseudocontrol 
	if((fatherGeno.first == 1 && fatherGeno.second == 0) || (motherGeno.first == 1 && motherGeno.second == 0)
		|| (fatherGeno2.first == 1 && fatherGeno2.second == 0) || (motherGeno2.first == 1 && motherGeno2.second == 0)
		|| (childGeno.first == 1 && childGeno.second == 0) || (childGeno2.first == 1 && childGeno2.second == 0)
		|| checkForMendelianError(fatherGeno, motherGeno, childGeno) || checkForMendelianError(fatherGeno2, motherGeno2, childGeno2)) 
	{			
		for(unsigned int i = 1; i <= 15; ++i)
		{				
			(*pc)->addGenotype(true, false, 1);
			(*pc)->addGenotype(true, false, 2);

			//write first SNP now and write second SNP later after all the data for the first SNP has been written
			writeBedGenotype(psBedFile, psBitCount, aBit, true, false);
				
			++pc;
		};	

		return;
	};
	
	//add genotype data to lists
	genotypes1.push_back(fatherGeno.first + motherGeno.first);		
	genotypes1.push_back(fatherGeno.first + motherGeno.second);			
	genotypes1.push_back(fatherGeno.second + motherGeno.first);			
	genotypes1.push_back(fatherGeno.second + motherGeno.second);
	
	//add genotype data to lists
	genotypes2.push_back(fatherGeno2.first + motherGeno2.first);		
	genotypes2.push_back(fatherGeno2.first + motherGeno2.second);			
	genotypes2.push_back(fatherGeno2.second + motherGeno2.first);			
	genotypes2.push_back(fatherGeno2.second + motherGeno2.second);
		
	unsigned int allele2CountChild1 = childGeno.first + childGeno.second;
	unsigned int allele2CountChild2 = childGeno2.first + childGeno2.second;	

	//create pseudocontrols for each pair of SNPs not transmitted
	bool skipped = false;
	unsigned int pseudoCount = 0;
	unsigned int noRepeats = cepg;

	for(unsigned int repeat = 0; repeat <= noRepeats; ++repeat) //repeat genotypes if assuming CEPG plus case
	{

		for(list<unsigned int >::const_iterator g1 = genotypes1.begin(); g1 != genotypes1.end(); ++g1)
		{
			for(list<unsigned int >::const_iterator g2 = genotypes2.begin(); g2 != genotypes2.end(); ++g2)
			{
				if(repeat == 0 && ((!skipped && (*g1 == allele2CountChild1) && (*g2 == allele2CountChild2)) ))						
				{
					skipped = true;
				}
				else
				{					
					if(pc == pseudoControls.end())
					{
						outErr("Failed to create pseudocontrol for subjects ");
						outErr(father->getIDName()); outErr(", ");
						outErr(mother->getIDName()); outErr(" and ");
						outErr(child->getIDName()); outErr("!\n");
						exit(1);
					};

					//write first SNP now and write second SNP later after all the data for the first SNP has been written
					if(*g1 == 2) {(*pc)->addGenotype(true, true, 1); writeBedGenotype(psBedFile, psBitCount, aBit, true, true);}
					else if(*g1 == 1) {(*pc)->addGenotype(false, true, 1); writeBedGenotype(psBedFile, psBitCount, aBit, false, true);}
					else {(*pc)->addGenotype(false, false, 1); writeBedGenotype(psBedFile, psBitCount, aBit, false, false);};					

					if(*g2 == 2) (*pc)->addGenotype(true, true, 2);
					else if(*g2 == 1) (*pc)->addGenotype(false, true, 2);
					else (*pc)->addGenotype(false, false, 2);					
										
					++pc;
					pseudoCount++;
				}; 

			};
		};

	};

};

//!
void Trio::writePseudoControlSecondSNPData(ofstream & psBedFile, unsigned int & psBitCount, int & aBit)
{
	//write second SNP case info first
	pair<bool, bool> childGeno2 = child->getGenotype(2);
	writeBedGenotype(psBedFile, psBitCount, aBit, childGeno2.first, childGeno2.second);

	pair<bool, bool> aGenotype;
	for(list<Subject * >::iterator p = pseudoControls.begin(); p != pseudoControls.end(); ++p)
	{
		aGenotype = (*p)->getGenotype(2);
		writeBedGenotype(psBedFile, psBitCount, aBit, aGenotype.first, aGenotype.second);	
	};
};

//! Unmark all subjects in pedigree for when adding extra trios.
void Pedigree::restoreAllSubjectsToPedigree()
{
	for(map<unsigned int, Subject *>::const_iterator s = subjects.begin(); s != subjects.end(); ++s)
	{
		s->second->restoreToPedigree();
	};
};

