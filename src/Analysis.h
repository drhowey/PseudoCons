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
 * along with PseudoCons. If not, see <http://www.gnu.org/licenses/>.
 ************************************************************************/


/*! \file Analysis.h
    \brief This file organises the analysis.
    
*/

#ifndef __ANALYSIS
#define __ANALYSIS

#include <string>
#include <fstream>
#include <map>

#include "Data.h"
#include "main.h"

using namespace std; // initiates the "std" or "standard" namespace


//! Organises the correct analysis to perform.
class Analysis
{
private:

	//parameters for the options of the analysis to perform	
	string inputFilename;		
	string outputFilename;
	string probandFilename;
	bool extraAffectedTrios;		
	unsigned int noOfSnps;
	unsigned int snpNo1, snpNo2;
	string snpName1, snpName2;
	bool interaction;
	bool cepg;
	string infoFilename;
	bool doInfo;
	bool infoFileFatherGeno;
	bool infoFileMotherGeno;
	unsigned int infoFileMaxSNPs;
	MapSnpAlleleIds mapSnpAlleleIds;

	ifstream inFileGeno;	
	ifstream readPedigree;

	//data used for processing the pedigrees
	MapIds subjectIds;
	MapIds pedigreeIds;	
	AllPedigreeData allPedigreeData;
	AllTrioData allTrioData;
	list<Subject *> orderListOfSubjects; //to use for binary SNP data
	map<unsigned int, set<unsigned int> > probandSubjectIds; //pedigree id, subject ids

	//variables for reading in binary data
	unsigned int bitCount;	
	int one;
	int aBit;	
	char buffer[1];
	int allele1, allele2;	

public:

	Analysis(string & i, string & of, string pro, bool & eat, unsigned int & npc, unsigned int & s1, unsigned int & s2, string & sn1, string & sn2, bool & c, string & ifn, bool & iffg, bool & ifmg, unsigned int & ifms)
		:  inputFilename(i), outputFilename(of), probandFilename(pro), extraAffectedTrios(eat), snpNo1(s1), snpNo2(s2), snpName1(sn1), snpName2(sn2), cepg(c),
		  infoFilename(ifn), infoFileFatherGeno(iffg), infoFileMotherGeno(ifmg), infoFileMaxSNPs(ifms)
	{	
		bitCount = 9;
		one = '\1';
		interaction = (snpNo2 != 0) || (snpName2 != "");		
		allTrioData.setNoPseudoConsAndCEPG(npc, cepg);
		doInfo = (infoFilename != "");
	};

	//! Delete analysis things
	~Analysis()
	{		
		
	};

	void runAnalysis();
	void setUpPedigreeData();
	void addPedigreeData();
	void extractTrioData();
	void setNoOfSnps();
	bool addProcessedPedigree(Pedigree * pedigree,	map<unsigned int, set<unsigned int> > & probandSubjectIds);
	bool addProcessedPedigreeProbandSubject(bool & found, const unsigned int & probandSubjectId, Pedigree * pedigree);
	void addGenotypeDataToSubject(Subject * subject, const unsigned int & geno = 1);
	void doNotAddGenotypeDataToSubject();
	void outputCasesAndPseudoControls();
};


#endif
