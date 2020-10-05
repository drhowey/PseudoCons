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



/*! \file Analysis.cpp
    \brief This file contains the methods the various analyse.
    
*/
#include <iostream>
#include <sstream>
#include <list>
#include <map>
#include <set>
#include <math.h>

using namespace std; // initiates the "std" or "standard" namespace

#include "Analysis.h"
#include "main.h"

//! Runs the chosen analysis.
void Analysis::runAnalysis()
{ 
	//setup the pedigree structure
	setUpPedigreeData();

	//get trios from the pedigrees from which psuedocontrols will be created
	extractTrioData();

	//loop thro' SNPs outputting pseudocontrol data with cases
	outputCasesAndPseudoControls();
};


//! Processes the pedigree file - determines pedigree file type for analysis.
void Analysis::setUpPedigreeData()
{	

	unsigned int length = (unsigned)inputFilename.length();
	string fileExtension = ".ped";
	if(length >= 4) fileExtension = inputFilename.substr(length-4,4);

	//determine the type of pedigree file
	if(fileExtension[0] == '.' &&
		(fileExtension[1] == 'b' || fileExtension[1] == 'B') &&
		(fileExtension[2] == 'e' || fileExtension[2] == 'E') &&
		(fileExtension[3] == 'd' || fileExtension[3] == 'D'))
	{		

		//open file in binary mode		
		inFileGeno.open(inputFilename.c_str(), ios::binary);
		if(!inFileGeno.is_open())
		{
			outErr("\nCannot read binary genotype file: "); outErr(inputFilename); outErr("!\n");
			exit(1);
		};
		
		char buffer[3];
		inFileGeno.read(buffer, 3);

		//check the plink magic numbers for the file type
		//3rd number indicates format of genotype data, 1 => subjects x SNPs, 0 => SNPs x subjects
		unsigned int magicNumber1 = buffer[0];
		unsigned int magicNumber2 = buffer[1];

		if(magicNumber1 != 108 || magicNumber2 != 27)
		{
			out("Detected an old version .bed file, please update this file using PLINK!\n");			
			inFileGeno.close();
			exit(1);
		};

		//determine binary file type
		unsigned int mode = buffer[2];
		if(mode == 0)
		{			
			out("Detected a .bed file in \"individual-major\" mode, please update this file using PLINK!\n");			
			inFileGeno.close();
			exit(1);
		};		
	
		bitCount = 9;

		//open corresponding pedigree file for the .bed file
		string famFileName = inputFilename.substr(0, length-4) + ".fam"; 
		readPedigree.open(famFileName.c_str());

		if(!readPedigree.is_open())
		{
			outErr("\nCannot read corresponding family file: "); outErr(famFileName); outErr("!\n");
			exit(1);
		};
	}
	else if(((fileExtension[0] == 'g' || fileExtension[0] == 'G') &&
		(fileExtension[1] == 'z' || fileExtension[1] == 'Z') &&
		(fileExtension[2] == 'i' || fileExtension[2] == 'I') &&
		(fileExtension[3] == 'p' || fileExtension[3] == 'P'))
		|| 
		(fileExtension[1] == '.' &&
		(fileExtension[2] == 'g' || fileExtension[1] == 'G') &&
		(fileExtension[3] == 'z' || fileExtension[2] == 'Z')))
	{
		outErr("gzipped files are not handled, please use a .bed file instead!\n");
		exit(1);
	};
	
	//count the number of SNPs and check the chromosome is in 1-22
	setNoOfSnps();

	//check not too many SNPs for info file
	if(doInfo && (infoFileFatherGeno || infoFileMotherGeno) && !interaction && noOfSnps > infoFileMaxSNPs)
	{
		outErr("\nNumber of SNPs ("); outErr(noOfSnps); outErr(") exceeds max allowed for info file ("); outErr(infoFileMaxSNPs); outErr(")!\n\n");
		outErr("Change max allowed with \"-info-maxsnps n\" at your own peril!\n\n");
		exit(1);
	};

	//output SNP column info for info file
	if(doInfo && (infoFileFatherGeno || infoFileMotherGeno))
	{
		unsigned int noSNPs = noOfSnps;
		if(interaction) noSNPs = 2;

		unsigned int startSNP1 = 5;
		unsigned int endSNP1 = 5 + 2*noSNPs - 1;
		unsigned int startSNP2 = 5 + 2*noSNPs;
		unsigned int endSNP2 = 5 + 4*noSNPs - 1;
		unsigned int startSNP3 = 5 + 4*noSNPs;
		unsigned int endSNP3 = 5 + 6*noSNPs - 1;
		string firstGeno, secondGeno, thirdGeno;
		if(infoFileFatherGeno && infoFileMotherGeno) {firstGeno = "father"; secondGeno = "mother"; thirdGeno = "child/pseudocontrol";}
		else if(infoFileMotherGeno) {firstGeno = "mother"; secondGeno = "child/pseudocontrol"; thirdGeno = "";}
		else {firstGeno = "father"; secondGeno = "child/pseudocontrol"; thirdGeno = "";};

		out("          Columns "); out(startSNP1); out("-"); out(endSNP1); out(": "); out(firstGeno); out(" genotype info\n");
		out("          Columns "); out(startSNP2); out("-"); out(endSNP2); out(": "); out(secondGeno); out(" genotype info\n");
		if(thirdGeno != "")
		{
			out("          Columns "); out(startSNP3); out("-"); out(endSNP3); out(": "); out(thirdGeno); out(" genotype info\n");
		};
	};

	//now add the pedigree structure
	addPedigreeData();
};

bool isChromosome1to22(string & chromosome)
{
	if(chromosome == "1" || chromosome == "2" || chromosome == "3" 
		|| chromosome == "4" || chromosome == "5" || chromosome == "6" 
		|| chromosome == "7" || chromosome == "8" || chromosome == "9" 
		|| chromosome == "10" || chromosome == "11" || chromosome == "12" 
		|| chromosome == "13" || chromosome == "14" || chromosome == "15" 
		|| chromosome == "16" || chromosome == "17" || chromosome == "18" 
		|| chromosome == "19" || chromosome == "20" || chromosome == "21" || chromosome == "22")
		return true;

	return false;
};

//! Determine how many SNPs there are from the .bim file and check chromosomes.
void Analysis::setNoOfSnps()
{
	bool bim = false;
	unsigned int snpCount = 0;
	string chromosome, snpIdentifier, geneticDistance, basePairPosition;
	string alleleName1, alleleName2;
	string prevSnpIdentifier = "";
	
	unsigned int length = (unsigned)inputFilename.length();	
	string mapFileName = inputFilename.substr(0, length-4) + ".bim";	

	
	ifstream readSnpFile(mapFileName.c_str());
	if(!readSnpFile.is_open())
	{
		outErr("\nCannot read SNP file: "); outErr(mapFileName); outErr("!\n");
		exit(1);
	};

	//output a copy of bim file for the pseudocontrols also
	unsigned int lengthOut = (unsigned)outputFilename.length();
	string fileExtension = "";
	if(lengthOut >= 4) fileExtension = outputFilename.substr(lengthOut-4, 4);

	string psBimFileName;
	if(fileExtension == ".bed" || fileExtension == ".BED") psBimFileName = outputFilename.substr(0, lengthOut-4) + ".bim";	
	else psBimFileName = outputFilename + ".bim";

	ofstream psBimFile(psBimFileName.c_str());	
	unsigned int noSNPsWrittenInter = 0;

	do{
		//read in .bim file lines
		readSnpFile >> chromosome >> snpIdentifier >> geneticDistance >> basePairPosition >> alleleName1 >> alleleName2;		
	
		if(snpIdentifier != prevSnpIdentifier)
		{
				snpCount++;
				if(interaction)
				{
					if(snpNo1 == snpCount || snpName1 == snpIdentifier)
					{
						snpNo1 = snpCount;
						psBimFile << chromosome << " " << snpIdentifier << " " << geneticDistance << " " << basePairPosition << " " << alleleName1 << " " <<  alleleName2 << "\n";			
						noSNPsWrittenInter++;
					};

					if(snpNo2 == snpCount || snpName2 == snpIdentifier)
					{
						snpNo2 = snpCount;
						psBimFile << chromosome << " " << snpIdentifier << " " << geneticDistance << " " << basePairPosition << " " << alleleName1 << " " <<  alleleName2 << "\n";			
						noSNPsWrittenInter++;
					};
				}
				else
				{
					//write for psuedo control
					psBimFile << chromosome << " " << snpIdentifier << " " << geneticDistance << " " << basePairPosition << " " << alleleName1 << " " <<  alleleName2 << "\n";			
				};

				//need to keep track of SNP names if outputting the info file which may include SNP data 
				if(infoFileFatherGeno || infoFileMotherGeno)
				{
					map<bool, string> alleleNames;
					alleleNames[false] =  alleleName1; 
					alleleNames[true] =  alleleName2;
					mapSnpAlleleIds.addSNPAlleleNames(snpCount, alleleNames);
				};
		};

		prevSnpIdentifier = snpIdentifier;
		
		if(!isChromosome1to22(chromosome))
		{
			outErr("PseudoCons is not designed for use with chromosome "); outErr(chromosome); outErr("!\n\n");
			outErr("Please remove all chromosome "); outErr(chromosome); outErr(" SNPs from your input files!\n\n");
			exit(1);
		};

	}while(!readSnpFile.eof());

	noOfSnps = snpCount;

	if(interaction && snpNo1 == 0)
	{
		outErr("SNP not found with name: "); outErr(snpName1); outErr("!\n");
		exit(1);
	};

	if(interaction && snpNo2 == 0)
	{
		outErr("SNP not found with name: "); outErr(snpName2); outErr("!\n");
		exit(1);
	};

	if(interaction && noSNPsWrittenInter != 2)
	{		
		outErr("SNP numbers "); outErr(snpNo1); outErr(" and "); outErr(snpNo2); outErr(" are not valid!\n");
		exit(1);
	};

	readSnpFile.close();
	psBimFile.close();
};

//add the IDs of the subjects in the proband file first so that they appear
//in the list of subjects first and are then given priority chosen for analysis
void addProbandSubjectIds(string & probandFileName, MapIds & pedigreeIds, MapIds & subjectIds, map<unsigned int, set<unsigned int> > & probandSubjectIds)
{
	ifstream probandFile;
	string pedigreeName, subjectIdName, subjectIdName2;
	unsigned int subjectId, pedigreeId;
	set<unsigned int> pedSubjectIds;

	probandFile.open(probandFileName.c_str());	
	
	if(!probandFile.is_open())
	{
		outErr("Proband file "); outErr(probandFileName); outErr(" not found!\n");
		exit(1);
	};		

	do{
		probandFile >> pedigreeName >> subjectIdName;

		if(probandFile.eof()) break;

		subjectIdName2 = pedigreeName + "-" + subjectIdName;
		subjectId = subjectIds.getId(subjectIdName2);
		pedigreeId = pedigreeIds.getId(pedigreeName);

		//add subject Id to the coressponding pedigree Id
		map<unsigned int, set<unsigned int> >::iterator psi = probandSubjectIds.find(pedigreeId);
		if(psi != probandSubjectIds.end())
		{
			psi->second.insert(subjectId);
		}
		else
		{
			pedSubjectIds.clear();
			pedSubjectIds.insert(subjectId);
			probandSubjectIds[pedigreeId] = pedSubjectIds;//pedigreeName + " " + subjectIdName; 
		};
	}while(!probandFile.eof());
	
	probandFile.close();
};

//! Returns the sex ID corresponding to the sex ID in file.
unsigned int getSexId(string & sexIdName)
{
	if(sexIdName == "1") return 1;
	else if(sexIdName == "2") return 2;

	return 0;
};

//! Returns the affect status corresponding to the affected value in file.
bool getAffected(string & affectedIdName)
{
	if(affectedIdName == "1") return false;
	else if(affectedIdName == "2") return true;

	return false;
};

//! Add pedigree data from the pedigree file.
void Analysis::addPedigreeData()
{	
	string pedigreeName, subjectIdName, subjectIdNameOrig, fatherIdName, motherIdName, sexIdName, affectedIdName;
	unsigned int subjectId, fatherId, motherId, pedigreeId, sexId;
	bool affected;	

	//open file to record the proband
	if(probandFilename != "") addProbandSubjectIds(probandFilename, pedigreeIds, subjectIds, probandSubjectIds);

	Subject * subject;
	
	unsigned int snpId = 1;
	unsigned int count = 0;
	unsigned int totalNoOfSubjects = 0;
	unsigned int totalNoOfMales = 0;
	unsigned int totalNoOfFemales = 0;
	unsigned int totalNoOfUnknownSex = 0;
	unsigned int totalNoOfSubjectsAffected = 0;
	
	//read in pedigree data
	do{

		//read in all of the details of a subject and map Ids of parents and the subject
		readPedigree >> pedigreeName >> subjectIdNameOrig >> fatherIdName >> motherIdName >> sexIdName >> affectedIdName;

		if(readPedigree.eof()) break; //check if the last row has been past
		
		pedigreeId = pedigreeIds.getId(pedigreeName);

		subjectIdName = pedigreeName + "-" + subjectIdNameOrig;
		subjectId = subjectIds.getId(subjectIdName);

		if(fatherIdName == "0") fatherId = 0;
		else
		{
			fatherIdName = pedigreeName + "-" + fatherIdName;
			fatherId = subjectIds.getId(fatherIdName);
		};

		if(motherIdName == "0") motherId = 0;
		else
		{
			motherIdName = pedigreeName + "-" + motherIdName;
			motherId = subjectIds.getId(motherIdName);
		};

		sexId = getSexId(sexIdName);

		affected = getAffected(affectedIdName);
	
		if(affected) totalNoOfSubjectsAffected++;
		
		//create a subject		
		if(interaction)	subject = new SubjectTwoSnps(subjectId, fatherId, motherId, sexId, affected, subjectIdNameOrig);
		else subject = new SubjectOneSnp(subjectId, fatherId, motherId, sexId, affected, subjectIdNameOrig);

		orderListOfSubjects.push_back(subject);		

		//add the subject to the pedigree to which they belong
		allPedigreeData.addSubjectToPedigree(subjectId, subject, pedigreeId, pedigreeName);

		totalNoOfSubjects++;
		if(sexId == 1) totalNoOfMales++;
		else if(sexId == 2) totalNoOfFemales++;
		else totalNoOfUnknownSex++;

		
	}while(!readPedigree.eof());

	double femalePercent, malePercent, unknownSexPercent, affectedPercent;
	malePercent = (((double)(totalNoOfMales))/((double)(totalNoOfSubjects)))*100;
	femalePercent = (((double)(totalNoOfFemales))/((double)(totalNoOfSubjects)))*100;
	unknownSexPercent = (((double)(totalNoOfUnknownSex))/((double)(totalNoOfSubjects)))*100;
	affectedPercent = (((double)(totalNoOfSubjectsAffected))/((double)(totalNoOfSubjects)))*100;

	out("\nNumber of subjects: "); out(totalNoOfSubjects);
	out("\n          Males: "); out(totalNoOfMales); out(" ("); out(malePercent); out("%)");
	out("\n          Females: "); out(totalNoOfFemales); out(" ("); out(femalePercent); out("%)");
	out("\n          Unknown sex: "); out(totalNoOfUnknownSex); out(" ("); out(unknownSexPercent); out("%)");
	out("\n          Affected: "); out(totalNoOfSubjectsAffected); out(" ("); out(affectedPercent); out("%)");
	out("\n          Unaffected: "); out(totalNoOfSubjects-totalNoOfSubjectsAffected); out(" ("); out(100-affectedPercent); out("%)");
	out("\nNumber of SNPs: "); out(noOfSnps); out("\n\n");
	
};


//! Adds a proband subject to counted data as a trio or something else searching for the best subset in the pedigree
bool Analysis::addProcessedPedigreeProbandSubject(bool & found, const unsigned int & probandSubjectId, Pedigree * pedigree)
{
	//find the proband subject in the pedigree 
	Subject * probandSubject = pedigree->getSubject(probandSubjectId);
	if(probandSubjectId == 0) return false;

	//try to find a case parent trio
	Trio * aTrio;
	aTrio = pedigree->findCaseParentTrio(found, probandSubject);

	if(found)
	{
		allTrioData.addTrio(aTrio);		

		if(extraAffectedTrios)
		{			
			aTrio->getChild()->removeFromPedigree();
		};
	};

	if(found) return true;
	
	return false;
};

//! Analyses a given pedigree to find the best trio to count.
bool Analysis::addProcessedPedigree(Pedigree * pedigree, map<unsigned int, set<unsigned int> > & probandSubjectIds)
{
	bool found = false;	
	bool foundParentCons = false;

	//add proband subject first
	map<unsigned int, set<unsigned int> >::const_iterator ps = probandSubjectIds.find(pedigree->getPedId());	

	if(ps != probandSubjectIds.end())
	{
		//add all proband subjects for one pedigree
		for(set<unsigned int>::const_iterator pro = ps->second.begin(); pro != ps->second.end(); ++pro)
		{			
			addProcessedPedigreeProbandSubject(found, *pro, pedigree);
		};
		
		//if proband subject(s) exist then we are done - unless we have chosen to find extras and there may be some to find
		if(!extraAffectedTrios) return found;
	};

	bool foundOne;
	Trio * aTrio;

	//try to find an affected child and parents, if extra trios are allowed, keep searching while removing found trios from search					
	do{
		foundOne = false;
		aTrio = pedigree->findCaseParentTrio(foundOne);
		if(foundOne)
		{
			allTrioData.addTrio(aTrio);
			found = true;

			if(extraAffectedTrios)
			{					
				aTrio->getChild()->removeFromPedigree();
			};
		};

	}while(extraAffectedTrios && foundOne);

	if(found) return true;
	
	return false;
};

void Analysis::extractTrioData()
{
	unsigned int pedigreesNotCounted = 0;
	map<unsigned int, unsigned int> noOfpedigreeSubjects; //pedigree id, subjects count
	
	//check all the proband subjects
	for(map<unsigned int, set<unsigned int> >::const_iterator psids = probandSubjectIds.begin(); psids != probandSubjectIds.end(); ++psids)
	{
		for(set<unsigned int>::const_iterator s = psids->second.begin(); s != psids->second.end(); ++s)
		{
			allPedigreeData.subjectExistsAndAffected(psids->first, *s, subjectIds);
		};
	};
			
	//loop through the pedigrees adding trios
	for(map<unsigned int, Pedigree *>::iterator ped = allPedigreeData.pedigrees.begin(); ped != allPedigreeData.pedigrees.end(); ++ped)
	{				
		if(!(addProcessedPedigree(ped->second, probandSubjectIds))) pedigreesNotCounted++;
	};

	//allPedigreeData.restoreAllSubjectsToAllPedigrees(); //only needed if doing SNP by SNP
	allPedigreeData.outputSummary();

	unsigned int totalNoOfPedigrees = (unsigned)allPedigreeData.pedigrees.size();
	unsigned int totalNoOfTrios = (unsigned)allTrioData.trios.size();
	
	out("Number of trios used to create pseudocontrols: "); out(totalNoOfTrios); out("\n");
	out("Number of pedigrees with no pseudocontrols: "); out(pedigreesNotCounted); out("\n\n");
};
	

//! Sets the genotype data for the given subject and SNP for a binary pedigree file.
void Analysis::addGenotypeDataToSubject(Subject * subject, const unsigned int & geno)
{
	
	//read in the next piece of data
	if(bitCount == 9)
	{
		
		inFileGeno.read(buffer, 1);
		if(inFileGeno.eof())
		{			
			outErr("Error: something is up with reading the binary (.bed) SNP file!\n");
			exit(1);
		};

		aBit = buffer[0];

		bitCount = 1;
	};

	allele1 = aBit & one; //read the least significant bit				
	aBit = aBit >> 1; //shift bits to the right
	allele2 = aBit & one; //read the new least significant bit				
	aBit = aBit >> 1; //shift bits to the right for next time

	bitCount += 2;	

	//add genotype data to the subject, let 0 = allele 1 in BIM file, 1 = allele 2 in BIM file
	subject->addGenotype((allele1 == 1), (allele2 == 1), geno);
	
};

//! Does not set the genotype data for the given subject but does read in unwanted SNP data.
void Analysis::doNotAddGenotypeDataToSubject()
{
	//read in the next piece of data
	if(bitCount == 9)
	{
		
		inFileGeno.read(buffer, 1);
		if(inFileGeno.eof())
		{			
			outErr("Error: something is up with reading the binary (.bed) SNP file!\n");
			exit(1);
		};

		aBit = buffer[0];

		bitCount = 1;
	};

	allele1 = aBit & one; //read the least significant bit				
	aBit = aBit >> 1; //shift bits to the right
	allele2 = aBit & one; //read the new least significant bit				
	aBit = aBit >> 1; //shift bits to the right for next time

	bitCount += 2;	
};

//! Creates and outputs data for the pseudocontrols.
void Analysis::outputCasesAndPseudoControls()
{
	unsigned int length = (unsigned)outputFilename.length();
	string fileExtension = "";
	if(length >= 4) fileExtension = outputFilename.substr(length-4, 4);
	string bedFileName, famFileName;

	if(fileExtension == ".bed" || fileExtension == ".BED")
	{
		bedFileName = outputFilename;		
		famFileName = outputFilename.substr(0, length-4) + ".fam";
	}
	else
	{
		bedFileName = outputFilename + ".bed"; 		
		famFileName = outputFilename + ".fam"; 
	};
	
	//create pseudocontrol subjects and output family file with cases
	ofstream psFamFile(famFileName.c_str());
	allTrioData.createPseudoControls(subjectIds, psFamFile, interaction);
	psFamFile.close();

	ofstream psBedFile(bedFileName.c_str(), ios::binary);
	unsigned int psBitCount = 0;
	int aBit = 0;

	//write out initial binary pedigree file bytes, first 2 byte are magic numbers the third to indicate SNP major (subjects x SNPs)
	char buffer[3];	
	buffer[0] = 108;
	buffer[1] = 27;
	buffer[2] = 1;
	psBedFile.write(buffer, 3);

	if(interaction)
	{
		unsigned int genotype = 1;
		//loop thro' SNPs and add the pseudocontrol data
		for(unsigned int snpID = 1; snpID <= noOfSnps; ++snpID)
		{
			if(snpID == snpNo1 || snpID == snpNo2)
			{
				if(snpID == snpNo2) genotype++;
				bitCount = 9;
				//add SNP data for all subjects
				for(list<Subject *>::iterator sub = orderListOfSubjects.begin(); sub != orderListOfSubjects.end(); ++sub)
				{			
					addGenotypeDataToSubject(*sub, genotype);
				};			
			}
			else
			{
				bitCount = 9;
				//read in unwanted SNP data for all subjects
				for(list<Subject *>::iterator sub = orderListOfSubjects.begin(); sub != orderListOfSubjects.end(); ++sub)
				{
					doNotAddGenotypeDataToSubject();
				};
			};
			
		};

		//add corresponding data to all pseudocontrols for the two SNPs
		allTrioData.updatePsuedoControlJointSNPData(psBedFile, psBitCount, aBit);

		//store SNP data for info file
		if(doInfo)
		{
			allTrioData.storeSNPDataForInfoFile(1);
			allTrioData.storeSNPDataForInfoFile(2);
		};
	}
	else
	{
		//loop thro' SNPs and add the pseudocontrol data
		for(unsigned int snpID = 1; snpID <= noOfSnps; ++snpID)
		{

			bitCount = 9;
			//add SNP data for all subjects
			for(list<Subject *>::iterator sub = orderListOfSubjects.begin(); sub != orderListOfSubjects.end(); ++sub)
			{			
				addGenotypeDataToSubject(*sub);
			};
	
			//add corresponding data to all pseudocontrols for this SNP
			allTrioData.updatePsuedoControlSNPData(psBedFile, psBitCount, aBit);

			//start new byte as starting a new SNP
			writeLastByteBeforeNextSNP(psBedFile, psBitCount, aBit);

			//store SNP data for info file
			if(doInfo) allTrioData.storeSNPDataForInfoFile();
		};
	};

	psBedFile.close();

	if(doInfo) allTrioData.outputInfoFile(infoFilename, mapSnpAlleleIds, infoFileFatherGeno, infoFileMotherGeno, interaction, cepg, snpNo1, snpNo2);
	
	allTrioData.outputPseudoControlInfo();
};

