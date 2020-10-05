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


/*! \file main.cpp
    \brief This file reads in the initial input files and options.
    
    This file also outputs usage instructions and program details.
*/

#include <iostream>
#include <ostream>
#include <fstream>
#include <string>
#include <time.h>

using namespace std; // initiates the "std" or "standard" namespace
 
#include "main.h"
#include "Analysis.h"

bool outputToScreen = true; 
ofstream logFile;


//! Converts an integer to a string
string toString(unsigned int & i)
{
	ostringstream aStringStream;
	aStringStream << i;

	return aStringStream.str();
};

//! Returns a string of the run time
string getTime(const double & t)
{
	double time = t;
	unsigned int days = 0;
	unsigned int hours = 0;
	unsigned int minutes = 0;
	unsigned int seconds = 0;

	string ans = "";
	days = (unsigned int) (time / 86400); time -= days*86400;
	hours = (unsigned int) (time / 3600); time -= hours*3600;
	minutes = (unsigned int) (time / 60); time -= minutes*60;
	seconds = (unsigned int) time;

	if(days == 1) ans += "1 day";
	else if(days > 0) ans += toString(days) + " days";

	if(hours > 0)
	{
		if(days != 0)
		{
			if(minutes == 0 && seconds == 0) ans += " and ";
			else ans += ", ";
		};

		if(hours == 1) ans += "1 hour";
		else ans += toString(hours) + " hours";
	};

	if(minutes > 0)
	{
		if(ans != "")
		{
			if(seconds == 0) ans += " and ";
			else ans += ", ";
		};

		if(minutes == 1) ans += "1 minute";
		else ans += toString(minutes) + " minutes";
	};

	if(seconds > 0)
	{
		if(ans != "")
		{
			ans += " and ";			
		};

		if(seconds == 1) ans += "1 second";
		else ans += toString(seconds) + " seconds";
	};

	if(ans == "") ans = "less than one second";

	return ans;
};

//! Output program title to screen
void header()
{
	out("\nPseudoCons: pseudocontrols from pedigree data v1.12\n");
	out("---------------------------------------------------\n");
	out("Copyright 2013-2014 Richard Howey, GNU General Public License, v3\n");
	out("Institute of Genetic Medicine, Newcastle University\n\n");
};

//! Output program usage to screen
void usage()
{
		header();
	 	
		out("Usage:\n\t ./pseudocons [options] inputFile\n");
		out(" or ./pseudocons -pf parameterfile\n\n");

		out("Options:\n");
		
		out("  -i filename         -- input filename\n");		
		out("  -o filename         -- output filename\n");
		out("  -pc1                -- one pseudocontrol per trio\n");
		out("  -pc3                -- three pseudocontrols per trio\n");
		out("  -pc15               -- 15 pseudocontrols per trio (2 SNP interaction only)\n");
		out("  -snpnos snp1 snp2   -- SNP numbers of pair to create 15 pseudocontrols per trio\n");
		out("  -snpnames snp1 snp2 -- as above using SNP names\n");
		out("  -cepg               -- use CEPG giving 3, 7 and 31 pseudos for -pc1, 3 and 15\n");
		out("  -info info.dat      -- output info file\n");
		out("  -info-fa info.dat   -- output info file with father+child genotypes\n");
		out("  -info-ma info.dat   -- output info file with mother+child genotypes\n");
		out("  -info-fama info.dat -- output info file with father+mother+child genotypes\n");
		out("  -info-maxsnps n     -- max number of SNPs allowed in info file (default=20)\n");
		out("  -pro                -- proband filename\n");
		out("  -xtrio              -- allow extra trios\n");
		out("  -log                -- log filename\n");		
		out("  -so                 -- suppress output to screen\n\n");
		out("Default options:\n");
		out("  -pc1\n\n");		
		
};


//! Get an option value either from the parameter file or the command line
void getOptionValue(unsigned int & anUnInt, bool & useParaFile, int & argcount, int & argc, char * argv[], ifstream & readParaFile)
{	
	if(useParaFile)
	{
		if(readParaFile.eof()) return;
		readParaFile >> anUnInt;		
	}
	else
	{
		argcount++; if(argcount >= argc) return;				
		anUnInt = atoi(argv[argcount]);		
	};
};

//! Get an option value either from the parameter file or the command line
void getOptionValue(double & aDouble, bool & useParaFile, int & argcount, int & argc, char * argv[], ifstream & readParaFile)
{
	if(useParaFile)
	{
		if(readParaFile.eof()) return;		
		readParaFile >> aDouble;		
	}
	else
	{
		argcount++; if(argcount >= argc) return;			
		aDouble = atof(argv[argcount]);		
	};
};

//! Get an option value either from the parameter file or the command line
void getOptionValue(string & aString, bool & useParaFile, int & argcount, int & argc, char * argv[], ifstream & readParaFile)
{
	if(useParaFile)
	{
		if(readParaFile.eof()) return;		
		readParaFile >> aString;
	}
	else
	{
		argcount++; if(argcount >= argc) return;		
		aString = argv[argcount];
	};
};

//! Gets the log filename from the results file name
string getDefaultLogFileName(string & outFileName)
{
	unsigned int filenameLength = outFileName.length();
	string logFileName;

	//find extension
	unsigned int a = filenameLength - 1;
	while(a > 0)
	{
		if(outFileName.substr(a, 1) == ".") break;
		a--;
	};

	if(a > 0) logFileName = outFileName.substr(0, a) + ".log";
	else  logFileName = outFileName + ".log";

	return logFileName;
};

//! The start of the program
int main(int argc, char * argv[])
{
	time_t start,end;
	double dif;
	time(&start);

	int argcount = 1;
	string option;

	string inputFilename = "";
	string probandFilename = "";
	string outputFilename = "";	
	string logFilename = "";
	string infoFilename = "";
	bool infoFileFatherGeno = false;
	bool infoFileMotherGeno = false;
	unsigned int infoFileMaxSNPs = 20;
	bool extraAffectedTrios = false;
	unsigned int noPseudoCons = 1;
	unsigned int snpNo1 = 0, snpNo2 = 0;
	string snpName1 = "", snpName2 = "";
	bool cepg = false;

	outputToScreen = true;	

	bool useParaFile = false;
	string paraFilename = "";
	ifstream readParaFile;
	
	if(argcount < argc) option = argv[argcount];

	//deal with parameter file
	if(option == "-pf")
	{
		argcount++; 
		if(argcount < argc) paraFilename = argv[argcount];

		//open parameter file		
		readParaFile.open(paraFilename.c_str());
		if(!readParaFile.is_open())
		{
			header();
			outErr("Cannot read parameter file: "); outErr(paraFilename); outErr("!\n");
			exit(1);
		};

		argcount++; 
		useParaFile = true;
	};

	
	//set given options
	while((!useParaFile && argcount < argc && argv[argcount][0] == '-') || (useParaFile && !readParaFile.eof()))
	{

		if(useParaFile)
		{
			//find the start of the next command
			do{
				readParaFile >> option;
				if(option.length() >= 2 && option.substr(0,1) == "-") break;				
			}while(!readParaFile.eof());
		}
		else
		{
			option = argv[argcount];
		};

		if(useParaFile && readParaFile.eof()) break;

		if(option ==  "-i") getOptionValue(inputFilename, useParaFile, argcount, argc, argv, readParaFile);		
		else if(option ==  "-o") getOptionValue(outputFilename, useParaFile, argcount, argc, argv, readParaFile);
		else if(option ==  "-log") getOptionValue(logFilename, useParaFile, argcount, argc, argv, readParaFile);
		else if(option ==  "-pro") getOptionValue(probandFilename, useParaFile, argcount, argc, argv, readParaFile);
		else if(option ==  "-snpnos")
		{
				getOptionValue(snpNo1, useParaFile, argcount, argc, argv, readParaFile);
				getOptionValue(snpNo2, useParaFile, argcount, argc, argv, readParaFile);
		}
		else if(option ==  "-snpnames")
		{
				getOptionValue(snpName1, useParaFile, argcount, argc, argv, readParaFile);
				getOptionValue(snpName2, useParaFile, argcount, argc, argv, readParaFile);
		}
		else if(option == "-so") outputToScreen = false;
		else if(option == "-pc1") noPseudoCons = 1;	
		else if(option == "-pc3") noPseudoCons = 3;	
		else if(option == "-pc15") noPseudoCons = 15;
		else if(option == "-cepg") cepg = true;
		else if(option == "-xtrio") extraAffectedTrios = true;
		else if(option == "-info")	getOptionValue(infoFilename, useParaFile, argcount, argc, argv, readParaFile);	
		else if(option == "-info-fa")
		{
				getOptionValue(infoFilename, useParaFile, argcount, argc, argv, readParaFile);		
				infoFileFatherGeno = true;
		}
		else if(option == "-info-ma")
		{
				getOptionValue(infoFilename, useParaFile, argcount, argc, argv, readParaFile);		
				infoFileMotherGeno = true;
		}
		else if(option == "-info-fama")
		{
				getOptionValue(infoFilename, useParaFile, argcount, argc, argv, readParaFile);		
				infoFileFatherGeno = true;
				infoFileMotherGeno = true;
		}
		else if(option == "-info-maxsnps")	getOptionValue(infoFileMaxSNPs, useParaFile, argcount, argc, argv, readParaFile);
		else
		{			
			header();
			if(useParaFile) {outErr("Unrecognised option: "); outErr(option); outErr("\n\n");}
			else {outErr("Unrecognised command line switch: "); outErr(option); outErr("\n\n");};			
    		exit(1);
		};

		if(!useParaFile) argcount++;
	};

	if(argc == 1)
	{
		usage();
		exit(1);
	};
	
	if(logFilename == "") logFilename = getDefaultLogFileName(outputFilename);

	logFile.open(logFilename.c_str());

	if(inputFilename == "")
	{
		outErr("\nInput file not set!\n\n");
		usage();
		exit(1);
	};	

	if(outputFilename == "")
	{
		outErr("\nOutput file not set!\n\n");
		usage();
		exit(1);
	};	

	//output options to screen	
	header();
	out("Parameters:\n");
	out("Input file: "); out(inputFilename); out("\n");		
	out("Output file: "); out(outputFilename); out("\n");
	out("Log file: "); out(logFilename); out("\n");
	if(snpNo2 != 0) {out("Interaction using SNP numbers "); out(snpNo1); out(" and "); out(snpNo2); out("\n");};
	if(snpName2 != "") {out("Interaction using SNP names "); out(snpName1); out(" and "); out(snpName2); out("\n");};
	out("Number of pseudocontrols per trio: "); out(noPseudoCons); 
	
	if(cepg)
	{
		if(noPseudoCons == 1) out(" + 2 = 3\n");
		else if(noPseudoCons == 3) out(" + 4 = 7\n");
		else if(noPseudoCons == 15) out(" + 16 = 31\n");

		out("Using Conditional on Exchangeable Parental Genotypes (CEPG)\n");
	}
	else out("\n");

	if(probandFilename != "") {out("Proband file: "); out(probandFilename); out("\n");};
	if(extraAffectedTrios) {out("Allowing extra trios\n");}; 

	if(noPseudoCons == 15 && (snpName2 == "" && snpNo2 == 0))
	{
		outErr("\nThe 15 pseudocontrol option is only available for 2 SNP interaction pseudocontrols!\n");
		exit(1);
	};

	if(noPseudoCons == 3 && (snpName2 != "" || snpNo2 != 0))
	{
		outErr("\nThe 3 pseudocontrol option is not available for 2 SNP interaction pseudocontrols!\n");
		exit(1);
	};

	if(noPseudoCons == 1 && (snpName2 != "" || snpNo2 != 0))
	{
		outErr("\nThe 1 pseudocontrol option is not available for 2 SNP interaction pseudocontrols!\n");
		exit(1);
	};

	if(snpName2 != "" && snpNo2 != 0)
	{
		outErr("\nCannot set SNPs for interaction using both numbers and names!\n");
		exit(1);
	};
	
	//run checks for infofile
	if(infoFilename != "")
	{
		out("\nInfo file: "); out(infoFilename); out("\n");
		out("          Column 1: family ID\n");
		out("          Column 2: individual ID\n");
		out("          Column 3: corresponding case ID\n");
		out("          Column 4: pseudocontrol set ID\n");
		//add SNP column info later
	};

	//create analysis option and run analysis
	Analysis anAnalysis(inputFilename, outputFilename, probandFilename, extraAffectedTrios, noPseudoCons, snpNo1, snpNo2, snpName1, snpName2, cepg, infoFilename, infoFileFatherGeno, infoFileMotherGeno, infoFileMaxSNPs);

	anAnalysis.runAnalysis();

	time(&end);
	dif = difftime(end, start);
	out("Run time: "); out(getTime(dif)); out("\n\n");

	logFile.close();
};

