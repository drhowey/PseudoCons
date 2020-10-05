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


/*! \file main.h
    \brief This file defines a global variable to suppress output to screen. 
    
*/

#ifndef __MAIN
#define __MAIN

#include <iostream>
#include <sstream>

extern bool outputToScreen; 
extern ofstream logFile; 

template<typename T>
//! Outputs message to screen and log file
void out(const T & text)
{	
	if(outputToScreen) cout << text;
	logFile << text;
};

template<typename T>
//! Outputs error message to screen and log file
void outErr(const T & text)
{
	cerr << text;
	logFile << text;
};

//! Converts an integer to a string
string toString(unsigned int & i);

//! Returns a string of the run time
string getTime(const double & t);

#endif
