// build_IL_tbl.cpp

// DRC 2025-04-23 (from 2023 PhD dissertation codebase)

// Based on C++11 standard. I compile on the NIH HPC using "g++ --std=c++11 -o [utility] [utility.cpp]

// call: ./build_IL_tbl input_file output_file

/*
This program takes as input a sequence file, replaces each instance of 'I' with 'L', and writes each modified sequence and the original sequence (tab-separated) to the output file.

The input file looks like this:
	ABCDEIL
	ABCDELL
	ABCIEDF

The output files looks like this:
	ABCDELL	ABCDEIL
	ABCDELL	ABCDELL
	ABCLEDF ABCIEDF

The output list is used to match the original sequences to results from a mass spectrometry analysis which returns 'L' for instances of 'L' or 'I' (LC-MS/MS does not distinguish between leucine (L) and isoleucine (I) because they have the same mass; MS programs typically return 'L' for this ambiguous amino acid).

The output file has no column headers, just two tab-separated columns. The output file extension is taken from the CL "output_file" parameter -- I typically use ".tbl" for headerless files with tab-separated columns.
*/

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <algorithm>

// Namespace 'using' declarations
using std::string;
using std::ifstream;
using std::ofstream;
using std::getline;
using std::find;
using std::cout;
using std::istringstream;
using std::endl;
using std::vector;
using std::to_string;
using std::unordered_map;
using std::array;
using std::stoi;
using std::replace;

// Declare functions
int open_check_ofstream(ofstream& the_stream, string the_label);

// Declare ofstreams
ofstream out_stream_IL;

int main( int argc, char **argv){
	
	// get command line arguments
	string input_file = argv[1];
	string output_file = argv[2];
	
	ifstream in_stream(input_file);
	if( !in_stream.is_open() ){
		cout << "Error -- could not open db file" << endl;
	}

	open_check_ofstream(out_stream_IL, output_file);

	string line = "";
	
	// Major loop
	
	while(getline(in_stream, line)){
		
		istringstream iss(line);
		string line_sequence;

		if( !(iss >> line_sequence) ) { break; } // add error call

		string line_copy = line_sequence;
		replace(line_copy.begin(), line_copy.end(), 'I', 'L');

		out_stream_IL << line_copy << "\t" << line_sequence << "\n";
		
	}
	in_stream.close();
	out_stream_IL.close();
	return(0);
}


// Open ofstream and check status
int open_check_ofstream(ofstream& the_stream, string the_label){

	the_stream.open(the_label);
	if( !the_stream.is_open() ){
		cout << "Error -- could not open file: " << the_label << endl;
	}
	
	return(0);
}
