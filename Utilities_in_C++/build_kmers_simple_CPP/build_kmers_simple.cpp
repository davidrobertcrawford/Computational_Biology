// build_kmers_simple.cpp

// DRC 2025-04-23 (from 2023 PhD dissertation codebase)

// Based on C++11 standard. I compile on the NIH HPC using "g++ --std=c++11 -o [utility] [utility.cpp]

// call: build_kmers_simple input_file output_file kmer_length

/*
This program takes as input a two-column (tab-separated) file and returns as output a kmer list.

A FASTA file consists of paired header (starts with ">") and sequence lines:
	>header_id
	NUCLEOTIDEORAMINOACIDSEQUENCE
	>header_id2
	ANOTHERSEQUENCE

The input file is a FASTA file that has been converted to a two-column (tab-separated) file:
	>header_id	NUCLEOTIDEORAMINOACIDSEQUENCE
	>header_id2	ANOTHERSEQUENCE

The output kmer list (for kmer_length==6) has the format:
	HERSEQ
	THERSE
	RSEQUE
	...
	EOTIDE
	EORAMI
	UCLEOT

The command line output for the example input file (above) is:
	Kmer size: 6
	Number of kmers: 31

FASTA files are typically used for nucleotide or amino acid sequences. The kmer list is simply a list of unique strings of length k (provided via the kmer_length argument) found in one or more of the FASTA sequences. The data in the header lines is ignored.

The output file has no column header, just an alphabetized list of kmers. The output file extension is taken from the CL "output_file" parameter -- I typically use ".tbl" for headerless files with tab-separated columns.
*/

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <unordered_map>

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

// Declare functions
int open_check_ofstream(ofstream& the_stream, string the_label);

// Declare ofstreams
ofstream out_stream_solo;


int main( int argc, char **argv){
	
	// get command line arguments
	string input_file = argv[1];
	string output_file = argv[2];
	string kmer_length_char = argv[3];
	
	const int kmer_length = stoi(kmer_length_char);
	
	ifstream in_stream(input_file);
	if( !in_stream.is_open() ){
		cout << "Error -- could not open db file" << endl;
	}

	open_check_ofstream(out_stream_solo, output_file);

	string line = "";
	
	
	// Build master array
	
	unordered_map<string, int> kmer_dict;
	
	// Major loop
	
	while(getline(in_stream, line)){
		
		istringstream iss(line);
		string line_label;
		string line_sequence;

		if( !(iss >> line_label >> line_sequence) ) { break; } // add error call

		auto line_size = line_sequence.size();
		
		int start = 0;
		int stop = line_size - kmer_length + 1;
		
		while(start < stop){
			string current_kmer = line_sequence.substr(start, kmer_length);
			++kmer_dict[current_kmer];
			start++;
		}
	}
	in_stream.close();

	// Save to disk
	for(const auto &w : kmer_dict){
		out_stream_solo << w.first  << "\n";
	}
	auto map_size = kmer_dict.size();
	cout << "Kmer size: " << kmer_length_char << "\n";
	cout << "Number of kmers: " << to_string(map_size) << endl;
	out_stream_solo.close();
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
