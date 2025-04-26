// "translate_FSud.cpp"

// Stop codons are assigned "B", "X", or "Z"
// call:  "input_file.tbl" "output_file_FSu.tbl" "output_file_FSd.tbl" 

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>

// Namespace using declarations
using std::string;
using std::cin;
using std::cout;
using std::ifstream;
using std::ofstream;
using std::endl;
using std::vector;
using std::find;
using std::getline;
using std::istringstream;
using std::to_string;
using std::replace;


// Smallest peptide kept in results
const int min_peptide_length = 6;
const int max_flank_size = 60; //nt

// Declare functions
void process_tbl_line(string str_label, string str_nt);
int open_check_ofstream(ofstream& the_stream, string the_label);
vector<string> convert_to_codons(string to_convert);
string translate_nt_to_aa(string to_translate);

// Declare ofstreams
ofstream output_aa_1p;
ofstream output_aa_1n;


int main(int argc, char **argv){
	
	// Build input & output file names
	string input_file = argv[1];
	string output_file_FSu = argv[2];
	string output_file_FSd = argv[3];
		
	open_check_ofstream(output_aa_1p, output_file_FSu);
	open_check_ofstream(output_aa_1n, output_file_FSd);	

	// Main input file & istream
	ifstream main_file;
	main_file.open(input_file);
	
	if ( !main_file.is_open() ){
		cout << "Error -- could not open file: " << input_file << endl;
	} else { cout << "Main file opened" << endl;}

	// Loop through each line, processing w/ major function
	string line = "";
	
	while (getline(main_file, line)){
   		istringstream iss(line);
	    string line_label;
		string line_sequence;
    
   		if( !(iss >> line_label >> line_sequence) ) { break; } // error

		auto line_seq_size = line_sequence.size();
		
		if(line_seq_size>=18){
			process_tbl_line(line_label, line_sequence);
		} 
//else{
//			output_excluded << line_label << "\t" << line_sequence << "\n";
//		}
	}
	

	// Close file streams
	main_file.close();
	output_aa_1p.close();
	output_aa_1n.close();

	return(0);
	
}
	

// Process label & nt sequence from tbl
void process_tbl_line(string str_label, string str_nt){
	
	// Loop through frames 
	for(int frame = 0; frame < 3; frame++){
	
		string frame_label = str_label + "__frame_" + to_string(frame);
		const string frame_nt = str_nt.substr(frame);
		const auto nt_size = frame_nt.size();
		
//		//Full read-through translation & write-out
//		string peptide = translate_nt_to_aa(frame_nt);			
//		string peptide_label = frame_label + "_noFS";
//		output_aa_0 << peptide_label << "\t" << peptide << "\n";
		
		string peptide_label = "";
		
		// Start position for indels
		int indel_start_position = 2;			
		
		// The +2< condition allows adding 1nt; the +4< allows subtracting
		while (indel_start_position+2 < nt_size){
			if(indel_start_position+4 < nt_size){
				// Build deletion strings -- skips nt at indel_start_position+1
				
				int head_start = 0;
			
				// Determine how much of the head to cut off
				if( (indel_start_position+1) > max_flank_size){
					head_start = (indel_start_position+1)-max_flank_size;
				}
			
				int tail_length = max_flank_size;
			
				// Determine how much of the tail to cut off
				if( (nt_size-(indel_start_position+2)) < max_flank_size){
					tail_length = nt_size-(indel_start_position+2);
				}
			
				string head = frame_nt.substr(head_start, (indel_start_position+1-head_start));
				string tail = frame_nt.substr(indel_start_position+2, tail_length);
				string deletion_string = head + tail;
				string peptide_1n = translate_nt_to_aa(deletion_string);
				auto peptide_length = peptide_1n.size();
						
				if( peptide_length >= min_peptide_length ){
				
					peptide_label = frame_label + "_deletion_" + to_string(indel_start_position+1);
						
					//	Write the translation to file
					output_aa_1n << peptide_label << "\t" << peptide_1n << "\n";
				
				}
			} // end if for the deletions

			// Build insertion strings -- add nt at indel_start_position
			
			int head_start = 0;
			
			// Determine how much of the head to cut off
			if( (indel_start_position+1) > max_flank_size){
				head_start = (indel_start_position+1)-max_flank_size;
			}
			
			int tail_length = max_flank_size+1;
			
			// Determine how much of the tail to cut off
			if( (nt_size-indel_start_position) < max_flank_size){
				tail_length = nt_size-indel_start_position;
			}
			string head = frame_nt.substr(head_start, (indel_start_position+1-head_start));
			string tail = frame_nt.substr(indel_start_position, tail_length);
			string insertion_string = head + tail;
			string peptide_1p = translate_nt_to_aa(insertion_string);
			auto peptide_length = peptide_1p.size();
				
			if (peptide_length >= min_peptide_length ){
					
				peptide_label = frame_label + "_insertion_" + to_string(indel_start_position+1);
					
				//	Write the translation to file
				output_aa_1p << peptide_label << "\t" << peptide_1p << "\n";
			}
			
			indel_start_position = indel_start_position+3;
			
		} // close while loop over indel positions
	} // close loop over frames
} // close function



// Open ofstream and check status
int open_check_ofstream(ofstream& the_stream, string the_label){

	the_stream.open(the_label);
	if( !the_stream.is_open() ){
		cout << "Error -- could not open file: " << the_label << endl;
	}
	return(0);
}



// Return vector of triplets from string
vector<string> convert_to_codons(string to_convert){
	vector<string> codon_vec;
	int start_position = 0;
	auto string_size = to_convert.size();
	
	int size_to_reserve = string_size/3;
	codon_vec.reserve(size_to_reserve);
	
	while( start_position + 2 < string_size ){
		string temp_codon = to_convert.substr(start_position, 3);
		codon_vec.push_back(temp_codon);
		start_position = start_position + 3;
	}
	
	return(codon_vec);
}
		
	
	
// Return amino acid string
string translate_nt_to_aa(string to_translate){

	vector<string> codon_vec = convert_to_codons(to_translate);
	
	string translated_aa = "";
	string current_aa = "";
	
	auto size_to_translate = codon_vec.size();
	int string_size_to_reserve = size_to_translate;
	translated_aa.reserve(string_size_to_reserve);
	
	for(auto vec_it = codon_vec.begin(); vec_it != codon_vec.end(); ++vec_it){
		
		string current_codon = *vec_it;
		
		if(current_codon.size()==3){
			string nt_0 = current_codon.substr(0,1);
			string nt_1 = current_codon.substr(1,1);
			string nt_2 = current_codon.substr(2,1);
			
			if(nt_0 == "G"){
				if(nt_1 == "G"){
					current_aa = "G";
				} else if (nt_1 == "C"){
					current_aa = "A";
				} else if (nt_1 == "T"){
					current_aa = "V";
				} else {//if (nt_1 == "A"){
					if( nt_2 == "A" || nt_2 == "G"){
						current_aa = "E";
					} else { //if (nt_2 == "C" || nt_2 == "T"){
						current_aa = "D";
					}
				}
			} else if (nt_0 == "A"){
				if(nt_1 == "C"){
					current_aa = "T";
				} else if (nt_1 == "G"){
					if( nt_2 == "G" || nt_2 == "A"){
						current_aa = "R";
					} else { //if (nt_2 == "C" || nt_2 == "T"){
						current_aa = "S";
					}
				} else if (nt_1 == "T"){
					if (nt_2 == "G"){
						current_aa = "M";
					} else { //if (nt_2 == "A" || nt_2 == "C" || nt_2 == "T"){
						current_aa = "I";
					}
				} else { //if (nt_1 == "A"){
					if (nt_2 == "G" || nt_2 == "A"){
						current_aa = "K";
					} else { //if (nt_2 == "C" || nt_2 == "T"){
						current_aa = "N";
					}
				}
			} else if (nt_0 == "T"){
				if(nt_1 == "C"){
					current_aa = "S";
				} else if(nt_1 == "A"){
					if(nt_2 == "T" || nt_2 == "C"){
						current_aa = "Y";
					} else if (nt_2 == "A"){
						current_aa = "B"; // also a stop codon
					} else { // if nt_2 == "G"
						current_aa = "Z"; //return(translated_aa);
					}
				} else if (nt_1=="G"){
					if(nt_2=="T" || nt_2=="C"){
						current_aa = "C";
					} else if (nt_2=="A"){
						current_aa = "X"; //return(translated_aa);
					} else { //if (nt_2=="G"){
						current_aa = "W";
					}
				} else {//if (nt_1 == "T"){
					if(nt_2 == "T" || nt_2 == "C"){
						current_aa = "F";
					} else { //(nt_2 == "A" || nt_2 == "G"){
						current_aa = "L";
					}
				}
			} else { //if (nt_0 == "C"){
				if(nt_1 == "T"){
					current_aa = "L";
				} else if (nt_1 == "C"){
					current_aa = "P";
				} else if (nt_1 == "G"){
					current_aa = "R";
				} else {//if (nt_1 == "A"){
					if(nt_2 == "T" || nt_2 == "C"){
						current_aa = "H";
					} else { //if (nt_2 == "A" || nt_2 == "G"){
						current_aa = "Q";
					}
				}
			}			
		} // end search through codon map										
			
		translated_aa = translated_aa + current_aa;
		current_aa = "";

	} // end loop through codons
	return(translated_aa);
}
