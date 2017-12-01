#ifndef SOAP_SNP_HH_
#define SOAP_SNP_HH_
#ifndef TIME_BREAK //[cuiyb]++
#define TIME_BREAK
#include <iostream>
#include <stdio.h> //[cuiyb]++ printf
#include <fstream>
#include <sstream>
#include <string.h>
#include <map>
#include <vector>
#include <cmath>
#include <math.h>
#include <iomanip>
#include <cassert>
#include <stdlib.h>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <sys/time.h> //[cuiyb]++ timeval, gettimeofday()
typedef unsigned long long ubit64_t;
typedef unsigned int ubit32_t;
typedef double rate_t;
typedef unsigned char small_int;
using namespace std;
typedef std::string Chr_name;
__attribute__((target(mic))) const size_t capacity = sizeof(ubit64_t)*8/4;
const char abbv[17]={'A','M','W','R','M','C','Y','S','W','Y','T','K','R','S','K','G','N'};
const ubit64_t glf_base_code[8]={1,2,8,4,15,15,15,15}; // A C T G
const ubit64_t glf_type_code[10]={0,5,15,10,1,3,2,7,6,11};// AA,CC,GG,TT,AC,AG,AT,CG,CT,GT

__attribute__((target(mic))) const int max_depth = 255; // max depth of sites
const ubit64_t read_len = 100; //[cuiyb] default read length
const ubit32_t mic_num = 3;
const ubit32_t mic_size = 336000; // size of window, processed by one MIC, 1120=224*5
const ubit32_t cpu_size = 0;
const int global_win_size = mic_size*mic_num + cpu_size;
const ubit64_t win_size = mic_size*mic_num + cpu_size;

//global variables for offload region, define them in normal_dis.cc
__attribute__((target(mic))) extern double * new_p_matrix; // a_adjust*coord*base*array*length_of(type_likely), 32M size
__attribute__((target(mic))) extern rate_t * mat_p_prior; // global for mat->p_prior
__attribute__((target(mic))) extern ubit64_t* bin_seq; // Sequence in binary format, for new_get_bin_base()

// Some global variables
class Files {
public:
	ifstream soap_result, ref_seq, dbsnp, region;
	ofstream consensus, baseinfo, o_region;
	fstream matrix_file;
	Files(){
		soap_result.close();
		ref_seq.close();
		dbsnp.close();
		consensus.close();
		baseinfo.close();
		matrix_file.close();
		region.close();
		o_region.close();
	};
};

class Parameter {
public:
	char q_min; // The char stands for 0 in fastq
	char q_max; // max quality score
	small_int read_length; // max read length
	bool is_monoploid; // Is it an monoploid? chrX,Y,M in man.
	bool is_snp_only;  // Only output possible SNP sites?
	bool refine_mode; // Refine prior probability using dbSNP
	bool rank_sum_mode; // Use rank sum test to refine HET quality
	bool binom_mode; // Use binomial test to refine HET quality
	bool transition_dominant; // Consider transition/transversion ratio?
	int glf_format; // Generate Output in GLF format: New File Definitions! Since May 1, 2009
	bool region_only; // Only report consensus in specified region
	std::string glf_header; // Header of GLF format
	rate_t althom_novel_r, het_novel_r; // Expected novel prior
	rate_t althom_val_r, het_val_r; // Expected Validated dbSNP prior
	rate_t althom_unval_r, het_unval_r; // Expected Unvalidated dbSNP prior
	rate_t global_dependency, pcr_dependency; // Error dependencies, 1 is NO dependency
	bool is_matrix_out; //[cuiyb]++ output matrix
	Chr_name chr_name; //[cuiyb]++
	ubit32_t ref_len; //[cuiyb] from Chr_info, length of reference sequence
	int cpu_thread_num; //[cuiyb]++
// Default onstruction
	Parameter(){
		q_min = 64;
		q_max = 64+40;
		read_length = 45;
		is_monoploid = is_snp_only = refine_mode = rank_sum_mode = binom_mode = transition_dominant = region_only = is_matrix_out =false;
		glf_format = 0;
		glf_header = "";
		althom_novel_r=0.0005, het_novel_r=0.0010;
		althom_val_r=0.05, het_val_r=0.10;
		althom_unval_r=0.01, het_unval_r=0.02;
		global_dependency= log10(0.9), pcr_dependency= log10(0.5); // In Log10 Scale
		chr_name = ""; //[cuiyb]++
		ref_len = 0; //[cuiyb]++
		cpu_thread_num = 1; //[cuiyb]++
	};
};

class Soap_format {
	// Soap alignment result
	std::string read_id, read, qual, chr_name;
	int hit, read_len, position, mismatch;
	char ab, strand;
public:
	Soap_format(){;};
	friend std::istringstream & operator>>(std::istringstream & alignment, Soap_format & soap) {
		alignment>>soap.read_id>>soap.read>>soap.qual>>soap.hit>>soap.ab>>soap.read_len>>soap.strand>>soap.chr_name>>soap.position>>soap.mismatch;
		//cerr<<soap<<endl;
		//exit(1);
		if(soap.mismatch>200) {
			int indel_pos,indel_len;
			string temp("");
			alignment>>indel_pos;
			indel_len = soap.mismatch-200;
			for(int i=0; i!=indel_len; i++) {
				temp = temp+'N';
			}
			soap.read = soap.read.substr(0,indel_pos)+temp+soap.read.substr(indel_pos,soap.read_len-indel_pos);
			soap.qual = soap.qual.substr(0,indel_pos)+temp+soap.qual.substr(indel_pos,soap.read_len-indel_pos);
			//cerr<<soap<<endl;
		}
		else if (soap.mismatch>100) {
			int indel_pos,indel_len;
			alignment>>indel_pos;
			indel_len = soap.mismatch-100;
			soap.read = soap.read.substr(0,indel_pos) + soap.read.substr(indel_pos+indel_len, soap.read_len-indel_pos-indel_len);
			soap.qual = soap.qual.substr(0,indel_pos) + soap.qual.substr(indel_pos+indel_len, soap.read_len-indel_pos-indel_len);
			//cerr<<soap<<endl;
		}
		//cerr<<soap.position<<'\t';
		soap.position -= 1;
		//cerr<<soap.position<<endl;
		return alignment;
	}
	friend std::ostream & operator<<(std::ostream & o, Soap_format & soap) {
		o<<soap.read_id<<'\t'<<soap.read<<'\t'<<soap.qual<<'\t'<<soap.hit<<'\t'<<soap.ab<<'\t'<<soap.read_len<<'\t'<<soap.strand<<'\t'<<soap.chr_name<<'\t'<<soap.position<<'\t'<<soap.mismatch;
		return o;
	}
	char get_base(std::string::size_type coord) {
		return read[coord];
	}
	char get_qual(std::string::size_type coord) {
		return qual[coord];
	}
	bool is_fwd(){
		return (strand=='+');
	}
	int get_read_len(){
		return read_len;
	}
	int get_pos(){
		return position;
	}
	std::string get_chr_name(){
		return chr_name;
	}
	int get_hit(){
		return hit;
	}
	bool is_unique(){
		return (hit==1);
	}
	bool is_N(int coord) {
		return (read[coord] == 'N');
	}
};

// dbSNP information
class Snp_info {
	bool validated;
	bool hapmap_site;
	bool indel_site;
	rate_t * freq; // Arrary of 4 elements recording frequency of ACTG
public:
	Snp_info(){
		validated=hapmap_site=indel_site=false;
		freq = new rate_t [4];
		memset(freq,0,sizeof(rate_t)*4);
	}
	Snp_info(const Snp_info & other) {
		validated = other.validated;
		hapmap_site = other.hapmap_site;
		indel_site = other.indel_site;
		freq = new rate_t [4];
		memcpy(freq, other.freq, sizeof(rate_t)*4);
	}
	~Snp_info(){
		delete [] freq;
	}
	friend std::istringstream& operator>>(std::istringstream & s, Snp_info & snp_form) {
		s>>snp_form.hapmap_site>>snp_form.validated>>snp_form.indel_site>>snp_form.freq[0]>>snp_form.freq[1]>>snp_form.freq[2]>>snp_form.freq[3];
		return s;
	}
	Snp_info & operator=(Snp_info& other) {
		this->validated = other.validated;
		this->hapmap_site = other.hapmap_site;
		this->indel_site = other.indel_site;
		this->freq = new rate_t [4];
		memcpy(this->freq, other.freq, sizeof(rate_t)*4);
		return *this;

	}
	bool is_validated(){
		return validated;
	}
	bool is_hapmap(){
		return hapmap_site;
	}
	bool is_indel(){
		return indel_site;
	}
	rate_t get_freq(char bin_base_2bit) {
		return freq[bin_base_2bit];
	}
};

// function re-written for offload
// fro class Call_win
int new_soap2cns(std::ifstream & alignment, std::ofstream & consensus, std::ofstream & baseinfo, Parameter * para);
int new_initialize(ubit64_t start);
int new_recycle();
__attribute__((target(mic))) void new_sort_base(ubit32_t array[],int left, int right);
__attribute__((target(mic))) ubit64_t new_get_bin_base(std::string::size_type pos);

void seek_pointer(std::ifstream & is);
int overlap_length(int coord, int win_size);
void add_overlap(int * array, int * length);
void shareBases (long long int start, int win_size);

void get_start_time(void);
void get_end_time(int &time_item);
void print_time_cost(int time_item, string process);

// from Prob_matix class
int p_matrix_gen(std::ifstream & alignment, Parameter * para, std::fstream &mat_out);
int new_p_matrix_gen(rate_t * p_matrix);
int p_matrix_read(std::fstream &mat_in, Parameter * para);
int p_matrix_write(std::fstream &mat_out, rate_t * p_matrix, Parameter * para);
int p_prior_gen(Parameter * para);

// from Genome class
int genome_read(std::ifstream &fasta, Parameter * para);
int genome_binarize(std::string & seq, Parameter * para);

#endif //TIME_BREAK in the head of file
#endif /*SOAP_SNP_HH_*/
