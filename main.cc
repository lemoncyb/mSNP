#include "soap_snp.h"
#include <getopt.h>
#include <mpi.h> // [cuiyb]++
#include <sstream>


using namespace std;

int usage() {
	cerr<<"mSNP"<<endl;
	cerr<<"Compulsory Parameters:"<<endl;
	cerr<<"-i <FILE> Input SORTED Soap Result"<<endl;
	cerr<<"-d <FILE> Reference Sequence in fasta format"<<endl;
	cerr<<"-o <FILE> Output consensus file"<<endl;
	cerr<<"Optional Parameters:(Default in [])"<<endl;
	cerr<<"-z <Char> ASCII chracter standing for quality==0 [@]"<<endl;
	cerr<<"-g <Double> Global Error Dependency Coefficient, 0.0(complete dependent)~1.0(complete independent)[0.9]"<<endl;
	cerr<<"-p <Double> PCR Error Dependency Coefficient, 0.0(complete dependent)~1.0(complete independent)[0.5]"<<endl;
	cerr<<"-r <Double> novel altHOM prior probability [0.0005]"<<endl;
	cerr<<"-e <Double> novel HET prior probability [0.0010]"<<endl;
	cerr<<"-t set transition/transversion ratio to 2:1 in prior probability"<<endl;
	cerr<<"-s <FILE> Pre-formated dbSNP information"<<endl;
	cerr<<"-2 specify this option will REFINE SNPs using dbSNPs information [Off]"<<endl;
	cerr<<"-a <Double> Validated HET prior, if no allele frequency known [0.1]"<<endl;
	cerr<<"-b <Double> Validated altHOM prior, if no allele frequency known[0.05]"<<endl;
	cerr<<"-j <Double> Unvalidated HET prior, if no allele frequency known [0.02]"<<endl;
	cerr<<"-k <Double> Unvalidated altHOM rate, if no allele frequency known[0.01]"<<endl;
	cerr<<"-u Enable rank sum test to give HET further penalty for better accuracy. [Off]"<<endl;
	//cerr<<"-n Enable binomial probability calculation to give HET for better accuracy. [Off]"<<endl;
	cerr<<"-m Enable monoploid calling mode, this will ensure all consensus as HOM and you probably should SPECIFY higher altHOM rate. [Off]"<<endl;
	cerr<<"-q Only output potential SNPs. Useful in Text output mode. [Off]"<<endl;
	cerr<<"-M <FILE> Output the quality calibration matrix; the matrix can be reused with -I if you rerun the program"<<endl;
	cerr<<"-I <FILE> Input previous quality calibration matrix. It cannot be used simutaneously with -M"<<endl;
	cerr<<"-L <short> maximum length of read [45]"<<endl;
	cerr<<"-Q <short> maximum FASTQ quality score [40]"<<endl;
	cerr<<"-F <int> Output format. 0: Text; 1: GLFv2; 2: GPFv2.[0]"<<endl;
	cerr<<"-E <String> Extra headers EXCEPT CHROMOSOME FIELD specified in GLFv2 output. Format is \"TypeName1:DataName1:TypeName2:DataName2\"[""]"<<endl;
	cerr<<"-T <int> number of threads on CPU"<<endl;
	//cerr<<"-T <FILE> Only call consensus on regions specified in FILE. Format: ChrName\\tStart\\tEnd."<<endl;
	//cerr<<"-S <FILE> Output summary of consensus"<<endl;
	cerr<<"-h Display this help"<<endl;
	exit(1);
	return 0;
}

int readme() {
	return usage();
}

int main ( int argc, char * argv[]) {
	// [cuiyb]++
	MPI_Init(&argc, &argv);
	
	// [cuiyb]++
	time_t start_main, end_main;
	struct timeval tv;
	gettimeofday(&tv, NULL);
	start_main = tv.tv_sec*1000 + tv.tv_usec/1000;
	
	// This part is the default values of all parameters
	Parameter * para = new Parameter;
	std::string alignment_name, consensus_name("");
	bool is_matrix_in = false; // Generate the matrix or just read it?
	int c;
	Files files;
	
	// Parse options
	while((c=getopt(argc,argv,"i:d:o:z:g:p:r:e:ts:2a:b:j:k:unmqM:I:L:Q:S:F:E:T:h")) != -1) {
		switch(c) {
			case 'i':
			{
				// Soap Alignment Result
				files.soap_result.clear();
				files.soap_result.open(optarg);
				if( ! files.soap_result) {
					cerr<<"No such file or directory:"<<optarg<<endl;
					//exit(1);
				}
				alignment_name = optarg;
				break;
			}
			case 'd':
			{
				// The reference genome in fasta format
				files.ref_seq.clear();
				files.ref_seq.open(optarg);
				if( ! files.ref_seq) {
					cerr<<"No such file or directory:"<<optarg<<endl;
					//exit(1);
				}
				files.ref_seq.clear();
				break;
			}
			case 'o':
			{
				int rank;
				MPI_Comm_rank(MPI_COMM_WORLD, &rank);
				stringstream ss;
				string s;
				ss << rank;
				ss >> s;
				consensus_name = optarg+s;
				break;
			}
			case 'z':
			{
				// The char stands for quality==0 in fastq format
				para->q_min = optarg[0];
				if(para->q_min == 33) {
					clog<<"Standard Fastq System Set"<<endl;
				}
				else if(para->q_min == 64) {
					clog<<"Illumina Fastq System Set"<<endl;
				}
				else {
					clog<<"Other types of Fastq files?? Are you sure?"<<endl;
				}
				para->q_max = para->q_min + 40;
				break;
			}
			case 'g':
			{
				para->global_dependency= log10(atof(optarg));
				break;
			}
			case 'p':
			{
				para->pcr_dependency= log10(atof(optarg));
				break;
			}
			case 'r':
			{
				para->althom_novel_r = atof(optarg);
				break;
			}
			case 'e':
			{
				para->het_novel_r=atof(optarg);
				break;
			}
			case 't':
			{
				para->transition_dominant=true;
				break;
			}
			case 's':
			{
				// Optional: A pre-formated dbSNP table
				files.dbsnp.clear();
				files.dbsnp.open(optarg);
				if( ! files.ref_seq) {
					cerr<<"No such file or directory:"<<optarg<<endl;
					//exit(1);
				}
				files.dbsnp.clear();
				break;
			}
			case '2':
			{
				// Refine prior probability based on dbSNP information
				para->refine_mode = true;
				break;
			}
			case 'a':
			{
				para->althom_val_r=atof(optarg);
				break;
			}
			case 'b':
			{
				para->het_val_r=atof(optarg);
				break;
			}
			case 'j':
			{
				para->althom_unval_r=atof(optarg);
				break;
			}
			case 'k':
			{
				para->het_unval_r=atof(optarg);
				break;
			}
			case 'u':
			{
				para->rank_sum_mode = true;
				break;
			}
			case 'n':
			{
				para->binom_mode = true;
				break;
			}
		 	case 'm':
			{
				para->is_monoploid=1;
				break;
			}
			case 'q':
			{
				para->is_snp_only=1;
				break;
			}
			case 'M':
			{
				files.matrix_file.close(); files.matrix_file.clear();
				// Output the calibration matrix
				files.matrix_file.open(optarg, fstream::out);
				if( ! files.matrix_file) {
					cerr<<"Cannot creat file :"<<optarg<<endl;
					//exit(1);
				}
				files.matrix_file.clear();
				para->is_matrix_out = true;
				break;
			}
			case 'I':
			{
				files.matrix_file.close(); files.matrix_file.clear();
				// Input the calibration matrix
				files.matrix_file.open(optarg, fstream::in);
				if( ! files.matrix_file) {
					cerr<<"No such file or directory:"<<optarg<<endl;
					//exit(1);
				}
				files.matrix_file.clear();
				is_matrix_in = true;
				break;
			}
			case 'S':
			{
				//files.summary.open(optarg);
				//// Output the summary of consensus
				//if( ! files.summary ) {
				//	cerr<<"No such file or directory: "<<optarg<<endl;
				//	exit(1);
				//}
				break;
			}
			case 'L':
			{
				para->read_length = atoi(optarg);
				break;
			}
			case 'Q':
			{
				para->q_max = optarg[0];
				if(para->q_max < para->q_min) {
					cerr<< "FASTQ quality character error: Q_MAX > Q_MIN" <<endl;
				}
				break;
			}
			case 'F': {
				para->glf_format = atoi(optarg);
				break;
			}
			case 'E': {
				para->glf_header = optarg;
				break;
			}
			case 'T': {
			/*	files.region.clear();
				files.region.open(optarg);
				files.region.clear();
				para->region_only = true;
			*/
				para->cpu_thread_num = atoi(optarg);
				break;
			}
			case 'h':readme();break;
			case '?':usage();break;
			default: cerr<<"Unknown error in command line parameters"<<endl;
		}
	}
	// [cuiyb]++ Parse options done
	
	if( consensus_name=="" || !files.ref_seq || !files.soap_result ) {
		// These are compulsory parameters
		usage();
	}
	
	// [cuiyb]++ Open input file
	if ( ! para->glf_format ) {
		// Normal SOAPsnp tab-delimited text format
		files.consensus.clear();
		files.consensus.open(consensus_name.c_str());
		if( ! files.consensus ) {
			cerr<<"Cannot creat file:" <<consensus_name<<endl;
			//exit(1);
		}
		files.consensus.clear();
	}
	else {
		// SOAPsnp-defined GLF and baseinfo format
		files.consensus.clear();
		files.consensus.open(consensus_name.c_str(), ios::binary);
		if(!files.consensus) {
			cerr<<"Cannot creat file:" <<consensus_name<<endl;
			//exit(1);
		}
		files.consensus.clear();

		files.baseinfo.clear();
		string baseinfo_name = consensus_name + ".baseinfo";
		files.baseinfo.open(baseinfo_name.c_str());
		if(!files.baseinfo) {
			cerr<<"Cannot creat file:" <<baseinfo_name<<endl;
			//exit(1);
		}
		files.baseinfo.clear();

		files.o_region.clear();
		string o_region_name = consensus_name + ".index";
		files.o_region.open(o_region_name.c_str());
		if(!files.o_region) {
			cerr<<"Cannot creat file:" <<o_region_name<<endl;
			//exit(1);
		}
		files.o_region.clear();
	}

	genome_read(files.ref_seq, para);
	files.ref_seq.close();
	files.dbsnp.close();
	clog<<"Reading Chromosome and dbSNP information Done."<<endl;

	if( ! is_matrix_in) {
		//Read the soap result and give the calibration matrix
		p_matrix_gen(files.soap_result, para, files.matrix_file);
	}
	else {
		p_matrix_read(files.matrix_file, para);
	}
	files.matrix_file.close();
	clog<<"Correction Matrix Done!"<<endl;

	p_prior_gen(para); //[cuiyb] this function costs some time
	
	//Call the consensus
	files.soap_result.close();
	files.soap_result.clear();
	files.soap_result.open(alignment_name.c_str());
	files.soap_result.clear();

#ifdef TIME_BREAK
	time_t start_time, end_time;
	gettimeofday(&tv, NULL);
	start_time = tv.tv_sec*1000 + tv.tv_usec/1000;
#endif //TIME_BREAK
	
	new_soap2cns(files.soap_result, files.consensus, files.baseinfo, para);

#ifdef TIME_BREAK
	gettimeofday(&tv, NULL);
	end_time = tv.tv_sec*1000 + tv.tv_usec/1000;
	fprintf(stderr, "Call soap2cns() costs %.2f sec.\n", (float)(end_time - start_time)/1000);
#endif //TIME_BREAK
	
	files.soap_result.close();
	files.consensus.close();
	cerr<<"Consensus Done!"<<endl;

	fprintf(stderr, "The whole program costs %.2f sec\n", (float)(end_time - start_main)/1000);
	
	MPI_Finalize(); // [cuiyb] finalize all MPI processes
	return 0;
}
