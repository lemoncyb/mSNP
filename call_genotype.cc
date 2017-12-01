#include "soap_snp.h"
#include <cassert>
#include "omp.h"
#include <mpi.h>


// calculated log10(q_score)
__attribute__((target(mic))) double * log10_q_score;

// member variables of Pos_info
__attribute__((target(mic))) unsigned char * ori_array;// = new unsigned char[win_size+read_len];
__attribute__((target(mic))) int * 	pos_array;//  = new int [win_size+read_len];
__attribute__((target(mic))) int * count_uni_array;//  = new int [(win_size+read_len)*4];
__attribute__((target(mic))) int * q_sum_array;//  = new int [(win_size+read_len)*4];
__attribute__((target(mic))) int * depth_array;//  = new int [win_size+read_len];
__attribute__((target(mic))) int * dep_uni_array;//  = new int [win_size+read_len];
__attribute__((target(mic))) int * repeat_time_array;//  = new int [win_size+read_len];
__attribute__((target(mic))) int * count_all_array;//  = new int [(win_size+read_len)*4];
__attribute__((target(mic))) ubit32_t * base_array;//  = new ubit32_t [(win_size+read_len)*max_depth];
__attribute__((target(mic))) int * base_array_length;//  = new int [win_size+read_len];

// array stores results
__attribute__((target(mic))) char * type1_array;//  = new char[win_size]; // assign -1 for N sites
__attribute__((target(mic))) char * base1_array;//  = new char[win_size];
__attribute__((target(mic))) char * base2_array;//  = new char[win_size];
__attribute__((target(mic))) int * q_cns_array;//  = new int[win_size];
__attribute__((target(mic))) double * rank_sum_array;//  = new double[win_size];

// pointer used to divide jobs for MIC and CPU
unsigned char * ori_ptr;
int * pos_ptr, *count_uni_ptr, *q_sum_ptr, *depth_ptr, *dep_uni_ptr, *repeat_time_ptr, *count_all_ptr, *base_array_length_ptr;
ubit32_t * base_array_ptr;
char * type1_ptr, *base1_ptr, *base2_ptr;
int * q_cns_ptr;
double * rank_sum_ptr;

#ifdef TIME_BREAK
	//[cuiyb]++ for time statistics
	time_t start_time, end_time;
	struct timeval tv;
	int time_callcns = 0; //time cost of likelihood
	int time_recycle = 0; //time cost of recycle()
#endif //TIME_BREAK

// function re-written for offload

// new_soap2cns(), add parameter read_length, contains all variables of Call_win
int new_soap2cns(std::ifstream & alignment, std::ofstream & consensus, std::ofstream & baseinfo, Parameter * para){
	
	//initialize variables
	Chr_name current_chr = ""; //[cuiyb]++
	unsigned int last_start = 0;
	bool recycled = false;
	int bin_seq_length = (para->ref_len%capacity==0)? (para->ref_len/capacity) : (1+(para->ref_len/capacity)) ;
	
	// calculate log10(q_score) before hand
	log10_q_score = new double[256];
	int q_score;
	#pragma ivdep
	for (q_score = 0; q_score < 256; q_score++)
		log10_q_score[q_score] = log10(q_score);

	// member variables of Pos_info
	ori_array = new unsigned char[win_size+read_len];
	pos_array = new int [win_size+read_len];
	count_uni_array = new int [(win_size+read_len)*4];
	q_sum_array = new int [(win_size+read_len)*4];
	depth_array = new int [win_size+read_len];
	dep_uni_array = new int [win_size+read_len];
	repeat_time_array = new int [win_size+read_len];
	count_all_array = new int [(win_size+read_len)*4];
	base_array = new ubit32_t [(win_size+read_len)*max_depth];
	base_array_length = new int [win_size+read_len];
	
	memset(count_uni_array, 0, sizeof(int)*(win_size+read_len)*4);
	memset(q_sum_array, 0, sizeof(int)*(win_size+read_len)*4);
	memset(depth_array, 0, sizeof(int)*(win_size+read_len));
	memset(dep_uni_array, 0, sizeof(int)*(win_size+read_len));
	memset(repeat_time_array, 0, sizeof(int)*(win_size+read_len));
	memset(count_all_array, 0, sizeof(int)*(win_size+read_len)*4);
	memset(base_array_length, 0, sizeof(int)*(win_size+read_len));
	
	type1_array = new char[win_size]; // assign -1 for N sites
	base1_array = new char[win_size];
	base2_array = new char[win_size];
	q_cns_array = new int[win_size];
	rank_sum_array = new double[win_size];
	
	// pointer used to divide jobs for MIC and CPU
	ori_ptr = ori_array;
	pos_ptr = pos_array;
	count_uni_ptr = count_uni_array;
	q_sum_ptr = q_sum_array;
	depth_ptr = depth_array;
	dep_uni_ptr = dep_uni_array;
	repeat_time_ptr = repeat_time_array;
	count_all_ptr = count_all_array;
	base_array_length_ptr = base_array_length;
	base_array_ptr = base_array;
	
	type1_ptr = type1_array;
	base1_ptr = base1_array;
	base2_ptr = base2_array;
	q_cns_ptr = q_cns_array;
	rank_sum_ptr = rank_sum_array;
	
	
	new_initialize(0);
	
	cerr << "soap2cns() malloc space on MIC start" << endl;
#ifdef TIME_BREAK
	gettimeofday(&tv, NULL);
	start_time = tv.tv_sec*1000 + tv.tv_usec/1000;
#endif //TIME_BREAK
	for (int mic_id = 0; mic_id < mic_num; mic_id++) {
	//upload new_p_matrix, mat->p_prior, para
	#pragma offload_transfer target(mic:mic_id) \
			in(new_p_matrix:length(256*256*4*16)					alloc_if(1) free_if(0) align(64)) \
			in(mat_p_prior:length(8*4*4)							alloc_if(1) free_if(0) align(64)) \
			in(bin_seq:length(bin_seq_length)						alloc_if(1) free_if(0) align(64)) \
			in(log10_q_score:length(256)							alloc_if(1) free_if(0) align(64))
				 
	//malloc space for sites on MIC
	#pragma offload_transfer target(mic:mic_id) \
			nocopy(ori_array:length(mic_size) 						alloc_if(1) free_if(0) align(64)) \
			nocopy(pos_array:length(mic_size) 						alloc_if(1) free_if(0) align(64)) \
			nocopy(count_uni_array:length(mic_size*4) 				alloc_if(1) free_if(0) align(64)) \
			nocopy(q_sum_array:length(mic_size*4) 					alloc_if(1) free_if(0) align(64)) \
			nocopy(depth_array:length(mic_size)						alloc_if(1) free_if(0) align(64)) \
			nocopy(dep_uni_array:length(mic_size) 					alloc_if(1) free_if(0) align(64)) \
			nocopy(repeat_time_array:length(mic_size) 				alloc_if(1) free_if(0) align(64)) \
			nocopy(count_all_array:length(mic_size*4) 				alloc_if(1) free_if(0) align(64)) \
			nocopy(base_array:length(mic_size*max_depth)			alloc_if(1) free_if(0) align(64)) \
			nocopy(base_array_length:length(mic_size) 				alloc_if(1) free_if(0) align(64))

	//malloc space for results arrays on MIC
	#pragma offload_transfer target(mic:mic_id) \
			nocopy(type1_array:length(mic_size)						alloc_if(1) free_if(0) align(64)) \
			nocopy(base1_array:length(mic_size)						alloc_if(1) free_if(0) align(64)) \
			nocopy(base2_array:length(mic_size)						alloc_if(1) free_if(0) align(64)) \
			nocopy(q_cns_array:length(mic_size)						alloc_if(1) free_if(0) align(64)) \
			nocopy(rank_sum_array:length(mic_size)					alloc_if(1) free_if(0) align(64))
	}
#ifdef TIME_BREAK
	gettimeofday(&tv, NULL);
	end_time = tv.tv_sec*1000 + tv.tv_usec/1000;
	fprintf(stderr, "Initialize MIC and malloc space costs %.2f sec.\n", (float)(end_time - start_time)/1000);
#endif //TIME_BREAK
	cerr << "soap2cns() malloc space on MIC done" << endl;

	Soap_format soap;
	long long int win_start = 0;
	//current_chr = genome->chromosomes.end();
//	current_chr = "";
	seek_pointer(alignment); // divide windows
	for(std::string line; getline(alignment, line);) {
		std::istringstream ss(line);
		if( ss >> soap ) {
			//new_deal_read(soap, consensus, baseinfo, para);
			int coord, sub;
			if(soap.get_pos() < 0) {
				return 0;
			}
			// get the coordinate of the first base of current window
			if(win_start == 0) {
				win_start = soap.get_pos();
			}
			
			if (current_chr == "" ) {
				current_chr = para->chr_name;
				new_initialize(0);
				last_start = 0;
				cerr<<"Processing "<<current_chr<<endl;
			}
			
			recycled = false;
			shareBases(win_start, win_size);
			while (soap.get_pos()/win_size > last_start/win_size ) {
				//new_pro_win(consensus, baseinfo, para);
				if(last_start > pos_array[win_size-1]) {
					new_initialize((last_start/win_size)*win_size);
				}
				#ifdef TIME_BREAK
				get_start_time(); //[cuiyb]++
				#endif //TIME_BREAK
				//new_call_cns(current_chr, win_size, para, consensus, baseinfo);
				////////////////////////////////////////////////////////////
				omp_set_nested(true);
				#pragma omp parallel for num_threads(mic_num+1) private(ori_array, pos_array, count_uni_array, q_sum_array, depth_array, dep_uni_array,\
				repeat_time_array, count_all_array, base_array, base_array_length, type1_array, base1_array, base2_array, q_cns_array, rank_sum_array)
				for (int mic_id = 0; mic_id < mic_num+1; mic_id++) {
				
				//divide jobs for MIC and CPU
				ori_array = ori_ptr + mic_id*mic_size;
				pos_array = pos_ptr + mic_id*mic_size;
				count_uni_array = count_uni_ptr + mic_id*mic_size*4;
				q_sum_array = q_sum_ptr + mic_id*mic_size*4;
				depth_array = depth_ptr + mic_id*mic_size;
				dep_uni_array = dep_uni_ptr + mic_id*mic_size;
				repeat_time_array = dep_uni_ptr + mic_id*mic_size;
				count_all_array = count_all_ptr + mic_id*mic_size*4;
				base_array = base_array_ptr + mic_id*mic_size*max_depth;
				base_array_length = base_array_length_ptr + mic_id*mic_size;
				
				type1_array = type1_ptr + mic_id*mic_size;
				base1_array = base1_ptr + mic_id*mic_size;
				base2_array = base2_ptr + mic_id*mic_size;
				q_cns_array = q_cns_ptr + mic_id*mic_size;
				rank_sum_array = rank_sum_ptr + mic_id*mic_size;
				
				//substitute member variables of Parameter with local variables, for offload
				small_int para_read_length = para->read_length;
				int para_glf_format = para->glf_format;
				bool para_is_snp_only = para->is_snp_only;
				rate_t para_global_dependency = para->global_dependency;
				rate_t para_pcr_dependency = para->pcr_dependency;
				bool para_is_monoploid = para->is_monoploid;
				int call_length = (mic_id < mic_num) ? mic_size : cpu_size;
				int num_threads = (mic_id < mic_num) ? 224 : para->cpu_thread_num;
				
				#pragma offload target(mic:mic_id) if(mic_id < mic_num)\
				nocopy(new_p_matrix:length(256*256*4*16)				alloc_if(0) free_if(0)) \
				nocopy(mat_p_prior:length(8*4*4)						alloc_if(0) free_if(0)) \
				nocopy(bin_seq:length(bin_seq_length)					alloc_if(0) free_if(0)) \
				nocopy(log10_q_score:length(256)						alloc_if(0) free_if(0)) \
					\
				nocopy(ori_array:length(mic_size) 						alloc_if(0) free_if(0)) \
				in(pos_array:length(mic_size) 							alloc_if(0) free_if(0)) \
				in(count_uni_array:length(mic_size*4) 					alloc_if(0) free_if(0)) \
				in(q_sum_array:length(mic_size*4) 						alloc_if(0) free_if(0)) \
				in(depth_array:length(mic_size)							alloc_if(0) free_if(0)) \
				in(dep_uni_array:length(mic_size) 						alloc_if(0) free_if(0)) \
				in(repeat_time_array:length(mic_size) 					alloc_if(0) free_if(0)) \
				in(count_all_array:length(mic_size*4) 					alloc_if(0) free_if(0)) \
				in(base_array:length(mic_size*max_depth) 				alloc_if(0) free_if(0)) \
				in(base_array_length:length(mic_size) 					alloc_if(0) free_if(0)) \
					\
				out(type1_array:length(mic_size)						alloc_if(0) free_if(0)) \
				out(base1_array:length(mic_size)						alloc_if(0) free_if(0)) \
				out(base2_array:length(mic_size)						alloc_if(0) free_if(0)) \
				out(q_cns_array:length(mic_size)						alloc_if(0) free_if(0)) \
				out(rank_sum_array:length(mic_size)						alloc_if(0) free_if(0))
				#pragma omp parallel for num_threads(num_threads)
					for(int j=0; j < call_length; j++){
						//	printf("new_call_cns() enter offload region \n");
						type1_array[j] = -1; // assign -1 for N sites as default
						//[cuiyb] all the variables are thread private, so declare them in the loop.
						std::string::size_type coord;
						small_int k;
						ubit64_t o_base, strand;
						char allele1, allele2, genotype, type, type1/*best genotype*/, type2/*suboptimal genotype*/, base1, base2, base3;
						int i, q_score, q_adjusted, qual1, qual2, qual3, q_cns, all_count1, all_count2, all_count3;
						int global_dep_count;//[cuiyb], *pcr_dep_count;
						//pcr_dep_count = new int [para->read_length*2];
						int pcr_dep_count[para_read_length*2];
						double  rank_sum_test_value, binomial_test_value;
						bool is_out;
						double real_p_prior[16];
						double likelihoods[10];
						memset(likelihoods, 0, sizeof(double)*10);
						//substitute mat->type_likely, mat->type_prob
						double new_type_likely[16];
						double new_type_prob[17];
						
						ori_array[j] = (new_get_bin_base(pos_array[j]))&0xF;
						if( ((ori_array[j]&4) !=0)/*an N*/ && depth_array[j] == 0) {
							continue;
						}
						base1 = 0, base2 = 0, base3 = 0;
						qual1 = -1, qual2= -2, qual3 = -3;
						all_count1 = 0, all_count2 =0, all_count3 =0;
						if(dep_uni_array[j]) {
							// This position is uniquely covered by at least one nucleotide
							for(i=0;i!=4;i++) {
								// i is four kind of alleles
								if(q_sum_array[j*4+i] >= qual1) {
									base3 = base2;
									qual3 = qual2;
									base2 = base1;
									qual2 = qual1;
									base1 = i;
									qual1 = q_sum_array[j*4+i];
								}
								else if (q_sum_array[j*4+i]>=qual2) {
									base3 = base2;
									qual3 = qual2;
									base2 = i;
									qual2  = q_sum_array[j*4+i];
								}
								else if (q_sum_array[j*4+i]>=qual3) {
									base3 = i;
									qual3  = q_sum_array[j*4+i];
								}
								else {
									;
								}
							}
							if(qual1 == 0) {
								// Adjust the best base so that things won't look ugly if the pos is not covered
								base1 = (ori_array[j]&7);
							}
							else if(qual2 ==0 && base1 != (ori_array[j]&7)) {
								base2 = (ori_array[j]&7);
							}
							else {
								;
							}
						}
						else {
							// This position is covered by all repeats
							for(i=0;i!=4;i++) {
								if(count_all_array[j*4+i] >= all_count1) {
									base3 = base2;
									all_count3 = all_count2;
									base2 = base1;
									all_count2 = all_count1;
									base1 = i;
									all_count1 = count_all_array[j*4+i];
								}
								else if (count_all_array[j*4+i]>=all_count2) {
									base3 = base2;
									all_count3 = all_count2;
									base2 = i;
									all_count2  = count_all_array[j*4+i];
								}
								else if (count_all_array[j*4+i]>=all_count3) {
									base3 = i;
									all_count3  = count_all_array[j*4+i];
								}
								else {
									;
								}
							}
							if(all_count1 ==0) {
								base1 = (ori_array[j]&7);
							}
							else if(all_count2 ==0 && base1 != (ori_array[j]&7)) {
								base2 = (ori_array[j]&7);
							}
							else {
								;
							}
						}
				
						new_sort_base(base_array+j*max_depth, 0, base_array_length[j]-1);  //std::sort has sort base_array in ascending order, base is sorted too.


						memset(new_type_likely, 0.0, sizeof(double)*16);
						
						ubit64_t pre_base = 0;
						global_dep_count = -1;
						memset(pcr_dep_count, 0, sizeof(int)*2*para_read_length);
						
						for (int ix = 0; ix != base_array_length[j]; ++ix){
							o_base = (base_array[j*max_depth+ix]&0x18000)>>15; /*1,1000,0000,0000,0000*/
							if (count_uni_array[j*4+o_base]==0) {continue;}
							if (o_base != pre_base){ /* check current o_base and reset global_dep_count & pcr_dep_count */
								pre_base = o_base;
								global_dep_count = -1;
								memset(pcr_dep_count, 0, sizeof(int)*2*para_read_length);
							}
							q_score = 63-((base_array[j*max_depth+ix]&0x7E00)>>9); /*0111,1110,0000,0000*/
							coord = (base_array[j*max_depth+ix]&0x1FE)>>1; /*0000,0001,1111,1110*/
							strand = base_array[j*max_depth+ix]&0x1; /*0000,0000,0000,0001*/
							
							if(pcr_dep_count[strand*para_read_length+coord]==0)
								global_dep_count += 1;
							
							pcr_dep_count[strand*para_read_length+coord] += 1;
							
							//q_adjusted = int( pow(10, (log10(q_score) +(pcr_dep_count[strand*para_read_length+coord]-1)*para_pcr_dependency +global_dep_count*para_global_dependency)) +0.5);
							q_adjusted = int( pow(10, (log10_q_score[q_score] +(pcr_dep_count[strand*para_read_length+coord]-1)*para_pcr_dependency +global_dep_count*para_global_dependency)) +0.5);
							if(q_adjusted < 1)
								q_adjusted = 1;
							
							const unsigned int new_idxBase = (((unsigned int)(q_adjusted)<< 10)|((coord)<< 2)|(o_base))*16; // offset=16, type_likely
							#pragma ivdep
							for(int l = 0; l < 16; l++){
								new_type_likely[l] += new_p_matrix[new_idxBase+l];
							}
						}
				
						// Calculate Posterior Probability
						memcpy(real_p_prior, &mat_p_prior[((ubit64_t)ori_array[j]&0x7)<<4], sizeof(double)*16);
				
						memset(new_type_prob,0,sizeof(rate_t)*17);
						type2=type1=16;
						// divide the for loop into two loop, one for calculation, one for juge
						#pragma ivdep
						for (genotype = 0; genotype < 16; genotype++) {
							new_type_prob[genotype] = new_type_likely[genotype] + log10(real_p_prior[genotype]) ;
						}	
						for (allele1=0; allele1!=4; allele1++) {
							for (allele2=allele1; allele2!=4; allele2++) {
								genotype = allele1<<2|allele2;
								if (new_type_prob[genotype] >= new_type_prob[type1] || type1 == 16) {
									// find the genotype which type_prob is biggest.
									type2 = type1;
									type1 = genotype;
								}
								else if (new_type_prob[genotype] >= new_type_prob[type2] || type2 ==16) {
									// find the genotype which type_prob is second biggest.
									type2 = genotype;
								}
								else {
									;
								}
							}
						}
						
						is_out = true; // Check if the position needs to be output, useful in snp-only mode
				
						rank_sum_test_value = 1.0;
				
						q_cns = (int)(10*(new_type_prob[type1]-new_type_prob[type2])+10*log10(rank_sum_test_value));
				
						if ( (type1&3) == ((type1>>2)&3)) { // Called Homozygous
							if (qual1>0 && base1 != (type1&3)) {
								//Wired: best base is not the consensus!
								q_cns = 0;
							}
							else if (/*qual2>0 &&*/ q_cns > qual1-qual2) {
								// Should not bigger than this
								q_cns = qual1-qual2;
							}
							else {
								;
							}
						}
						else {	// Called Heterozygous
							if(q_sum_array[j*4+base1]>0 && q_sum_array[j*4+base2]>0 && type1 == (base1<base2 ? (base1<<2|base2):(base2<<2|base1))) {
								// The best bases are in the heterozygote
								if (q_cns > qual2-qual3) {
									q_cns = qual2-qual3;
								}
							}
							else {	// Ok, wired things happened
								q_cns = 0;
							}
						}
						if(q_cns>99) {
							q_cns = 99;
						}
						if (q_cns<0) {
							q_cns = 0;
						}
						//store results to arrays
						type1_array[j] = type1;
						base1_array[j] = base1;
						base2_array[j] = base2;
						q_cns_array[j] = q_cns;
						rank_sum_array[j] = rank_sum_test_value;
					} //end of calculation for loop
				
					#pragma omp critical
					for(std::string::size_type j=0; j < call_length; j++){
						char base1 = base1_array[j], base2 = base2_array[j];
						if (type1_array[j] != -1){ // for non-N sites
							// ChrID\tPos\tRef\tCns\tQual\tBase1\tAvgQ1\tCountUni1\tCountAll1\tBase2\tAvgQ2\tCountUni2\tCountAll2\tDepth\tRank_sum\tCopyNum\tSNPstauts\n"
							if(! para->is_snp_only || (abbv[(type1_array[j])] != "ACTGNNNN"[(ori_array[j]&0x7)] && depth_array[j] > 0)) {
								if(base1<4 && base2<4) {
									consensus<<current_chr<<'\t'<<(pos_array[j]+1)<<'\t'<<("ACTGNNNN"[(ori_array[j]&0x7)])<<'\t'<<abbv[(type1_array[j])]<<'\t'<<q_cns_array[j]<<'\t'
											 <<("ACTGNNNN"[base1])<<'\t'<<(q_sum_array[j*4+base1]==0?0:q_sum_array[j*4+base1]/count_uni_array[j*4+base1])<<'\t'<<count_uni_array[j*4+base1]<<'\t'<<count_all_array[j*4+base1]<<'\t'
											 <<("ACTGNNNN"[base2])<<'\t'<<(q_sum_array[j*4+base2]==0?0:q_sum_array[j*4+base2]/count_uni_array[j*4+base2])<<'\t'<<count_uni_array[j*4+base2]<<'\t'<<count_all_array[j*4+base2]<<'\t'
											 <<depth_array[j]<<'\t'<<showpoint<<rank_sum_array[j]<<'\t'<<(depth_array[j]==0?255:(double)(repeat_time_array[j])/depth_array[j])<<'\t'<<((ori_array[j]&8)?1:0)<<endl;
								}
								else if(base1<4) {
									consensus<<current_chr<<'\t'<<(pos_array[j]+1)<<'\t'<<("ACTGNNNN"[(ori_array[j]&0x7)])<<'\t'<<abbv[(type1_array[j])]<<'\t'<<q_cns_array[j]<<'\t'
											 <<("ACTGNNNN"[base1])<<'\t'<<(q_sum_array[j*4+base1]==0?0:q_sum_array[j*4+base1]/count_uni_array[j*4+base1])<<'\t'<<count_uni_array[j*4+base1]<<'\t'<<count_all_array[j*4+base1]<<'\t'
											 <<"N\t0\t0\t0\t"
											 <<depth_array[j]<<'\t'<<showpoint<<rank_sum_array[j]<<'\t'<<(depth_array[j]==0?255:(double)(repeat_time_array[j])/depth_array[j])<<'\t'<<((ori_array[j]&8)?1:0)<<endl;
								}
								else {
									consensus<<current_chr<<'\t'<<(pos_array[j]+1)<<"\tN\tN\t0\tN\t0\t0\t0\tN\t0\t0\t0\t0\t1.000\t255.000\t0"<<endl;
								}
							}
						}
						else
							consensus<<current_chr<<'\t'<<(pos_array[j]+1)<<"\tN\tN\t0\tN\t0\t0\t0\tN\t0\t0\t0\t0\t1.000\t255.000\t0"<<endl;
					}
		
				} // end of multiple mics loop
				////////////////////////////////////////////////////////////
#ifdef TIME_BREAK
				get_end_time(time_callcns); //[cuiyb]++
#endif //TIME_BREAK

				last_start = pos_array[win_size-1];
				recycled = false;
				
				new_recycle();
				recycled = true;
				last_start = pos_array[win_size-1];
			}
			last_start=soap.get_pos();
			for(coord=0; coord<soap.get_read_len(); coord++){
				if( (soap.get_pos()+coord)/win_size == soap.get_pos()/win_size ) {
					// In the same sliding window
					sub = (soap.get_pos()+coord) % win_size;
				}
				else {
					sub = (soap.get_pos()+coord) % win_size + win_size; // Use the tail to store the info so that it won't intervene the uncalled bases
				}
				
				if (sub >= (win_size + para->read_length))
				{
					cerr << "\tERROR: the parameter '-L' was wrong, you should set the longest read's length!!!" 
						 << endl;
					//exit(0);
				}
				
				depth_array[sub] += 1;
				repeat_time_array[sub] += soap.get_hit();
				if((soap.is_N(coord)) || soap.get_qual(coord)<para->q_min || dep_uni_array[sub] >= 0xFF) {
					// An N, low quality or meaningless huge depth
					continue;
				}
				if(soap.get_hit() == 1) {
					//exit((fprintf(stderr, "Wo Cao!\n")));
					dep_uni_array[sub] += 1;
					// Update the covering info: 4x2x64x64 matrix, base x strand x q_score x read_pos, 2-1-6-6 bits for each
					if(soap.is_fwd()) {
						// Binary strand: 0 for plus and 1 for minus
						/* [cuiyb] new order in base_array for sorting, 2*6*8*1 array,base*q_score*read_pos*strand, 2-6-8-1 bits for each */
						base_array[sub*max_depth+base_array_length[sub]++]=(((ubit32_t)((soap.get_base(coord)&0x6)>>1))<<15 | ((ubit32_t)(63-soap.get_qual(coord)+para->q_min))<<9 | coord<<1 | 0);
					}
					else {
						base_array[sub*max_depth+base_array_length[sub]++]=(((ubit32_t)((soap.get_base(coord)&0x6)>>1))<<15 | ((ubit32_t)(63-soap.get_qual(coord)+para->q_min))<<9 | (soap.get_read_len()-1-coord)<<1 | 1 );
					}
					count_uni_array[sub*4+((soap.get_base(coord)>>1)&3)] += 1;
					q_sum_array[sub*4+((soap.get_base(coord)>>1)&3)] += (soap.get_qual(coord)-para->q_min);
				}
				count_all_array[sub*4+((soap.get_base(coord)>>1)&3)] += 1;
			}
		}
	}
	
#ifdef TIME_BREAK
	get_start_time(); //[cuiyb]++
#endif //TIME_BREAK

	omp_set_nested(true);
	#pragma omp parallel for num_threads(mic_num+1) private(ori_array, pos_array, count_uni_array, q_sum_array, depth_array, dep_uni_array,\
	repeat_time_array, count_all_array, base_array, base_array_length, type1_array, base1_array, base2_array, q_cns_array, rank_sum_array)
	for (int mic_id = 0; mic_id < mic_num+1; mic_id++) {
	
	//divide jobs for MIC and CPU
	ori_array = ori_ptr + mic_id*mic_size;
	pos_array = pos_ptr + mic_id*mic_size;
	count_uni_array = count_uni_ptr + mic_id*mic_size*4;
	q_sum_array = q_sum_ptr + mic_id*mic_size*4;
	depth_array = depth_ptr + mic_id*mic_size;
	dep_uni_array = dep_uni_ptr + mic_id*mic_size;
	repeat_time_array = dep_uni_ptr + mic_id*mic_size;
	count_all_array = count_all_ptr + mic_id*mic_size*4;
	base_array = base_array_ptr + mic_id*mic_size*max_depth;
	base_array_length = base_array_length_ptr + mic_id*mic_size;
	
	type1_array = type1_ptr + mic_id*mic_size;
	base1_array = base1_ptr + mic_id*mic_size;
	base2_array = base2_ptr + mic_id*mic_size;
	q_cns_array = q_cns_ptr + mic_id*mic_size;
	rank_sum_array = rank_sum_ptr + mic_id*mic_size;
	
	//substitute member variables of Parameter with local variables, for offload
	small_int para_read_length = para->read_length;
	int para_glf_format = para->glf_format;
	bool para_is_snp_only = para->is_snp_only;
	rate_t para_global_dependency = para->global_dependency;
	rate_t para_pcr_dependency = para->pcr_dependency;
	bool para_is_monoploid = para->is_monoploid;
	int call_length = (mic_id < mic_num) ? mic_size : cpu_size;
	int num_threads = (mic_id < mic_num) ? 224 : para->cpu_thread_num;
	
	#pragma offload target(mic:mic_id) if(mic_id < mic_num)\
	nocopy(new_p_matrix:length(256*256*4*16)					alloc_if(0) free_if(0)) \
	nocopy(mat_p_prior:length(8*4*4)							alloc_if(0) free_if(0)) \
	nocopy(bin_seq:length(bin_seq_length)						alloc_if(0) free_if(0)) \
	nocopy(log10_q_score:length(256)							alloc_if(0) free_if(0)) \
		\
	nocopy(ori_array:length(mic_size) 						alloc_if(0) free_if(0)) \
	in(pos_array:length(mic_size) 							alloc_if(0) free_if(0)) \
	in(count_uni_array:length(mic_size*4) 					alloc_if(0) free_if(0)) \
	in(q_sum_array:length(mic_size*4) 						alloc_if(0) free_if(0)) \
	in(depth_array:length(mic_size)							alloc_if(0) free_if(0)) \
	in(dep_uni_array:length(mic_size) 						alloc_if(0) free_if(0)) \
	in(repeat_time_array:length(mic_size) 					alloc_if(0) free_if(0)) \
	in(count_all_array:length(mic_size*4) 					alloc_if(0) free_if(0)) \
	in(base_array:length(mic_size*max_depth) 				alloc_if(0) free_if(0)) \
	in(base_array_length:length(mic_size) 					alloc_if(0) free_if(0)) \
		\
	out(type1_array:length(mic_size)						alloc_if(0) free_if(0)) \
	out(base1_array:length(mic_size)						alloc_if(0) free_if(0)) \
	out(base2_array:length(mic_size)						alloc_if(0) free_if(0)) \
	out(q_cns_array:length(mic_size)						alloc_if(0) free_if(0)) \
	out(rank_sum_array:length(mic_size)						alloc_if(0) free_if(0))
   #pragma omp parallel for num_threads(num_threads)
		for(int j=0; j < call_length; j++){
			//	printf("new_call_cns() enter offload region \n");
			type1_array[j] = -1; // assign -1 for N sites as default
			
			std::string::size_type coord;
			small_int k;
			ubit64_t o_base, strand;
			char allele1, allele2, genotype, type, type1/*best genotype*/, type2/*suboptimal genotype*/, base1, base2, base3;
			int i, q_score, q_adjusted, qual1, qual2, qual3, q_cns, all_count1, all_count2, all_count3;
			int global_dep_count;//[cuiyb], *pcr_dep_count;
			//pcr_dep_count = new int [para->read_length*2];
			int pcr_dep_count[para_read_length*2];
			double  rank_sum_test_value, binomial_test_value;
			bool is_out;
			double real_p_prior[16];
			double likelihoods[10];
			memset(likelihoods, 0, sizeof(double)*10);
			//substitute mat->type_likely, mat->type_prob
			double new_type_likely[16];
			double new_type_prob[17];
			
			ori_array[j] = (new_get_bin_base(pos_array[j]))&0xF;
			if( ((ori_array[j]&4) !=0)/*an N*/ && depth_array[j] == 0) {
				continue;
			}
			base1 = 0, base2 = 0, base3 = 0;
			qual1 = -1, qual2= -2, qual3 = -3;
			all_count1 = 0, all_count2 =0, all_count3 =0;
			if(dep_uni_array[j]) {
				// This position is uniquely covered by at least one nucleotide
				for(i=0;i!=4;i++) {
					// i is four kind of alleles
					if(q_sum_array[j*4+i] >= qual1) {
						base3 = base2;
						qual3 = qual2;
						base2 = base1;
						qual2 = qual1;
						base1 = i;
						qual1 = q_sum_array[j*4+i];
					}
					else if (q_sum_array[j*4+i]>=qual2) {
						base3 = base2;
						qual3 = qual2;
						base2 = i;
						qual2  = q_sum_array[j*4+i];
					}
					else if (q_sum_array[j*4+i]>=qual3) {
						base3 = i;
						qual3  = q_sum_array[j*4+i];
					}
					else {
						;
					}
				}
				if(qual1 == 0) {
					// Adjust the best base so that things won't look ugly if the pos is not covered
					base1 = (ori_array[j]&7);
				}
				else if(qual2 ==0 && base1 != (ori_array[j]&7)) {
					base2 = (ori_array[j]&7);
				}
				else {
					;
				}
			}
			else {
				// This position is covered by all repeats
				for(i=0;i!=4;i++) {
					if(count_all_array[j*4+i] >= all_count1) {
						base3 = base2;
						all_count3 = all_count2;
						base2 = base1;
						all_count2 = all_count1;
						base1 = i;
						all_count1 = count_all_array[j*4+i];
					}
					else if (count_all_array[j*4+i]>=all_count2) {
						base3 = base2;
						all_count3 = all_count2;
						base2 = i;
						all_count2  = count_all_array[j*4+i];
					}
					else if (count_all_array[j*4+i]>=all_count3) {
						base3 = i;
						all_count3  = count_all_array[j*4+i];
					}
					else {
						;
					}
				}
				if(all_count1 ==0) {
					base1 = (ori_array[j]&7);
				}
				else if(all_count2 ==0 && base1 != (ori_array[j]&7)) {
					base2 = (ori_array[j]&7);
				}
				else {
					;
				}
			}
	
			new_sort_base(base_array+j*max_depth, 0, base_array_length[j]-1);  //std::sort has sort base_array in ascending order, base is sorted too.
			
			for(genotype=0;genotype!=16;genotype++){
				new_type_likely[genotype] = 0.0;
			}
			
			ubit64_t pre_base = 0;
			global_dep_count = -1;
			memset(pcr_dep_count, 0, sizeof(int)*2*para_read_length);
			
			for (int ix = 0; ix != base_array_length[j]; ++ix){
				o_base = (base_array[j*max_depth+ix]&0x18000)>>15; /*1,1000,0000,0000,0000*/
				if (count_uni_array[j*4+o_base]==0) {continue;}
				if (o_base != pre_base){ /* check current o_base and reset global_dep_count & pcr_dep_count */
					pre_base = o_base;
					global_dep_count = -1;
					memset(pcr_dep_count, 0, sizeof(int)*2*para_read_length);
				}
				q_score = 63-((base_array[j*max_depth+ix]&0x7E00)>>9); /*0111,1110,0000,0000*/
				coord = (base_array[j*max_depth+ix]&0x1FE)>>1; /*0000,0001,1111,1110*/
				strand = base_array[j*max_depth+ix]&0x1; /*0000,0000,0000,0001*/
				
				if(pcr_dep_count[strand*para_read_length+coord]==0)
					global_dep_count += 1;
				
				pcr_dep_count[strand*para_read_length+coord] += 1;
				
				
				q_adjusted = int( pow(10, (log10_q_score[q_score] +(pcr_dep_count[strand*para_read_length+coord]-1)*para_pcr_dependency +global_dep_count*para_global_dependency)) +0.5);
				if(q_adjusted < 1)
					q_adjusted = 1;
				
				const unsigned int new_idxBase = (((unsigned int)(q_adjusted)<< 10)|((coord)<< 2)|(o_base))*16; // offset=16, type_likely
				#pragma ivdep
				for(int l = 0; l < 16; l++){
					new_type_likely[l] += new_p_matrix[new_idxBase+l];
				}
			}
	
			// Calculate Posterior Probability
			memcpy(real_p_prior, &mat_p_prior[((ubit64_t)ori_array[j]&0x7)<<4], sizeof(double)*16);
	
			memset(new_type_prob,0,sizeof(rate_t)*17);
			type2=type1=16;
			// divide the for loop into two loop, one for calculation, one for juge
			#pragma ivdep
			for (genotype = 0; genotype < 16; genotype++) {
				new_type_prob[genotype] = new_type_likely[genotype] + log10(real_p_prior[genotype]) ;
			}	
			for (allele1=0; allele1!=4; allele1++) {
				for (allele2=allele1; allele2!=4; allele2++) {
					genotype = allele1<<2|allele2;
					if (new_type_prob[genotype] >= new_type_prob[type1] || type1 == 16) {
						// find the genotype which type_prob is biggest.
						type2 = type1;
						type1 = genotype;
					}
					else if (new_type_prob[genotype] >= new_type_prob[type2] || type2 ==16) {
						// find the genotype which type_prob is second biggest.
						type2 = genotype;
					}
					else {
						;
					}
				}
			}
						
			is_out = true; // Check if the position needs to be output, useful in snp-only mode
	
			rank_sum_test_value = 1.0;
	
			q_cns = (int)(10*(new_type_prob[type1]-new_type_prob[type2])+10*log10(rank_sum_test_value));
	
			if ( (type1&3) == ((type1>>2)&3)) { // Called Homozygous
				if (qual1>0 && base1 != (type1&3)) {
					//Wired: best base is not the consensus!
					q_cns = 0;
				}
				else if (/*qual2>0 &&*/ q_cns > qual1-qual2) {
					// Should not bigger than this
					q_cns = qual1-qual2;
				}
				else {
					;
				}
			}
			else {	// Called Heterozygous
				if(q_sum_array[j*4+base1]>0 && q_sum_array[j*4+base2]>0 && type1 == (base1<base2 ? (base1<<2|base2):(base2<<2|base1))) {
					// The best bases are in the heterozygote
					if (q_cns > qual2-qual3) {
						q_cns = qual2-qual3;
					}
				}
				else {	// Ok, wired things happened
					q_cns = 0;
				}
			}
			if(q_cns>99) {
				q_cns = 99;
			}
			if (q_cns<0) {
				q_cns = 0;
			}
			//store results to arrays
			type1_array[j] = type1;
			base1_array[j] = base1;
			base2_array[j] = base2;
			q_cns_array[j] = q_cns;
			rank_sum_array[j] = rank_sum_test_value;
		} //end of calculation for loop
		
		#pragma omp critical
		for(std::string::size_type j=0; j < call_length; j++){
			char base1 = base1_array[j], base2 = base2_array[j];
			if (type1_array[j] != -1){ // for non-N sites
				// ChrID\tPos\tRef\tCns\tQual\tBase1\tAvgQ1\tCountUni1\tCountAll1\tBase2\tAvgQ2\tCountUni2\tCountAll2\tDepth\tRank_sum\tCopyNum\tSNPstauts\n"
				if(! para->is_snp_only || (abbv[(type1_array[j])] != "ACTGNNNN"[(ori_array[j]&0x7)] && depth_array[j] > 0)) {
					if(base1<4 && base2<4) {
						consensus<<current_chr<<'\t'<<(pos_array[j]+1)<<'\t'<<("ACTGNNNN"[(ori_array[j]&0x7)])<<'\t'<<abbv[(type1_array[j])]<<'\t'<<q_cns_array[j]<<'\t'
								 <<("ACTGNNNN"[base1])<<'\t'<<(q_sum_array[j*4+base1]==0?0:q_sum_array[j*4+base1]/count_uni_array[j*4+base1])<<'\t'<<count_uni_array[j*4+base1]<<'\t'<<count_all_array[j*4+base1]<<'\t'
								 <<("ACTGNNNN"[base2])<<'\t'<<(q_sum_array[j*4+base2]==0?0:q_sum_array[j*4+base2]/count_uni_array[j*4+base2])<<'\t'<<count_uni_array[j*4+base2]<<'\t'<<count_all_array[j*4+base2]<<'\t'
								 <<depth_array[j]<<'\t'<<showpoint<<rank_sum_array[j]<<'\t'<<(depth_array[j]==0?255:(double)(repeat_time_array[j])/depth_array[j])<<'\t'<<((ori_array[j]&8)?1:0)<<endl;
					}
					else if(base1<4) {
						consensus<<current_chr<<'\t'<<(pos_array[j]+1)<<'\t'<<("ACTGNNNN"[(ori_array[j]&0x7)])<<'\t'<<abbv[(type1_array[j])]<<'\t'<<q_cns_array[j]<<'\t'
								 <<("ACTGNNNN"[base1])<<'\t'<<(q_sum_array[j*4+base1]==0?0:q_sum_array[j*4+base1]/count_uni_array[j*4+base1])<<'\t'<<count_uni_array[j*4+base1]<<'\t'<<count_all_array[j*4+base1]<<'\t'
								 <<"N\t0\t0\t0\t"
								 <<depth_array[j]<<'\t'<<showpoint<<rank_sum_array[j]<<'\t'<<(depth_array[j]==0?255:(double)(repeat_time_array[j])/depth_array[j])<<'\t'<<((ori_array[j]&8)?1:0)<<endl;
					}
					else {
						consensus<<current_chr<<'\t'<<(pos_array[j]+1)<<"\tN\tN\t0\tN\t0\t0\t0\tN\t0\t0\t0\t0\t1.000\t255.000\t0"<<endl;
					}
				}
			}
			else
				consensus<<current_chr<<'\t'<<(pos_array[j]+1)<<"\tN\tN\t0\tN\t0\t0\t0\tN\t0\t0\t0\t0\t1.000\t255.000\t0"<<endl;
		}
	
	} // end of multiple mics loop
	//////////////////////////////////////////////////////
	#ifdef TIME_BREAK
	get_end_time(time_callcns); //[cuiyb]++
	#endif //TIME_BREAK

	for (int mic_id = 0; mic_id < mic_num; mic_id++) {
	// free sapce for new_p_matrix
	#pragma offload_transfer target(mic:mic_id) \
			nocopy(new_p_matrix:length(256*256*4*16)				alloc_if(0) free_if(1)) \
			nocopy(mat_p_prior:length(8*4*4)						alloc_if(0) free_if(1)) \
			nocopy(bin_seq:length(bin_seq_length)					alloc_if(0) free_if(1)) \
			nocopy(log10_q_score:length(256)						alloc_if(0) free_if(1))

	// free space for sites on MIC
	#pragma offload_transfer target(mic:mic_id) \
			nocopy(ori_array:length(mic_size) 						alloc_if(0) free_if(1)) \
			nocopy(pos_array:length(mic_size) 						alloc_if(0) free_if(1)) \
			nocopy(count_uni_array:length(mic_size*4) 				alloc_if(0) free_if(1)) \
			nocopy(q_sum_array:length(mic_size*4) 					alloc_if(0) free_if(1)) \
			nocopy(depth_array:length(mic_size)						alloc_if(0) free_if(1)) \
			nocopy(dep_uni_array:length(mic_size) 					alloc_if(0) free_if(1)) \
			nocopy(repeat_time_array:length(mic_size) 				alloc_if(0) free_if(1)) \
			nocopy(count_all_array:length(mic_size*4) 				alloc_if(0) free_if(1)) \
			nocopy(base_array:length(mic_size*max_depth)			alloc_if(0) free_if(1)) \
			nocopy(base_array_length:length(mic_size) 				alloc_if(0) free_if(1))
	
		//malloc space for results arrays on MIC
	#pragma offload_transfer target(mic:mic_id) \
			nocopy(type1_array:length(mic_size)						alloc_if(0) free_if(1)) \
			nocopy(base1_array:length(mic_size)						alloc_if(0) free_if(1)) \
			nocopy(base2_array:length(mic_size)						alloc_if(0) free_if(1)) \
			nocopy(q_cns_array:length(mic_size)						alloc_if(0) free_if(1)) \
			nocopy(rank_sum_array:length(mic_size)					alloc_if(0) free_if(1))
	}
			
	//free variables
	delete [] new_p_matrix;
	
	delete [] type1_array;
	delete [] base1_array;
	delete [] base2_array;
	delete [] q_cns_array;
	delete [] rank_sum_array;
	
	delete [] ori_array;
	delete [] pos_array;
	delete [] count_uni_array;
	delete [] q_sum_array;
	delete [] depth_array;
	delete [] dep_uni_array;
	delete [] repeat_time_array;
	delete [] count_all_array;
	delete [] base_array;
	delete [] base_array_length;
	
	#ifdef TIME_BREAK
	print_time_cost(time_callcns, "call_cns");
	print_time_cost(time_recycle, "recycle()");
	#endif //TIME_BREAK
	
	alignment.close();
	consensus.close();
	baseinfo.close();
	
	return 1;
}

// initialize sites array
int new_initialize(ubit64_t start) {
	std::string::size_type i;
	#pragma ivdep
	for(i=0; i != read_len + win_size ; i++) {
		pos_array[i] = i+start;
	}
	return 1;
}

int new_recycle() {
#ifdef TIME_BREAK
	get_start_time(); //[cuiyb]++
#endif //TIME_BREAK
	std::string::size_type i;
	for(i=0; i != read_len ; i++) {
		pos_array[i] = pos_array[i+win_size];
		ori_array[i] = ori_array[i+win_size];
		depth_array[i] = depth_array[i+win_size];
		repeat_time_array[i] = repeat_time_array[i+win_size];
		dep_uni_array[i] = dep_uni_array[i+win_size];
		base_array_length[i] = base_array_length[i+win_size];
		memcpy(base_array+i*max_depth, base_array+(i+win_size)*max_depth, sizeof(ubit32_t)*base_array_length[i+win_size]);
		memcpy(count_uni_array+i*4, count_uni_array+(i+win_size)*4, sizeof(int)*4);
		memcpy(q_sum_array+i*4, q_sum_array+(i+win_size)*4, sizeof(int)*4);
		memcpy(count_all_array+i*4, count_all_array+(i+win_size)*4, sizeof(int)*4);
	}

	// divide the for loop for vectorization
	#pragma ivdep
	for(i=read_len; i != read_len+win_size; i++) {
		ori_array[i] = 0xFF;
		depth_array[i] = 0;
		repeat_time_array[i] = 0;
		dep_uni_array[i] = 0;
		memset(base_array+i*max_depth,0,sizeof(ubit32_t)*base_array_length[i]); //[cuiyb]++
		base_array_length[i] = 0;
		memset(count_uni_array+i*4,0,sizeof(int)*4);
		memset(q_sum_array+i*4,0,sizeof(int)*4);
		memset(count_all_array+i*4,0,sizeof(int)*4);
	}
	for(i=read_len; i != read_len+win_size; i++) {
		pos_array[i] = pos_array[i-1]+1;
	}
#ifdef TIME_BREAK
	get_end_time(time_recycle); //[cuiyb]++
#endif //TIME_BREAK
	return 1;
}

/* [cuiyb]++ to sort base_array, fast sort, increasing */
void new_sort_base(ubit32_t array[],int left,int right)
{	
	if (right <= 0)
		return;

	int i = left, j = right;
	ubit32_t tmp;
	ubit32_t key = array[(left+right)/2];
	
	/* partition */
	while (i <= j){
		while ( array[i] < key)
			i++;
		while ( array[j] > key)
			j--;
		if (i <= j) {
			//std::swap(array[i],array[j]);
			tmp = array[i];
			array[i] = array[j];
			array[j] = tmp;
			i++;
			j--;
		}
	}

	/* recursion */
	if (left < j)
		new_sort_base(array, left, j);
	if (i < right)
		new_sort_base(array, i, right);
}

ubit64_t new_get_bin_base(std::string::size_type pos) {
	return (bin_seq[pos/capacity]>>(pos%capacity*4))&0xF; // All 4 bits
}

//[cuiyb]++ Read-based windows division strategy
void seek_pointer(std::ifstream & is)
{
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	is.seekg (0, ios::end);
	long long int fileSize = is.tellg(); // return the file size in bytes
	is.seekg (0, ios::beg);
	long long int offset = (fileSize/size) * rank;
	is.seekg (offset, ios::cur);
	
	// the file pointer may be in the middle of current line, skip to the head of next line
	if (rank != 0) {
		std::string line;
		getline(is, line);
	}
}

// calculate the long of overlapped bases
int overlap_length(int coord, int win_size) {
	int length = 0;
	for (int i = coord; i < win_size; i++)
		length += base_array_length[coord];
	
	return length;
}

void add_overlap(int * array, int * length) {
	int i = 0;
	while (length[i] > 0) {
		base_array[base_array_length[i]] = array[i];
	}
}

void shareBases (long long int start, int win_size) {
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	long long int next_start; 

	MPI_Request request;
	if (rank != 0) {
		MPI_Isend(&start, 1, MPI_LONG_LONG_INT, rank-1, 0, MPI_COMM_WORLD, &request);
		int * temp_array = (int*)calloc(sizeof(int), max_depth*read_len);
		int * temp_length = (int*)calloc(sizeof(int), read_len);
		MPI_Recv(temp_array, max_depth*read_len, MPI_UNSIGNED, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(temp_length, read_len, MPI_INT, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		add_overlap(temp_array, temp_length);
	}
	if (rank != (size-1)) {
		MPI_Recv(&next_start, 1, MPI_LONG_LONG_INT, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Send(&base_array[next_start-start], overlap_length(next_start-start, win_size), MPI_UNSIGNED, rank+1, 0, MPI_COMM_WORLD);
		MPI_Send(&base_array_length[next_start-start], win_size-(next_start-start)+1, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
	}
    if (rank != 0)
		MPI_Wait(&request, MPI_STATUS_IGNORE);
}

//[cuiyb]++ for time statistics
#ifdef TIME_BREAK
void get_start_time(void)
{
	gettimeofday(&tv, NULL);	
	start_time = tv.tv_sec*1000 + tv.tv_usec/1000;
}

void get_end_time(int &time_item)
{
	gettimeofday(&tv, NULL);
	end_time = tv.tv_sec*1000 + tv.tv_usec/1000;
	time_item += (end_time - start_time);
}

void print_time_cost(int time_item, string process)
{
	cerr << "Process " << process << " costs " << ((float)time_item/1000) << " sec.\n";
}
#endif // TIME_BREAK
