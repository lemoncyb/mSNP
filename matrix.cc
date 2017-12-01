#include "soap_snp.h"

int p_matrix_gen(std::ifstream & alignment, Parameter * para, std::fstream &mat_out) {
	// Read Alignment files
	Soap_format soap;
	rate_t * p_matrix = new rate_t [256*256*4*4]; // 8bit: q_max, 8bit: read_len, 4bit: number of types of all mismatch/match 4x4
	ubit64_t * count_matrix = new ubit64_t [256*256*4*4];

	Chr_name current_chr;

	current_chr = "";
	ubit64_t ref(0);
	std::string::size_type coord;
	for(std::string line; getline(alignment, line);) {
		std::istringstream s(line);
		if( s>> soap ) {
			if(soap.get_pos() < 0) {
				continue;
			}
			//cerr<<soap<<endl;
			// In the overloaded "+" above, soap.position will be substracted by 1 so that coordiates start from 0
		//	if (current_chr == genome->chromosomes.end() || current_chr->first != soap.get_chr_name()) {
			if (current_chr == "" || current_chr != soap.get_chr_name()) {
			//	current_chr = genome->chromosomes.find(soap.get_chr_name());
				current_chr = soap.get_chr_name();
				
				}
			else {
				;
			}
			if (soap.get_pos()+soap.get_read_len()>=para->ref_len) {
				continue;
			}
			if (soap.is_unique()) {
				for(coord = 0; coord != soap.get_read_len(); coord++) {
					if (soap.is_N(coord)) {
						;
					}
					else {
						if(! (soap.get_pos()+coord<para->ref_len)) {
							cerr<<soap<<endl;
							cerr<<"The program found the above read has exceed the reference length:\n";
							cerr<<"The read is aligned to postion: "<<soap.get_pos()<<" with read length: "<<soap.get_read_len()<<endl;
							cerr<<"Reference: "<<current_chr<<" FASTA Length: "<<para->ref_len<<endl;
							//exit(255);
						}
						//ref = current_chr->second->get_bin_base(soap.get_pos()+coord);
						ref = new_get_bin_base(soap.get_pos()+coord);
						if ( (ref&12) !=0 ) {
							// This is an N on reference or a dbSNP which should be excluded from calibration
							;
						}
						else {
							if(soap.is_fwd()) {
								// forward strand
								count_matrix[(((ubit64_t)soap.get_qual(coord))<<12) | (coord<<4) | ((ref&0x3)<<2) | (soap.get_base(coord)>>1)&3] += 1;
							}
							else {
								// reverse strand
								count_matrix[(((ubit64_t)soap.get_qual(coord))<<12) | ((soap.get_read_len()-1-coord)<<4) | ((ref&0x3)<<2) | (soap.get_base(coord)>>1)&3] += 1;
							}
						}
					}
				}
			}
		}
	}
	ubit64_t o_base/*o_based base*/, t_base/*theorecical(supposed) base*/, type, sum[4], same_qual_count_by_type[16], same_qual_count_by_t_base[4], same_qual_count_total, same_qual_count_mismatch;
	char q_char/*fastq quality char*/;

	const ubit64_t sta_pow=10; // minimum number to say statistically powerful
	for(q_char=para->q_min; q_char<=para->q_max ;q_char++) {
		memset(same_qual_count_by_type, 0, sizeof(ubit64_t)*16);
		memset(same_qual_count_by_t_base, 0, sizeof(ubit64_t)*4);
		same_qual_count_total = 0;
		same_qual_count_mismatch = 0;
		for(coord=0; coord != para->read_length ; coord++) {
			for(type=0;type!=16;type++) {
				// If the sample is small, then we will not consider the effect of read cycle.
				same_qual_count_by_type[type] += count_matrix[ ((ubit64_t)q_char<<12) | coord <<4 | type];
				same_qual_count_by_t_base[(type>>2)&3] += count_matrix[ ((ubit64_t)q_char<<12) | coord <<4 | type];
				same_qual_count_total += count_matrix[ ((ubit64_t)q_char<<12) | coord <<4 | type];
				//cerr<<(int)type<<'\t'<<same_qual_count_by_type[type]<<'\t'<<same_qual_count_by_t_base[0]<<endl;
				if(type % 5 != 0) {
					// Mismatches
					same_qual_count_mismatch += count_matrix[ ((ubit64_t)q_char<<12) | coord <<4 | type];
				}
			}
		}
		for(coord=0; coord != para->read_length ; coord++) {
			//cerr<<(q_char)<<'\t'<<coord;
			memset(sum, (ubit64_t)0, sizeof(ubit64_t)*4);
			// Count of all ref base at certain coord and quality
			for(type=0;type!=16;type++) {
				sum[(type>>2)&3] += count_matrix[ ((ubit64_t)q_char<<12) | (coord <<4) | type]; // (type>>2)&3: the ref base
				//cerr<<sum[type&3]<<endl;
			}
			for(t_base=0; t_base!=4; t_base++) {
				for(o_base=0; o_base!=4; o_base++) {
					//cerr<<'\t'<<count_matrix[ ((ubit64_t)q_char<<12) | (coord <<4) | (t_base<<2) | o_base];
					if (count_matrix[ ((ubit64_t)q_char<<12) | (coord <<4) | (t_base<<2) | o_base] > sta_pow) {
						// Statistically powerful
						p_matrix [ ((ubit64_t)(q_char-para->q_min)<<12) | (coord <<4) | (t_base<<2) | o_base] = ((double)count_matrix[ ((ubit64_t)q_char<<12) | (coord <<4) | (t_base<<2) | o_base]) / sum[t_base];
					}
					else if (same_qual_count_by_type[t_base<<2|o_base] > sta_pow) {
						// Smaller sample, given up effect from read cycle
						p_matrix [ ((ubit64_t)(q_char-para->q_min)<<12) | (coord <<4) | (t_base<<2) | o_base] =  ((double)same_qual_count_by_type[t_base<<2|o_base]) / same_qual_count_by_t_base[t_base];
					}
					else if (same_qual_count_total > 0){
						// Too small sample, given up effect of mismatch types
						if (o_base == t_base) {
							p_matrix [ ((ubit64_t)(q_char-para->q_min)<<12) | (coord <<4) | (t_base<<2) | o_base] = ((double)(same_qual_count_total-same_qual_count_mismatch))/same_qual_count_total;
						}
						else {
							p_matrix [ ((ubit64_t)(q_char-para->q_min)<<12) | (coord <<4) | (t_base<<2) | o_base] = ((double)same_qual_count_mismatch)/same_qual_count_total;
						}
					}
					else {
						;
					}

					// For these cases like:
					// Ref: G o_base: G x10 Ax5. When calculate the probability of this allele to be A,
					// If there's no A in reference gives observation of G, then the probability will be zero,
					// And therefore exclude the possibility of this pos to have an A
					// These cases should be avoid when the dataset is large enough
					// If no base with certain quality is o_based, it also doesn't matter
					if( (p_matrix [ ((ubit64_t)(q_char-para->q_min)<<12) | (coord <<4) | (t_base<<2) | o_base]==0) || p_matrix [ ((ubit64_t)(q_char-para->q_min)<<12) | (coord <<4) | (t_base<<2) | o_base] ==1) {
						if (o_base == t_base) {
							p_matrix [ ((ubit64_t)(q_char-para->q_min)<<12) | (coord <<4) | (t_base<<2) | o_base] = (1-pow(10, -((q_char-para->q_min)/10.0)));
							if(p_matrix [ ((ubit64_t)(q_char-para->q_min)<<12) | (coord <<4) | (t_base<<2) | o_base]<0.25) {
								p_matrix [ ((ubit64_t)(q_char-para->q_min)<<12) | (coord <<4) | (t_base<<2) | o_base] = 0.25;
							}
						}
						else {
							p_matrix [ ((ubit64_t)(q_char-para->q_min)<<12) | (coord <<4) | (t_base<<2) | o_base] = (pow(10, -((q_char-para->q_min)/10.0))/3);
							if(p_matrix [ ((ubit64_t)(q_char-para->q_min)<<12) | (coord <<4) | (t_base<<2) | o_base]>0.25) {
								p_matrix [ ((ubit64_t)(q_char-para->q_min)<<12) | (coord <<4) | (t_base<<2) | o_base] = 0.25;
							}
						}
					}
				}
			}
		}
		//cerr<<'\t'<<same_qual_count_by_type[0]<<'\t'<<same_qual_count_by_t_base[0]<<endl;
	}
	
	//[cuiyb] generate new_p_matrix
	new_p_matrix_gen(p_matrix);
	
	if (para->is_matrix_out)
		p_matrix_write(mat_out, p_matrix, para);
	
	delete [] p_matrix;
	delete [] count_matrix;

	//memset(&p_matrix[((ubit64_t)para->q_max-para->q_min+1)<<12], 0, 256*256*4*4 - (para->q_max-para->q_min+1)*256*4*4);
	// Note: from now on, the first 8 bit of p_matrix is its quality score, not the FASTQ char
	return 1;
}

//[cuiyb]++ generate new_p_matrix
int new_p_matrix_gen(rate_t * p_matrix)
{
	for(int q_adj = 0; q_adj < 256; q_adj++) {
		for(int coord = 0; coord < 256; coord++) {
			for(int o_base = 0; o_base < 4; o_base++) {
				const unsigned int new_idxBase = (((unsigned int)(q_adj)<< 10)|((coord)<< 2)|(o_base))*16; // offset=16, type_likely
				
				//[cuiyb] non-zero elements of type_lilely
				for(int allele1 = 0; allele1 < 4; allele1++) {
					for(int allele2 = allele1; allele2 < 4; allele2++) {
						new_p_matrix[new_idxBase + (allele1<<2|allele2)] =
							log10(
								0.5*p_matrix[(((ubit64_t)(q_adj))<< 12)|((coord)<< 4)|((allele1)<< 2)|(o_base)]
								+
								0.5*p_matrix[(((ubit64_t)(q_adj))<< 12)|((coord)<< 4)|((allele2)<< 2)|(o_base)]);
					}
				}
				//[cuiyb] set other elements to be zero, index=4,8,9,12,13,14
				new_p_matrix[new_idxBase + 4] = 0;
				new_p_matrix[new_idxBase + 8] = 0;
				new_p_matrix[new_idxBase + 9] = 0;
				new_p_matrix[new_idxBase + 12] = 0;
				new_p_matrix[new_idxBase + 13] = 0;
				new_p_matrix[new_idxBase + 14] = 0;
			}
		}
	}
	return 1;
}

int p_matrix_read(std::fstream &mat_in, Parameter * para) {
	int q_char, type;
	std::string::size_type coord;
	rate_t * p_matrix = new rate_t [256*256*4*4]; // 8bit: q_max, 8bit: read_len, 4bit: number of types of all mismatch/match 4x4
	for(std::string line; getline(mat_in, line);) {
		std::istringstream s(line);
		s>>q_char>>coord;
		for(type=0;type!=16;type++) {
			s>>p_matrix [ ((ubit64_t)q_char<<12) | (coord <<4) | type];
			//cerr<<q_char<<"|"<<coord<<"|"<<p_matrix [ ((ubit64_t)q_char<<12) | (coord <<4) | type]<<endl;
			//exit(0);
		}
	}
	
	//[cuiyb] generate new_p_matrix
	new_p_matrix_gen(p_matrix);
	
	delete [] p_matrix;
	
	return 1;
}

int p_matrix_write(std::fstream &mat_out, rate_t * p_matrix, Parameter * para) {
	for( char q_char = para->q_min; q_char <= para->q_max; q_char++ ) {
		for( std::string::size_type coord=0; coord != para->read_length; coord++) {
			mat_out<<((ubit64_t)q_char-para->q_min)<<'\t'<<coord;
			for(char type=0;type!=16;type++) {
				mat_out<<'\t'<<scientific<<showpoint<<setprecision(16)<<p_matrix [ ((ubit64_t)(q_char-para->q_min)<<12) | (coord <<4) | type];
			}
			mat_out<<endl;
		}
	}
	return 1;
}
