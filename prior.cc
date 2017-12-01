#include "soap_snp.h"

int p_prior_gen(Parameter * para) {
	char t_base, allele1, allele2;
	// Note, the above parameter should be changed to a more reasonable one
	for(t_base=0;t_base!=4;t_base++) {
		for(allele1=0;allele1!=4;allele1++) {
			for(allele2=allele1;allele2!=4;allele2++) {
				if(allele1 == t_base && allele2 == t_base) {
					// refHOM
					mat_p_prior[t_base<<4|allele1<<2|allele2] = 1;
				}
				else if (allele1 == t_base || allele2 == t_base) {
					// refHET: 1 ref 1 alt
					mat_p_prior[t_base<<4|allele1<<2|allele2] = para->het_novel_r;
				}
				else if (allele1 == allele2) {
					// altHOM
					mat_p_prior[t_base<<4|allele1<<2|allele2] = para->althom_novel_r;
				}
				else {
					// altHET: 2 diff alt base
					mat_p_prior[t_base<<4|allele1<<2|allele2] = para->het_novel_r * para->althom_novel_r;
				}
				if( para->transition_dominant && ((allele1^t_base) == 0x3 || (allele2^t_base) == 0x3)) {
					// transition
					mat_p_prior[t_base<<4|allele1<<2|allele2] *= 4;
				}
				//std::cerr<<"ACTG"[t_base]<<"\t"<<"ACTG"[allele1]<<"ACTG"[allele2]<<"\t"<<p_prior[t_base<<4|allele1<<2|allele2]<<endl;
			}
		}
	}
	for(allele1=0;allele1!=4;allele1++) {
		for(allele2=allele1;allele2!=4;allele2++) {
			// Deal with N
			mat_p_prior[0x4<<4|allele1<<2|allele2] = (allele1==allele2? 1: (2*para->het_novel_r)) * 0.25 *0.25;
			mat_p_prior[0x5<<4|allele1<<2|allele2] = (allele1==allele2? 1: (2*para->het_novel_r)) * 0.25 *0.25;
			mat_p_prior[0x6<<4|allele1<<2|allele2] = (allele1==allele2? 1: (2*para->het_novel_r)) * 0.25 *0.25;
			mat_p_prior[0x7<<4|allele1<<2|allele2] = (allele1==allele2? 1: (2*para->het_novel_r)) * 0.25 *0.25;
		}
	}
	return 1;
}

