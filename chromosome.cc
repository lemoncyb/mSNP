#include "soap_snp.h"


int genome_binarize(std::string & seq, Parameter * para) {
	para->ref_len = seq.length();
	//cerr<<len<<endl;
	// 4bit for each base
	// Allocate memory
	if (para->ref_len%capacity==0) {
		bin_seq = new ubit64_t [(para->ref_len/capacity)];
		memset(bin_seq,0,sizeof(ubit64_t)*(para->ref_len/capacity));
	}
	else {
		bin_seq = new ubit64_t [1+(para->ref_len/capacity)];
		memset(bin_seq,0,sizeof(ubit64_t)*(1+(para->ref_len/capacity)));
	}

	// Add each base, 7 is 0b111
	for(std::string::size_type i=0;i!=seq.length();i++) {
		bin_seq[i/capacity] |= ((((ubit64_t)seq[i]>>1)&7)<<(i%capacity*4));
	}
	return 1;
}

int genome_read(std::ifstream &fasta, Parameter * para) {
	std::string seq(""); //[cuiyb] store reference sequence
//	Chr_name current_name("");
	//map<Chr_name, Chr_info*>::iterator chr_iter;
	for(std::string buff;getline(fasta,buff);) {
		if('>' == buff[0]) {
			// Fasta id
			// Deal with previous chromosome
			
			std::string::size_type i;
			for(i=1;!isspace(buff[i]) && i != buff.length();i++) {
				;
			}
			Chr_name new_chr_name(buff,1,i-1);
		
			para->chr_name = new_chr_name;
			seq = "";
		}
		else {
			seq += buff;
		}
	}

	genome_binarize(seq, para); //[cuiyb] adopt binarize() as member function of Genome

	return 1;
}
