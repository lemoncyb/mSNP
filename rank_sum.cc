#include "soap_snp.h"


int rank_table_gen() {
	// When N <= 63, (so that n1<=31), use this table to test
	ubit64_t i, n1, N, T1;
	rate_t p_left, p_right;

	// Calculate the factorials
	double * fact = new double [64];
	fact[0]=(double)1.0;
	for(i=1;i!=64;i++) {
		fact[i] = fact[i-1]*i;
	}

	ubit64_t * rank_sum= new ubit64_t [64*64*2048]; // 6bit: N; 5bit: n1; 11bit; T1
	ubit64_t * mat_p_rank = new ubit64_t [64*64];
	memset(rank_sum, 0, sizeof(ubit64_t)*64*64*2048);
	rank_sum[0]=1;
	for(N=1;N!=64;N++) {
		//cerr<<N<<endl;
		for(n1=0;n1<=N;n1++) {
			//cerr<<n1<<endl;
			for(T1=(1+n1)*n1/2;T1<=(N+N-n1+1)*n1/2;T1++) {
				//cerr<<T1<<endl;
				// Dynamic programming to generate the table
				rank_sum[N<<17|n1<<11|T1] = rank_sum[((N-1)<<17)|(n1<<11)|T1] + ((T1>=N && n1>0) ? rank_sum[((N-1)<<17)|((n1-1)<<11)|(T1-N)]:0);
				// Here, the p_rank is not cumulative
				mat_p_rank[(N<<17)|(n1<<11)|T1] = rank_sum[N<<17|n1<<11|T1] / (fact[N]/(fact[n1]*fact[N-n1]));
				//cerr<<p_rank[(N<<17)|(n1<<11)|T1]<<endl;
			}
			p_left = 0.0, p_right =1.0;
			for(T1=(1+n1)*n1/2;T1<=(N+N-n1+1)*n1/2;T1++) {
				p_right = 1.0 - p_left;
				p_left += mat_p_rank[(N<<17)|(n1<<11)|T1];
				mat_p_rank[N<<17|n1<<11|T1] = (p_left<p_right?p_left:p_right);
				//std::cerr<<N<<"\t"<<n1<<"\t"<<T1<<"\t"<<p_rank[(N<<17)|(n1<<11)|T1]<<endl;
			}
		}
	}
	delete [] rank_sum;
	delete [] fact;
	return 1;
}


double normal_test(int n1, int n2, double T1, double T2) {
	double u1, u2;
	u1 = (T1 - n1*(n1+n2+1)/2) / sqrt(n1*n2*(n1+n2+1)/(double)12);
	u2 = (T2 - n2*(n1+n2+1)/2) / sqrt(n1*n2*(n1+n2+1)/(double)12);
	return (fabs(u1)>fabs(u2)?u1:u2);
}

double table_test(double *mat_p_rank, int n1, int n2, double T1, double T2) {
	if(n1<=n2) {
		return mat_p_rank[(n1+n2)<<17|n1<<11|(int)(T1)]+(T1-(int)T1)*(mat_p_rank[(n1+n2)<<16|n1<<11|(int)(T1+1)]-mat_p_rank[(n1+n2)<<17|n1<<11|(int)(T1)]);
	}
	else {
		return mat_p_rank[(n1+n2)<<17|n2<<11|(int)(T2)]+(T2-(int)T2)*(mat_p_rank[(n1+n2)<<16|n2<<11|(int)(T2+1)]-mat_p_rank[(n1+n2)<<17|n2<<11|(int)(T2)]);
	}
}
