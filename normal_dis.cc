#include "soap_snp.h"

//global variables for offload region
double * new_p_matrix = new double[256*256*4*16]; // a_adjust*coord*base*array*length_of(type_likely), 32M size
rate_t * mat_p_prior = new rate_t [8*4*4]; // global for mat->p_prior

ubit64_t* bin_seq = NULL; // Sequence in binary format, for new_get_bin_base()
