/** Computes likelihoods for each pair of indices, given matching
 * probabilities for each and a matrix of prior probabilities.  Return
 * the index pair that yields the maximum likelihood and the total
 * likelihood.  (The length of p5_ must be a multiple of 32 to make
 * vectorization easier.)
 *
 * @param v_  matrix of dimension (n7,n5_*32) containing the prior
 * @param p7_ vector of length n7 containing matching probabilities for
 *            the first index
 * @param n7  length of vector p7_
 * @param p5_ vector of length (n5_*32) containing matching
 *            probabilities for the second index
 * @param n5  length of vector p5_ divided by 32
 * @param pi7 pointer to location that receives index of the first index
 *            that yields the maximum likelihood (ignored if null)
 * @param pi5 pointer to location that receives index of the second index
 *            that yields the maximum likelihood (ignored if null)
 * @return the total likelihood
 */
double c_unmix_total( const double *restrict v_
                    , const double *restrict p7_, unsigned n7
                    , const double *restrict p5_, unsigned n5_
                    , unsigned *pi7, unsigned *pi5 )
{
    unsigned n5 = n5_ * 32 ;
    const double *restrict  v = v_ ; // __builtin_assume_aligned(  v_, 16 ) ;
    const double *restrict p5 = p5_ ; // __builtin_assume_aligned( p5_, 16 ) ;
    const double *restrict p7 = p7_ ; // __builtin_assume_aligned( p7_, 16 ) ;

    double acc = 0 ;
    double max = 0 ;
    unsigned mi7 = 0 ;
    unsigned mi5 = 0 ;
    for( unsigned i = 0, k = 0 ; i != n7 ; ++i, k += n5 ) {
        double p7i = p7[i] ;
        for( unsigned j = 0 ; j != n5 ; ++j ) {
            double p = v[k+j] * p7i * p5[j] ;
            acc += p ;
            if( p > max ) {
                max = p ;
                mi7 = i ;
                mi5 = j ;
            }
        }
    }
    if( pi7 ) *pi7 = mi7 ;
    if( pi5 ) *pi5 = mi5 ;
    return acc ;
}

/** Computes posterior probabilities for each pair of indices, given
 * matching probabilities for each and a matrix of prior probabilities,
 * the total likelihood and the index pair that yields the maximum
 * likelihood.  The posterior is added to an accumulator, and a quality
 * score is returned.  (The length of p5_ must be a multiple of 32 to
 * make vectorization easier.)
 *
 * @param w_ matrix of dimension (n7,n5_*32) to which the posterior is added (ignored if null)
 * @param v_ matrix of dimension (n7,n5_*32) containing the prior
 * @param p7_ vector of length n7 containing matching probabilities for the first index
 * @param n7 length of vector p7_
 * @param p5_ vector of length (n5_*32) containing matching probabilities for the second index
 * @param n5 length of vector p5_ divided by 32
 * @param total the total likelihood
 * @param mi7 index of the first index that yields the maximum likelihood
 * @param mi5 index of the second index that yields the maximum likelihood
 * @return the posterior probability for any other than the most likely assignment
 */
double c_unmix_qual( double *restrict w_
                   , const double *restrict v_
                   , const double *restrict p7_, unsigned n7
                   , const double *restrict p5_, unsigned n5_
                   , double total, unsigned mi7, unsigned mi5 )
{
    unsigned n5 = n5_ * 32 ;
    double        *restrict w = w_ ; // __builtin_assume_aligned(  w_, 16 ) ;
    const double *restrict  v = v_ ; // __builtin_assume_aligned(  v_, 16 ) ;
    const double *restrict p5 = p5_ ; // __builtin_assume_aligned( p5_, 16 ) ;
    const double *restrict p7 = p7_ ; // __builtin_assume_aligned( p7_, 16 ) ;
    double acc = 0 ;

    total = 1.0 / total ;
    for( unsigned i = 0, k = 0 ; i != n7 ; ++i ) {
        double p7i = p7[i] ;
        for( unsigned j = 0 ; j != n5 ; ++j, ++k ) {
            double p = total * v[k] * p7i * p5[j] ;
            if( w ) w[k] += p ;
            if( mi7 != i || mi5 != j ) acc += p ;
        }
    }
    return acc ;
}

