void MatrixCopy_( int m, int n, double *a, double *u )
{
    int max;
    
    max = m * n ;
    
    for ( int i = 0 ; i < max ; i++ )
	    u[ i ] = a[ i ];

}
