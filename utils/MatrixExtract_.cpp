void MatrixExtract_( int m, int n, int i_s, int i_f, int j_s, int j_f, double *A, double *ans )
{
    int row, col;

    row = i_f - (i_s-1);
    col = j_f - (j_s-1);

    for(int i = (i_s-1); i < i_f; i++)
    {
        for(int j = (j_s-1); j < j_f; j++ )
        {
            ans[col*(i-(i_s-1))+(j-(j_s-1))] = A[n*i+j];
        }
    }
}
