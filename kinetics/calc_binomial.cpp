#define _USE_MATH_DEFINES 
int calc_binomial(int n, int k)
{
    int dp[20][20];
    
	if (k == 0 || k == n)
        return 1;
    
    else
	    return dp[n][k] = calc_binomial(n - 1, k - 1) + calc_binomial(n - 1, k);
    
}


/*
int main()
{
	int a = calc_binomial( 5, 3);
	
	cout<< a << endl;
	
}

*/
