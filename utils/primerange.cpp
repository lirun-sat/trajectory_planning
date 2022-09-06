int primerange(int start, int end, int* result, int count)
{
    // 求 [start, end] 之间的 素数， 返回值 count 表示 素数总个数 
    int* result_temp = new int[1000];
    int result_count = 0;
    bool prime;  //定义bool变量
            
    
    for(int i = start; i < (end+1); i++)
    {
        prime = true;    //先令 prime 为真
        
        if((i == 0) || (i == 1))
        {
            // cout << "i=0 or i=1" <<endl;
            continue;
        }
       
        for(int j = 2; j < i; j++)    //对 2 到 m 进行循环
        {
            if(i % j == 0)    //若 i 整除 j 为 0，令 prime 为假，循环终止
            {
                prime = false;
                break;
            }
        }
        if(prime)    //若 prime 为真，输出 n
        {
            result_temp[result_count] = i;
            result_count += 1;
        }
        
    }
    for (int i = 0; i < result_count; i++)
    {
        result[i] = result_temp[i];
    }
        
    count = result_count;
    
    delete[] result_temp;
    
    return count;
    
}


/*

int main()
{
    int* result = new int[300];
    int count;
    
    count = primerange(0, 200, result, count);
    
    for(int i = 0; i< count; i++)
        cout<< result[i] << endl; 
        
}

*/
























