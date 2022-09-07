#include <stdio.h>
#include <stdlib.h>

#define fscanf_s fscanf


int readinput(const char *inputfile, double ***pts, int *out) 
{
    int npoints = 0;

    int idx = 0;

    FILE *fp;

    /* Open file. */
    /*
        C 库函数 :
        FILE *fopen(const char *filename, const char *mode) 
        使用给定的模式 mode 打开 filename 所指向的文件。
        下面是 fopen() 函数的声明:
        FILE *fopen(const char *filename, const char *mode)

        filename -- 字符串，表示要打开的文件名称。
        mode -- 字符串，表示文件的访问模式，可以是以下表格中的值：
        "r"	打开一个用于读取的文件。该文件必须存在。
        "w"	创建一个用于写入的空文件。如果文件名称与已存在的文件相同，则会删除已有文件的内容，文件被视为一个新的空文件。
        "a"	追加到一个文件。写操作向文件末尾追加数据。如果文件不存在，则创建文件。
        "r+"	打开一个用于更新的文件，可读取也可写入。该文件必须存在。
        "w+"	创建一个用于读写的空文件。
        "a+"	打开一个用于读取和追加的文件。

        该函数返回一个 FILE 指针。否则返回 NULL，且设置全局变量 errno 来标识错误。

    */
    
    fp = fopen(inputfile, "r");

    if ((fp) == NULL) 
    {
        fprintf(stdout, "ERROR: input file %s not found!\n", inputfile);
        fprintf(stdout, "  -> The file must be in the folder from which this " "program is launched\n\n");
        return 1;
    }

    /* Read number of input vertices. */
    if (fscanf_s(fp, "%d", &npoints) != 1)
        return 1;

    /* Allocate memory. */
    double **arr = (double **)malloc(npoints * sizeof(double *));
    for (int i = 0; i < npoints; i++)
        arr[i] = (double *)malloc(3 * sizeof(double));

    /* Read and store vertices' coordinates. */
    for (idx = 0; idx < npoints; idx++) 
    {
      if (fscanf_s(fp, "%lf %lf %lf\n", &arr[idx][0], &arr[idx][1], &arr[idx][2]) != 3)
          return 1;
    }

    fclose(fp);

    *pts = arr;
    *out = idx;

    return (0);
}


/*
int main() 
{
    int nvrtx1, nvrtx2;    

    char inputfileA[40] = "userP.dat", inputfileB[40] = "userQ.dat";

    double(**vrtx1) = NULL, (**vrtx2) = NULL;   

    if (readinput(inputfileA, &vrtx1, &nvrtx1))
        return (1);  

    if (readinput(inputfileB, &vrtx2, &nvrtx2))
        return (1);    

    for (int i = 0; i < nvrtx1; ++i) 
    {
        printf("%.2f ", vrtx1[i][0]);
        printf("%.2f ", vrtx1[i][1]);
        printf("%.2f\n", vrtx1[i][2]);
    }    

    for (int i = 0; i < nvrtx2; ++i) 
    {
        printf("%.2f ", vrtx2[i][0]);
        printf("%.2f ", vrtx2[i][1]);
        printf("%.2f\n", vrtx2[i][2]);
    }    

    // Free memory 
    for (int i = 0; i < nvrtx1; i++)
      free(vrtx1[i]);
    free(vrtx1);
    
    for (int i = 0; i < nvrtx2; i++)
      free(vrtx2[i]);
    free(vrtx2);  
       
    return (0);
}
*/