#include "stdio.h"
#include "math.h"
#include "stdlib.h"

/////////////////////////////////////矩阵转置////////////////////////////////////
void matrix_transpose(double *A,double *C,int hA,int wA)
//矩阵转置 原矩阵A=hA*wA;转置后矩阵为C=wA*hA;
{
	int i,j;
	for (i=0;i<hA;i++)
		for(j=0;j<wA;j++)
		{
			C[j*hA+i]=A[i*wA+j];
		}
}
/////////////////////////////////////矩阵相乘///////////////////////////////////
void matrix_multiply(double *A,double *B,double *C ,int hA,int wA,int wB)
//矩阵相乘 C=A*B，其中A矩阵是hA*wA，B矩阵是wA*wB 的，则C矩阵是hA*wB;
{
	int i,j,k;
	for (i=0;i<hA;i++)
		for(j=0;j<wB;j++)
		{
			double sum=0;
			for (k=0;k<wA;k++)
			{
				double a=A[i*wA+k];
				double b=B[k*wB+j];
				sum+=a*b;
			}
			C[i*wB+j]=sum;		
		}
}
/////////////////////////////////矩阵求逆////////////////////////////////
//利用的是AX=B，X=A'B,这里B=E；进行初等行变换求解，把左边化为单位阵，
//右边就是A矩阵的逆矩阵；

void swap(double *a,int i,int line,int n) 
// exchange line//交换行位置，i控制行号，line也是行号，//n是矩阵列数
{
	int j;
	double temp;
	for(j=0;j<n;j++)
	{
		temp=a[i*n+j];
		a[i*n+j]=a[line*n+j];//将第line行元素换到第i行
		a[line*n+j]=temp;
	}
}
void unitmatrix(double *q,int n)   //形成单位矩阵，矩阵q是一个单位阵
{
	int i,j;
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			if(i==j)
			{
				q[i*n+j]=1;
			}
			else
				q[i*n+j]=0;
		}
	}
}

void calculate(double *p,double *q,int n) //形成上三角阵
{
	int i,j,k,m,line;
	double max,temp,mmul;
	for(i=0;i<n;i++)
	{
		max=fabs(p[i*n+i]);
		temp=p[i*n+i];
		line=i;
		for(j=i+1;j<n;j++)       
			//选出每一列中最大值，并且用line记录行号,
           //用temp记录该最大值，用j来控制行号
		{
			if(fabs(p[j*n+i])>max)
			{
				max=fabs(p[j*n+i]); 
				temp=p[j*n+i]; 
				line=j;
			}
		}  
		if(max<=1e-10)  
		{
			printf("no inverse array\n");
			return;
		}
		if(line!=i)
		{ 
			swap(p,i,line,n);//将每一列中最大行换到i行
			swap(q,i,line,n);
		}
		for(k=0;k<n;k++)
		{
			p[i*n+k]/=temp;//将i行的每个数都除以i行i列的值，将i行i列化为1
			q[i*n+k]/=temp;
		}
		for(k=i+1;k<n;k++)//将i列i行下面其它行的值都化为0
		{
			mmul=p[k*n+i];
			for(m=0;m<n;m++)
			{
				p[k*n+m]-=p[i*n+m]*mmul;//每一行都减去行首元素*i行对应列的值
				q[k*n+m]-=q[i*n+m]*mmul;
			}
		}
	} 
	
}

void backcalculate(double*p,double*q,int n)//形成单位矩阵  
{
	int i,j,k;
	double mmul;
	for(i=n-1;i>0;i--)
	{
		for(j=i-1;j>=0;j--)//从下往上每一行进行计算
		{
			mmul=p[j*n+i];
			p[j*n+i]-=p[i*n+i]*mmul;//化为0
			for(k=0;k<n;k++)//q每列也如此处理
			{
				q[j*n+k]-=q[i*n+k]*mmul;//每一行都减去行末元素*i行对应列的值
			}
		}
	}
	
}


void matrix_invert(double *A,double *C,int n)//矩阵求逆，A为要求的矩阵，C为求得的逆矩阵
{
	unitmatrix(C,n);//将C初始化为单位阵
    calculate(A,C,n);//将A化为了上三角阵
    backcalculate(A,C,n);//将A化为了单位阵，这时C即为所求
}
