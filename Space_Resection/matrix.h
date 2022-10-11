#include "stdio.h"
#include "math.h"
#include "stdlib.h"

/////////////////////////////////////����ת��////////////////////////////////////
void matrix_transpose(double *A,double *C,int hA,int wA)
//����ת�� ԭ����A=hA*wA;ת�ú����ΪC=wA*hA;
{
	int i,j;
	for (i=0;i<hA;i++)
		for(j=0;j<wA;j++)
		{
			C[j*hA+i]=A[i*wA+j];
		}
}
/////////////////////////////////////�������///////////////////////////////////
void matrix_multiply(double *A,double *B,double *C ,int hA,int wA,int wB)
//������� C=A*B������A������hA*wA��B������wA*wB �ģ���C������hA*wB;
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
/////////////////////////////////��������////////////////////////////////
//���õ���AX=B��X=A'B,����B=E�����г����б任��⣬����߻�Ϊ��λ��
//�ұ߾���A����������

void swap(double *a,int i,int line,int n) 
// exchange line//������λ�ã�i�����кţ�lineҲ���кţ�//n�Ǿ�������
{
	int j;
	double temp;
	for(j=0;j<n;j++)
	{
		temp=a[i*n+j];
		a[i*n+j]=a[line*n+j];//����line��Ԫ�ػ�����i��
		a[line*n+j]=temp;
	}
}
void unitmatrix(double *q,int n)   //�γɵ�λ���󣬾���q��һ����λ��
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

void calculate(double *p,double *q,int n) //�γ���������
{
	int i,j,k,m,line;
	double max,temp,mmul;
	for(i=0;i<n;i++)
	{
		max=fabs(p[i*n+i]);
		temp=p[i*n+i];
		line=i;
		for(j=i+1;j<n;j++)       
			//ѡ��ÿһ�������ֵ��������line��¼�к�,
           //��temp��¼�����ֵ����j�������к�
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
			swap(p,i,line,n);//��ÿһ��������л���i��
			swap(q,i,line,n);
		}
		for(k=0;k<n;k++)
		{
			p[i*n+k]/=temp;//��i�е�ÿ����������i��i�е�ֵ����i��i�л�Ϊ1
			q[i*n+k]/=temp;
		}
		for(k=i+1;k<n;k++)//��i��i�����������е�ֵ����Ϊ0
		{
			mmul=p[k*n+i];
			for(m=0;m<n;m++)
			{
				p[k*n+m]-=p[i*n+m]*mmul;//ÿһ�ж���ȥ����Ԫ��*i�ж�Ӧ�е�ֵ
				q[k*n+m]-=q[i*n+m]*mmul;
			}
		}
	} 
	
}

void backcalculate(double*p,double*q,int n)//�γɵ�λ����  
{
	int i,j,k;
	double mmul;
	for(i=n-1;i>0;i--)
	{
		for(j=i-1;j>=0;j--)//��������ÿһ�н��м���
		{
			mmul=p[j*n+i];
			p[j*n+i]-=p[i*n+i]*mmul;//��Ϊ0
			for(k=0;k<n;k++)//qÿ��Ҳ��˴���
			{
				q[j*n+k]-=q[i*n+k]*mmul;//ÿһ�ж���ȥ��ĩԪ��*i�ж�Ӧ�е�ֵ
			}
		}
	}
	
}


void matrix_invert(double *A,double *C,int n)//�������棬AΪҪ��ľ���CΪ��õ������
{
	unitmatrix(C,n);//��C��ʼ��Ϊ��λ��
    calculate(A,C,n);//��A��Ϊ����������
    backcalculate(A,C,n);//��A��Ϊ�˵�λ����ʱC��Ϊ����
}
