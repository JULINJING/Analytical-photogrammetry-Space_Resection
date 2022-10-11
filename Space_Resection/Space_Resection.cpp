// Space_Resection.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//单像空间后方交会

#include <iostream>
#include <fstream>
#include <iomanip>//输入输出
#include "matrix.h"//矩阵计算
using namespace std;

int N = 4;//控制点总数
#define PI acos(-1)//计算PI

/*
数据文件data.txt内容如下
影像坐标xy和地面坐标XYZ
-86.15       -68.99    36589.41       25273.32       2195.17
-53.40       82.21     37631.08       31324.51       728.69
-14.78       -76.63    39100.97       24934.98       2386.50
10.46        64.43     40426.54       30319.81       757.31
*/

//读取数据文件(文件名，影像x，影像y，地面X，地面Y，地面Z)
void get_ini_data(const char* filename, double* x, double* y, double* X, double* Y, double* Z);

int main()
{

	//已知内方位元素
	double f = 153.24;//摄影机主距(mm)
	double x0 = 0;
	double y0 = 0;//内方位元素xy

	double M = 50000;//摄影比例尺

	//未知外方位元素
	double Xs, Ys, Zs, phi, omega, kappa;//6个外方位元素

	//旋转矩阵R
	double a1, a2, a3, b1, b2, b3, c1, c2, c3;//9个元素

	//声明动态数组
	double* A = new double[2 * N * 6];//存放外方位元素系数A的矩阵
	double* AT = new double[2 * N * 6];//存储外方位元素系数转置矩阵AT
	double delta[6];//6个外方位元素改正值矩阵

	//存储坐标
	double* X = new double[N]; double* Y = new double[N]; double* Z = new double[N];//存储地面坐标
	double* x = new double[N]; double* y = new double[N];//存储影像坐标
	double* xj = new double[N]; double* yj = new double[N];//存储由共线方程解算的近似像点坐标

	//中间过程运算矩阵
	double* temp1 = new double[36];
	double* temp2 = new double[6];
	double* temp3 = new double[36];

	//平差声明
	double Zup;//普遍分母
	double* L = new double[2 * N];
	double Q[6];//存储6个权倒数
	double* V = new double[2 * N];//存储像点观测值改正数
	double Vsum = 0, m0;//单位权中误差
	double m[6];//存储六个外方位元素的中误差

	//获得已知数据
	get_ini_data("data.txt", x, y, X, Y, Z);

	//确定观测数据初值
	f /= 1000;
	for (int i = 0; i < N; i++)
	{
		x[i] /= 1000;
		y[i] /= 1000;
	}//单位转换mm到m
	phi = omega = kappa = 0.0;

	double sum_Xs = 0, sum_Ys = 0, sum_Zs = 0;
	for (int i = 0; i < N; i++) {
		sum_Xs += X[i];
		sum_Ys += Y[i];
		sum_Zs += Z[i];
	}
	Xs = sum_Xs / N;
	Ys = sum_Ys / N;
	Zs = M * f + sum_Zs / N;

	//迭代
	//6.检查计算是否收敛
	int iteration = 0;
	while (fabs(delta[0]) >= 0.001 || fabs(delta[1]) >= 0.001 || fabs(delta[2]) >= 0.001 
	|| fabs(delta[3]) >= 1.0e-6 || fabs(delta[4]) >= 1.0e-6 || fabs(delta[5]) >= 1.0e-6)//阈值10-3和10-6
	{
		//1.计算旋转矩阵R
		a1 = cos(phi) * cos(kappa) - sin(phi) * sin(omega) * sin(kappa);
		a2 = -cos(phi) * sin(kappa) - sin(phi) * sin(omega) * cos(kappa);
		a3 = -sin(phi) * cos(omega);
		b1 = cos(omega) * sin(kappa);
		b2 = cos(omega) * cos(kappa);
		b3 = -sin(omega);
		c1 = sin(phi) * cos(kappa) + cos(phi) * sin(omega) * sin(kappa);
		c2 = -sin(phi) * sin(kappa) + cos(phi) * sin(omega) * cos(kappa);
		c3 = cos(phi) * cos(omega);

		for (int i = 0; i < N; i++)
		{
			//2.共线方程求近似像点坐标
			xj[i] = x0 - f * (a1 * (X[i] - Xs) + b1 * (Y[i] - Ys) + c1 * (Z[i] - Zs)) / (a3 * (X[i] - Xs) + b3 * (Y[i] - Ys) + c3 * (Z[i] - Zs));
			yj[i] = y0 - f * (a2 * (X[i] - Xs) + b2 * (Y[i] - Ys) + c2 * (Z[i] - Zs)) / (a3 * (X[i] - Xs) + b3 * (Y[i] - Ys) + c3 * (Z[i] - Zs));

			//3.计算组成误差方程式
			//计算Zup，即普遍分母
			Zup = a3 * (X[i] - Xs) + b3 * (Y[i] - Ys) + c3 * (Z[i] - Zs);

			//逐点计算误差方程的系数矩阵A
			A[12 * i] = (a1 * f + a3 * (x[i] - x0)) / Zup;
			A[12 * i + 1] = (b1 * f + b3 * (x[i] - x0)) / Zup;
			A[12 * i + 2] = (c1 * f + c3 * (x[i] - x0)) / Zup;
			A[12 * i + 3] = (y[i] - y0) * sin(omega) - cos(omega) * ((x[i] - x0) * ((x[i] - x0) * cos(kappa) - (y[i] - y0) * sin(kappa)) / f + f * cos(kappa));
			A[12 * i + 4] = -f * sin(kappa) - (x[i] - x0) * ((x[i] - x0) * sin(kappa) + (y[i] - y0) * cos(kappa)) / f;
			A[12 * i + 5] = (y[i] - y0);
			A[12 * i + 6] = (a2 * f + a3 * (y[i] - y0)) / Zup;
			A[12 * i + 7] = (b2 * f + b3 * (y[i] - y0)) / Zup;
			A[12 * i + 8] = (c2 * f + c3 * (y[i] - y0)) / Zup;
			A[12 * i + 9] = -(x[i] - x0) * sin(omega) - cos(omega) * ((y[i] - y0) * ((x[i] - x0) * cos(kappa) - (y[i] - y0) * sin(kappa)) / f - f * sin(kappa));
			A[12 * i + 10] = -f * cos(kappa) - (y[i] - y0) * (x[i] * sin(kappa) + (y[i] - y0) * cos(kappa)) / f;
			A[12 * i + 11] = -(x[i] - x0);
			//计算常数项矩阵L
			L[2 * i] = x[i] - xj[i];
			L[2 * i + 1] = y[i] - yj[i];
		}

		//4.计算组成法方程式
		//计算法方程系数矩阵ATA和常数项矩阵ATL
		matrix_transpose(A, AT, 2 * N, 6);//计算AT
		matrix_multiply(AT, A, temp1, 6, 2 * N, 6);//计算AT*A
		matrix_multiply(AT, L, temp2, 6, 2 * N, 1);//计算AT*L

		//5.解求外方位元素
		matrix_invert(temp1, temp3, 6);//计算(AT*A)逆存入temp3
		matrix_multiply(temp3, temp2, delta, 6, 6, 1);//计算(AT*A)-1*(AT*L)求得外方位元素delta[6]	
		//外方位元素改正后的近似值
		Xs += delta[0];
		Ys += delta[1];
		Zs += delta[2];
		phi += delta[3];
		omega += delta[4];
		kappa += delta[5];

		iteration += 1;
	}//迭代结束

	//精度评定
	for (int i = 0; i < 6; i++)
	{
		int j = i;
		Q[i] = temp3[i * 6 + j];//将ATA逆的对角线元素附给Q[i]
	}
	for (int i = 0; i < N; i++)//V=AX-L
	{
		V[2 * i] = A[12 * i] * delta[0] + A[12 * i + 1] * delta[1] + A[12 * i + 2] * delta[2] + A[12 * i + 3] * delta[3] + A[12 * i + 4] * delta[4] + A[12 * i + 5] * delta[5] - L[2 * i];
		V[2 * i + 1] = A[12 * i + 6] * delta[0] + A[12 * i + 7] * delta[1] + A[12 * i + 8] * delta[2] + A[12 * i + 9] * delta[3] + A[12 * i + 10] * delta[4] + A[12 * i + 11] * delta[5] - L[2 * i + 1];
		Vsum += V[2 * i] * V[2 * i] + V[2 * i + 1] * V[2 * i + 1];
	}
	Vsum /= 2 * N - 6;
	m0 = sqrt(Vsum);//中误差m0

	//输入结果入txt文件
	ofstream fileout("计算结果.txt", ios::out);
	if (!fileout)
	{
		cout << "存储文件失败!" << endl;
		exit(1);
	}
	fileout << "1.影像坐标和地面坐标的信息如下所示:" << endl;
	fileout << endl;
	for (int i = 0; i < N; i++) {
		fileout << "第" << i + 1 << "个坐标对为：x=" << x[i] << "米" << " y=" << y[i] << "米" << " X=" << X[i] << "米" << " Y=" << Y[i] << "米" << " Z=" << Z[i] << "米" << endl;
	}
	fileout << endl;
	fileout << "2.在解算中一共经过了" << iteration << "次迭代计算,解算得到的六个外方位元素的结果为：" << endl;
	fileout << endl;
	//fileout << setiosflags(ios::fixed) << setprecision(3);//以小数位固定为3位输出
	fileout << "3个外方位线元素为：" << endl;
	fileout << "Xs=" << Xs << "米" << endl;
	fileout << "Ys=" << Ys << "米" << endl;
	fileout << "Zs=" << Zs << "米" << endl;
	fileout << endl;
	fileout << "3个外方位角元素为：" << endl;
	fileout << "phi=" << phi * 180 * 3600 / PI << "秒" << " or " << phi << "rad" << endl;
	fileout << "omega=" << omega * 180 * 3600 / PI << "秒" << " or " << omega << "rad" << endl;
	fileout << "kappa=" << kappa * 180 * 3600 / PI << "秒" << " or " << kappa << "rad" << endl;
	fileout << endl;
	//fileout << setiosflags(ios::fixed) << setprecision(5);//以小数位固定为5位输出
	fileout << "3.由解算得到的外方位元素计算得到的旋转矩阵R为：" << endl;
	fileout << endl;
	fileout << setw(14) << a1 << setw(14) << a2 << setw(14) << a3 << endl;
	fileout << setw(14) << b1 << setw(14) << b2 << setw(14) << b3 << endl;
	fileout << setw(14) << c1 << setw(14) << c2 << setw(14) << c3 << endl;
	fileout << endl;
	fileout << "4.解算得到的单位权中误差和外方位元素精度为：" << endl;
	fileout << endl;
	fileout << "单位权中误差为：" << m0 << endl;
	fileout << "六个外方元素精度依次是：Xs,Ys,Zs,phi,omega,kappa" << endl;
	for (int i = 0; i < 6; i++)
	{
		m[i] = m0 * sqrt(Q[i]);
		if (i < 3)
		{
			fileout << m[i] << "米" << endl;
		}
		else
		{
			fileout << m[i] * 180 * 3600 / PI << "秒" << " or " << m[i] << "rad" << endl;
		}
	}
	fileout.close();

	//运行框显示结果
	cout << "**********这里是单像空间后方交会程序的结果显示界面**********" << endl;
	cout << endl;
	cout << "1.影像坐标和地面坐标的信息如下所示:" << endl;
	cout << endl;
	for (int i = 0; i < N; i++) {
		cout << "第" << i + 1 << "个坐标对为：x=" << x[i] << "米" << " y=" << y[i] << "米" << " X=" << X[i] << "米" << " Y=" << Y[i] << "米" << " Z=" << Z[i] << "米" << endl;
	}
	cout << endl;
	//cout << setiosflags(ios::fixed) << setprecision(3);//以小数位固定为3位输出
	cout << "2.在解算中一共经过了" << iteration << "次迭代计算,解算得到的六个外方位元素的结果为：" << endl;
	cout << endl;
	cout << "3个外方位线元素为：" << endl;
	cout << "Xs=" << Xs << "米" << endl;
	cout << "Ys=" << Ys << "米" << endl;
	cout << "Zs=" << Zs << "米" << endl;
	cout << endl;
	cout << "3个外方位角元素为：" << endl;
	cout << "phi=" << phi * 180 * 3600 / PI << "秒" << " or " << phi << "rad" << endl;
	cout << "omega=" << omega * 180 * 3600 / PI << "秒" << " or " << omega << "rad" << endl;
	cout << "kappa=" << kappa * 180 * 3600 / PI << "秒" << " or " << kappa << "rad" << endl;
	cout << endl;
	//cout << setiosflags(ios::fixed) << setprecision(5);//以小数位固定为5位输出
	cout << "3.由解算得到的外方位元素计算得到的旋转矩阵R为：" << endl;
	cout << endl;
	cout << setw(14) << a1 << setw(14) << a2 << setw(14) << a3 << endl;
	cout << setw(14) << b1 << setw(14) << b2 << setw(14) << b3 << endl;
	cout << setw(14) << c1 << setw(14) << c2 << setw(14) << c3 << endl;
	cout << endl;
	cout << "4.解算得到的单位权中误差和外方位元素精度为：" << endl;
	cout << endl;
	cout << "单位权中误差为：" << m0 << endl;
	cout << "六个外方元素精度依次是：Xs,Ys,Zs,phi,omega,kappa" << endl;
 	for (int i=0;i<6;i++)
 	{
		m[i] = m0 * sqrt(Q[i]);
 		if (i<3)
 		{			
			cout << m[i] << "米" << endl;
 		}
 		else 
 		{
			cout << m[i] * 180 * 3600 / PI << "秒" << " or " << m[i] << "rad" << endl;
 		}		
 	}
	cout << endl << "**********结果显示完毕，感谢使用！**********" << endl;

	//删除数组申请的内存空间
	delete[]A;
	delete[] X, Y, Z;
	delete[] x, y;
	delete[] xj, yj;
	delete[] AT;
	delete[] temp1, temp2, temp3;
	delete[] V;
	return 0;
}

//函数具体代码
void get_ini_data(const char* filename, double* x, double* y, double* X, double* Y, double* Z) {
	ifstream filein(filename, ios::in);
	// 以输入方式打开一个文件
	if (!filein)
	{
		cout << "打开文件失败!" << endl;
		exit(1);
	}
	for (int i = 0; i < N; i++)
		filein >> x[i] >> y[i] >> X[i] >> Y[i] >> Z[i];//读入每行点对的像方和物方坐标
	filein.close();
}

// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门使用技巧: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
