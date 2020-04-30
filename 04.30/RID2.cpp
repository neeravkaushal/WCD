#include<iostream>
#include<cmath>
#include<cstdlib>
#include<ctime>
using namespace std;

double getRandFloat(double min, double max)
	{
	return min + (double) rand() / (double) (RAND_MAX/(max-min));
	}

int getRandInt(long min, long max)
	{
	return min + ( (long)   rand() / (long)   (RAND_MAX/(max-min+1)) );
	}


double getDistance(double a1[], double a2[], long n)
	{
	double distance=0.0;
	for (unsigned int i=0; i<n; i++)
		distance += pow(a2[i]-a1[i],2.0);
	return pow(distance,0.5);
	}


double getAngle(double b[], double a[], double c[], long n)
	{
	double ab[n], ac[n];
	double origin[3]={0.0,0.0,0.0};
	double dotproduct=0.0;
	for (unsigned int i=0; i<n; i++){
		ab[i] = b[i]-a[i];
		ac[i] = c[i]-a[i]; }
	for (unsigned int i=0; i<n; i++)
		dotproduct += (ab[i]*ac[i]);
	return  acos (dotproduct/(getDistance(ab,origin,n) * getDistance(ac,origin,n)));
	}


void setCoordinates(long choice, double point[])
	{
	double x,y,z,Y;
	double H=4.5, R=3.65;
	x = -R + ( (double)   rand() / (double)   (RAND_MAX/(R+R)) );
	if (choice==1){
		Y = pow(R*R-x*x,0.5);
		y = -Y + ( (double)   rand() / (double)   (RAND_MAX/(Y+Y)) );
		z = H; }
	else if (choice==2){
		z = 0 + ( (double)   rand() / (double)   (RAND_MAX/(H)) );
		y = pow(R*R-x*x,0.5);
		int dir = 1 + ( (long)   rand() / (long)   (RAND_MAX/(2-1+1)) );
		if (dir==2)
			y = -y; }
	else{
		Y = pow(R*R-x*x,0.5);
		y = -Y + ( (double)   rand() / (double)   (RAND_MAX/(Y+Y)) );
		z = 0;}
	point[0]=x;
	point[1]=y;
	point[2]=z;
	}


int main(int argc , char* argv[])
	{
	srand(time(0));
   const long   all_cases=atol(argv[1]);
   const long   pi=3.141597;
   const double c=299792458.0/1.33, v=1.33*c, H=4.5;// R=3.65;
	const double c1[3] = {0.0, 0.0, 0.0};
   const double c2[3] = {1.85*cos(2*pi/3) , 1.85*sin(2*pi/3), 0.0};
   const double c3[3] = {1.85*cos(4*pi/3) , 1.85*sin(4*pi/3), 0.0};
   const double c4[4] = {1.85*cos(0)      , 1.85*sin(0)     , 0.0};
	long   counts=0, iter=0;
	long   choice_A, choice_B, flag;
   double A[3], B[3], C[3];
	double theta, alpha, nAB, L, xc, zc;

	while (iter<all_cases)
		{
		choice_A = getRandInt(1,2);
		choice_B = getRandInt(2,3);
		setCoordinates(choice_A, A);
		setCoordinates(choice_B, B);
		while (B[2] > A[2])
			{
			choice_A = getRandInt(1,2);
			choice_B = getRandInt(2,3);
			setCoordinates(choice_A, A);
			setCoordinates(choice_B, B);
			}
		C[0]=A[0];
		C[1]=A[1];
		C[2]=0.0;
		theta = getAngle(B,A,C,3);
		nAB = getDistance(A,B,3);
		double detectors[1][4][3]  = {{{c1[0],c1[1],c1[2]}, {c2[0],c2[1],c2[2]}, {c3[0],c3[1],c3[2]}, {c4[0],c4[1],c4[2]}}};
		flag = 0;
		for (int index=0; index<4; index++)
			{
			alpha = getAngle(B,A,detectors[0][index],3);
			L     = getDistance(A,detectors[0][index],3);
			xc    = L*cos(alpha) - (c*L*sin(alpha)) / (pow(v*v-c*c,0.5));
			zc    = (nAB-xc)*cos(theta) + B[2];
			if ((zc>0.0) && (zc<H)) {
				counts+=1;
				flag = 1; }
			if (flag==1)
				break;
			}
		iter+=1;
		}

	cout<<"Total Particles  : "<<all_cases<<endl;
   cout<<"Particles RIDed  : "<<counts<<endl;
   cout<<"Fraction         : "<<(double)counts/(double)all_cases<<endl;
   return 0;
}
