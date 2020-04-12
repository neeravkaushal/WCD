#include<iostream>
#include<vector>
#include<cmath>
#include<cstdlib>
#include<ctime>
using namespace std;


float get_random(float min, float max)
	{
	float number = min + (float) rand() / (float) (RAND_MAX/(max-min));
	return number;
	}


float get_distance(vector<float> v1, vector<float> v2)
	{
	float distance=0.0;
	for (int i=0; i<v1.size(); i++)
		distance += ( pow(v2[i]-v1[i],2.0) );
	distance = pow(distance,0.5);
	return distance;
	}

int main()
{
    const long all_cases=400000000;
	const long ntimes=300000;
    const long pi=3.141597;
    long cases=0, counts=0;
    const double c  = 299792458.0/1.33;
    const double H  = 4.5;
    const float  R  = 3.65;
    const double v  = 1.33 * c;
	srand(time(0));
    vector<float> c1 = {0, 0};
    vector<float> c2 = {1.85*cos(2*pi/3) , 1.85*sin(2*pi/3)};
    vector<float> c3 = {1.85*cos(4*pi/3) , 1.85*sin(4*pi/3)};
    vector<float> c4 = {1.85*cos(0     ) , 1.85*sin(0     )};
    vector<double> times;

    for (long i=0; i<=ntimes; i++)
        times.push_back(1e-11 + i*(1e-6-1e-11)/(float)ntimes);

	for (int i=1; i<=all_cases; i++)
		{
		float h = get_random(0.0, 4*H);
		float theta2 = get_random(pi/2, pi);
		float phi = get_random(0.0, 2*pi);
		float theta = pi - theta2;

		float x1,y1;
		if (h<=H)
			{
			x1 = get_random(-R,R);
			y1 = pow(R*R-x1*x1,0.5);
			}
		else
			{
			h = H;
			x1 = get_random(-R,R);
			float ymax = pow(R*R-x1*x1,0.5);
			y1 = get_random(-ymax, ymax);
			}

		float x2 = x1 + h*tan(theta)*cos(phi);
		float y2 = x2 + h*tan(theta)*sin(phi);
		vector<float> cp = {x2,y2};

		if (x2*x2 + y2*y2 <= R*R)
			{
			cases+=1;
			float d1 = get_distance(cp,c1);
			float d2 = get_distance(cp,c2);
			float d3 = get_distance(cp,c3);
			float d4 = get_distance(cp,c4);
			vector<float> Ls = {d1,d2,d3,d4};

			for (int L : Ls)
				{
				float num_zc = -pow( (-pow(c,4.0)*L*L*cos(theta)*cos(theta) + pow(c*L*v*cos(theta),2.0)),0.5)
							   -c*c*L*sin(theta) + L*v*v*sin(theta);
				float den = (c*c-v*v);
				float zc  = (num_zc/den)*cos(theta);
				if ((zc>0.0) && (zc<H))
					{
					counts += 1;
					break;
					}
				}
			}

	}
	cout<<"Total Particles  : "<<all_cases<<endl;
	cout<<"Particles entered: "<<cases<<endl;
	cout<<"Particles RIDed  : "<<counts<<endl;
	cout<<"Fraction         : "<<(float)counts/(float)cases<<endl;
    return 0;
}
