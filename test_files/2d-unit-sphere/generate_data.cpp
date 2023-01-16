/* 
Author: Vinaik Chhetri
*/
#include <iostream>
#include <vector>
#include <math.h>
#include "../../Box.h"
#include "../../Interval.h"
#include "../../Point.h"
#include "../../derivatives.h"
#include <random>



#include <fstream>
#include <cstdlib>

using namespace std;


//print vector of points
inline
std::ostream&
operator<< (std::ostream& os, vector < Point >& v)
{
    cout << "vector of points" << "\n";
    for (int i = 0; i < v.size(); i++)
    {

        os << v[i] << " ";



    } os << "\n";

    return os;
}



//2d unit sphere

int main(){

const int nrolls=20;  // number of experiments

std::default_random_engine generator;
std::normal_distribution<double> distribution(0.0,1.0);

vector<double> points;
double deno;

ofstream input ("2d-unit-circle.txt" , ios::app);

for (int i=0; i<nrolls; i++) {

double number1 = distribution(generator);
double number2 = distribution(generator);

if ((number1 == 0.0) && (number2 == 0.0) ) i=i-1;
deno = (number1*number1) + (number2*number2);
deno = sqrt(deno);
deno = 1/deno;
number1 = number1*deno;
number2 = number2*deno;

input << number1 << " ";
input << number2<< "\n";

points.push_back(number1);
points.push_back(number2);
}

input.close();

vector<Point> data;
vector<double> point;

for(int i =0; i<points.size(); i+=2){
   point.push_back(points[i]);
   point.push_back(points[i+1]);
   data.push_back(Point(point));
   point.erase(point.begin(),point.begin()+point.size());
}
cout<<" data " << data;
cout<<" data.size "<<data.size();

return 0;
}

