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


//100d unit cube

int main(){
  

std::random_device rd;  //Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
std::uniform_real_distribution<> dis(-1.0, 1.0);
vector<double> points;

ofstream input ("100d-unit-cube.txt" , ios::app);
double number;
for(int n = 0; n < 20; n++)
{
for (int d = 0; d < 100; d++) {
    // Use dis to transform the random unsigned int generated by gen into a 
    // double in [1, 2). Each call to dis(gen) generates a new random double
    number = dis(gen);
    points.push_back(number);
    input << number<< " ";
}
input <<"\n";
}
input.close();


vector<Point> data;
vector<double> point;

for(int i =0; i<20; i++){
   for(int j=0;j<100;j++){
     point.push_back(points[i*100+j]);
   }
   data.push_back(Point(point));
   point.erase(point.begin(),point.begin()+point.size()); 
}



cout<<" data " << data;
cout<<" data.size "<<data.size(); 

return 0;
}

