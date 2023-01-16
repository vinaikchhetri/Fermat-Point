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

   
inline int all0(Point &p){
    int count=0;
    for(int i =0;i<p.point.size();i++){
        if(p.point[i]==0) count++;
    }
    //cout<<"\n"<<"counter"<<count;
    return count;
}
   
   
   

//100d unit sphere

int main(){
  

 //4.100-dim unit sphere
  
const int nrolls=20;  // number of experiments

std::default_random_engine generator;
std::normal_distribution<double> distribution(0.0,1.0);

vector<Point> data;
Point p(100);
ofstream input ("100d-unit-sphere.txt" , ios::app);
double number;
for (int i=0; i<nrolls; i++) {
  for(int j =0;j<100;j++){
      number = distribution(generator);
      p.point[j]=number;
      if(j==99){
           if (all0(p)==p.point.size()) j--;
      }

  }
  number = normP(p); 
  p = p/number; data.push_back(p); 
}

for(int i = 0; i < data.size(); i++){
for(int j = 0; j < 100; j++){
    input << data[i].point[j]<< " ";
}
input << "\n";
}

input.close();


cout<<" data " << data;
cout<<" data.size "<<data.size();

return 0;
}


