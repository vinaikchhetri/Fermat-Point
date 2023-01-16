/* 
Author: Vinaik Chhetri
*/
#include <iostream>
#include <vector>
#include <math.h>
#include "Box.h"
#include "Interval.h"
#include "Point.h"
#include "derivatives.h"
#include <random>

#include <fstream>
#include <cstdlib>




using namespace std;


//print vector of points
inline std::ostream& operator<< (std::ostream& os, vector < Point >& v)
{
    cout << "vector of points" << "\n";
    
    for (int i = 0; i < v.size(); i++){    
        os << v[i] << " ";
    } 
    os << "\n";
    return os;
}


//Compute the barycenter of the focis'.
inline Point foci_mid(vector<Point>& foci) {
    Point ret(foci[0].point.size());

    for (int i = 0; i < foci.size(); i++) {
        ret = ret + foci[i];
    }
    double den = (double)foci.size();
    ret = ret / den;
    return ret;
}

//Weizfield Iteration and return the next approximation for the fermat point.
inline Point weizfield(vector<Point>& foci, Point& p) {
    Point a(p.point.size());
    double norm;
    Point num = (p.point.size());
    double den = 0;

    Point diff(p.point.size());
    Point div(p.point.size());
    Point ret(p.point.size());
    
    for (int i = 0; i < foci.size(); i++) {
        a = foci[i];
        diff = p - a;
        norm = normP(diff);
        div = (a / norm);
        num = num + div;
        den = den + (1 / norm);
    }
    ret = num / den;
    return ret;
}

//Count the number of dimensions of the Jaccobian(Fermat point) where the entry<10^(-|number|)
inline void how_many(Point &p){
    int count=0;

    for(int i =0;i<p.point.size();i++){
        if(p.point[i]<0.0005) count++;
    }
    cout<<"\n"<<"counter"<<count;
}

//Count the number of dimensions of the Jaccobian(Fermat point) where the entry == 0
inline int all0(Point &p){
    int count=0;

    for(int i =0;i<p.point.size();i++){
        if(p.point[i]==0) count++;
    }
    //cout<<"\n"<<"counter"<<count;
    return count;
}

// Main loop of the algorithm to compute the best time to stop iterating the wiezfield iteration.
int main()
{  
  //Read the data and construct a vector of points.
  
  vector<double> points;
  ifstream inFile;
  clock_t begin = clock();

  //Change file name!!!!!!

  //inFile.open("test_files/test1.txt");
  //inFile.open("test_files/2d-unit-sphere/2d-unit-circle.txt");
  //inFile.open("test_files/3d-unit-sphere/3d-unit-circle.txt");
  //inFile.open("test_files/100d-unit-cube/100d-unit-cube.txt");
  inFile.open("test_files/100d-unit-sphere/100d-unit-sphere.txt");

  //Change dimension!!!!!
  int dim_para = 100;

  double x;
  //Read the points from txt file.
  while (inFile >> x) {
  points.push_back(x);
}

inFile.close();
 
vector<Point> data;
vector<double> point;
//Put the data inside a vector of points.
for(int i =0; i<(points.size()/dim_para); i++){
    for(int j = 0; j< dim_para; j++){
        point.push_back(points[i*dim_para+j]);
    }
   data.push_back(Point(point));
   point.erase(point.begin(),point.begin()+point.size());
}

cout<<data;


int dim = data[0].point.size();
Point p(dim);
//clock_t begin = clock();
p = foci_mid(data);

double eps = 0.0005;//0.005; //epsilon
double l = eps; //Width of the box.
bool flag = true; //Flag to stop the wiezfield iteration.
Box B(dim);
Box NB(dim);
Box B10(dim);
Box NB10(dim);
Point ret(dim);
double sc = 0.1;
int i = 0;
while (flag) {
    cout<<i<<"\n";
    B = center2box(p, l);
    NB = Newton(B, data);
   
    B10 = scaleBox(B, sc);
    NB10 = Newton(B10, data);
   
    if (inclusion(NB, B) == true) {
    ret = B.center();
    cout << "Fermat point : " <<ret<<"\n"<<"iteration number "<<i<<"\n";
    cout<<"width "<<l<<"\n";
    flag = false;
    }

    else if (exclusion(NB10, B10) == true) {
        //cout << " B10 " << B10 << " NB10 " << NB10;
        //cout << "2.1 " << l<<"\n";
        l = min(10 * l, eps);
        //cout << "2.2 " << l<< "\n";
    }

    else {
        //cout << "3.1 " << l<< "\n";
        l = l / 10;
        //cout << "3.2 " << l<< "\n";
    }

    p = weizfield(data, p);
    i++;
}
print_time(begin);
p = jaccobian(p,data);
cout <<"Jaccobian of the approximated fermat point : " <<p;

}
//print vector of points
/*
inline std::ostream& operator<< (std::ostream& os, vector < Point >& v)
{
    cout << "vector of points" << "\n";
    
    for (int i = 0; i < v.size(); i++){    
        os << v[i] << " ";
    } 
    os << "\n";
    return os;
}


//Compute the barycenter of the focis'.
inline Point foci_mid(vector<Point>& foci) {
    Point ret(foci[0].point.size());

    for (int i = 0; i < foci.size(); i++) {
        ret = ret + foci[i];
    }
    double den = (double)foci.size();
    ret = ret / den;
    return ret;
}

inline Point weizfield(vector<Point>& foci, Point& p) {
    Point a(p.point.size());
    double norm;
    Point num = (p.point.size());
    double den = 0;

    Point diff(p.point.size());
    Point div(p.point.size());
    Point ret(p.point.size());
    
    for (int i = 0; i < foci.size(); i++) {
        a = foci[i];
        diff = p - a;
        norm = normP(diff);
        div = (a / norm);
        num = num + div;
        den = den + (1 / norm);
    }
    ret = num / den;
    return ret;
}
int main(){
    
    try{
        vector<double> points;
        ifstream inFile;
        

        //Change file name!!!!!!

        inFile.open("test_files/test1.txt");
        //inFile.open("test_files/2d-unit-sphere/2d-unit-circle.txt");
        //inFile.open("test_files/3d-unit-sphere/3d-unit-circle.txt");
        //inFile.open("test_files/100d-unit-cube/100d-unit-cube.txt");
        //inFile.open("test_files/100d-unit-sphere/100d-unit-sphere.txt");

        //Change dimension!!!!!
        int dim_para = 2;

        double x;
        //Read the points from txt file.
        while (inFile >> x) {
        points.push_back(x);
        }

        inFile.close();

        vector<Point> data;
        vector<double> point;
        //Put the data inside a vector of points.
        for(int i =0; i<(points.size()/dim_para); i++){
        for(int j = 0; j< dim_para; j++){
        point.push_back(points[i*dim_para+j]);
        }
        data.push_back(Point(point));
        point.erase(point.begin(),point.begin()+point.size());
        }

        //cout<<data;


        int dim = data[0].point.size();
        Point p(dim);
        Point prev(dim);
        double diff;
        //clock_t begin = clock();
        p = foci_mid(data);
        cout<<p;
        diff=3;
        int c=0;
        //while(diff>=pow(10,-16)){
        while(diff!=0){ 
            c++;
            p = weizfield(data, p);
            diff = normP(p-prev);
            //cout<<p;
            cout<<diff<<"\n";
            prev=p;
            if(c>10)
                break;
        }
        
    


    
     

    }
    catch(exception &e){
        cout<<"exception caught:"<<e.what();
    }
    
    return 0;
}*/

