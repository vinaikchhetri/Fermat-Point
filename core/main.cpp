/*
Author: Vinaik Chhetri
*/

#include <iostream>
#include <vector>
#include <math.h>
#include <random>
#include <fstream>
#include <cstdlib>
#include "Point.h"
#include "Box.h"
#include "Interval.h"
#include "derivatives.h"

#define CORE_LEVEL 3

#include "CORE.h"



using namespace std;

/* original
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
*/

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

inline Point weizfield(vector<Point>& foci, Point& p) {
    Point a(p.point.size());
    double norm;
    Point num = (p.point.size());
    double den = 0;

    Point diff(p.point.size());
    Point div(p.point.size());
    Point ret(p.point.size());
    double sub;//Used BigFloat before.
    
    for (int i = 0; i < foci.size(); i++) {
        a = foci[i];
        diff = p - a;
        norm = normP(diff);
        //cout<<"hererere"<<norm.approx(120,120);
        
        sub = norm.approx(120,120);
        //sub = norm;

        div = a/sub;
        num = num + div;
        den = den + (Expr(1) / sub);
    }
    
    sub = den.approx(120,120);
    //sub = den;


    ret = num / sub;

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
  setPositionalFormat();
  setDefaultInputDigits(6);
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
cout<<points.size()<<"\n";//
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

cout<<data.size()<<"\n";//


int dim = data[0].point.size();
Point p(dim);
//clock_t begin = clock();
p = foci_mid(data);


double eps = 0.005;//0.005; //epsilon
//double eps = 0.000000000000000000005;
double l = eps; //Width of the box.
bool flag = true; //Flag to stop the wiezfield iteration.
Box B(dim);
Box NB(dim);
Box B10(dim);
Box NB10(dim);
Point ret(dim);
double sc = 0.1;
int i = 0;
/*B = center2box(p, l);
NB = Newton(B, data);
cout<<NB;
*/
while (flag) {
    //cout<<i;
    B = center2box(p, l);
    NB = Newton(B, data);
    
    //cout<<NB;
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
        l = min(Expr(10) * l, eps);
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
/*
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

inline Point weizfield(vector<Point>& foci, Point& p) {
    Point a(p.point.size());
    double norm;
    Point num = (p.point.size());
    double den = 0;

    Point diff(p.point.size());
    Point div(p.point.size());
    Point ret(p.point.size());
    double sub;//Used BigFloat before.
    
    for (int i = 0; i < foci.size(); i++) {
        a = foci[i];
        diff = p - a;
        norm = normP(diff);
        //cout<<"hererere"<<norm.approx(120,120);
        sub = norm.approx(960,960);
        div = a/sub;
        num = num + div;
        den = den + (Expr(1) / sub);
    }
    sub = den.approx(960,960);
    ret = num / sub;

    return ret;
}





//Trial

int main(){
    setPositionalFormat();
    setDefaultInputDigits(6);
    try{
        vector<double> points;
        ifstream inFile;
       

      
        //Change file name!!!!!!

        //inFile.open("test_files/test1.txt");
        //inFile.open("test_files/2d-unit-sphere/2d-unit-circle.txt");
        //inFile.open("test_files/3d-unit-sphere/3d-unit-circle.txt");
        inFile.open("test_files/100d-unit-cube/100d-unit-cube.txt");
        //inFile.open("test_files/100d-unit-sphere/100d-unit-sphere.txt");

        //Change dimension!!!!!
        int dim_para = 100;

        double x;
        //Read the points from txt file.
        //BigFloat y = 2.0;
        //readFromFile(y,inFile);
        //setDefaultRelPrecision(100);
        //setDefaultAbsPrecision(100);
        //setDefaultInputDigits(CORE_INFTY);
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
        cout << setprecision(5);
        cout<<data;


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
            if(c>150)
                break;
        }
        cout<<c<<"\n";
        cout<<"diff"<<diff<<"\n";
        //Expr e = "1/3";
        //cout<<setprecision(20)<<e.approx(120,120)<<"\n";
        //cout<<setprecision(40)<<e.approx(80,80)<<"\n";
       
        /*
        #include "Box.h"
#include "Interval.h"
#include "derivatives.h"
    */


    
    /* 

    }
    catch(exception &e){
        cout<<"exception caught:"<<e.what();
    }
    
    return 0;
}
*/

/*
int main(){
    setPositionalFormat();
    setDefaultInputDigits(6);
    /*
    vector<vector<double>> I = init(5);
    cout<<I;
    vector<vector<double>> I2 = init2(5);
    cout<<I2;
    */
/*
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
    //BigFloat y = 2.0;
    //readFromFile(y,inFile);
    //setDefaultRelPrecision(100);
    //setDefaultAbsPrecision(100);
    //setDefaultInputDigits(CORE_INFTY);
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
    cout << setprecision(5);
    cout<<data;
    vector<vector<Interval>> I;
    I.push_back({Interval(1,1),Interval(2,2)});
    I.push_back({Interval(3,3),Interval(4,4)});
    cout<<I;
    vector<vector<Interval>> J;
    J.push_back({Interval(1,1),Interval(3,3)});
    J.push_back({Interval(2,2),Interval(4,4)});
    cout<<J;
    //vector < vector < Interval >> K = I*J;
    J = I*J;
    cout<<J;
    /*
    vector<Interval> I = {Interval(-3,-2),Interval(-50,-20)};
    Box b(I);
    vector < vector < Interval >> h = hessian(b,data);
    cout<<h;
    vector<double>aa = {1,2};
    Point a(aa);
    a = jaccobian(a,data);
    cout<<a;
    vector<vector<double>> I2 = init2(5);
    I2[0][2] = 2;
    I2[0][4] = 2;
    I2[4][2] = 2;
    I2[4][4] = 2;
    cout<<I2;
    I2= inverse(I2);
    I2 = transpose(I2);
    cout<<I2;
    I2 = real_hessian(data,a);
    cout<<I2;*/
    /*vector<Interval> I = {Interval(-3,-2),Interval(-50,-20)};
    Box b(I);
    vector<Interval> J = {Interval(-4,-2),Interval(-100,-10)};
    Box c(J);
    double sc = 2;
    Box newsc = scaleBox(b,sc);
    //cout<<newsc;
    bool val = inclusion(b,c);
    cout<<val;
    */
    //Interval II = Interval(1,2);
    //Interval JJ = Interval(-4,-2);
    //Interval intt = intersection(II,JJ);
    //cout<<intt<<"\n";
    //cout<<intt.isempty || false;
    
    /*
    Interval I(1,2);
    cout<<"I"<<I;
    Interval J(-CORE_INFTY,CORE_INFTY);
    Interval K = I-1;
    cout<<"J"<<J;
    cout<<"K"<<K;
    */
    /*Expr e = "-1/3";
    cout<<e<<"\n";
    cout<<setprecision(100)<<e.approx(240,240)<<"\n";
    */
    //(e.approx(120,120)>-CORE_INFTY)?cout<<e.approx(120,120):cout<<"nope";

//}
    

