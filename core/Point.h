/* 
Author: Vinaik Chhetri
*/

#ifndef POINT_H
#define POINT_H

#include<iostream>
#include<cmath>
#include<fstream>
#include<sstream>
#include<vector>
#include <time.h>


#define CORE_LEVEL 3

#include "CORE.h"
using namespace std;


class Point {
public:
    Point(int dimension);
    Point(vector<double>& input_point);
    vector<double> point;
};

//Point constructor to initialise point vector with 0's.

inline Point::Point(int dimension){
    for(int i=0;i<dimension;i++){
        point.push_back(0);
    }
}

//Point constructor to initilise point vector variable to a vector of doubles.

inline Point::Point(vector<double>& input_point){
    point = input_point;
}

//Check if point p == point q.

inline bool operator==(const Point& p, const Point& q) {
    bool ret = true;
    for(int i=0; i<p.point.size();i++){
        ret = ret && (p.point[i] == q.point[i]);
    }
    return ret;
}

//?

inline bool operator!=(const Point& p, const Point& q) {
    return !(p == q);
}

//Return sum of 2 points.

inline Point operator+(const Point& p, const Point& q) {
    vector<double> sumP;
    for(int i=0; i<p.point.size();i++){
        sumP.push_back(p.point[i]+q.point[i]);
    }
    Point ret(sumP);
    return ret;
}

//Return substraction of 2 points.

inline Point operator-(const Point& p, const Point& q) {
    vector<double> subP;
    for(int i=0; i<p.point.size();i++){
        subP.push_back(p.point[i]-q.point[i]);
    }
    Point ret(subP);
    return ret;
}

//Return point divided by scalar.
inline Point operator/(const Point& p, double& a) {
    if (a==0) 
        throw runtime_error(" at Point.h->Point divided by 0");

    vector<double> divP;
     double sub = a.approx(120,120);
     //double sub = a;
    double ssub;
    for(int i=0; i<p.point.size();i++){
        ssub = p.point[i];
        ssub = ssub.approx(120,120);
        divP.push_back(ssub/sub);
    }
    Point ret(divP);
    return ret;
}

//Return point * scalar

inline Point operator*(const Point& p, double& a) {
    vector<double> mulP;
    double ssub;
    for(int i=0; i<p.point.size();i++){
        ssub = p.point[i];
        ssub = ssub.approx(120,120);
        mulP.push_back(ssub*a);
    }
    Point ret(mulP);
    return ret;
}

//Print point.

inline ostream& operator<<(std::ostream& os, Point& p) {
    double sub;
    os<<'('<<" ";
    for(int i=0; i<p.point.size();i++){
        sub = p.point[i];
        sub =  sub.approx(120,120);
        os << setprecision(10) << sub << " ";
    }
    os<< ')'<<"\n";
    return os;
}
/*
//*****need to be modified!!!
inline istream& operator>>(std::istream& is, Point& p) {
    double x, y;
    char c1, c2, c3;
    is >> c1 >> x >> c2 >> y >> c3;
    if (!is) return is;
    if (c1 != '(' || c2 != ',' || c3 != ')') {
        throw std::runtime_error("Input file has wrong format.");
        return is;
    }
   
    //p = Point(x, y);
    return is;
}
*/

//Return norm of a point.

inline double normP(const Point& p) {
    double ret = 0;
    double sub;
    double ssub;
    for(int i=0; i<p.point.size();i++){
        ssub = p.point[i];
        ssub = ssub.approx(120,120);
        ret = ret + (ssub * ssub) ;
    }
     sub = ret;
    sub = ret.approx(120,120);
    return sqrt(sub);
}

//Return distance between 2 points.

inline double dist(const Point &p1, const  Point &p2) {
    // return sqrt(pow(p2.x - p1.x, 2) + pow(p2.y - p1.y, 2));
    return normP(p1 - p2);
}

/*
inline int precision(double epsilon, double x){
    int dig_after  = ceil( log10( 1.0/ std::min( 0.1, epsilon) ) );
    int dig_before = ceil(log10( std::max( 1.0,std::abs(x) ) ) ) ;
    return dig_after+dig_before;
}
*/

//clock operations.

inline void print_time(clock_t begin){
    /* The function print_time(begin) shows the time since begin, defined as:
     clock_t begin = clock(); */
    clock_t end = clock();
    double elapsed_secs = (machine_double)(end - begin) / CLOCKS_PER_SEC;
    std::cout << "Elapsed time : " << elapsed_secs<< " sec"<<std::endl;
}

#endif /* POINT_H */


/*
int main(){
    Point p(2);
    vector<double> d = {1,1};
    Point newp(d);
    vector<double> e = {1,2};
    Point newpe(e);
   
    double c;
    c= dist(d,e);
   
    cout<<c<<"\n";
   
    return 0;
}
*/
