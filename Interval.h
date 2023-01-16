/*
 * File:   Interval.h
 * Author: Vinaik Chhetri
 *
 *
 * For more details on the algorithms see
 * "Certified Approximation Algorithms for the Fermat Point" by
 * Kolja Junginger, Ioannis Mantas, Evanthia Papadopoulou, Martin Suderland, and Chee Yap:
 * http://www1.pub.informatik.uni-wuerzburg.de/eurocg2020/data/uploads/eurocg20_proceedings.pdf
 */

#ifndef INTERVAL_H
#define INTERVAL_H

#include<cmath>
#include<vector>
#include<iostream>
#include <algorithm>

using namespace std;

class Interval {
public:
    Interval(double a, double b);
    Interval();
    double min;
    double max;
    bool isempty = false; //if b>a then invalid interval = empty interval.
    

    double length() const {
        double l;
        if (!isempty)
            l = max - min;
        else
            l = std::nan("1");
        return l;
    }
};

//Constructor to initialize start and end points of an interval. 
inline Interval::Interval(double a, double b)
{
    min = a;
    max = b;
    /*if (max < min){
        throw runtime_error(" at Interval.h: constructor(a,b) -> invalid interval");
    }
    */
    if (max < min)
        isempty = true;
}

//Constructor to create an interval only containing 0.
inline Interval::Interval()
{
    min  = 0;
    max = 0;
}

//Print Interval
inline ostream& operator<<(ostream& os, Interval& I) {
    return os << '[' << I.min << ',' << I.max << "]\n";
}

//2 Intervals equal?
inline bool operator==(const Interval& I, const Interval& J) {
    return (I.min == J.min && I.max == J.max);
}

//?
inline bool operator!=(const Interval& I, const Interval& J) {
    return !(I == J);
}

//Interval addition
inline Interval operator+(const Interval& I, const Interval& J) {
    Interval sumIJ = Interval(I.min + J.min, I.max + J.max);
    return sumIJ;
}

//Interval substraction
inline Interval operator-(const Interval& I, const Interval& J) {
    Interval subIJ = Interval(I.min - J.max, I.max - J.min);
    return subIJ;
}

//Interval - scalar
inline Interval operator-(const Interval& I, const double a_cord_i) {
    Interval subIJ = Interval(I.min - a_cord_i, I.max - a_cord_i);
    return subIJ;
}

//Interval * Interval 
inline Interval operator*(const Interval& I, const Interval& J) {
    vector<double> options = {I.min*J.min, I.min*J.max, I.max*J.min, I.max*J.max};
    double  min = * std::min_element(options.begin(),options.end());
    double  max = * std::max_element(options.begin(),options.end());
    Interval mul = Interval(min, max);
    return mul;
}

//Interval divided by Interval
inline Interval operator/(const Interval& I, const Interval& J) {
   
    double inf = std::numeric_limits<double>::infinity();
   
    /*
     double Jmin_inverse;
     if (J.min==0){
     //        Jmin_inverse = numeric_limits<double>::max();
     Jmin_inverse = inf;
     }
     else { Jmin_inverse =   1.0/J.min;}
     
     double Jmax_inverse;
     if (J.max==0){
     //        Jmax_inverse = numeric_limits<double>::max();
     Jmax_inverse = inf;
     }
     else { Jmax_inverse = 1.0/J.max; }
     */
    Interval Jinverse;
    if (0 < J.min || J.max < 0){
        Jinverse = Interval(1.0/J.max,1.0/J.min);
    }
    else if (J.min == 0){
        Jinverse = Interval(1.0/J.max, inf);
    }
    else if (J.max == 0){
        Jinverse = Interval(-inf,1.0/J.min);
    }
    else
        Jinverse = Interval(-inf,inf);
   
    Interval div = I*Jinverse;
    return div;
}

// checks whether Interval I is a subset of J.
inline bool subset(const Interval& I,const Interval& J) {
    
    //    double eps = pow(10,-12); //std::numeric_limits<double>::min();
    bool answer = ( (J.min <= I.min) && (I.max <= J.max) );
    return answer;
}

// returns the intersection of Intervals I and J.
inline Interval intersection(const Interval& I,const Interval& J) {
    
    Interval Int = Interval( max(I.min,J.min), min(I.max,J.max) );
   
    if (Int.min  > Int.max){
        Int.isempty = true;
        //throw runtime_error(" at Interval.h: funtion intersection-> invalid interval");
    }
    //std::cout << "Int " << Int.isempty;
    return Int;
}

//returns union of 2 intervals.
inline Interval convexhull(const Interval& I,const Interval& J) {
    return Interval(std::min(I.min,J.min), std::max(I.max,J.max) );
}

#endif /* INTERVAL_H */

