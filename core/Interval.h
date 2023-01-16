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

#define CORE_LEVEL 3

#include "CORE.h"

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
    double inf = CORE_INFTY;
    double  min;
    double  max;
    bool rflag = false;
    bool lflag = false;
    /*
    bool rflag = false;
    bool lflag = false;
    if (I.min<=-inf || J.min<=-inf){
        min = -inf;
        lflag = true;
    }
        
    if (I.max<=-inf || J.max<=-inf){
        max = -inf;
        rflag = true;
    }
      
    if (I.min>=inf || J.min>=inf){
        min = inf;
        lflag = true;
    }
       
    if (I.max>=inf || J.max>=inf){
        max = inf;
        rflag = true;
    }
    if lflag==true && rflag==true
        return Interval(min, max);
    if lflag==true 
        return Interval(min, I.max + J.max);
    if rflag==true
        return Interval(I.min + J.min,, max);
    */
    //Under the assumption that we do not have interrvals [-inf,-inf] and[inf,inf]
    if (I.min<=-inf || J.min<=-inf){
        min = -inf;
        lflag = true;
    }
    if (I.max>=inf || J.max>=inf){
        max = inf;
        rflag = true;
    }
    if (lflag==true && rflag==true)
        return Interval(min, max);
    if (lflag==true) 
        return Interval(min, I.max + J.max);
    if (rflag==true)
        return Interval(I.min + J.min, max);

    Interval sumIJ = Interval(I.min + J.min, I.max + J.max);
    return sumIJ;
}

//Interval substraction
inline Interval operator-(const Interval& I, const Interval& J) {
    //Under the assumption that we do not have interrvals [-inf,-inf] and[inf,inf]
    double inf = CORE_INFTY;
    double  min;
    double  max;
    bool rflag = false;
    bool lflag = false;
    if (I.min<=-inf || J.max>=inf){
        min = -inf;
        lflag = true;
    }        
    if (I.max>=inf || J.min<=-inf){
        max = inf;
        rflag = true;
    }
    if (lflag==true && rflag==true)
        return Interval(min, max);
    if (lflag==true) 
        return Interval(min, I.max - J.min);
    if (rflag==true)
        return Interval(I.min - J.max, max);
    Interval subIJ = Interval(I.min - J.max, I.max - J.min);
    return subIJ;
}

//Interval - scalar
inline Interval operator-(const Interval& I, const double a_cord_i) {
    //Under the assumption that we do not have interrvals [-inf,-inf] and[inf,inf]
    double inf = CORE_INFTY;
    double  min;
    double  max;
    bool rflag = false;
    bool lflag = false;
    if (I.min<=-inf){
        min = -inf;
        lflag = true;
    }
    if (I.max>=inf){
        max = inf;
        rflag = true;
    }
    if (lflag==true && rflag==true)
        return Interval(min, max);
    if (lflag==true) 
        return Interval(min, I.max - a_cord_i);
    if (rflag==true)
        return Interval(I.min - a_cord_i, max);
    Interval subIJ = Interval(I.min - a_cord_i, I.max - a_cord_i);
    return subIJ;
}

//Interval * scalar 
inline Interval operator*(const Interval& I, const double a) {
   

    double i = I.min*a;
    double j = I.max*a;
    Interval mul = Interval(i, j);
    return mul;
}

//Interval * Interval 
inline Interval operator*(const Interval& I, const Interval& J) {
    double inf = CORE_INFTY;
    double sub1 = I.min;
    sub1 = sub1.approx(120,120);
    double sub2 = I.max;
    sub2 = sub2.approx(120,120);
    double sub3 = J.min;
    sub3 = sub3.approx(120,120);
    double sub4 = J.max;
    sub4 = sub4.approx(120,120);
    vector<double> options = {sub1*sub3, sub1*sub4, sub2*sub3, sub2*sub4};
    double  min = * std::min_element(options.begin(),options.end());
    double  max = * std::max_element(options.begin(),options.end());
    if (min<=-inf)
        min = -inf;
    if (max<=-inf)
        max = -inf;
    if (min>=inf)
        min = inf;
    if (max>=inf)
        max = inf;
    Interval mul = Interval(min, max);
    return mul;
}

//Interval divided by Interval
inline Interval operator/(const Interval& I, const Interval& J) {
   
    //double inf = std::numeric_limits<double>::infinity();
    double inf  = CORE_INFTY;
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
    double sub1;
    double sub2;
    Interval Jinverse;
    if (0 < J.min || J.max < 0){
        sub1 = J.min;
        sub1=sub1.approx(120,120);
        sub2 = J.max;
        sub2 = sub2.approx(120,120);
        Jinverse = Interval(Expr(1)/sub2,Expr(1)/sub1);
    }
    else if (J.min == 0){
        sub2 = J.max;
        sub2 = sub2.approx(120,120);
        Jinverse = Interval(Expr(1)/sub2, inf);
    }
    else if (J.max == 0){
        sub1 = J.min;
        sub1=sub1.approx(120,120);
        Jinverse = Interval(-inf,Expr(1)/sub1);
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

