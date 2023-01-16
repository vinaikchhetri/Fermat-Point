/*
 * File:   Box.h
 * Author: Vinaik Chhetri
 *
 *
 * For more details on the algorithms see
 * "Certified Approximation Algorithms for the Fermat Point" by
 * Kolja Junginger, Ioannis Mantas, Evanthia Papadopoulou, Martin Suderland, and Chee Yap:
 * http://www1.pub.informatik.uni-wuerzburg.de/eurocg2020/data/uploads/eurocg20_proceedings.pdf
 */

#ifndef BOX_H
#define BOX_H

#include "Interval.h"
#include "Point.h"
#include<vector>
using namespace std;

class Box {
public:
    Box(vector<Interval>& I);
    Box(int dimension);
   
    vector<Interval> inte;
   
    double width() const{
        return inte[0].max-inte[0].min;
    }
   
    double radius() const;
   
    Point center() const;
};

//Initialise box with vector of intervals.
inline Box::Box(vector<Interval>& I){
    inte = I;
}

//initialise a box = origin.
inline Box::Box(int dimension){
    for (int i = 0; i < dimension; i++){
          inte.push_back(Interval(0,0));
        }
}

//Radius of a box
double Box:: radius() const{
    double ret = 0;
    for(int i =0; i<inte.size(); i++){
        ret+= pow(width(),2);
    }
    ret = sqrt(ret)/2;
   
    return ret;
}

//Center of box
Point Box::center() const{
    vector<double> c;
    for(int i=0;i<inte.size();i++){
        c.push_back((inte[i].max+inte[i].min)/2);
    }
    Point C(c);
    return C;
}

//Print box.
inline std::ostream& operator<<(std::ostream& os, Box& B) {
    cout<<"Box intervals"<<"\n";
    for(int i=0; i<B.inte.size();i++){
        os << B.inte[i];
    }
   
    return os;
}

// Checks whether Box B1 is a subset of B2,
// !!!Caution!!! Different from Martin's matlab implementation!!! Inputs swapped!
inline bool inclusion(Box& B1, Box& B2){
    bool answer = true;
    for(int i=0;i<B1.inte.size();i++){
        answer = answer && subset(B1.inte[i] , B2.inte[i]);
    }
    return answer;
}

// Returns true if intersection of Boxes B1 and B2 is empty
inline bool exclusion(Box& B1, Box& B2){
    Interval Int;
    Int = intersection(B1.inte[0] , B2.inte[0]);
    bool answer  = Int.isempty;
    
   
    for(int i=1;i<B1.inte.size();i++){
        Int = intersection(B1.inte[i] , B2.inte[i]);
        answer = answer || Int.isempty;
    }
    //cout << "excluded "<<answer;
    return answer;
}

//Returns the intersected box between B1 and B2.
//Only run this function if exclusion evaluates to True.
inline Box intersection(Box& B1 , Box& B2){
    Interval Int;
    vector<Interval> intes;
    for(int i=0; i<B1.inte.size(); i++){
        Int = intersection(B1.inte[i] , B2.inte[i]);
        intes.push_back(Int);
    }
    return Box(intes);
}

//Returns true if point P in box B.
inline bool contained(const Point& P , const Box& B){
    bool answer = true;
    for(int i=0;i<B.inte.size();i++){
        answer = answer && ( (P.point[i]>=B.inte[i].min) && (P.point[i]<=B.inte[i].max) );
    }
    return answer;
}

//Return box centered at point center.
inline Box center2box(const Point& center, double width){
    vector<Interval> intes;
    vector<double> I;
    for(int i=0;i<center.point.size();i++){
        intes.push_back(Interval(center.point[i] - width/2, center.point[i] + width/2));
    }
    return Box(intes);
}

//Computes Box cB = c*B where the dilated box has its center at the midpoint of the box B
inline Box scaleBox(const Box& B, double c){
    Point Center = B.center();
    double new_width = c * B.width();
    Box cB = center2box(Center,new_width);
    return cB;
}

#endif /* BOX_H */

