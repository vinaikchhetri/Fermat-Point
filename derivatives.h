/*
 * File:   derivatives.h
 * Author: Vinaik Chhetri
 *
 * Created on March 3, 2020, 9:16 AM
 *
 * * For more details on the algorithms see
 * "Certified Approximation Algorithms for the Fermat Point" by
 * Kolja Junginger, Ioannis Mantas, Evanthia Papadopoulou, Martin Suderland, and Chee Yap:
 * http://www1.pub.informatik.uni-wuerzburg.de/eurocg2020/data/uploads/eurocg20_proceedings.pdf
 */

#ifndef DERIVATIVES_H
#define DERIVATIVES_H

#include<vector>
#include "Point.h"
#include "Interval.h"
#include "Box.h"

using namespace std;

//Custom data_type used inplace of a simple matrix while computing the inverse of a matrix
//to track permutations.

struct permuted_mat{
  bool flag;
  vector < vector < double >>v;
  vector < int >index; 
  int nz; //index of first row with non-zero entry as a specific column index.
};

Box Newton (const Box &B, const vector < Point > &foci);

ostream & operator<< (std::ostream & os, vector < vector < Interval >> &v);
ostream & operator<< (std::ostream & os, vector < int >&v);
ostream & operator<< (std::ostream & os, vector < Interval > &v);
ostream & operator<< (std::ostream & os, vector < vector < double >>&v);

vector < Interval > point2Interval (const Point & p);
vector < vector < Interval >> mat2Interval (const vector < vector <double >>&mat);

Point term_i (const Point & p, const Point & a_i);
Point jaccobian (const Point & p, const vector < Point > &foci);

vector < vector < Interval >> hes_eval (const Box & B, const Point & a);
vector < vector < Interval >> hessian (const Box & B, const vector < Point > &foci);

vector < vector < double >>init (const int dim);
vector < vector < double >>init2 (const int dim);

bool check_a0 (const vector < double >&row, const int ind);
permuted_mat permute (const permuted_mat & permuted_act, const int ind);
bool sign_comp (const double a, const double b);
vector < vector < double >>permute_b (const permuted_mat & permute);
vector < vector < double >>change_lower (const vector < vector <double >>&lower, const int from, const int to);
vector < vector < double >>transpose (const vector < vector < double >>&mat);
vector < vector < Interval >>transposeInt (const vector < vector < Interval >>&mat);
vector < vector < double >>inverse (const vector < vector < double >>&hes);

vector < vector < double >>real_hessian (const vector < Point > &foci, const Point & p);

vector < Interval > jacobian_box (const Box & B, const vector < Point > &foci);

vector < Interval > operator- (const Box & B, const Point & mB);
vector < Interval > operator/ (const vector < Interval > &v, const Interval & I);
vector < Interval > operator+ (const vector < Interval > &v1, const vector < Interval > &v2);
Interval operator* (const vector < Interval > &row, const vector < Interval > &col);
vector < Interval > operator- (const vector < Interval > &row1, const vector < Interval > &row2);


vector < vector < Interval >> operator- (const vector < vector <Interval >> &m1, const vector < vector <Interval >> &m2);
vector < Interval > operator* (const vector < vector < Interval >> &mat, const vector < Interval > &col);
vector < vector < Interval >> operator* (const vector < vector <Interval >> &inv, const vector < vector <Interval >> &hes);
vector < vector < Interval >> operator+ (const vector < vector <Interval >> &M1, const vector < vector <Interval >> &M2);



//Returns the Newton operator of the box.
inline Box Newton (const Box &B, const vector < Point > &foci){
  /* This function Newton maps a Box B to a new box N(B) using the Newton formula by Krawzcyk: enter reference...
     N(B) = mB - K * f(mB) + (I - K * J_f(B)) * (B - mB)
          = final - partB + (C1 * C2)
          = final - partB + C
     
     N(B)  : Newton operator for the box B.
     mB    : center of the box.
     K     : Inverse of Real Hessian matrix evaluated at the center of the box i.e, mB.
     f(mB) : Jaccobian of then fermat function evaluated at mB.
     I     : Identity matrix.
     J_f(B): Interval Hessian matrix = Jaccobian (Jaccobian of the box).
   */

  //mB
  Point center (foci[0].point.size ()); 
  center = B.center ();

  //f(mB)
  Point gradient_phi = jaccobian (center, foci);

  //J_f(B)
  vector < vector < Interval >> box_hessian = hessian (B, foci);

  //Real Hessian Matrix.
  vector < vector < double >>real_hessian_mat;

  //K transposed.
  vector < vector < double >>inv_transposed;

  //K
  vector < vector < double >>inv;

  //K * J_f(B)
  vector < vector < Interval >> K_Jf;

  //(I - K * J_f(B))
  vector < vector < Interval >> C1;

  //Identity matrix.
  vector < vector < double >>I;

  //(B - mB)
  vector < Interval > C2;

  //(I - K * J_f(B)) * (B - mB)
  vector < Interval > C;

  //Real Hessian Matrix.
  real_hessian_mat = real_hessian (foci, center);

  //K transposed.
  inv_transposed = inverse (real_hessian_mat);

  //K * J_f(B)
  //K_Jf =  mat2Interval (inv_transposed) * box_hessian;
  K_Jf =  mat2Interval (transpose(inv_transposed)) * transposeInt(box_hessian);

  //Identity matrix.
  I = init2 (inv_transposed.size ());

  //(I - K * J_f(B))
  C1 = mat2Interval (I) - K_Jf;

  //(B - mB)
  C2 = B - center;

  //(I - K * J_f(B)) * (B - mB)
  C = C1 * C2;

  
  vector < Interval > f_mb;
  vector < Interval > partB;
  
  //f(mB) to intervals.
  f_mb = point2Interval (gradient_phi);

  //Inverse i.e. K
  inv = transpose (inv_transposed);

  //Turn inv to intervals and store it in K.
  vector < vector < Interval >> K;
  K = mat2Interval (inv);

  //K * f(mB)
  partB = K * f_mb;

  //turn center to intervals.
  vector < Interval > final = point2Interval (center);

  //mB - K * f(mB) + (I - K * J_f(B)) * (B - mB)
  final = final - partB;
  final = final + C;
   
  Box NB(final);
  return NB;
}

/*Start of Printing matrices and vector*/

//print matrix of intervals.
inline std::ostream & operator<< (std::ostream & os, vector < vector < Interval >> &v)
{
  cout << "dxd Interval matrix" << "\n";
  
  for (int i = 0; i < v.size (); i++){
      for (int j = 0; j < v[0].size (); j++){
        os << v[i][j] << " ";
      }
      os << "\n";
    }
  return os;
}

//print matrix of doubles
inline std::ostream & operator<< (std::ostream & os, vector < vector < double >>&v){
  cout << "dxd real number matrix" << "\n";
  
  for (int i = 0; i < v.size (); i++){
    for (int j = 0; j < v[0].size (); j++){
      os << v[i][j] << " ";
    }
    os << "\n";
  }
  return os;
}

//print vector of integers.
inline std::ostream & operator<< (std::ostream & os, vector < int >&v){
  cout << "indices" << "\n";
  
  for (int i = 0; i < v.size (); i++){
    os << v[i] << " ";
  } 
  os << "\n";
  return os;
}

//print vector of intervals
inline std::ostream & operator<< (std::ostream & os, vector < Interval > &v){
  cout << "vector of intervals" << "\n";
  
  for (int i = 0; i < v.size (); i++){
    os << v[i] << " ";
  } 
  os << "\n";
  return os;
}

/*End of Printing matrices and vector*/

/*Start of converting point and matrices to intervals*/

//Convert from point to interval.
inline vector < Interval > point2Interval (const Point & p){
  vector < Interval > ret;
  Interval temp;

  for (int i = 0; i < p.point.size (); i++){
    temp = Interval (p.point[i], p.point[i]);
    ret.push_back (temp);
    }
    return ret;
}

//Convert from matrix of doubles to matrix of intervals.
inline vector < vector < Interval >> mat2Interval (const vector < vector < double >>&mat){
  vector < vector < Interval >> ret;
  vector < Interval > row;
  Interval temp;

  for (int i = 0; i < mat.size (); i++){
      for (int j = 0; j < mat[0].size (); j++){
        temp = Interval (mat[i][j], mat[i][j]);
        row.push_back (temp);
      }
      ret.push_back (row);
      row.erase (row.begin (), row.begin () + row.size ());
    }
    return ret;
}

/*End of converting point and matrices to intervals*/

/*Start of computing Jaccobian of p = sum over all a of {(p-a)/norm(p-a)}*/

//Computing (p-a)/norm(p-a) for a specific a.
inline Point term_i (const Point & p, const Point & a_i){
  Point ret (p.point.size ());
  double norm;
  ret = p - a_i;
  norm = normP (ret);
  ret = ret / norm;
  return ret;
}

//sum over all a of {(p-a)/norm(p-a)}.
inline Point jaccobian (const Point & p, const vector < Point > &foci){
  Point sum (p.point.size ());

  for (int i = 0; i < foci.size (); i++){
    sum = sum + term_i (p, foci[i]);
  }
  return sum;
}

/*End of computing Jaccobian of p.*/

/*start of computing box hessian*/

//Hessian box matrix = sum over all focis of {dxd matrix correspoding to a specific foci}
//Compute {dxd matrix correspoding to a specific foci}.
inline vector < vector < Interval >>hes_eval (const Box & B, const Point & a){
  vector < vector < Interval >> ret;
  vector < Interval > row;
  int d = B.inte.size ();
  Point mB = B.center ();
  double n = normP (mB - a);
  Interval den = Interval (n - B.radius (), n + B.radius ()); //This is an approximation for ||B-a||
  
  //Computing ||B-a||^2
  Interval den_sq;
  if (den.min < 0 && den.max > 0){
    den_sq.min = 0;
  }
  den_sq.min = den.min * den.min;
  den_sq.max = den.max * den.max;
  
  //Computing ||B-a||^3 for the denominator for each of the entries.
  Interval den_cub (pow (den.min, 3), pow (den.max, 3));

  Interval f;
  Interval s;

  for (int i = 0; i < d; i++){
    for (int j = 0; j < d; j++){
    //if i==j then compute ||B-a||^2 - (B_i-a_i)^2
      if (i == j){
        f = B.inte[i] - a.point[i];
        if (f.min < 0 && f.max > 0){
          den_sq.min = 0;
        }
        f.min = f.min * f.min;
        f.max = f.max * f.max;
        f = den_sq - f;
        f = f / den_cub;
        row.push_back (f);
      }
      else{
        f = (B.inte[i] - a.point[i]);
        s = (B.inte[j] - a.point[j]);
        f = f * s;
        f = f * Interval (-1.0, -1.0);
        f = f / den_cub;
        row.push_back (f);
      }
    }
    ret.push_back (row);
    row.erase (row.begin (), row.begin () + row.size ());
  }
  return ret;
}

// Compute the box approximation for Hessian of fermat distance function.
inline vector < vector < Interval >>hessian (const Box & B, const vector < Point > &foci){
  int k = foci.size ();
  vector < vector < Interval >> hes_pt;
  vector < vector < Interval >> ret;
  vector < Interval > row;
  //Hessian box matrix will be computed as follows sum over all a of {dxd matrix}
  ret = hes_eval (B, foci[0]); //Compute the dxd matrix corresponding to the 1st foci.

  //sum over all a of {dxd matrix}
  for (int i = 1; i < foci.size (); i++){
    hes_pt = hes_eval (B, foci[i]);
    ret = ret + hes_pt;
  }
  return ret;
}

/*end of computing box hessian*/

/*Start of Initialisation of matrices.*/
//return a zero filled matrix.
inline vector < vector < double >>init (const int dim){
  vector < vector < double >>upper;
  vector < double >row;
  
  for (int j = 0; j < dim; j++){
    row.push_back (0);
  }
  for (int i = 0; i < dim; i++){
    upper.push_back (row);
  }
  return upper;
}

//return an identity matrix.
inline vector < vector < double >>init2 (const int dim){
  vector < vector < double >>upper;
  vector < double >row;
  
  for (int j = 0; j < dim; j++){
    row.push_back (0);
  }
  for (int i = 0; i < dim; i++){
    row[i] = 1;
    upper.push_back (row);
    row[i] = 0;
  }
  return upper;
}

/*end of Initialisation of matrices.*/

/*Start of PLU decomposition of real matrices.*/

//Check if the there is a zero in the ind_th index of the vector.
inline bool check_a0 (const vector < double >&row, const int ind){
  if (row[ind] == 0){
    return true;
  }
  return false;
}

//Swap the ind_th row of the matrix with a following row containing non-zero double in the ind_th column.
inline permuted_mat permute (const permuted_mat & permuted_act, const int ind){ 
//parameters : struct containing current state of the matrix, index of the row to be swapped.
int i = ind + 1;
permuted_mat permuted;
permuted = permuted_act;

vector < double >t0;
vector < double >ti;
vector < vector < double >>v = permuted.v;

//Loop until row containing non-zero double in the ith column is found.
while (i < v.size () && v[i][ind] == 0){
  i++;
}

permuted.nz = i; //keep track of the index of the row we just found. Will be used in computing the L matrix.

//Swapping the rows.
//erasing the rows first and storing in temporary variables.
t0 = v[ind];
ti = v[i];
v.erase (v.begin () + ind, v.begin () + ind + 1);
v.erase (v.begin () + i - 1, v.begin () + i);

//inserting rows in the new locations.
v.insert (v.begin () + ind, ti);
v.insert (v.begin () + i, t0);
permuted.v = v;
permuted.flag = true;

//permute the indices to reflect the permutation that took place.
int temp = permuted.index[i];
permuted.index[i] = permuted.index[ind];
permuted.index[ind] = temp;

return permuted;
}

//Checking if signs of 2 doubles are the same.
inline bool sign_comp (const double a, const double b){
  bool val;
  if (a <= 0 && b <= 0){
    val = true;
  }
  else if (a > 0 && b > 0){
    val = true;
  }
  else
    val = false;
  return val;
}

//permute the identity matrix(I) in (P)(LU)(D) = (P)(I) based on permutation P.
inline vector < vector < double >> permute_b (const permuted_mat & permute){
vector < int >indices = permute.index; //Indices which store the row indices in permuted order.
vector < vector < double >>I;
vector < vector < double >>ret; //store the permuted identity matrix.
I = init2 (indices.size ()); //identity matrix.

for (int i = 0; i < indices.size (); i++){
ret.push_back (I[indices[i]]); //Push into ret rows of identity matrix in the order of indices.
}
return ret;
}

//As matrix A is permuted, the lower matrix being generated must also be permuted.
//returns permuted lower matrix.
inline vector < vector < double >> change_lower (const vector < vector < double >>&lower, const int from , const int to){
  double temp;
  vector < vector < double >>ret;
  ret = lower;

  for (int i = 0; i < from; i++){
    temp = ret[from][i];
    ret[from][i] = ret[to][i];
    ret[to][i] = temp;
  }
  return ret;
}

//transpose matrix of doubles.
inline vector < vector < double >> transpose (const vector < vector < double >>&mat){
  vector < vector < double >>ret;
  vector < double >vec;

  for (int j = 0; j < mat[0].size (); j++){
    for (int i = 0; i < mat.size (); i++){
      vec.push_back (mat[i][j]);
    }
    ret.push_back (vec);
    vec.erase (vec.begin (), vec.begin () + vec.size ());
  }
  return ret;
}

//transpose matrix of Intervals.
inline vector < vector < Interval >> transposeInt (const vector < vector < Interval >>&mat){
  vector < vector < Interval >>ret;
  vector < Interval >vec;

  for (int j = 0; j < mat[0].size (); j++){
    for (int i = 0; i < mat.size (); i++){
      vec.push_back (mat[i][j]);
    }
    ret.push_back (vec);
    vec.erase (vec.begin (), vec.begin () + vec.size ());
  }
  return ret;
}

//Compute the inverse transposed!!! of a dxd real matrix using PLU decomposition.
inline vector < vector < double >> inverse (const vector < vector < double >>&hes){
  vector < double >vec;
  //Store dxd real hessian matrix in hes.
  vector < vector < double >>hesc;
  hesc = hes;

  //Initialise the L and U matrices.
  int d = hesc.size ();
  vector < vector < double >>lower;
  lower = init2 (d); //Initialise it as an identity matrix.
  vector < vector < double >>upper;
  upper = init (d); //Initilise it as a 0-filled matrix.

  /*Create and initialise a custom data type
     i.e. permuted to keep track of the orders of rows(permutations),
     the current matrix and
     flag to determine if any permutations have occured. */
  
  permuted_mat permuted;
  bool flag_permute = false;
  for (int i = 0; i < d; i++){
    permuted.index.push_back (i);
  }
  permuted.v = hesc;
  permuted.flag = false;

  //Insert the first row of the matrix to be decomposed into upper matrix
  //upper.insert (upper.begin (), permuted.v[0]);
  upper.erase(upper.begin (),upper.begin ()+1);
  upper.insert (upper.begin (), permuted.v[0]);

  double s;
  double m;
  double value;

  //Begin decomposition
  for (int i = 0; i < d - 1; i++){
  //check if the ith row, ith column of the matrix is nonzero.
    if (check_a0 (permuted.v[i], i) == true){
    //Permute current state of the matrix, upper and lower matrices.
      permuted = permute (permuted, i);
      upper = permuted.v;
      lower = change_lower (lower, i, permuted.nz);
    }

    //Reduce ith cols of rows below row (i+1)
    for (int j = i + 1; j < d; j++){
    //Skip the row which has zero in the ith col.
      if (permuted.v[j][i] == 0){
        vec = permuted.v[j];
      }
      //If row does not have 0 in the ith col then...
      else
      {
      //Compare the signs of the element at location i,j and location j,j.
        if (sign_comp (permuted.v[j][i], permuted.v[i][i]) == true){
          s = -1.0;
        }
        else{
          s = 1.0;
        }
        
        //Compute A[i][j]/A[j][j]
        m = s * (permuted.v[j][i] / permuted.v[i][i]);
        //Update the lower matrix.
        lower[j][i] = -1 * m;
        
        //Push back 0's i number of times as first i cols for rows below row i are 0s..
        for (int ind = 0; ind <= i; ind++){
          vec.push_back (0.0);
        }
        
        //As for d-i cols which do not necessarily turn into 0s, compute A[j][k] + (m * A[i][k]) for k = i+1...d
        for (int k = i + 1; k < d; k++){
          value = permuted.v[j][k] + (m * permuted.v[i][k]);
          vec.push_back (value);
        }
      }
      //Replace the existing row to the reduced version.
       upper.erase (upper.begin () + j, upper.begin () + j + 1);
       upper.insert (upper.begin () + j, vec);
           
      //erase vec for the next row.
      vec.erase (vec.begin (), vec.begin () + vec.size ());
    }
    permuted.v = upper; //Change the matrix to the current reduced version.
  }
  //L and U matrices have been computed above.
   
  /*Compute the inverse of the matrix using backward and forward substitution*/
   
  double sum = 0.0;
  vector < vector < double >>permuted_I;
  vector < vector < double >>matrix;
  vector < vector < double >>inv_mat;
  permuted_I = permute_b (permuted);
  vector < double >col;

  //Forward substitution.
  for (int i = 0; i < d; i++){
    for (int j = 0; j < d; j++){
      sum = permuted_I[j][i];
      if (j == 0){
        sum = sum / lower[j][j];
        col.push_back (sum);
      }
      else
      {
        for (int k = 0; k < j; k++){
          sum = sum - lower[j][k] * col[k];
        }
        sum = sum / lower[j][j];
        col.push_back (sum);
      }
    }
    matrix.push_back (col);
    col.erase (col.begin (), col.begin () + col.size ());
  }

  for (int i = 0; i < d; i++){
    col.push_back (0);
  }
   
  //Backward substitution.
  for (int rowe = 0; rowe < d; rowe++){
    for (int rowu = d - 1; rowu >= 0; rowu--){
      sum = matrix[rowe][rowu];

      if (rowu == d - 1){
        sum = sum / upper[rowu][rowu];
        col[rowu] = sum;
      }
      else{
        for (int k = d - 1; k > rowu; k--){
          sum = sum - (upper[rowu][k] * col[k]);
        }
        sum = sum / upper[rowu][rowu];
        col[rowu] = sum;
      }
    }
    inv_mat.push_back (col);
    for (int i = 0; i < d; i++){
      col.push_back (0);
    }
  }
  //return the inverse transposed matrix.
  return inv_mat;
}

/*End of PLU decomposition of real matrices.*/


/*Start of computing the real hessian matrix*/
inline vector < vector < double >>real_hessian (const vector < Point > &foci, const Point & p){
  int d = p.point.size ();
  vector < vector < double >>mat = init (d);
  double den;
  double norm_sq;
  Point diff (d);
  double computation;

  for (int a = 0; a < foci.size (); a++){
    diff = p - foci[a];
    den = normP (diff);
    norm_sq = pow (den, 2);
    den = pow (den, 3);

    for (int i = 0; i < d; i++){
      for (int j = 0; j < d; j++){
        if (i == j){
          computation = norm_sq - pow ((p.point[i] - foci[a].point[i]), 2);
          mat[i][j] += computation / den;
        }
        else
        {
          computation = -1 * (p.point[i] - foci[a].point[i]) * (p.point[j] - foci[a].point[j]);
          mat[i][j] += computation / den;
        }
      }
    }
  }
  return mat;
}
/*End of computing the real hessian matrix*/

/*Start of computing the Jacobian of the box*/
inline vector < Interval > jacobian_box (const Box & B, const vector < Point > &foci){
  Point a (foci[0].point.size ());
  Point center = B.center ();
  Point t (a.point.size ());
  double r = B.radius ();
  double n;
  Interval den;
  double norm1;
  double norm2;
  vector < Interval > ret;
  vector < Interval > retr;

  for (int i = 0; i < a.point.size (); i++){
    retr.push_back (Interval ());
  }
  
  for (int i = 0; i < foci.size (); i++){
    a = foci[i];
    t = center - a;
    n = normP (t);
    norm1 = n - r;
    norm2 = n + r;

    den = Interval (norm1, norm2);

    ret = B - a;
    ret = ret / den;
    retr = retr + ret;
  }
  return retr;
}
/*End of computing the Jacobian of the box*/


/*Operations of vector of intervals.*/
inline vector < Interval > operator- (const Box & B, const Point & mB){
  vector < Interval > ret;
  Interval interval;
  
  for (int i = 0; i < B.inte.size (); i++){
    interval = B.inte[i] - mB.point[i];
    ret.push_back (interval);
  }
  return ret;
}

inline vector < Interval > operator/ (const vector < Interval > &v, const Interval & I){
  vector < Interval > ret;
  Interval interval;
  
  for (int i = 0; i < v.size (); i++){
    interval = v[i] / I;
    ret.push_back (interval);
  }
  return ret;
}

inline vector < Interval > operator+ (const vector < Interval > &v1, const vector < Interval > &v2){
  vector < Interval > ret;
  Interval interval;
  
  for (int i = 0; i < v1.size (); i++){
    interval = v1[i] + v2[i];
    ret.push_back (interval);
  }
  return ret;
}

inline Interval operator* (const vector < Interval > &row, const vector < Interval > &col){
  Interval interval;

  for (int i = 0; i < row.size (); i++){
    interval = interval + (row[i] * col[i]);
  }
  return interval;
}

inline vector < Interval > operator- (const vector < Interval > &row1, const vector < Interval > &row2){
  Interval interval;
  vector < Interval > ret;

  for (int i = 0; i < row1.size (); i++){
    interval = row1[i] - row2[i];
    ret.push_back (interval);
  }
  return ret;
}



/*
/* Matrix operations */

inline vector < vector < Interval >> operator- (const vector < vector < Interval >> &m1, const vector < vector < Interval >> &m2){
  vector < Interval > row;
  vector < vector < Interval >> ret;
  for (int i = 0; i < m1.size (); i++){
    row = m1[i] - m2[i];
    ret.push_back (row);
  }
  return ret;
}

inline vector < Interval > operator* (const vector < vector < Interval >> &mat, const vector < Interval > &col){
  Interval interval;
  vector < Interval > ret;

  for (int i = 0; i < mat.size (); i++){
    interval = mat[i] * col;
    ret.push_back (interval);
  }
  return ret;
  }

//multiplication of matrices is defined as row to row. not row to col.
inline vector < vector < Interval >> operator* (const vector < vector < Interval >> &inv, const vector < vector < Interval >> &hes){
  Interval interval;
  vector < Interval > row;
  vector < vector < Interval >> ret;


  for (int i = 0; i < inv.size (); i++){
    for (int j = 0; j < hes.size (); j++){
      interval = inv[i] * hes[j];
      row.push_back (interval);
    }
    ret.push_back (row);
    row.erase (row.begin (), row.begin () + row.size ());
  }
  return ret;
}

inline vector < vector < Interval >> operator+ (const vector < vector < Interval >> &M1, const vector < vector < Interval >> &M2){
//Computes Bout = B1-B2 for two intervals
  int rows1 = M1.size ();
  int rows2 = M2.size ();
  if (rows1 != rows2){
    throw runtime_error("ERROR: Matrices do not have the fitting size for addition.");
  }
  vector < vector < Interval >> ret;
  vector < Interval > row;
  for (int i = 0; i < rows1; i++){
   row = M1[i]+M2[i];
   ret.push_back (row);
  }
  /*for (int i = 0; i < rows1; i++){
    for (int j = 0; j < M1.size (); j++){
      Interval M1i = M1[i][j];
      Interval M2i = M2[i][j];
      M1i = M1i + M2i;
      row.push_back (M1i);
    }
    ret.push_back (row);
    row.erase (row.begin (), row.begin () + row.size ());
  }*/
  return ret;
}

#endif /* DERIVATIVES_H */

