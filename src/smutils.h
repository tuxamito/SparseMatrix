/************************************************************************/
/* Copyright (C) 2017 Daniel Rubio Bonilla                              */
/* email: <danielrubiob_at_gmail_ dot _com>                             */
/* link: https://github.com/tuxamito/SparseMatrix                       */
/*                                                                      */
/* This file is free software: you can redistribute it and/or modify    */
/* it under the terms of the GNU LesserGeneral Public License as        */
/* published by the Free Software Foundation, either version 3 of the   */
/* License, or (at your option) any later version.                      */
/*                                                                      */
/* This file is distributed in the hope that it will be useful,         */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        */
/* GNU General Public License for more details.                         */
/*                                                                      */
/* You should have received a copy of the GNU Lesser General Public     */
/* License along with this program.  If not,                            */
/* see <http://www.gnu.org/licenses/>.                                  */
/************************************************************************/

// Date: 20171208
// Version: 0.1

#ifndef SMUTILS_H_
#define SMUTILS_H_

#include <vector>
#include "sparsematrixst.h"

using namespace std;

// generate a random number between min and max
double drand(double min, double max) {
  double diff = max - min;
  double r = drand48() * diff;
  return min + r;
}

// generate a ST Matrix of size r,c and with maximum er non Zero elements
// per row with values between min and max
template<typename T>
SparseMatrixST<T> generateSM(IT r, IT c, IT er, double min=-10, double max=10) {
  SparseMatrixST<T> a(r, c);
  vector<int> cols(er);

  for(int i=0; i<r; i++) {
    // select columns, including diagonal
    for(int j=0; j<er; j++) {
      cols[j] = rand() % c;
    }
    if(i<c)
      cols[0] = i;
    sort(cols.begin(), cols.end());

    // generate values
    for(int j=0; j<er; j++) {
      T val = (T)drand(min, max);
      a.set(val, i, cols[j]);
    }
  }

  return a;
}

// generate a vector of size, with values between min and max
template<typename T>
vector<T> generateVector(IT size, double min=-10, double max=10) {
  vector<T> v(size);
  for(IT i=0; i<size; i++) {
    v[i] = drand(min, max);
  }
  return v;
}

#endif // SMUTILS_H_
