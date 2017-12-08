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

#ifndef SPARSEMATRIXST_H_
#define SPARSEMATRIXST_H_

#include <vector>
#include <utility>
#include <algorithm>
#include <functional>
#include <ostream>
#include <string>
#include <sstream>

// Integer definition, it should be signed!
#define IT int32_t

using namespace std;

template<typename T>
class SparseMatrixMT;

template<typename T>
class SparseMatrixST
{
  template<typename T2>
  friend class SparseMatrixST;
  friend class SparseMatrixMT<T>;

public:
  /**********************************************************************/
  /**********************************************************************/
  /*                              CREATION                              */
  /**********************************************************************/

  // Empty Matrix with no size
  SparseMatrixST() {create(0, 0);}
  // Squared Matrix with size NxN
  SparseMatrixST(IT n) {create(n, n);}
  // Matrix with N rows and M columns
  SparseMatrixST(IT n, IT m) {create(n, m);}
  // Copy a Matrix
  SparseMatrixST(const SparseMatrixST& matrix) {destroy(); deepCopy(matrix);}
  // Destructor
  ~SparseMatrixST() {destroy();}

  // Reset an existing Matrix object to a square one of NxN
  void configure(IT n)   {
    configure(n, n);
  }

  // Reset an existing Matrix object to one of N rows and M columns
  void configure(IT n, IT m) {
    destroy();
    create(n, m);
  }

  /**********************************************************************/
  /**********************************************************************/
  /*                         Matrix Information                         */
  /**********************************************************************/

  // Number of rows of the Matrix
  IT rowNumber() const {
    return _n;
  }

  // Number of columns of the Matrix
  IT colNumber() const {
    return _m;
  }

  // Number of non Zero elements in the Matrix
  IT nonZeros() const {
    return rowptr[_n];
  }

  // Get the Nth non Zero element of the Matrix
  // the row and column is stored in r and c
  T elementN(IT n, IT* r, IT* c) {
    IT _r = 0;
    while(n >= rowptr[_r+1]) {
      _r++;
    }
    *r = _r;
    *c = colind[n];
    return vals[n];
  }

  // Get the value of the element at the given row and column
  T get(IT row, IT col) const {
    IT i1 = rowptr[row];
    IT i2 = rowptr[row+1];

    if(row==col)
      return diag[row];

    for(IT i=i1; i<i2; i++) {
      IT _c = colind[i];
      if(col == _c)
        return vals[i];
      else if (col < _c)
        break;
    }

    return T();
  }

  // Get the diagonal of the Matrix
  vector<T> diagonal() const {
    return diag;
  }

  // Get the Nth element of the diagonal of the Matrix
  T diagonalElement(IT i) {
    return diag[i];
  }

  // Get a vector with all the elements of a column
  vector<T> column(const IT col) const {
    vector<T> res(_m, T());

    IT index = 0;
    IT nnz = rowptr[_n];
    for(IT i=0; i<nnz; i++) {
      if(colind[i] == col) {

        while(rowptr[index+1] <= i) {
          index++;
        }
        res[index] = vals[i];
        index++;
      }
    }

    return res;
  }

  // Get a vector of position/value pairs
  // with all the non Zero elements of a column
  vector<pair<IT,T>> columnCr(const IT col) const {
    vector<pair<IT,T>> res(_m);
    IT n = 0;
    IT index = 0;
    IT nnz = rowptr[_n];

    for(IT i=0; i<nnz; i++) {
      if(colind[i] == col) {
        while(rowptr[index+1] <= i) {
          index++;
        }
        res[n] = make_pair(index, vals[i]);
        index++;
        n++;
      }
    }
    res.resize(n);
    return res;
  }

  /**********************************************************************/
  /**********************************************************************/
  /*                         Matrix Operations                          */
  /**********************************************************************/

  // Set a value in the Matrix at the given position
  // if the value is Zero it is not inserted, or the
  // previous value is removed
  void set(T val, IT row, IT col) {
    if(row == col)
      diag[row] = val;

    // erase element if null
    if(val == T()) {
      if(get(row, col) != T()) {
        IT istart = rowptr[row];
        IT iend = rowptr[row+1];

        for(IT i=istart; i<iend; i++) {
          IT _c = colind[i];
          if(col == _c) {
            vals.erase(vals.begin()+i);
            colind.erase(colind.begin()+i);
            for(IT i=row+1; i<_n+1; i++) {
              rowptr[i]--;
            }
            return;
          }
        }
      }
      return;
    }

    IT i1 = rowptr[row];
    IT i2 = rowptr[row+1];
    IT _c = 0;
    for(; i1<i2; i1++) {
      _c = colind[i1];
      if(col <= _c) {
        break;
      }
    }

    if(col < _c) {
      vals.insert(vals.begin()+i1, val);
      colind.insert(colind.begin()+i1, col);
      for(IT i=row+1; i<_n+1; i++) {
        rowptr[i]++;
      }
    }
    else if(col == _c) {
      if(i1 == i2) {
        vals.insert(vals.begin()+i1, val);
        colind.insert(colind.begin()+i1, col);
        for(IT i=row+1; i<_n+1; i++) {
          rowptr[i]++;
        }
      }
      else {
        vals[i1] = val;
      }
    }
    else if(col > _c) {
      vals.insert(vals.begin()+i1, val);
      colind.insert(colind.begin()+i1, col);
      for(IT i=row+1; i<_n+1; i++) {
        rowptr[i]++;
      }
    }
  }

  //Compute transpose of the current Matrix
  SparseMatrixST<T> transpose() const {
    SparseMatrixST<T> r(_m, _n);
    IT nnz = rowptr[_n];
    r.colind = vector<IT>(nnz);
    r.vals = vector<T>(nnz);

    // generate new rowptr
    vector<IT> colcount = vector<IT>(_m, 0);
    for(IT i=0; i<nnz; i++) {
      colcount[colind[i]]++;
    }
    r.rowptr[0] = 0;
    for(IT i=0; i<_m; i++) {
      r.rowptr[i+1] = r.rowptr[i] + colcount[i];
    }

    // generate vals and colind
    vector<IT> colsc = vector<IT>(_m, 0);
    for(IT i=0; i<_n; i++) {
      IT jstart = rowptr[i];
      IT jend   = rowptr[i+1];
      for(IT j = jstart; j < jend; j++) {
        IT _c = colind[j];
        T a = vals[j];
        IT _p = r.rowptr[_c] + colsc[_c];
        r.vals[_p] = a;
        r.colind[_p] = i;
        if(i == _c) {
          r.diag[i] = a;
        }
        colsc[_c]++;
      }
    }

    return r;
  }

  //Compute r = A*v for CSR matrix A and dense vector v
  vector<T> multiply(const vector<T> &v) const {
    vector<T> r(_n);

    for(IT i=0; i<_n; ++i) {
      T sum = T();
      IT jstart = rowptr[i];
      IT jend   = rowptr[i+1];
      for(IT j=jstart; j<jend; j++) {
        sum += vals[j] * v[colind[j]];
      }
      r[i] = sum;
    }

    return r;
  }

  //Compute R = A*M for CSR matrices A, M and R
  SparseMatrixST<T> multiply(const SparseMatrixST<T> &m) const {
    SparseMatrixST<T> r(_n, m._m);

    multiplyPhaseOne(m, r);
    IT size = r.rowptr[r._n];
    r.vals = vector<T>(size, T());
    r.colind = vector<IT>(size, 0);
    multiplyPhaseTwo(m, r);
    r.sortColumns();

    return r;
  }

  //Compute r = A*v for CSR matrix A and dense vector v
  vector<T> operator*(const vector<T> &v) const {
    return multiply(v);
  }

  // Matrix of the result of multiplication of the Matrix with
  // another Matrix
  SparseMatrixST<T> operator*(const SparseMatrixST<T> &m) const {
    return multiply(m);
  }

  // Compare two matrices
  friend bool operator==(const SparseMatrixST& l, const SparseMatrixST& r) {
    // Check sizes
    if(l._n != r._n)
      return false;
    if(l._m != r._m)
      return false;

    // Compare containers' content
    if(l.vals != r.vals)
      return false;
    if(l.colind != r.colind)
      return false;
    if(l.rowptr != r.rowptr)
      return false;

    return true;
  }

  // Compare two matrices
  friend bool operator!=(const SparseMatrixST& l, const SparseMatrixST& r) {
    return !(l == r);
  }

  // Copy a matrix
  void operator=(const SparseMatrixST& r) const {
    deepCopy(r);
  }

  // Less comparator
  SparseMatrixST<bool> binOpL(const SparseMatrixST<T> &m) const {
    return binOp<bool>(m, std::less<T>());
  }

  // Less comparator
  SparseMatrixST<bool> operator<(const SparseMatrixST<T> &m) const {
    return binOpL(m);
  }

  // Less-Equal comparator
  SparseMatrixST<bool> binOpLE(const SparseMatrixST<T> &m) const {
    return binOp<bool>(m, std::less_equal<T>());
  }

  // Less-Equal comparator
  SparseMatrixST<bool> operator<=(const SparseMatrixST<T> &m) const {
    return binOpLE(m);
  }

  // Greater comparator
  SparseMatrixST<bool> binOpG(const SparseMatrixST<T> &m) const {
    return binOp<bool>(m, std::greater<T>());
  }

  // Greater comparator
  SparseMatrixST<bool> operator>(const SparseMatrixST<T> &m) const {
    return binOpG(m);
  }

  // Greater-Equal comparator
  SparseMatrixST<bool> binOpGE(const SparseMatrixST<T> &m) const {
    return binOp<bool>(m, std::greater_equal<T>());
  }

  // Greater-Equal comparator
  SparseMatrixST<bool> operator>=(const SparseMatrixST<T> &m) const {
    return binOpGE(m);
  }

  // Matrix with the minimum values of both Matrices
  SparseMatrixST<T> minValues(const SparseMatrixST<T> &m) const {
    return binOp<T>(m, std::min<T>());
  }

  // Matrix with the maximum values of both Matrices
  SparseMatrixST<T> maxValues(const SparseMatrixST<T> &m) const {
    return binOp<T>(m, std::max<T>());
  }

  // Division of the elements (if divisor is 0, then element is set to 0)
  SparseMatrixST<T> binOpDiv(const SparseMatrixST<T> &m) const {
    return binOp<T>(m, sdiv());
  }

  // Division of the elements (if divisor is 0, then element is set to 0)
  SparseMatrixST<T> operator/(const SparseMatrixST<T> &m) const {
    return binOpDiv(m);
  }

  // Multiplication of value of the different elements
  SparseMatrixST<T> binOpMul(const SparseMatrixST<T> &m) const {
    return binOp<T>(m, std::multiplies<T>());
  }

  // Addition of values of the different elements
  SparseMatrixST<T> binOpAdd(const SparseMatrixST<T> &m) const {
    return binOp<T>(m, std::plus<T>());
  }

  // Addition of values of the different elements
  SparseMatrixST<T> operator+(const SparseMatrixST<T> &m) const {
    return binOpAdd(m);
  }

  // Subtraction of values of the different elements
  SparseMatrixST<T> binOpSub(const SparseMatrixST<T> &m) const {
    return binOp<T>(m, std::minus<T>());
  }

  // Subtraction of values of the different elements
  SparseMatrixST<T> operator-(const SparseMatrixST<T> &m) const {
    return binOpSub(m);
  }

  // Matrix with true values for all non-different elements in a Matrix
  SparseMatrixST<bool> binOpNE(const SparseMatrixST<T> &m) const {
    return binOp<bool>(m, std::not_equal_to<T>());
  }

  /**********************************************************************/
  /**********************************************************************/
  /*                        Matrix I/O operations                       */
  /**********************************************************************/

  // Convert the Matrix into a string
  string toStringList() const {
    string s;

    s = to_string(_n) + " " + to_string(_m) + " " + to_string(rowptr[_n]) +"\n";
    for(IT i=0; i<_n; i++) {
      IT jstart = rowptr[i];
      IT jend   = rowptr[i+1];
      for(IT j = jstart; j < jend; j++) {
        T v = vals[j];
        IT c = colind[j];
        s += to_string(i) + " " + to_string(c) + " " + to_string(v) + "\n";
      }
    }

    return s;
  }

  // Load a Matrix from a string
  // TODO: only working with float/double types
  void fromStringList(string s) {
    string line;
    bool first = true;
    istringstream f(s);
    destroy();
    while(std::getline(f, line)) {
      vector<std::string> v = split(line, ' ');
      if(first) {
        first = false;
        create(stoi(v[0]), stoi(v[1]));
      } else {
        stringstream convert(v[2]);
        T value;
        convert >> value;
        set(value, stoi(v[0]), stoi(v[1]));
      }
    }
  }

  // Stream of all the values of the Matrix
  friend ostream& operator<< (ostream& stream, const SparseMatrixST<T> &m) {
    IT nrow = m._n;
    IT ncol = m._m;
    for(IT i = 0; i<nrow; i++) {
      for(IT j = 0; j<ncol; j++) {
        stream << m.get(i, j) << " ";
      }
      if(i!=nrow-1)
        stream << endl;
    }
    return stream;
  }

private:
  template <typename T2, typename binary_op>
  SparseMatrixST<T2> binOp(const SparseMatrixST<T> &m, const binary_op& op) const {
    SparseMatrixST<T2> r(_n, _m);
    IT nrow = _n;

    // Allocate maximum amount of needed memory,
    // later resize the vectors
    IT size = rowptr[_n] + m.rowptr[_n];
    r.vals = vector<T2>(size, T2());
    r.colind = vector<IT>(size, 0);

    r.rowptr[0] = 0;
    IT nnz = 0;

    for(IT i = 0; i < nrow; i++) {
      IT A_pos = rowptr[i];
      IT B_pos = m.rowptr[i];
      IT A_end = rowptr[i+1];
      IT B_end = m.rowptr[i+1];

      //while not finished with either row
      while(A_pos < A_end && B_pos < B_end){
        IT A_j = colind[A_pos];
        IT B_j = m.colind[B_pos];

        if(A_j == B_j) {
          T result = op(vals[A_pos], m.vals[B_pos]);
          if(result != 0) {
            r.colind[nnz] = A_j;
            r.vals[nnz] = result;
            nnz++;
            if(i == A_j)
              r.diag[i] = result;
          }
          A_pos++;
          B_pos++;
        } else if (A_j < B_j) {
          T result = op(vals[A_pos], 0);
          if (result != 0) {
            r.colind[nnz] = A_j;
            r.vals[nnz] = result;
            nnz++;
            if(i == A_j)
              r.diag[i] = result;
          }
          A_pos++;
        } else { //B_j < A_j
          T result = op(0, m.vals[B_pos]);
          if (result != 0){
            r.colind[nnz] = B_j;
            r.vals[nnz] = result;
            nnz++;
            if(i == B_j)
              r.diag[i] = result;
          }
          B_pos++;
        }
      }

      //tail
      while(A_pos < A_end) {
        T result = op(vals[A_pos], 0);
        if (result != 0){
          r.colind[nnz] = colind[A_pos];
          r.vals[nnz] = result;
          nnz++;
          if(i == colind[A_pos])
            r.diag[i] = result;
        }
        A_pos++;
      }
      while(B_pos < B_end) {
        T result = op(0, m.vals[B_pos]);
        if (result != 0) {
          r.colind[nnz] = m.colind[B_pos];
          r.vals[nnz] = result;
          nnz++;
          if(i == m.colind[B_pos])
            r.diag[i] = result;
        }
        B_pos++;
      }

      r.rowptr[i+1] = nnz;
    }

    r.vals.resize(r.rowptr[r._n]);
    r.colind.resize(r.rowptr[r._n]);

    return r;
  }

  // Calculte the RowPtr values
  void multiplyPhaseOne(const SparseMatrixST<T> &m, SparseMatrixST<T> &r) const {
    IT ncol = r._m;
    IT nrow = r._n;
    vector<int> mask(ncol, -1);
    r.rowptr[0] = 0;
    IT nnz = 0;

    for(IT i = 0; i< nrow; i++) {
      IT row_nnz = 0;
      IT jstart = rowptr[i];
      IT jend   = rowptr[i+1];
      for(IT j=jstart; j<jend; j++) {
        IT c1 = colind[j];

        IT kstart = m.rowptr[c1];
        IT kend   = m.rowptr[c1+1];
        for(IT k = kstart; k < kend; k++) {
          IT c2 = m.colind[k];
          if(mask[c2] != i) {
            mask[c2] = i;
            row_nnz++;
          }
        }
      }

      IT next_nnz = nnz + row_nnz;
      nnz = next_nnz;
      r.rowptr[i+1] = nnz;
    }
  }

  // Calculte the RowPtr values
  void multiplyPhaseTwo(const SparseMatrixST<T> &m, SparseMatrixST<T> &r) const {
    IT ncol = r._m;
    IT nrow = r._n;
    vector<IT> next(ncol, -1);
    vector<T> sums(ncol, T());

    IT nnz = 0;
    r.rowptr[0] = 0;

    for(IT i = 0; i < nrow; i++) {
      IT head   = -2;
      IT length =  0;

      IT jstart = rowptr[i];
      IT jend   = rowptr[i+1];
      for(IT j = jstart; j < jend; j++) {
        IT c1 = colind[j];
        T v = vals[j];

        IT kstart = m.rowptr[c1];
        IT kend   = m.rowptr[c1+1];
        for(IT k = kstart; k < kend; k++) {
          IT c2 = m.colind[k];
          sums[c2] += v * m.vals[k];
          if(next[c2] == -1) {
            next[c2] = head;
            head = c2;
            length++;
          }
        }
      }

      for(IT j = 0; j < length; j++) {
        if(sums[head] != 0) {
          r.colind[nnz] = head;
          r.vals[nnz] = sums[head];
          nnz++;
          if(i == head)
            r.diag[i] = sums[head];
        }

        IT temp = head;
        head = next[head];
        next[temp] = -1;
        sums[temp] =  0;
      }

      r.rowptr[i+1] = nnz;
    }
  }

  void sortColumns() {
    IT nrow = _n;
    vector<pair<IT,T>> temp;

    for(IT i = 0; i < nrow; i++) {
      IT rstart = rowptr[i];
      IT rend   = rowptr[i+1];

      temp.resize(rend - rstart);
      for(IT j = rstart, n = 0; j < rend; j++, n++) {
        temp[n].first  = colind[j];
        temp[n].second = vals[j];
      }

      sort(temp.begin(), temp.end(),
           [](const pair<IT,T>& a, const pair<IT,T>& b)
           {return a.first < b.first;});

      for(IT j = rstart, n = 0; j < rend; j++, n++) {
        colind[j] = temp[n].first;
        vals[j]   = temp[n].second;
      }
    }
  }

  // full copy of all the parameters of the Matrix
  void deepCopy(const SparseMatrixST& matrix) {
    _m = matrix._m;
    _n = matrix._n;
    vals = matrix.vals;
    colind = matrix.colind;
    rowptr = matrix.rowptr;
    diag = matrix.diag;
  }

  // basic initiallization of the Matrix
  void create(IT n, IT m) {
    _n = n;
    _m = m;
    rowptr = vector<IT>(n+1, 0);

    IT d = n;
    if(m < n)
      d = m;
    diag = vector<T>(d, T());
  }

  // destroy Matrix
  void destroy() {
    _n = 0;
    _m = 0;

    vals.clear();
    rowptr.clear();
    colind.clear();
    diag.clear();
  }

  // returns 0 in case of division by 0
  // used for the division binary operation
  struct sdiv {
    T operator() (const T& x, const T& y) const {
      if (y == T())
        return T();
      else
        return x/y;
    }
  };

  // split a string into vector of string divided
  // by the given character
  const vector<string> split(const string& s, const char& c) {
    string buff{""};
    vector<string> v;

    for(auto n:s) {
      if(n != c)
        buff+=n;
      else if(n == c && buff != "") {
        v.push_back(buff);
        buff = "";
      }
    }
    if(buff != "")
      v.push_back(buff);

    return v;
  }

  IT _n, _m; // rows, columns
  vector<T> vals; // values
  vector<T> diag; // diagonal
  vector<IT> colind; // column indexes
  vector<IT> rowptr; // row start/end pointers
};

#endif // SPARSEMATRIXST_H_
