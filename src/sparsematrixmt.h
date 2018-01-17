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

// Date: 20180117
// Version: 0.2

#ifndef SPARSEMATRIXMT_H_
#define SPARSEMATRIXMT_H_

#include <omp.h>
#include <vector>
#include <utility>
#include <algorithm>
#include <functional>
#include <ostream>
#include <string>
#include <sstream>

#ifndef SM_NO_ST
#include "sparsematrixst.h"
#endif

// Integer definition, it should be signed!
#define IT int32_t

using namespace std;

template<typename T>
class SparseMatrixMT
{
  template<typename T2>
  friend class SparseMatrixMT;

public:
  /**********************************************************************/
  /**********************************************************************/
  /*                              CREATION                              */
  /**********************************************************************/

  // Empty Matrix with no size
  SparseMatrixMT() {create(0, 0, 0);}
  // Squared Matrix with size NxN divided in chunks
  SparseMatrixMT(IT n, int chunks=0) {create(n, n, chunks);}
  // Matrix with N rows and M columns divided in chunks
  SparseMatrixMT(IT n, IT m, int chunks=0) {create(n, m, chunks);}
  // Copy a Matrix
  SparseMatrixMT(const SparseMatrixMT& matrix) {destroy(); deepCopy(matrix);}

  // Copy a Single-Threaded Matrix into a Multi-Threaded Matrix
  // divided in chunks
#ifndef SM_NO_ST
  SparseMatrixMT(SparseMatrixST<T> m, IT chunks=0) {
    IT r, c;
    IT nnz = m.nonZeros();
    create(m.rowNumber(), m.colNumber(), chunks);
    for(IT i = 0; i<nnz; i++) {
      T v = m.elementN(i, &r, &c);
      set(v, r, c);
    }
  }
#endif

  // Copy a Multi-Threaded Matrix into a Single-Threaded Matrix
#ifndef SM_NO_ST
  SparseMatrixST<T> toST() const {
    SparseMatrixST<T> r(_n, _m);
    IT row, col;
    IT nnz = nonZeros();
    r.colind = vector<IT>(0);
    r.vals   = vector<T>(T());
    for(IT i=0; i<nnz; i++) {
      T v = elementN(i, &row, &col);
      r.set(v, row, col);
    }
    return r;
  }
#endif

  // Destructor
  ~SparseMatrixMT() {destroy();}

  // Reset an existing Matrix object to a square one of NxN
  // divided in chunks
  void configure(IT n, IT chunks)   {
    configure(n, n, chunks);
  }

  // Reset an existing Matrix object to one of N rows and M columns
  // divided in chunks
  void configure(IT n, IT m, IT chunks) {
    destroy();
    create(n, m, chunks);
  }

  /**********************************************************************/
  /**********************************************************************/
  /*                         Matrix Information                         */
  /**********************************************************************/

  // Number of rows of the Matrix
  IT rowNumber() const {
    return _n;
  }

  // Number of rows of the Matrix of the specified chunk
  IT rowNumberChunk(IT chunk) const {
    return _chunkRows[chunk];
  }

  // Number of columns of the Matrix
  IT colNumber() const {
    return _m;
  }

  // Number of division chunks
  IT chunksNumber() const {
    return _chunks;
  }

  // Number of non Zero elements in the Matrix
  IT nonZeros() const {
    IT nnz = 0;
    for(IT i=0; i<_chunks; i++) {
      nnz += rowptr[i][_chunkRows[i]];
    }
    return nnz;
  }

  // Number of non Zero elements in the Matrix in the specified chunk
  IT nonZerosChunk(IT chunk) const {
    return rowptr[chunk][_chunkRows[chunk]];
  }

  // Get the Nth non Zero element of the Matrix
  // the row and column is stored in r and c
  T elementN(IT n, IT* r, IT* c) const {
    IT chunk = 0;
    IT ns = rowptr[chunk][_chunkRows[chunk]];

    while(n >= ns) {
      chunk++;
      ns += rowptr[chunk][_chunkRows[chunk]];
    }

    IT localn = n - ns + rowptr[chunk][_chunkRows[chunk]];
    T v = elementNChunk(chunk, localn, r, c);
    *r = *r + _chunkRowsStart[chunk];

    return v;
  }

  // Get the Nth non Zero element of the specified chunk
  // the row and column is stored in r and c
  T elementNChunk(IT chunk, IT n, IT* r, IT* c) const {
    IT _r = 0;
    while(n >= rowptr[chunk][_r+1]) {
      _r++;
    }
    *r = _r;
    *c = colind[chunk][n];
    return vals[chunk][n];
  }

  // Get the value of the element at the given row and column
  T get(IT row, IT col) const {
    pair<IT,IT> chunklrow = _chunkOfRow[row];
    return getLocal(chunklrow.first, chunklrow.second, col);
  }

  // Get the value of the element at the given local row and column
  // of the specified chunk
  T getLocal(IT chunk, IT lrow, IT col) const {
    IT i1 = rowptr[chunk][lrow];
    IT i2 = rowptr[chunk][lrow+1];

    for(IT i=i1; i<i2; i++) {
      IT _c = colind[chunk][i];
      if(col == _c)
        return vals[chunk][i];
      else if (col < _c)
        break;
    }

    return T();
  }

  // Get a vector with all the elements of a column
  vector<T> column(const IT col) const {
    vector<T> res(_n, T());

#pragma omp parallel for
    for(IT cchunk = 0; cchunk < _chunks; cchunk++) {
      IT index = 0;
      IT nnz = rowptr[cchunk][_chunkRows[cchunk]];

      for(IT i=0; i<nnz; i++) {
        if(colind[cchunk][i] == col) {

          while(rowptr[cchunk][index+1] <= i) {
            index++;
          }
          res[index+_chunkRowsStart[cchunk]] = vals[cchunk][i];
          index++;
        }
      }
    }

    return res;
  }

  // Get a vector of position/value pairs
  // with all the non Zero elements of a column
  vector<pair<IT,T>> columnCr(const IT col) const {
    vector<pair<IT,T>> res(_n);
    IT n = 0;

    for(IT cchunk = 0; cchunk < _chunks; cchunk++) {
      IT index = 0;
      IT nnz = rowptr[cchunk][_chunkRows[cchunk]];

      for(IT i=0; i<nnz; i++) {
        if(colind[cchunk][i] == col) {
          while(rowptr[cchunk][index+1] <= i) {
            index++;
          }
          res[n] = make_pair(index+_chunkRowsStart[cchunk], vals[cchunk][i]);
          index++;
          n++;
        }
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
    pair<IT,IT> chunklrow = _chunkOfRow[row];
    setLocal(chunklrow.first, val, chunklrow.second, col);
  }

  // Set a value in the selected Matrix chunk at the given local
  // position if the value is Zero it is not inserted, or the
  // previous value is removed
  void setLocal(IT chunk, T val, IT lrow, IT col) {
    // erase element if null
    if(val == T() && getLocal(chunk, lrow, col) != T()) {
      IT istart = rowptr[chunk][lrow];
      IT iend = rowptr[chunk][lrow+1];

      for(IT i=istart; i<iend; i++) {
        IT _c = colind[chunk][i];
        if(col == _c) {
          vals[chunk].erase(vals[chunk].begin()+i);
          colind[chunk].erase(colind[chunk].begin()+i);
          for(IT i=lrow+1; i<_chunkRows[chunk]+1; i++) {
            rowptr[chunk][i]--;
          }
          return;
        }
      }
    }

    IT i1 = rowptr[chunk][lrow];
    IT i2 = rowptr[chunk][lrow+1];
    IT _c = 0;
    for(; i1<i2; i1++) {
      _c = colind[chunk][i1];
      if(col <= _c) {
        break;
      }
    }

    if(col < _c) {
      vals[chunk].insert(vals[chunk].begin()+i1, val);
      colind[chunk].insert(colind[chunk].begin()+i1, col);
      for(IT i=lrow+1; i<_chunkRows[chunk]+1; i++) {
        rowptr[chunk][i]++;
      }
    }
    else if(col == _c) {
      if(i1 == i2) {
        vals[chunk].insert(vals[chunk].begin()+i1, val);
        colind[chunk].insert(colind[chunk].begin()+i1, col);
        for(IT i=lrow+1; i<_chunkRows[chunk]+1; i++) {
          rowptr[chunk][i]++;
        }
      }
      else {
        vals[chunk][i1] = val;
      }
    }
    else if(col > _c) {
      vals[chunk].insert(vals[chunk].begin()+i1, val);
      colind[chunk].insert(colind[chunk].begin()+i1, col);
      for(IT i=lrow+1; i<_chunkRows[chunk]+1; i++) {
        rowptr[chunk][i]++;
      }
    }
  }

  //Compute the transpose of the current Matrix
  SparseMatrixMT<T> transpose(IT chunks=0) const {
    if(chunks==0) {
      chunks = _chunks;
    }
    SparseMatrixMT<T> r(_m, _n, chunks);

    vector<IT> colcount = vector<IT>(_m, 0);
    for(IT c=0; c<_chunks; c++) {
      IT nnz = rowptr[c][_chunkRows[c]];
      for(IT i=0; i<nnz; i++) {
        colcount[colind[c][i]]++;
      }
    }
    for(IT c=0; c<r._chunks; c++) {
      r.rowptr[c][0] = 0;
      IT rowschunk = r._chunkRows[c];
      for(IT i=0; i<rowschunk; i++) {
        r.rowptr[c][i+1] = r.rowptr[c][i] + colcount[r._chunkRowsStart[c]+i];
      }
      IT nnz = r.rowptr[c][r._chunkRows[c]];
      r.colind[c] = vector<IT>(nnz);
      r.vals[c]   = vector<T>(nnz);
    }

    // generate vals and colind
    vector<IT> colsc = vector<IT>(_m, 0);
    for(IT c=0; c<_chunks; c++) {
      for(IT i=0; i<_chunkRows[c]; i++) {
        IT jstart = rowptr[c][i];
        IT jend   = rowptr[c][i+1];
        for(IT j = jstart; j < jend; j++) {
          IT _c = colind[c][j];
          T a = vals[c][j];
          pair<IT,IT> cr = r._chunkOfRow[_c];
          IT _p = r.rowptr[cr.first][cr.second];
          _p += colsc[_c];
          r.vals[cr.first][_p] = a;
          r.colind[cr.first][_p] = i + _chunkRowsStart[c];
          colsc[_c]++;
        }
      }
    }

    return r;
  }

  //Compute r = A*v for CSR matrix A and dense vectors v
  vector<T> multiply(const vector<T> &v) const {
    vector<T> r(_n);

#pragma omp parallel for
    for(IT cchunk = 0; cchunk < _chunks; cchunk++) {
      IT nrow = _chunkRows[cchunk];
      for(IT i = 0; i < nrow; i++) {
        T sum = T();
        IT jstart = rowptr[cchunk][i];
        IT jend   = rowptr[cchunk][i+1];
        for(IT j=jstart; j<jend; j++) {
          sum += vals[cchunk][j] * v[colind[cchunk][j]];
        }

        r[_chunkRowsStart[cchunk] + i] = sum;
      }
    }

    return r;
  }

  //Compute R = A*M for CSR matrices A, M and R
  SparseMatrixMT<T> multiply(const SparseMatrixMT<T> &m) const {
    SparseMatrixMT<T> r(_n, m._m, _chunks);

#pragma omp parallel for
    for(IT cchunk = 0; cchunk < _chunks; cchunk++) {
      multiplyFastPhaseOne(cchunk, m, r);
      IT size = r.rowptr[cchunk][r._chunkRows[cchunk]];
      r.vals[cchunk] = vector<T>(size, T());
      r.colind[cchunk] = vector<IT>(size, 0);
      multiplyFastPhaseTwo(cchunk, m, r);
      r.sortColumns(cchunk);
    }

    return r;
  }

  //Compute r = A*v for CSR matrix A and dense vector v
  vector<T> operator*(const vector<T> &v) const {
    return multiply(v);
  }

  // Matrix of the result of multiplication of the Matrix with
  // another Matrix
  SparseMatrixMT<T> operator*(const SparseMatrixMT<T> &m) const {
    return multiply(m);
  }

  // Compare two matrices, true if they are equal (including partitioning)
  friend bool operator==(const SparseMatrixMT& l, const SparseMatrixMT& r) {
    // Check sizes
    if(l._n != r._n)
      return false;
    if(l._m != r._m)
      return false;
    if(l._chunks != r._chunks)
      return false;

    // Check data
    if(l.vals != r.vals)
      return false;
    if(l.colind != r.colind)
      return false;
    if(l.rowptr != r.rowptr)
      return false;

    // All is equal
    return true;
  }

  // Compare two matrices, true if they are different
  // or have different partitioning
  friend bool operator!=(const SparseMatrixMT& l, const SparseMatrixMT& r) {
    return !(l == r);
  }

  // Copy a matrix
  void operator=(const SparseMatrixMT& r) const {
    deepCopy(r);
  }

  // Less comparator
  SparseMatrixMT<bool> binOpL(const SparseMatrixMT<T> &m) const {
    return binOp<bool>(m, std::less<T>());
  }

  // Less comparator
  SparseMatrixMT<bool> operator<(const SparseMatrixMT<T> &m) const {
    return this->binOpL(m);
  }

  // Less-Equal comparator
  SparseMatrixMT<bool> binOpLE(const SparseMatrixMT<T> &m) const {
    return binOp<bool>(m, std::less_equal<T>());
  }

  // Less-Equal comparator
  SparseMatrixMT<bool> operator<=(const SparseMatrixMT<T> &m) const {
    return binOpLE(m);
  }

  // Greater comparator
  SparseMatrixMT<bool> binOpG(const SparseMatrixMT<T> &m) const {
    return binOp<bool>(m, std::greater<T>());
  }

  // Greater comparator
  SparseMatrixMT<bool> operator>(const SparseMatrixMT<T> &m) const {
    return binOpG(m);
  }

  // Greater-Equal comparator
  SparseMatrixMT<bool> binOpGE(const SparseMatrixMT<T> &m) const {
    return binOp<bool>(m, std::greater_equal<T>());
  }

  // Greater-Equal comparator
  SparseMatrixMT<bool> operator>=(const SparseMatrixMT<T> &m) const {
    return binOpGE(m);
  }

  // Matrix with the minimum values of both Matrices
  SparseMatrixMT<T> minValues(const SparseMatrixMT<T> &m) const {
    return binOp<T>(m, std::min<T>());
  }

  // Matrix with the maximum values of both Matrices
  SparseMatrixMT<T> maxValues(const SparseMatrixMT<T> &m) const {
    return binOp<T>(m, std::max<T>());
  }

  // Division of the elements (if divisor is 0, then element is set to 0)
  SparseMatrixMT<T> binOpDiv(const SparseMatrixMT<T> &m) const {
    return binOp<T>(m, sdiv());
  }

  // Division of the elements (if divisor is 0, then element is set to 0)
  SparseMatrixMT<T> operator/(const SparseMatrixMT<T> &m) const {
    return binOpDiv(m);
  }

  // Multiplication of value of the different elements
  SparseMatrixMT<T> binOpMul(const SparseMatrixMT<T> &m) const {
    return binOp<T>(m, std::multiplies<T>());
  }

  // Addition of values of the different elements
  SparseMatrixMT<T> binOpAdd(const SparseMatrixMT<T> &m) const {
    return binOp<T>(m, std::plus<T>());
  }

  // Addition of values of the different elements
  SparseMatrixMT<T> operator+(const SparseMatrixMT<T> &m) const {
    return binOpAdd(m);
  }

  // Subtraction of values of the different elements
  SparseMatrixMT<T> binOpSub(const SparseMatrixMT<T> &m) const {
    return binOp<T>(m, std::minus<T>());
  }

  // Subtraction of values of the different elements
  SparseMatrixMT<T> operator-(const SparseMatrixMT<T> &m) const {
    return binOpSub(m);
  }

  // Matrix with true values for all non-different elements in a Matrix
  SparseMatrixMT<bool> binOpNE(const SparseMatrixMT<T> &m) const {
    return binOp<bool>(m, std::not_equal_to<T>());
  }

  /**********************************************************************/
  /**********************************************************************/
  /*                        Matrix I/O operations                       */
  /**********************************************************************/

  // Stream of all the values of the Matrix
  friend ostream& operator<< (ostream& stream, const SparseMatrixMT<T> &m) {
    for(IT cchunk = 0; cchunk < m._chunks; cchunk++) {
      IT nrow = m._chunkRows[cchunk];
      IT ncol = m._m;
      for(IT i = 0; i<nrow; i++) {
        for(IT j = 0; j<ncol; j++) {
          stream << m.get(m._chunkRowsStart[cchunk]+i, j) << " ";
        }
        if(i!=nrow-1 || cchunk !=m._chunks-1)
          stream << endl;
      }
    }
    return stream;
  }

private:
  template <typename T2, typename binary_op>
  SparseMatrixMT<T2> binOp(const SparseMatrixMT<T> &m, const binary_op& op) const {
    SparseMatrixMT<T2> r(_n, _m, _chunks);

#pragma omp parallel for
    for(IT cchunk = 0; cchunk < _chunks; cchunk++) {
      IT nrow = _chunkRows[cchunk];
      IT size = rowptr[cchunk][nrow] + m.rowptr[cchunk][nrow];
      r.vals[cchunk] = vector<T2>(size, T2());
      r.colind[cchunk] = vector<IT>(size, 0);
      r.rowptr[cchunk][0] = 0;
      IT nnz = 0;

      for(IT i = 0; i < nrow; i++) {
        IT A_pos = rowptr[cchunk][i];
        IT B_pos = m.rowptr[cchunk][i];
        IT A_end = rowptr[cchunk][i+1];
        IT B_end = m.rowptr[cchunk][i+1];

        //while not finished with either row
        while(A_pos < A_end && B_pos < B_end){
          IT A_j = colind[cchunk][A_pos];
          IT B_j = m.colind[cchunk][B_pos];

          if(A_j == B_j) {
            T result = op(vals[cchunk][A_pos], m.vals[cchunk][B_pos]);
            if(result != 0) {
              r.colind[cchunk][nnz] = A_j;
              r.vals[cchunk][nnz] = result;
              nnz++;
            }
            A_pos++;
            B_pos++;
          } else if (A_j < B_j) {
            T result = op(vals[cchunk][A_pos], 0);
            if (result != 0) {
              r.colind[cchunk][nnz] = A_j;
              r.vals[cchunk][nnz] = result;
              nnz++;
            }
            A_pos++;
          } else { //B_j < A_j
            T result = op(0, m.vals[cchunk][B_pos]);
            if (result != 0){
              r.colind[cchunk][nnz] = B_j;
              r.vals[cchunk][nnz] = result;
              nnz++;
            }
            B_pos++;
          }
        }

        //tail
        while(A_pos < A_end) {
          T result = op(vals[cchunk][A_pos], 0);
          if (result != 0){
            r.colind[cchunk][nnz] = colind[cchunk][A_pos];
            r.vals[cchunk][nnz] = result;
            nnz++;
          }
          A_pos++;
        }
        while(B_pos < B_end) {
          T result = op(0, m.vals[cchunk][B_pos]);
          if (result != 0) {
            r.colind[cchunk][nnz] = m.colind[cchunk][B_pos];
            r.vals[cchunk][nnz] = result;
            nnz++;
          }
          B_pos++;
        }

        r.rowptr[cchunk][i+1] = nnz;
      }

      r.vals[cchunk].resize(r.rowptr[cchunk][r._chunkRows[cchunk]]);
      r.colind[cchunk].resize(r.rowptr[cchunk][r._chunkRows[cchunk]]);
    }
    return r;
  }

  // Calculte the RowPtr values
  void multiplyFastPhaseOne(IT cchunk, const SparseMatrixMT<T> &m, SparseMatrixMT<T> &r) const {
    IT ncol = r._m;
    IT nrow = r._chunkRows[cchunk];
    vector<int> mask(ncol, -1);
    r.rowptr[cchunk][0] = 0;
    IT nnz = 0;

    for(IT i = 0; i< nrow; i++) {
      IT row_nnz = 0;
      IT jstart = rowptr[cchunk][i];
      IT jend   = rowptr[cchunk][i+1];
      for(IT j=jstart; j<jend; j++) {
        IT c1 = colind[cchunk][j];

        IT kstart = m.getGenRPat(c1);
        IT kend   = m.getGenRPat(c1+1);
        for(IT k = kstart; k < kend; k++) {
          IT c2 = m.getGenCIat(k);
          if(mask[c2] != i) {
            mask[c2] = i;
            row_nnz++;
          }
        }
      }

      IT next_nnz = nnz + row_nnz;
      nnz = next_nnz;
      r.rowptr[cchunk][i+1] = nnz;
    }
  }

  void multiplyFastPhaseTwo(IT cchunk, const SparseMatrixMT<T> &m,SparseMatrixMT<T> &r) const {
    IT ncol = r._m;
    IT nrow = r._chunkRows[cchunk];
    vector<IT> next(ncol, -1);
    vector<T> sums(ncol, T());

    IT nnz = 0;
    r.rowptr[cchunk][0] = 0;

    for(IT i = 0; i < nrow; i++) {
      IT head   = -2;
      IT length =  0;

      IT jstart = rowptr[cchunk][i];
      IT jend   = rowptr[cchunk][i+1];
      for(IT j = jstart; j < jend; j++) {
        IT c1 = colind[cchunk][j];
        T v = vals[cchunk][j];

        IT kstart = m.getGenRPat(c1);
        IT kend   = m.getGenRPat(c1+1);
        for(IT k = kstart; k < kend; k++) {
          IT c2 = m.getGenCIat(k);
          sums[c2] += v * m.getGenValat(k);

          if(next[c2] == -1) {
            next[c2] = head;
            head = c2;
            length++;
          }
        }
      }

      for(IT j = 0; j < length; j++) {
        if(sums[head] != 0) {
          r.colind[cchunk][nnz] = head;
          r.vals[cchunk][nnz] = sums[head];
          nnz++;
        }

        IT temp = head;
        head = next[head];
        next[temp] = -1;
        sums[temp] =  0;
      }

      r.rowptr[cchunk][i+1] = nnz;
    }
  }

  IT getGenCIat(IT pos) const {
    IT v = 0;
    IT ind = 0;

    for(IT i=0; i<_chunks; i++) {
      if(pos < ind+(IT)colind[i].size()) {
        IT p = pos - ind;
        return colind[i][p];
      }
      ind += colind[i].size();
    }

    return v;
  }

  T getGenValat(IT pos) const {
    T v = 0;
    IT ind = 0;

    for(IT i=0; i<_chunks; i++) {
      if(pos < ind+(IT)vals[i].size()) {
        IT p = pos - ind;
        T v = vals[i][p];
        return v;
      }
      ind += vals[i].size();
    }

    return v;
  }

  IT getGenRPat(IT pos) const {
    IT v = 0;

    if(pos==0)
      return v;

    for(IT i=0; i<_chunks; i++) {
      for(IT j=1; j<=_chunkRows[i]; j++) {
        v += rowptr[i][j] - rowptr[i][j-1];
        pos--;
        if(pos==0) {
          return v;
        }
      }
    }

    return v;
  }

  void sortColumns(IT cchunk) {
    IT nrow = _chunkRows[cchunk];
    vector<pair<IT,T>> temp;

    for(IT i = 0; i < nrow; i++) {
      IT rstart = rowptr[cchunk][i];
      IT rend   = rowptr[cchunk][i+1];

      temp.resize(rend - rstart);
      for(IT j = rstart, n = 0; j < rend; j++, n++) {
        temp[n].first  = colind[cchunk][j];
        temp[n].second = vals[cchunk][j];
      }

      sort(temp.begin(), temp.end(),
           [](const pair<IT,T>& a, const pair<IT,T>& b)
           {return a.first < b.first;});

      for(IT j = rstart, n = 0; j < rend; j++, n++) {
        colind[cchunk][j] = temp[n].first;
        vals[cchunk][j]   = temp[n].second;
      }
    }
  }

  // full copy of all the parameters of the Matrix
  void deepCopy(const SparseMatrixMT& matrix) {
    _m = matrix._m;
    _n = matrix._n;
    _chunks = matrix._chunks;

    _chunkRows = matrix._chunkRows;
    _chunkRowsStart = matrix._chunkRowsStart;
    _chunkRowsEnd = matrix._chunkRowsEnd;
    _chunkOfRow = matrix._chunkOfRow;

    vals = matrix.vals;
    colind = matrix.colind;
    rowptr = matrix.rowptr;
  }

  // basic initiallization of the Matrix
  void create(IT n, IT m, IT chunks) {
    _n = n;
    _m = m;
    if(chunks<=0)
      chunks = omp_get_max_threads();
    _chunks = chunks;

    IT rows = _n / chunks;
    _chunkRows = vector<IT>(chunks, rows);
    _chunkRowsStart = vector<IT>(chunks);
    _chunkRowsEnd = vector<IT>(chunks);
    _chunkOfRow = vector<pair <int,int>>(_n);

    IT left = _n - (chunks*rows);
    for(IT i=0; i<left; i++) {
      _chunkRows[i]++;
    }
    IT rc = 0;
    for(IT i=0; i<_chunks; i++) {
      _chunkRowsStart[i] = rc;
      rc += _chunkRows[i];
      _chunkRowsEnd[i] = rc-1;

      IT nlr=0;
      for(IT j=_chunkRowsStart[i]; j<=_chunkRowsEnd[i]; j++) {
        _chunkOfRow[j] = make_pair(i,nlr);
        nlr++;
      }
    }

    rowptr = vector<vector<IT>>(chunks);
    for(IT i=0; i<_chunks; i++) {
      rowptr[i] = vector<IT>(_chunkRows[i]+1, 0);
    }

    vals = vector<vector<T>>(chunks);
    colind = vector<vector<IT>>(chunks);
  }

  // destroy Matrix
  void destroy() {
    _n = 0;
    _m = 0;
    _chunks = 0;

    _chunkRows.clear();
    _chunkRowsStart.clear();
    _chunkRowsEnd.clear();
    _chunkOfRow.clear();

    vals.clear();
    colind.clear();
    rowptr.clear();
  }

  // returns 0 in case of division by 0
  struct sdiv {
    T operator() (const T& x, const T& y) const {
      if (y == T())
        return T();
      else
        return x/y;
    }
  };

  IT _n, _m; //rows, columns
  IT _chunks; // parallel divisions
  vector<IT> _chunkRows;
  vector<IT> _chunkRowsStart;
  vector<IT> _chunkRowsEnd;
  vector<pair<IT,IT>> _chunkOfRow;

  vector<vector<T>> vals;
  vector<vector<IT>> colind;
  vector<vector<IT>> rowptr;
};

#endif // SPARSEMATRIXMT_H_
