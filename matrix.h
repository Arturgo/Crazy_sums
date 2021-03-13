#pragma once
#include <vector>
using namespace std;

// We work with row vectors

// TODO : adapt the algorithm to sparse matrix

template<typename T>
vector<T> operator + (const vector<T>& a, const vector<T>& b) {
   vector<T> res(b.size(), T(0));
   for(size_t iCoeff = 0;iCoeff < b.size();iCoeff++) {
      res[iCoeff] = a[iCoeff] + b[iCoeff];
   }
   return res;
}

template<typename T>
vector<T> operator - (const vector<T>& a, const vector<T>& b) {
   vector<T> res(b.size(), T(0));
   for(size_t iCoeff = 0;iCoeff < b.size();iCoeff++) {
      res[iCoeff] = a[iCoeff] - b[iCoeff];
   }
   return res;
}

template<typename T>
vector<T> operator * (T a, const vector<T>& b) {
   vector<T> res(b.size(), T(0));
   for(size_t iCoeff = 0;iCoeff < b.size();iCoeff++) {
      res[iCoeff] = a * b[iCoeff];
   }
   return res;
}

template<typename T>
T operator * (const vector<T>& a, const vector<T>& b) {
   T res(0);
   
   for(size_t iCoeff = 0;iCoeff < b.size();iCoeff++) {
      res = res + a[iCoeff] * b[iCoeff];
   }
   
   return res;
}

template<typename T>
T squared_norm(const vector<T>& a) {
   return a * a;
}

template<typename T>
class Matrix {
public:
   size_t nbRows() const;
   size_t nbCols() const;
   Matrix(size_t nbRows, size_t nbCols, size_t value = 0);
   Matrix(vector<vector<T>> _coeffs);
   vector<vector<T>> coeffs;
};

template<typename T>
size_t Matrix<T>::nbRows() const {
   return coeffs.size();
}

template<typename T>
size_t Matrix<T>::nbCols() const {
   if (coeffs.size() == 0) {
      return 0;
   }
   return coeffs[0].size();
}

template<typename T>
Matrix<T>::Matrix(size_t nbRows, size_t nbCols, size_t value) {
   coeffs = vector<vector<T>>(nbRows, vector<T>(nbCols, value));
}

template<typename T>
Matrix<T>::Matrix(vector<vector<T>> _coeffs) {
   coeffs = _coeffs;
}

template<typename T>
Matrix<T> identity(size_t size) {
   Matrix<T> res = Matrix<T>(size, size);
   for(size_t i = 0;i < size;i++) {
      res.coeffs[i][i] = T(1);
   }
   return res;
}

template<typename T>
Matrix<T> operator + (const Matrix<T>& a, const Matrix<T>& b) {
   Matrix<T> res(a.nbCols(), a.nbRows());
   
   for(size_t iRow = 0;iRow < a.nbRows();iRow++) {
      for(size_t iCol = 0;iCol < a.nbCols();iCol++) {
         res.coeffs[iCol][iRow] = 
         	a.coeffs[iRow][iCol]
         	+ b.coeffs[iRow][iCol];
      }
   }
   
   return res;
}

template<typename T>
Matrix<T> operator - (const Matrix<T>& a, const Matrix<T>& b) {
   Matrix<T> res(a.nbCols(), a.nbRows());
   
   for(size_t iRow = 0;iRow < a.nbRows();iRow++) {
      for(size_t iCol = 0;iCol < a.nbCols();iCol++) {
         res.coeffs[iCol][iRow] = 
         	a.coeffs[iRow][iCol]
         	- b.coeffs[iRow][iCol];
      }
   }
   
   return res;
}

template<typename T>
Matrix<T> transpose(Matrix<T> mat) {
   Matrix<T> res(mat.nbCols(), mat.nbRows());
   for(size_t iRow = 0;iRow < mat.nbRows();iRow++) {
      for(size_t iCol = 0;iCol < mat.nbCols();iCol++) {
         res.coeffs[iCol][iRow] = mat.coeffs[iRow][iCol];
      }
   }
   return res;
}

template<typename T>
Matrix<T> row_echelon_form(Matrix<T> mat) {
   for(size_t iCol = 0;iCol < mat.nbCols();iCol++) {
      size_t non_zero = iCol;
      for(size_t iRow = iCol;iRow < mat.nbRows();iRow++) {
         if(!(mat.coeffs[iRow][iCol] == T(0))) {
            non_zero = iRow;
         }
      }
      
      if(iCol >= mat.nbRows()) break;
      swap(mat.coeffs[iCol], mat.coeffs[non_zero]);
      if(mat.coeffs[iCol][iCol] == T(0)) continue;
      
      mat.coeffs[iCol] = (T(1) / mat.coeffs[iCol][iCol]) * mat.coeffs[iCol];
      
      for(size_t iRow = 0;iRow < mat.nbRows();iRow++) {
         if(iRow == iCol) continue;
         if(!(mat.coeffs[iRow][iCol] == T(0))) {
            mat.coeffs[iRow] = mat.coeffs[iRow] - mat.coeffs[iRow][iCol] * mat.coeffs[iCol];
         }
      }
   }
   return mat;
}

template<typename T>
Matrix<T> kernel_basis(Matrix<T> mat) {
   Matrix<T> id = identity<T>(mat.nbRows());
   
   for(size_t iCol = 0;iCol < mat.nbCols();iCol++) {
      size_t non_zero = iCol;
      for(size_t iRow = iCol;iRow < mat.nbRows();iRow++) {
         if(!(mat.coeffs[iRow][iCol] == T(0))) {
            non_zero = iRow;
         }
      }
      
      if(iCol >= mat.nbRows()) break;
      
      swap(mat.coeffs[iCol], mat.coeffs[non_zero]);
      swap(id.coeffs[iCol], id.coeffs[non_zero]);
      
      if(mat.coeffs[iCol][iCol] == T(0)) continue;
      
      id.coeffs[iCol] = (T(1) / mat.coeffs[iCol][iCol]) * id.coeffs[iCol];
      mat.coeffs[iCol] = (T(1) / mat.coeffs[iCol][iCol]) * mat.coeffs[iCol];
      
      for(size_t iRow = 0;iRow < mat.nbRows();iRow++) {
         if(iRow == iCol) continue;
         if(!(mat.coeffs[iRow][iCol] == T(0))) {
            id.coeffs[iRow] = id.coeffs[iRow] - mat.coeffs[iRow][iCol] * id.coeffs[iCol];
            mat.coeffs[iRow] = mat.coeffs[iRow] - mat.coeffs[iRow][iCol] * mat.coeffs[iCol];
         }
      }
   }
   
   Matrix<T> basis(0, 0);
   
   for(size_t iRow = 0;iRow < mat.nbRows();iRow++) {
      bool is_zero = true;
      
      for(size_t iCol = 0;iCol < mat.nbCols();iCol++) {
         is_zero &= mat.coeffs[iRow][iCol] == T(0);
      }
      
      if(is_zero) {
         basis.coeffs.push_back(id.coeffs[iRow]);
      }
   }
   
   return basis;
}

template<typename T>
void debug(const Matrix<T>& mat) {
   for(size_t iRow = 0;iRow < mat.nbRows();iRow++) {
      for(size_t iCol = 0;iCol < mat.nbCols();iCol++) {
         cout << mat.coeffs[iRow][iCol] << " ";
      }
      cout << endl;
   }
}

template<typename T>
void debug(const vector<T>& v) {
   for(size_t iCol = 0;iCol < v.size();iCol++) {
      cout << v[iCol] << " ";
   }
   cout << endl;
}

template<typename T>
vector<T> projection(size_t i, const Matrix<T>& mat, vector<T> v) {
   vector<T> proj(v.size(), T(0));
   
   for(size_t iRow = i;iRow < mat.nbRows();iRow++) {
      proj = proj + ((v * mat.coeffs[iRow]) / squared_norm(mat.coeffs[iRow])) 
       * mat.coeffs[iRow];
   }
   
   return proj;
}

template<typename T>
Matrix<T> gram_schmidt(Matrix<T> mat) {
   for(size_t iRow = 0;iRow < mat.nbRows();iRow++) {
      for(size_t iProj = 0;iProj < iRow;iProj++) {
         mat.coeffs[iRow] = mat.coeffs[iRow]
          - ((mat.coeffs[iRow] * mat.coeffs[iProj]) / squared_norm(mat.coeffs[iProj])) 
          * mat.coeffs[iProj];
      }
   }
   
   return mat;
}

template<typename T>
vector<T> nearest_plane(Matrix<T> mat, vector<T> v) {
   if(mat.nbRows() == 0) return vector<T>(v.size(), T(0));
   
   Matrix<T> ortho = gram_schmidt(mat);
   T coeff = round((v * ortho.coeffs.back()) / squared_norm(ortho.coeffs.back()));
   
   vector<T> bn = mat.coeffs.back();
   mat.coeffs.pop_back();
   
   return coeff * bn + nearest_plane(mat, v - coeff * bn);
}

template<typename T>
Matrix<T> LLL(Matrix<T> mat, T delta) {
   Matrix<T> ortho = gram_schmidt(mat);
   
   for(size_t iRow = 1;iRow < mat.nbRows();iRow++) {
      mat.coeffs[iRow] = mat.coeffs[iRow] + nearest_plane(mat, ortho.coeffs[iRow] - mat.coeffs[iRow]);
   }
   
   ortho = gram_schmidt(mat);
   for(size_t iRow = 1;iRow < mat.nbRows();iRow++) {
      vector<T> projA = ortho.coeffs[iRow - 1];
      vector<T> projB = projection(iRow - 1, ortho, mat.coeffs[iRow]);
      
      T diff = delta * squared_norm(projA) - squared_norm(projB);
      if(diff.getNumerator() > 0) {
         swap(mat.coeffs[iRow], mat.coeffs[iRow - 1]);
         return LLL(mat, delta);
      }
   }
   
   return mat;
}
