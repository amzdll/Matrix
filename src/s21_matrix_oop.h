#ifndef CPP1_S21_MATRIXPLUS_S21_MATRIX_OOP_H_
#define CPP1_S21_MATRIXPLUS_S21_MATRIX_OOP_H_

#include <cmath>
#include <cstdlib>
#include <iostream>

class S21Matrix {
 private:
  int rows_, cols_;
  double** matrix_;

 public:
  S21Matrix();
  S21Matrix(int rows, int cols);
  S21Matrix(const S21Matrix& other);
  S21Matrix(S21Matrix&& other) noexcept;
  ~S21Matrix();

  // temp additional functions
  void display_matrix();

  // accessor && mutator
  int get_rows() const { return rows_; }
  int get_cols() const { return cols_; }
  void set_rows(int rows) { this->rows_ = rows; }
  void set_cols(int cols) { this->cols_ = cols; }

  //  core
  bool eq_matrix(const S21Matrix& other);
  void sum_matrix(const S21Matrix& other);
  void sub_matrix(const S21Matrix& other);
  void mul_number(double num);
  void mul_matrix(const S21Matrix& other);
  S21Matrix transpose();
  S21Matrix calc_complements();
  double determinant();
  S21Matrix inverse_matrix();

  //  overloading
  S21Matrix operator+(const S21Matrix& other);
  S21Matrix operator-(const S21Matrix& other);
  S21Matrix operator*(double num);
  S21Matrix operator*(const S21Matrix& other);
  bool operator==(const S21Matrix& other);
  S21Matrix& operator=(const S21Matrix& other);
  S21Matrix operator+=(const S21Matrix& other);
  S21Matrix operator-=(const S21Matrix& other);
  S21Matrix operator*=(double num);
  S21Matrix operator*=(const S21Matrix& other);
  double& operator()(int i, int j);

  // additional
  void filling_matrix();
  bool size_comparison(const S21Matrix& other) const;
  bool matrix_is_existing() const;
  void delete_matrix();
  void zeroing_matrix();
  void shortened_copy(const S21Matrix& other, int row, int column);
};

#endif  // CPP1_S21_MATRIXPLUS_S21_MATRIX_OOP_H_
