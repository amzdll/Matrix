#ifndef CPP1_S21_MATRIXPLUS_S21_MATRIX_OOP_H_
#define CPP1_S21_MATRIXPLUS_S21_MATRIX_OOP_H_

#include <cmath>
#include <cstdlib>
#include <iostream>

class S21Matrix {
public:
  S21Matrix();
  explicit S21Matrix(int rows, int cols);
  S21Matrix(const S21Matrix &other);
  S21Matrix(S21Matrix &&other) noexcept;
  ~S21Matrix();

  int get_rows() const;
  int get_cols() const;
  void set_rows(int rows);
  void set_cols(int cols) { this->cols_ = cols; }

  bool EqMatrix(const S21Matrix &other);
  void SumMatrix(const S21Matrix &other);
  void SubMatrix(const S21Matrix &other);
  void MulNumber(double num);
  void MulMatrix(const S21Matrix &other);
  S21Matrix Transpose();
  S21Matrix CalcComplements();
  double Determinant();
  S21Matrix InverseMatrix();

  S21Matrix operator+(const S21Matrix &other);
  S21Matrix operator-(const S21Matrix &other);
  S21Matrix operator*(double num);
  S21Matrix operator*(const S21Matrix &other);
  bool operator==(const S21Matrix &other);
  S21Matrix &operator=(const S21Matrix &other);
  S21Matrix operator+=(const S21Matrix &other);
  S21Matrix operator-=(const S21Matrix &other);
  S21Matrix operator*=(double num);
  S21Matrix operator*=(const S21Matrix &other);
  double &operator()(int i, int j);

  // additional to test
  void FillingMatrix();
  void ZeroingMatrix();

  void display_matrix();

private:
  int rows_, cols_;
  double **matrix_;

  bool MatrixIsExtend() const;
  bool SizeComparison(const S21Matrix &other) const;
  void ShortenedCopy(const S21Matrix &other, int row, int column);
  void DeleteMatrix();
};

#endif // CPP1_S21_MATRIXPLUS_S21_MATRIX_OOP_H_
