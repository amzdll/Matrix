#include "s21_matrix_oop.h"

#include <cstdio>


// TODO free throw ???

//  CONSTRUCTORS
S21Matrix::S21Matrix() {
  rows_ = 0;
  cols_ = 0;#include "s21_matrix_oop.h"
    
#include <cstdio>

S21Matrix::S21Matrix() {
  rows_ = 0;
  cols_ = 0;
  matrix_ = nullptr;
}

S21Matrix::S21Matrix(int rows, int cols) {
  if (rows > 0 && cols > 0) {
    rows_ = rows;
    cols_ = cols;
    matrix_ = new double *[rows_];
    for (int i = 0; i < rows_; ++i) {
      matrix_[i] = new double[cols_];
    }
    this->ZeroingMatrix();
  }
}

S21Matrix::S21Matrix(const S21Matrix &other) {
  if (other.MatrixIsExtend()) {
    rows_ = other.rows_;
    cols_ = other.cols_;
    matrix_ = new double *[other.rows_];
    for (int i = 0; i < other.rows_; ++i) {
      matrix_[i] = new double[cols_];
    }
    for (int i = 0; i < rows_; ++i) {
      for (int j = 0; j < cols_; ++j) {
        matrix_[i][j] = other.matrix_[i][j];
      }
    }
  }
}

S21Matrix::S21Matrix(S21Matrix &&other) noexcept {
  this->rows_ = other.rows_;
  this->cols_ = other.cols_;
  this->matrix_ = other.matrix_;
  other.rows_ = 0;
  other.cols_ = 0;
  other.matrix_ = nullptr;
}

S21Matrix::~S21Matrix() { DeleteMatrix(); }

int S21Matrix::get_rows() const { return rows_; }

int S21Matrix::get_cols() const { return cols_; }

void S21Matrix::set_rows(int rows) {
  if (this->matrix_) {
    S21Matrix temp_matrix(rows, cols_);
    for (int i = 0; i < rows_; ++i) {
      for (int j = 0; j < cols_; ++j) {
        temp_matrix.matrix_[i][j] = matrix_[i][j];
        temp_matrix.display_matrix();
        printf("-------------------------------------------------------\n");
      }
    }
    matrix_ = temp_matrix.matrix_;
    rows_ = rows;
  }
}

bool S21Matrix::EqMatrix(const S21Matrix &other) {
  bool result = false;
  if (this->SizeComparison(other) && this->MatrixIsExtend() &&
      other.MatrixIsExtend()) {
    result = true;
    for (int i = 0; i < this->rows_ && result; ++i) {
      for (int j = 0; j < this->cols_ && result; ++j) {
        if (fabs(this->matrix_[i][j] - other.matrix_[i][j]) > 1e-07) {
          result = false;
        }
      }
    }
  }
  return result;
}

void S21Matrix::SumMatrix(const S21Matrix &other) {
  if (this->SizeComparison(other)) {
    if (this->MatrixIsExtend() && other.MatrixIsExtend()) {
      for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
          matrix_[i][j] = matrix_[i][j] + other.matrix_[i][j];
        }
      }
    }
  } else {
    throw std::out_of_range("different dimensionality of matrices");
  }
}

void S21Matrix::SubMatrix(const S21Matrix &other) {
  if (this->SizeComparison(other)) {
    if (this->MatrixIsExtend() && other.MatrixIsExtend()) {
      for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
          matrix_[i][j] = matrix_[i][j] - other.matrix_[i][j];
        }
      }
    }
  } else {
    throw std::out_of_range("different dimensionality of matrices");
  }
}

void S21Matrix::MulNumber(const double num) {
  if (this->MatrixIsExtend()) {
    for (int i = 0; i < rows_; ++i) {
      for (int j = 0; j < cols_; ++j) {
        this->matrix_[i][j] *= num;
      }
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix &other) {
  if (cols_ == other.rows_) {
    if (this->MatrixIsExtend() && other.MatrixIsExtend()) {
      S21Matrix temp_result(rows_, other.cols_);
      for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < other.cols_; ++j) {
          for (int k = 0; k < cols_; ++k) {
            temp_result.matrix_[i][j] += matrix_[i][k] * other.matrix_[k][j];
          }
        }
      }
      *this = temp_result;
    }
  } else {
    throw std::out_of_range(
        "the number of columns of the first matrix does not equal the number "
        "of rows of the second matrix");
  }
}

S21Matrix S21Matrix::Transpose() {
  S21Matrix Transpose(this->cols_, this->rows_);
  if (this->MatrixIsExtend()) {
    for (int i = 0; i < this->rows_; ++i) {
      for (int j = 0; j < this->cols_; ++j) {
        Transpose.matrix_[j][i] = this->matrix_[i][j];
      }
    }
  }
  return Transpose;
}

S21Matrix S21Matrix::CalcComplements() {
  S21Matrix result(this->rows_, this->cols_);
  if (this->rows_ == this->cols_) {
    if (this->MatrixIsExtend()) {
      S21Matrix temp_matrix(this->rows_ - 1, this->cols_ - 1);
      double temp_result;
      for (int i = 0; i < this->rows_; ++i) {
        for (int j = 0; j < this->cols_; ++j) {
          temp_matrix.ShortenedCopy(*this, i, j);
          temp_result = temp_matrix.Determinant();
          result.matrix_[i][j] += pow(-1, i + j) * temp_result;
        }
      }
    }
  } else {
    throw std::out_of_range("the matrix is not square");
  }
  return result;
}

double S21Matrix::Determinant() {
  double result = 0;
  double temp_result = 0.0;
  if (this->rows_ == this->cols_) {
    if (this->MatrixIsExtend()) {
      if (this->rows_ == 1) {
        temp_result = this->matrix_[0][0];
      } else if (this->rows_ == 2) {
        temp_result = (this->matrix_[0][0] * this->matrix_[1][1]) -
                      (this->matrix_[1][0] * this->matrix_[0][1]);
      } else {
        S21Matrix temp_matrix(this->rows_ - 1, this->cols_ - 1);
        for (int i = 0; i < this->cols_; ++i) {
          temp_matrix.ShortenedCopy(*this, 0, i);
          result = temp_matrix.Determinant();
          temp_result += this->matrix_[0][i] * pow(-1, i) * result;
        }
      }
      result = temp_result;
    }
  } else {
    throw std::out_of_range("the matrix is not square");
  }
  return result;
}

S21Matrix S21Matrix::InverseMatrix() {
  S21Matrix temp_matrix(this->rows_, this->cols_);
  double deter = this->Determinant();
  if (deter != 0) {
    if (this->MatrixIsExtend() && (this->rows_ == this->cols_)) {
      temp_matrix = this->CalcComplements();
      temp_matrix = temp_matrix.Transpose();
      for (int i = 0; i < this->rows_; ++i) {
        for (int j = 0; j < this->cols_; ++j) {
          temp_matrix.matrix_[i][j] = 1 / deter * temp_matrix.matrix_[i][j];
        }
      }
    }
  } else {
    throw std::invalid_argument("the Determinant of the matrix is 0");
  }
  return temp_matrix;
}

// OVERLOADING
S21Matrix S21Matrix::operator+(const S21Matrix &other) {
  S21Matrix new_matrix(*this);
  new_matrix.SumMatrix(other);
  return new_matrix;
}

S21Matrix S21Matrix::operator-(const S21Matrix &other) {
  S21Matrix new_matrix(*this);
  new_matrix.SubMatrix(other);
  return new_matrix;
}

S21Matrix S21Matrix::operator*(const double num) {
  S21Matrix new_matrix(*this);
  new_matrix.MulNumber(num);
  return new_matrix;
}

S21Matrix S21Matrix::operator*(const S21Matrix &other) {
  S21Matrix new_matrix(*this);
  new_matrix.MulMatrix(other);
  return new_matrix;
}

bool S21Matrix::operator==(const S21Matrix &other) { return EqMatrix(other); }

S21Matrix &S21Matrix::operator=(const S21Matrix &other) {
  if (this != &other) {
    this->DeleteMatrix();
    rows_ = other.rows_;
    cols_ = other.cols_;
    matrix_ = new double *[rows_];
    for (int i = 0; i < rows_; ++i) {
      matrix_[i] = new double[cols_];
    }
    for (int i = 0; i < rows_; ++i) {
      for (int j = 0; j < cols_; ++j) {
        matrix_[i][j] = other.matrix_[i][j];
      }
    }
  }
  return *this;
}

S21Matrix S21Matrix::operator+=(const S21Matrix &other) {
  this->SumMatrix(other);
  return *this;
}

S21Matrix S21Matrix::operator-=(const S21Matrix &other) {
  this->SubMatrix(other);
  return *this;
}

S21Matrix S21Matrix::operator*=(const double num) {
  this->MulNumber(num);
  return *this;
}

S21Matrix S21Matrix::operator*=(const S21Matrix &other) {
  this->MulMatrix(other);
  return *this;
}

double &S21Matrix::operator()(int i, int j) {
  if (!(this->rows_ >= i && i >= 0 && this->cols_ >= j && j >= 0)) {
    throw std::out_of_range("different dimensionality of matrices");
  }
  return this->matrix_[i][j];
}

void S21Matrix::FillingMatrix() {
  double temp = 0;
  for (int i = 0; i < this->get_rows(); ++i) {
    for (int j = 0; j < this->get_cols(); ++j) {
      this->matrix_[i][j] = temp;
      temp++;
    }
  }
}

//  ADDITIONAL
bool S21Matrix::SizeComparison(const S21Matrix &other) const {
  bool result = false;
  if (cols_ == other.cols_ && rows_ == other.rows_) {
    result = true;
  }
  return result;
}

bool S21Matrix::MatrixIsExtend() const {
  bool result = false;
  if (cols_ > 0 && rows_ > 0 && matrix_ != nullptr) {
    result = true;
  }
  return result;
}

void S21Matrix::DeleteMatrix() {
  if (this->matrix_) {
    for (int i = 0; i < rows_; ++i) {
      delete this->matrix_[i];
    }
    delete this->matrix_;
  }
}

void S21Matrix::ZeroingMatrix() {
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      this->matrix_[i][j] = 0;
    }
  }
}

void S21Matrix::ShortenedCopy(const S21Matrix &other, int row, int column) {
  for (int i = 0, x = 0; i < other.rows_; ++i) {
    if (i != row) {
      for (int j = 0, y = 0; j < other.cols_; ++j) {
        if (j != column) {
          this->matrix_[x][y] = other.matrix_[i][j];
          y++;
        }
      }
      x++;
    }
  }
}

void S21Matrix::display_matrix() {
  for (int i = 0; i < this->get_rows(); ++i) {
    for (int j = 0; j < this->get_cols(); ++j) {
      printf("%lf\t", this->matrix_[i][j]);
    }
    printf("\n");
  }
}

int main() {
  S21Matrix a(3, 3);
  //  S21Matrix b(3, 3);
  a.FillingMatrix();
  //  b.FillingMatrix();
  a.display_matrix();
  a.set_rows(4);
//  printf("\n");
  a.display_matrix();
  return 0;
}
  matrix_ = nullptr;
}

S21Matrix::S21Matrix(int rows, int cols) {
  if (rows > 0 && cols > 0) {
    rows_ = rows;
    cols_ = cols;
    matrix_ = new double *[rows_];
    for (int i = 0; i < rows_; ++i) {
      matrix_[i] = new double[cols_];
    }
    this->zeroing_matrix();
  }
}

S21Matrix::S21Matrix(const S21Matrix &other) {
  if (other.matrix_is_existing()) {
    rows_ = other.rows_;
    cols_ = other.cols_;
    matrix_ = new double *[other.rows_];
    for (int i = 0; i < other.rows_; ++i) {
      matrix_[i] = new double[cols_];
    }
    for (int i = 0; i < rows_; ++i) {
      for (int j = 0; j < cols_; ++j) {
        matrix_[i][j] = other.matrix_[i][j];
      }
    }
  }
}

S21Matrix::S21Matrix(S21Matrix &&other) noexcept {
  this->rows_ = other.rows_;
  this->cols_ = other.cols_;
  this->matrix_ = other.matrix_;
  other.rows_ = 0;
  other.cols_ = 0;
  other.matrix_ = nullptr;
}

S21Matrix::~S21Matrix() { delete_matrix(); }

// CORE
bool S21Matrix::eq_matrix(const S21Matrix &other) {
  bool result = false;
  if (this->size_comparison(other) && this->matrix_is_existing() &&
      other.matrix_is_existing()) {
    result = true;
    for (int i = 0; i < this->rows_ && result; ++i) {
      for (int j = 0; j < this->cols_ && result; ++j) {
        if (fabs(this->matrix_[i][j] - other.matrix_[i][j]) > 1e-07) {
          result = false;
        }
      }
    }
  }
  return result;
}

void S21Matrix::sum_matrix(const S21Matrix &other) {
  if (this->size_comparison(other)) {
    if (this->matrix_is_existing() && other.matrix_is_existing()) {
      for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
          matrix_[i][j] = matrix_[i][j] + other.matrix_[i][j];
        }
      }
    }
  } else {
    throw std::out_of_range("different dimensionality of matrices");
  }
}

void S21Matrix::sub_matrix(const S21Matrix &other) {
  if (this->size_comparison(other)) {
    if (this->matrix_is_existing() && other.matrix_is_existing()) {
      for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
          matrix_[i][j] = matrix_[i][j] - other.matrix_[i][j];
        }
      }
    }
  } else {
    throw std::out_of_range("different dimensionality of matrices");
  }
}

void S21Matrix::mul_number(const double num) {
  if (this->matrix_is_existing()) {
    for (int i = 0; i < rows_; ++i) {
      for (int j = 0; j < cols_; ++j) {
        this->matrix_[i][j] *= num;
      }
    }
  }
}

void S21Matrix::mul_matrix(const S21Matrix &other) {
  if (cols_ == other.rows_) {
    if (this->matrix_is_existing() && other.matrix_is_existing()) {
      S21Matrix temp_result(rows_, other.cols_);
      for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < other.cols_; ++j) {
          for (int k = 0; k < cols_; ++k) {
            temp_result.matrix_[i][j] += matrix_[i][k] * other.matrix_[k][j];
          }
        }
      }
      *this = temp_result;
    }
  } else {
    throw std::out_of_range(
        "the number of columns of the first matrix does not equal the number "
        "of rows of the second matrix");
  }
}

S21Matrix S21Matrix::transpose() {
  S21Matrix transposed_matrix(this->cols_, this->rows_);
  if (this->matrix_is_existing()) {
    for (int i = 0; i < this->rows_; ++i) {
      for (int j = 0; j < this->cols_; ++j) {
        transposed_matrix.matrix_[j][i] = this->matrix_[i][j];
      }
    }
  }
  return transposed_matrix;
}

S21Matrix S21Matrix::calc_complements() {
  S21Matrix result(this->rows_, this->cols_);
  if (this->rows_ == this->cols_) {
    if (this->matrix_is_existing()) {
      S21Matrix temp_matrix(this->rows_ - 1, this->cols_ - 1);
      double temp_result;
      for (int i = 0; i < this->rows_; ++i) {
        for (int j = 0; j < this->cols_; ++j) {
          temp_matrix.shortened_copy(*this, i, j);
          temp_result = temp_matrix.determinant();
          result.matrix_[i][j] += pow(-1, i + j) * temp_result;
        }
      }
    }
  } else {
    throw std::out_of_range("the matrix is not square");
  }
  return result;
}

double S21Matrix::determinant() {
  double result = 0;
  double temp_result = 0.0;
  if (this->rows_ == this->cols_) {
    if (this->matrix_is_existing()) {
      if (this->rows_ == 1) {
        temp_result = this->matrix_[0][0];
      } else if (this->rows_ == 2) {
        temp_result = (this->matrix_[0][0] * this->matrix_[1][1]) -
                      (this->matrix_[1][0] * this->matrix_[0][1]);
      } else {
        S21Matrix temp_matrix(this->rows_ - 1, this->cols_ - 1);
        for (int i = 0; i < this->cols_; ++i) {
          temp_matrix.shortened_copy(*this, 0, i);
          result = temp_matrix.determinant();
          temp_result += this->matrix_[0][i] * pow(-1, i) * result;
        }
      }
      result = temp_result;
    }
  } else {
    throw std::out_of_range("the matrix is not square");
  }
  return result;
}

S21Matrix S21Matrix::inverse_matrix() {
  S21Matrix temp_matrix(this->rows_, this->cols_);
  double deter = this->determinant();
  if (deter != 0) {
    if (this->matrix_is_existing() && (this->rows_ == this->cols_)) {
      temp_matrix = this->calc_complements();
      temp_matrix = temp_matrix.transpose();
      for (int i = 0; i < this->rows_; ++i) {
        for (int j = 0; j < this->cols_; ++j) {
          temp_matrix.matrix_[i][j] = 1 / deter * temp_matrix.matrix_[i][j];
        }
      }
    }
  } else {
    throw std::invalid_argument("the determinant of the matrix is 0");
  }
  return temp_matrix;
}

// OVERLOADING
S21Matrix S21Matrix::operator+(const S21Matrix &other) {
  S21Matrix new_matrix(*this);
  new_matrix.sum_matrix(other);
  return new_matrix;
}

S21Matrix S21Matrix::operator-(const S21Matrix &other) {
  S21Matrix new_matrix(*this);
  new_matrix.sub_matrix(other);
  return new_matrix;
}

S21Matrix S21Matrix::operator*(const double num) {
  S21Matrix new_matrix(*this);
  new_matrix.mul_number(num);
  return new_matrix;
}

S21Matrix S21Matrix::operator*(const S21Matrix &other) {
  S21Matrix new_matrix(*this);
  new_matrix.mul_matrix(other);
  return new_matrix;
}

bool S21Matrix::operator==(const S21Matrix &other) { return eq_matrix(other); }

S21Matrix &S21Matrix::operator=(const S21Matrix &other) {
  if (this != &other) {
    this->delete_matrix();
    rows_ = other.rows_;
    cols_ = other.cols_;
    matrix_ = new double *[rows_];
    for (int i = 0; i < rows_; ++i) {
      matrix_[i] = new double[cols_];
    }
    for (int i = 0; i < rows_; ++i) {
      for (int j = 0; j < cols_; ++j) {
        matrix_[i][j] = other.matrix_[i][j];
      }
    }
  }
  return *this;
}

S21Matrix S21Matrix::operator+=(const S21Matrix &other) {
  this->sum_matrix(other);
  return *this;
}

S21Matrix S21Matrix::operator-=(const S21Matrix &other) {
  this->sub_matrix(other);
  return *this;
}

S21Matrix S21Matrix::operator*=(const double num) {
  this->mul_number(num);
  return *this;
}

S21Matrix S21Matrix::operator*=(const S21Matrix &other) {
  this->mul_matrix(other);
  return *this;
}

double &S21Matrix::operator()(int i, int j) {
  if (!(this->rows_ >= i && i >= 0 && this->cols_ >= j && j >= 0)) {
    throw std::out_of_range("different dimensionality of matrices");
  }
  return this->matrix_[i][j];
}

// ADDITIONAL TEMP
void S21Matrix::display_matrix() {
  for (int i = 0; i < this->get_rows(); ++i) {
    for (int j = 0; j < this->get_cols(); ++j) {
      printf("%lf\t", this->matrix_[i][j]);
    }
    printf("\n");
  }
}

void S21Matrix::filling_matrix() {
  int temp = 0;
  for (int i = 0; i < this->get_rows(); ++i) {
    for (int j = 0; j < this->get_cols(); ++j) {
      this->matrix_[i][j] = temp;
      temp++;
    }
  }
}

//  ADDITIONAL
bool S21Matrix::size_comparison(const S21Matrix &other) const {
  bool result = false;
  if (cols_ == other.cols_ && rows_ == other.rows_) {
    result = true;
  }
  return result;
}

bool S21Matrix::matrix_is_existing() const {
  bool result = false;
  if (cols_ > 0 && rows_ > 0 && matrix_ != nullptr) {
    result = true;
  }
  return result;
}

void S21Matrix::delete_matrix() {
  if (this->matrix_) {
    for (int i = 0; i < rows_; ++i) {
      delete this->matrix_[i];
    }
    delete this->matrix_;
  }
}

void S21Matrix::zeroing_matrix() {
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      this->matrix_[i][j] = 0;
    }
  }
}

void S21Matrix::shortened_copy(const S21Matrix &other, int row, int column) {
  for (int i = 0, x = 0; i < other.rows_; ++i) {
    if (i != row) {
      for (int j = 0, y = 0; j < other.cols_; ++j) {
        if (j != column) {
          this->matrix_[x][y] = other.matrix_[i][j];
          y++;
        }
      }
      x++;
    }
  }
}
