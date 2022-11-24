#include "s21_matrix_oop.h"
#include "gtest/gtest.h"

TEST(S21Matrix_constructor_suite, true_test) {
  S21Matrix a;
  EXPECT_EQ(a.get_rows(), 0);
  EXPECT_EQ(a.get_cols(), 0);
}

TEST(S21Matrix_move_constructor_suite, true_test) {
  S21Matrix a(3, 3);
  S21Matrix b(std::move(a));
  EXPECT_EQ(a.get_rows(), 0);
  EXPECT_EQ(a.get_cols(), 0);
  EXPECT_EQ(b.get_rows(), 3);
  EXPECT_EQ(b.get_cols(), 3);
}

TEST(set_rows_suite, extend_test) {
  S21Matrix a(3, 3);
  S21Matrix b(4, 3);
  a.FillingMatrix();
  b.FillingMatrix();
  a.set_rows(4);
  b(3, 0) = 0;
  b(3, 1) = 0;
  b(3, 2) = 0;
  EXPECT_TRUE(a == b);
}

TEST(set_rows_suite, reduce_test) {
  S21Matrix a(3, 3);
  S21Matrix b(4, 3);
  a.FillingMatrix();
  b.FillingMatrix();
  b.set_rows(3);
  EXPECT_TRUE(a == b);
}

TEST(set_cols_suite, extend_test) {
  S21Matrix a(3, 3);
  a.FillingMatrix();
  a.set_cols(4);

  EXPECT_EQ(a(0, 0), 0);
  EXPECT_EQ(a(0, 1), 1.0);
  EXPECT_EQ(a(0, 2), 2.0);
  EXPECT_EQ(a(1, 0), 3.0);
  EXPECT_EQ(a(1, 1), 4.0);
  EXPECT_EQ(a(1, 2), 5.0);
  EXPECT_EQ(a(2, 0), 6.0);
  EXPECT_EQ(a(2, 1), 7.0);
  EXPECT_EQ(a(2, 2), 8.0);
  EXPECT_TRUE(0 - a(0, 3) < 1e-7);
  EXPECT_TRUE(0 - a(1, 3) < 1e-7);
  EXPECT_TRUE(0 - a(2, 3) < 1e-7);
}

TEST(set_cols_suite, reduce_test) {
  S21Matrix a(3, 3);
  a.FillingMatrix();
  a.set_cols(2);
  EXPECT_EQ(a(0, 0), 0);
  EXPECT_EQ(a(0, 1), 1.0);
  EXPECT_EQ(a(1, 0), 3.0);
  EXPECT_EQ(a(1, 1), 4.0);
  EXPECT_EQ(a(2, 0), 6.0);
  EXPECT_EQ(a(2, 1), 7.0);
}

TEST(EqMatrix_suite, true_test) {
  S21Matrix a(3, 3);
  S21Matrix b(3, 3);
  a.FillingMatrix();
  b.FillingMatrix();

  EXPECT_TRUE(a.EqMatrix(b));
  EXPECT_TRUE(b.EqMatrix(a));

  EXPECT_TRUE(a == b);
  EXPECT_TRUE(b == a);
}

TEST(EqMatrix_suite, false_test) {
  S21Matrix a(2, 3);
  S21Matrix b(3, 3);
  a.FillingMatrix();
  b.FillingMatrix();

  EXPECT_FALSE(a.EqMatrix(b));
  EXPECT_FALSE(b.EqMatrix(a));

  EXPECT_FALSE(a == b);
  EXPECT_FALSE(b == a);
}

TEST(EqMatrix_suite, accuracy_true_test) {
  S21Matrix a(3, 3);
  S21Matrix b(3, 3);
  a.FillingMatrix();
  b.FillingMatrix();
  for (int i = 0; i < a.get_rows(); i++) {
    for (int j = 0; j < a.get_cols(); j++) {
      a(i, j) += 0.00000001;
    }
  }

  EXPECT_TRUE(a.EqMatrix(b));
  EXPECT_TRUE(a == b);
}

TEST(EqMatrix_suite, accuracy_false_test) {
  S21Matrix a(3, 3);
  S21Matrix b(3, 3);
  a.FillingMatrix();
  b.FillingMatrix();
  for (int i = 0; i < a.get_rows(); i++) {
    for (int j = 0; j < a.get_cols(); j++) {
      a(i, j) += 0.000001;
    }
  }
  EXPECT_FALSE(a.EqMatrix(b));
  EXPECT_FALSE(a == b);
}

TEST(SumMatrix_suite, true_test) {
  S21Matrix a(3, 3);
  S21Matrix b(3, 3);
  S21Matrix c(3, 3);
  a.FillingMatrix();
  b.FillingMatrix();
  c.FillingMatrix();
  a.SumMatrix(c);
  b.SumMatrix(c);
  EXPECT_TRUE(a.EqMatrix(b));
}

TEST(SumMatrix_suite, exceptional_test) {
  S21Matrix a(3, 3);
  S21Matrix b(3, 4);
  a.FillingMatrix();
  b.FillingMatrix();
  ASSERT_THROW(a.SumMatrix(b), std::out_of_range);
}

TEST(SubMatrix_suite, true_test) {
  S21Matrix a(3, 3);
  S21Matrix b(3, 3);
  S21Matrix c(3, 3);
  a.FillingMatrix();
  b.FillingMatrix();
  c.FillingMatrix();

  a.SubMatrix(c);
  b.SubMatrix(c);
  EXPECT_TRUE(a.EqMatrix(b));
}

TEST(SubMatrix_suite, exceptional_test) {
  S21Matrix a(3, 3);
  S21Matrix b(3, 4);
  a.FillingMatrix();
  b.FillingMatrix();
  ASSERT_THROW(a.SubMatrix(b), std::out_of_range);
}

TEST(MulNumber_suite, true_test) {
  S21Matrix a(3, 3);
  S21Matrix b(3, 3);
  a.FillingMatrix();
  b.FillingMatrix();
  for (int i = 0; i < 4; ++i) {
    b.SumMatrix(a);
  }
  a.MulNumber(5);
  EXPECT_TRUE(a.EqMatrix(b));
}

TEST(MulMatrix_suite, three_three_true_test) {
  S21Matrix a(3, 3);
  S21Matrix b(3, 3);
  S21Matrix expected_result(3, 3);
  a.FillingMatrix();
  b.FillingMatrix();
  expected_result(0, 0) = 15;
  expected_result(0, 1) = 18;
  expected_result(0, 2) = 21;
  expected_result(1, 0) = 42;
  expected_result(1, 1) = 54;
  expected_result(1, 2) = 66;
  expected_result(2, 0) = 69;
  expected_result(2, 1) = 90;
  expected_result(2, 2) = 111;
  a.MulMatrix(b);
  EXPECT_TRUE(a.EqMatrix(expected_result));
}

TEST(MulMatrix_suite, four_four_true_test) {
  S21Matrix a(4, 4);
  S21Matrix b(4, 4);
  S21Matrix expected_result(4, 4);
  a.FillingMatrix();
  b.FillingMatrix();
  expected_result(0, 0) = 56;
  expected_result(0, 1) = 62;
  expected_result(0, 2) = 68;
  expected_result(0, 3) = 74;
  expected_result(1, 0) = 152;
  expected_result(1, 1) = 174;
  expected_result(1, 2) = 196;
  expected_result(1, 3) = 218;
  expected_result(2, 0) = 248;
  expected_result(2, 1) = 286;
  expected_result(2, 2) = 324;
  expected_result(2, 3) = 362;
  expected_result(3, 0) = 344;
  expected_result(3, 1) = 398;
  expected_result(3, 2) = 452;
  expected_result(3, 3) = 506;
  a.MulMatrix(b);
  EXPECT_TRUE(a.EqMatrix(expected_result));
}

TEST(MulMatrix_suite, five_two_true_test) {
  S21Matrix a(2, 5);
  S21Matrix b(5, 2);
  a.FillingMatrix();
  b.FillingMatrix();
  S21Matrix expected_result(2, 2);
  expected_result(0, 0) = 60;
  expected_result(0, 1) = 70;
  expected_result(1, 0) = 160;
  expected_result(1, 1) = 195;
  a.MulMatrix(b);
  EXPECT_TRUE(a.EqMatrix(expected_result));
}

TEST(MulMatrix_suite, exceptional_test) {
  S21Matrix a(3, 3);
  S21Matrix b(2, 3);
  a.FillingMatrix();
  b.FillingMatrix();
  ASSERT_THROW(a.MulMatrix(b), std::out_of_range);
}

TEST(Transpose_matrix_suite, true_test) {
  S21Matrix a(2, 3);
  S21Matrix expected_result(3, 2);
  a.FillingMatrix();
  a = a.Transpose();
  expected_result(0, 0) = 0;
  expected_result(0, 1) = 3;
  expected_result(1, 0) = 1;
  expected_result(1, 1) = 4;
  expected_result(2, 0) = 2;
  expected_result(2, 1) = 5;
  EXPECT_TRUE(a.EqMatrix(expected_result));
}

TEST(CalcComplement_suite, four_four_test) {
  S21Matrix b(4, 4);
  b.FillingMatrix();
  b(0, 0) = 50;
  S21Matrix expected_result(4, 4);
  b = b.CalcComplements();
  expected_result(0, 0) = 0.000000;
  expected_result(0, 1) = 0.000000;
  expected_result(0, 2) = 0.000000;
  expected_result(0, 3) = 0.000000;
  expected_result(1, 0) = 0.000000;
  expected_result(1, 1) = -200.000000;
  expected_result(1, 2) = 400.000000;
  expected_result(1, 3) = -200.000000;
  expected_result(2, 0) = 0.000000;
  expected_result(2, 1) = 400.000000;
  expected_result(2, 2) = -800.000000;
  expected_result(2, 3) = 400.000000;
  expected_result(3, 0) = 0.000000;
  expected_result(3, 1) = -200.000000;
  expected_result(3, 2) = 400.000000;
  expected_result(3, 3) = -200.000000;
  EXPECT_TRUE(b.EqMatrix(expected_result));
}

TEST(CalcComplement_suite, three_three_test) {
  S21Matrix a(3, 3);
  a.FillingMatrix();
  S21Matrix expected_result(3, 3);
  a = a.CalcComplements();
  expected_result(0, 0) = -3.000000;
  expected_result(0, 1) = 6.000000;
  expected_result(0, 2) = -3.000000;
  expected_result(1, 0) = 6.000000;
  expected_result(1, 1) = -12.000000;
  expected_result(1, 2) = 6.000000;
  expected_result(2, 0) = -3.000000;
  expected_result(2, 1) = 6.000000;
  expected_result(2, 2) = -3.000000;
  EXPECT_TRUE(a.EqMatrix(expected_result));
}

TEST(CalcComplement_suite, exceptional_test) {
  S21Matrix a(3, 4);
  a.FillingMatrix();
  ASSERT_THROW(a.CalcComplements(), std::out_of_range);
}

TEST(Determinant_suite, three_three_test) {
  S21Matrix a(3, 3);
  S21Matrix b(3, 3);
  a.FillingMatrix();
  b.FillingMatrix();
  b(2, 2) = 56;
  S21Matrix c(b);
  c(0, 0) = 1000;
  EXPECT_TRUE(a.Determinant() == 0);
  EXPECT_TRUE(b.Determinant() == -144);
  EXPECT_TRUE(c.Determinant() == 188856);
}

TEST(Determinant_suite, one_one_test) {
  S21Matrix a(1, 1);
  a(0, 0) = 1;
  EXPECT_TRUE(a.Determinant() == 1);
}

TEST(Determinant_suite, four_four_test) {
  S21Matrix a(4, 4);
  S21Matrix b(4, 4);
  a.FillingMatrix();
  b.FillingMatrix();
  b(0, 0) = 100;
  b(3, 3) = 56;
  S21Matrix c(b);
  c(0, 0) = 1000;
  EXPECT_TRUE(a.Determinant() == 0);
  EXPECT_TRUE(b.Determinant() == -16400);
  EXPECT_TRUE(c.Determinant() == -164000);
}

TEST(Determinant_suite, exceptional_test) {
  S21Matrix a(3, 4);
  a.FillingMatrix();
  ASSERT_THROW(a.Determinant(), std::out_of_range);
}

TEST(InverseMatrix_suite, three_three_test) {
  S21Matrix a(3, 3);
  S21Matrix expected_result_a(3, 3);
  a.FillingMatrix();
  a(0, 0) = 50;
  a = a.InverseMatrix();
  expected_result_a(0, 0) = 0.02;
  expected_result_a(0, 1) = -0.04;
  expected_result_a(0, 2) = 0.02;
  expected_result_a(1, 0) = -0.04;
  expected_result_a(1, 1) = -2.5866667;
  expected_result_a(1, 2) = 1.6266667;
  expected_result_a(2, 0) = 0.02;
  expected_result_a(2, 1) = 2.2933333;
  expected_result_a(2, 2) = -1.3133333;
  EXPECT_TRUE(a.EqMatrix(expected_result_a));
}

TEST(InverseMatrix_suite, exceptional_test) {
  S21Matrix a(3, 3);
  a.FillingMatrix();
  ASSERT_THROW(a.InverseMatrix(), std::invalid_argument);
}

TEST(index_operator_suite, true_test) {
  S21Matrix a(3, 3);
  a.FillingMatrix();
  EXPECT_TRUE(a(0, 0) == 0);
  EXPECT_TRUE(a(0, 1) == 1);
  EXPECT_TRUE(a(0, 2) == 2);
  EXPECT_TRUE(a(1, 0) == 3);
  EXPECT_TRUE(a(1, 1) == 4);
  EXPECT_TRUE(a(1, 2) == 5);
  EXPECT_TRUE(a(2, 0) == 6);
  EXPECT_TRUE(a(2, 1) == 7);
  EXPECT_TRUE(a(2, 2) == 8);
}

TEST(index_operator_suite, exceptional_test) {
  S21Matrix a(3, 3);
  a.FillingMatrix();
  ASSERT_THROW(a(0, -1), std::out_of_range);
  ASSERT_THROW(a(-1, 0), std::out_of_range);
  ASSERT_THROW(a(-1, -1), std::out_of_range);
  ASSERT_THROW(a(10, 0), std::out_of_range);
  ASSERT_THROW(a(0, 10), std::out_of_range);
  ASSERT_THROW(a(10, 10), std::out_of_range);
}

TEST(SumMatrix_operator_suite, true_test) {
  S21Matrix a(3, 3);
  S21Matrix b(3, 3);
  S21Matrix c(3, 3);
  a.FillingMatrix();
  b.FillingMatrix();
  c.FillingMatrix();
  a = a + c;
  b = b + c;
  c.MulNumber(2);
  EXPECT_TRUE(a == c);
  EXPECT_TRUE(b == c);
  a += c;
  b += c;
  c.MulNumber(2);
  EXPECT_TRUE(a == c);
  EXPECT_TRUE(b == c);
}

TEST(SumMatrix_operator_suite, exceptional_test) {
  S21Matrix a(3, 3);
  S21Matrix b(3, 4);
  a.FillingMatrix();
  b.FillingMatrix();
  ASSERT_THROW(a += b, std::out_of_range);
  ASSERT_THROW(a = a + b, std::out_of_range);
}

TEST(SubMatrix_operator_suite, true_test) {
  S21Matrix a(3, 3);
  S21Matrix b(3, 3);
  S21Matrix c(3, 3);
  a.FillingMatrix();
  b.FillingMatrix();
  c.FillingMatrix();
  c.MulNumber(2);
  a = a - c;
  b = b - c;
  c.MulNumber(-0.5);
  EXPECT_TRUE(a == c);
  EXPECT_TRUE(b == c);
  a -= c;
  b -= c;
  c.ZeroingMatrix();
  EXPECT_TRUE(a == c);
  EXPECT_TRUE(b == c);
}

TEST(SubMatrix_operator_suite, exceptional_test) {
  S21Matrix a(3, 3);
  S21Matrix b(3, 4);
  a.FillingMatrix();
  b.FillingMatrix();
  ASSERT_THROW(a -= b, std::out_of_range);
  ASSERT_THROW(a = a - b, std::out_of_range);
}

TEST(MulMatrix_operator_suite, four_four_true_test) {
  S21Matrix a(4, 4);
  S21Matrix b(4, 4);
  S21Matrix c(4, 4);
  S21Matrix expected_result(4, 4);
  a.FillingMatrix();
  b.FillingMatrix();
  c.FillingMatrix();
  expected_result(0, 0) = 56;
  expected_result(0, 1) = 62;
  expected_result(0, 2) = 68;
  expected_result(0, 3) = 74;
  expected_result(1, 0) = 152;
  expected_result(1, 1) = 174;
  expected_result(1, 2) = 196;
  expected_result(1, 3) = 218;
  expected_result(2, 0) = 248;
  expected_result(2, 1) = 286;
  expected_result(2, 2) = 324;
  expected_result(2, 3) = 362;
  expected_result(3, 0) = 344;
  expected_result(3, 1) = 398;
  expected_result(3, 2) = 452;
  expected_result(3, 3) = 506;
  a = a * b;
  b *= c;
  EXPECT_TRUE(a == expected_result);
  EXPECT_TRUE(b == expected_result);
}

TEST(MulMatrix_operator_suite, five_two_true_test) {
  S21Matrix a(2, 5);
  S21Matrix b(5, 2);
  S21Matrix c(2, 5);
  S21Matrix expected_result(2, 2);
  a.FillingMatrix();
  b.FillingMatrix();
  c.FillingMatrix();
  expected_result(0, 0) = 60;
  expected_result(0, 1) = 70;
  expected_result(1, 0) = 160;
  expected_result(1, 1) = 195;
  a = a * b;
  c *= b;
  EXPECT_TRUE(a == expected_result);
  EXPECT_TRUE(c == expected_result);
}

TEST(MulMatrix_operator_suite, exceptional_test) {
  S21Matrix a(3, 3);
  S21Matrix b(2, 3);
  a.FillingMatrix();
  b.FillingMatrix();
  ASSERT_THROW(a = a * b, std::out_of_range);
  ASSERT_THROW(a *= b, std::out_of_range);
}

TEST(MulNumber_operator_suite, num_test) {
  S21Matrix a(3, 3);
  S21Matrix b(3, 3);
  S21Matrix c(3, 3);
  a.FillingMatrix();
  b.FillingMatrix();
  c.FillingMatrix();
  for (int i = 0; i < 4; ++i) {
    c.SumMatrix(a);
  }
  a = a * 5;
  b *= 5;
  EXPECT_TRUE(a.EqMatrix(c));
  EXPECT_TRUE(b.EqMatrix(c));
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}