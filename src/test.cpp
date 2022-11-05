#include <gmock/gmock.h>

#include "gtest/gtest.h"
#include "s21_matrix_oop.h"

TEST(sum_matrix, true_test) {
  S21Matrix a(3, 3);
  S21Matrix b(3, 3);
  S21Matrix c(3, 3);
  a.filling_matrix();
  b.filling_matrix();
  c.filling_matrix();
  a.sum_matrix(c);
  b.sum_matrix(c);
  EXPECT_TRUE(a.eq_matrix(b));
}

TEST(sum_matrix, exceptional_test) {
  S21Matrix a(3, 3);
  S21Matrix b(3, 4);
  a.filling_matrix();
  b.filling_matrix();
  ASSERT_THROW(a.sum_matrix(b), std::out_of_range);
}

TEST(sum_matrix_operator, true_test) {
  S21Matrix a(3, 3);
  S21Matrix b(3, 3);
  S21Matrix c(3, 3);
  a.filling_matrix();
  b.filling_matrix();
  c.filling_matrix();
  a = a + c;
  b = b + c;
  c.mul_number(2);
  EXPECT_TRUE(a == c);
  EXPECT_TRUE(b == c);
  a += c;
  b += c;
  c.mul_number(2);
  EXPECT_TRUE(a == c);
  EXPECT_TRUE(b == c);
}

TEST(sum_matrix_operator, exceptional_test) {
  S21Matrix a(3, 3);
  S21Matrix b(3, 4);
  a.filling_matrix();
  b.filling_matrix();
  ASSERT_THROW(a += b, std::out_of_range);
  ASSERT_THROW(a = a + b, std::out_of_range);
}

TEST(sub_matrix, true_test) {
  S21Matrix a(3, 3);
  S21Matrix b(3, 3);
  S21Matrix c(3, 3);
  a.filling_matrix();
  b.filling_matrix();
  c.filling_matrix();

  a.sub_matrix(c);
  b.sub_matrix(c);
  EXPECT_TRUE(a.eq_matrix(b));
}

TEST(sub_matrix, exceptional_test) {
  S21Matrix a(3, 3);
  S21Matrix b(3, 4);
  a.filling_matrix();
  b.filling_matrix();
  ASSERT_THROW(a.sub_matrix(b), std::out_of_range);
}

TEST(sub_matrix_operator, true_test) {
  S21Matrix a(3, 3);
  S21Matrix b(3, 3);
  S21Matrix c(3, 3);
  a.filling_matrix();
  b.filling_matrix();
  c.filling_matrix();
  c.mul_number(2);
  a = a - c;
  b = b - c;
  c.mul_number(-0.5);
  EXPECT_TRUE(a == c);
  EXPECT_TRUE(b == c);
  a -= c;
  b -= c;
  c.zeroing_matrix();
  EXPECT_TRUE(a == c);
  EXPECT_TRUE(b == c);
}

TEST(sub_matrix_operator, exceptional_test) {
  S21Matrix a(3, 3);
  S21Matrix b(3, 4);
  a.filling_matrix();
  b.filling_matrix();
  ASSERT_THROW(a -= b, std::out_of_range);
  ASSERT_THROW(a = a - b, std::out_of_range);
}

TEST(mul_number, true_test) {
  S21Matrix a(3, 3);
  S21Matrix b(3, 3);
  a.filling_matrix();
  b.filling_matrix();
  for (int i = 0; i < 4; ++i) {
    b.sum_matrix(a);
  }
  a.mul_number(5);
  EXPECT_TRUE(a.eq_matrix(b));
}

TEST(mul_matrix, three_three_true_test) {
  S21Matrix a(3, 3);
  S21Matrix b(3, 3);
  S21Matrix expected_result(3, 3);
  a.filling_matrix();
  b.filling_matrix();
  expected_result(0, 0) = 15;
  expected_result(0, 1) = 18;
  expected_result(0, 2) = 21;
  expected_result(1, 0) = 42;
  expected_result(1, 1) = 54;
  expected_result(1, 2) = 66;
  expected_result(2, 0) = 69;
  expected_result(2, 1) = 90;
  expected_result(2, 2) = 111;
  a.mul_matrix(b);
  EXPECT_TRUE(a.eq_matrix(expected_result));
}

TEST(mul_matrix, four_four_true_test) {
  S21Matrix a(4, 4);
  S21Matrix b(4, 4);
  S21Matrix expected_result(4, 4);
  a.filling_matrix();
  b.filling_matrix();
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
  a.mul_matrix(b);
  EXPECT_TRUE(a.eq_matrix(expected_result));
}

TEST(mul_matrix, five_two_true_test) {
  S21Matrix a(2, 5);
  S21Matrix b(5, 2);
  a.filling_matrix();
  b.filling_matrix();
  S21Matrix expected_result(2, 2);
  expected_result(0, 0) = 60;
  expected_result(0, 1) = 70;
  expected_result(1, 0) = 160;
  expected_result(1, 1) = 195;
  a.mul_matrix(b);
  EXPECT_TRUE(a.eq_matrix(expected_result));
}

TEST(mul_matrix, exceptional_test) {
  S21Matrix a(3, 3);
  S21Matrix b(2, 3);
  a.filling_matrix();
  b.filling_matrix();
  ASSERT_THROW(a.mul_matrix(b), std::out_of_range);
}

TEST(mul_matrix_operator, four_four_true_test) {
  S21Matrix a(4, 4);
  S21Matrix b(4, 4);
  S21Matrix c(4, 4);
  S21Matrix expected_result(4, 4);
  a.filling_matrix();
  b.filling_matrix();
  c.filling_matrix();
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

TEST(mul_matrix_operator, five_two_true_test) {
  S21Matrix a(2, 5);
  S21Matrix b(5, 2);
  S21Matrix c(2, 5);
  S21Matrix expected_result(2, 2);
  a.filling_matrix();
  b.filling_matrix();
  c.filling_matrix();
  expected_result(0, 0) = 60;
  expected_result(0, 1) = 70;
  expected_result(1, 0) = 160;
  expected_result(1, 1) = 195;
  a = a * b;
  c *= b;
  EXPECT_TRUE(a == expected_result);
  EXPECT_TRUE(c == expected_result);
}

TEST(mul_matrix_operator, exceptional_test) {
  S21Matrix a(3, 3);
  S21Matrix b(2, 3);
  a.filling_matrix();
  b.filling_matrix();
  ASSERT_THROW(a = a * b, std::out_of_range);
  ASSERT_THROW(a *= b, std::out_of_range);
}

TEST(transpose_matrix, true_test) {
  S21Matrix a(2, 3);
  S21Matrix expected_result(3, 2);
  a.filling_matrix();
  a = a.transpose();
  expected_result(0, 0) = 0;
  expected_result(0, 1) = 3;
  expected_result(1, 0) = 1;
  expected_result(1, 1) = 4;
  expected_result(2, 0) = 2;
  expected_result(2, 1) = 5;
  EXPECT_TRUE(a.eq_matrix(expected_result));
}

TEST(calc_complemet, four_four_test) {
  S21Matrix b(4, 4);
  b.filling_matrix();
  b(0, 0) = 50;
  S21Matrix expected_result(4, 4);
  b = b.calc_complements();
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
  EXPECT_TRUE(b.eq_matrix(expected_result));
}

TEST(calc_complemet, three_three_test) {
  S21Matrix a(3, 3);
  a.filling_matrix();
  S21Matrix expected_result(3, 3);
  a = a.calc_complements();
  expected_result(0, 0) = -3.000000;
  expected_result(0, 1) = 6.000000;
  expected_result(0, 2) = -3.000000;
  expected_result(1, 0) = 6.000000;
  expected_result(1, 1) = -12.000000;
  expected_result(1, 2) = 6.000000;
  expected_result(2, 0) = -3.000000;
  expected_result(2, 1) = 6.000000;
  expected_result(2, 2) = -3.000000;
  EXPECT_TRUE(a.eq_matrix(expected_result));
}

TEST(calc_complemet, exceptional_test) {
  S21Matrix a(3, 4);
  a.filling_matrix();
  ASSERT_THROW(a.calc_complements(), std::out_of_range);
}

TEST(determinant, three_three_test) {
  S21Matrix a(3, 3);
  S21Matrix b(3, 3);
  a.filling_matrix();
  b.filling_matrix();
  b(2, 2) = 56;
  S21Matrix c(b);
  c(0, 0) = 1000;
  EXPECT_TRUE(a.determinant() == 0);
  EXPECT_TRUE(b.determinant() == -144);
  EXPECT_TRUE(c.determinant() == 188856);
}

TEST(determinant, four_four_test) {
  S21Matrix a(4, 4);
  S21Matrix b(4, 4);
  a.filling_matrix();
  b.filling_matrix();
  b(0, 0) = 100;
  b(3, 3) = 56;
  S21Matrix c(b);
  c(0, 0) = 1000;
  EXPECT_TRUE(a.determinant() == 0);
  EXPECT_TRUE(b.determinant() == -16400);
  EXPECT_TRUE(c.determinant() == -164000);
}

TEST(determinant, exceptional_test) {
  S21Matrix a(3, 4);
  a.filling_matrix();
  ASSERT_THROW(a.determinant(), std::out_of_range);
}

TEST(inverse_matrix, three_three_test) {
  S21Matrix a(3, 3);
  S21Matrix expected_result_a(3, 3);
  a.filling_matrix();
  a(0, 0) = 50;
  a = a.inverse_matrix();
  expected_result_a(0, 0) = 0.02;
  expected_result_a(0, 1) = -0.04;
  expected_result_a(0, 2) = 0.02;
  expected_result_a(1, 0) = -0.04;
  expected_result_a(1, 1) = -2.5866667;
  expected_result_a(1, 2) = 1.6266667;
  expected_result_a(2, 0) = 0.02;
  expected_result_a(2, 1) = 2.2933333;
  expected_result_a(2, 2) = -1.3133333;
  EXPECT_TRUE(a.eq_matrix(expected_result_a));
}

TEST(inverse_matrix, exceptional_test) {
  S21Matrix a(3, 3);
  a.filling_matrix();
  ASSERT_THROW(a.inverse_matrix(), std::invalid_argument);
}

TEST(index_operator, true_test) {
  S21Matrix a(3, 3);
  a.filling_matrix();
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

TEST(index_operator, exceptional_test) {
  S21Matrix a(3, 3);
  a.filling_matrix();
  ASSERT_THROW(a(0, -1), std::out_of_range);
  ASSERT_THROW(a(-1, 0), std::out_of_range);
  ASSERT_THROW(a(-1, -1), std::out_of_range);
  ASSERT_THROW(a(10, 0), std::out_of_range);
  ASSERT_THROW(a(0, 10), std::out_of_range);
  ASSERT_THROW(a(10, 10), std::out_of_range);
}

TEST(eq_matrix, true_test) {
  S21Matrix a(3, 3);
  S21Matrix b(3, 3);
  a.filling_matrix();
  b.filling_matrix();

  EXPECT_TRUE(a.eq_matrix(b));
  EXPECT_TRUE(b.eq_matrix(a));

  EXPECT_TRUE(a == b);
  EXPECT_TRUE(b == a);
}

TEST(eq_matrix, false_test) {
  S21Matrix a(2, 3);
  S21Matrix b(3, 3);
  a.filling_matrix();
  b.filling_matrix();

  EXPECT_FALSE(a.eq_matrix(b));
  EXPECT_FALSE(b.eq_matrix(a));

  EXPECT_FALSE(a == b);
  EXPECT_FALSE(b == a);
}

TEST(eq_matrix, accuracy_true_test) {
  S21Matrix a(3, 3);
  S21Matrix b(3, 3);
  a.filling_matrix();
  b.filling_matrix();
  for (int i = 0; i < a.get_rows(); i++) {
    for (int j = 0; j < a.get_cols(); j++) {
      a(i, j) += 0.00000001;
    }
  }

  EXPECT_TRUE(a.eq_matrix(b));
  EXPECT_TRUE(a == b);
}

TEST(eq_matrix, accuracy_false_test) {
  S21Matrix a(3, 3);
  S21Matrix b(3, 3);
  a.filling_matrix();
  b.filling_matrix();
  for (int i = 0; i < a.get_rows(); i++) {
    for (int j = 0; j < a.get_cols(); j++) {
      a(i, j) += 0.000001;
    }
  }
  EXPECT_FALSE(a.eq_matrix(b));
  EXPECT_FALSE(a == b);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}