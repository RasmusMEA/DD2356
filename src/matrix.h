#ifndef MATRIX_H
#define MATRIX_H

#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <utility>

class Matrix {
 public:
  double *m_data;
  size_t m_rows;
  size_t m_cols;

  // Constructors
  Matrix(size_t rows, size_t cols);
  ~Matrix();
  Matrix(const Matrix &other);
  Matrix(Matrix &&other) noexcept;
  Matrix &operator=(const Matrix &other);
  Matrix &operator=(Matrix &&other) noexcept;

  // Accessors
  size_t rows() const;
  size_t cols() const;

  double &operator()(size_t row, size_t col);
  const double &operator()(size_t row, size_t col) const;

  // Misc
  void print() const;

  Matrix rollDown() const;
  Matrix rollUp() const;
  Matrix rollLeft() const;
  Matrix rollRight() const;

  void writeToFile(FILE *stream) const;

  // If element in array is zero, then add constant d to it
  Matrix zeroCheck(double d) const;

  static Matrix abs(Matrix a);
  static Matrix sqrt(Matrix a);
  static Matrix maximum(const Matrix &a, const Matrix &b);
  static Matrix maximum(Matrix a, double b);

  static Matrix minimum(const Matrix &a, const Matrix &b);
  static Matrix minimum(Matrix a, double b);

  static Matrix clamp(Matrix a, double min, double max);

  static Matrix sin(Matrix a);
  static Matrix exp(Matrix a);

  // return a matrix of 1's and 0's if a(i,j) < b
  static Matrix filter_lt(Matrix a, double b);

  static Matrix ones(size_t rows, size_t cols);
  static double min(const Matrix &a);
  
  Matrix& operator+=(const Matrix & o);
  Matrix& operator+=(const double o);

  Matrix& operator-=(const Matrix & o);
  Matrix& operator-=(const double o);

  Matrix& operator*=(const Matrix & o);
  Matrix& operator*=(const double o);

  Matrix& operator/=(const Matrix & o);
  Matrix& operator/=(const double o);
};

// Math operations
Matrix operator+(const Matrix &a, const Matrix &b);
Matrix operator+(const Matrix &a, double b);
Matrix operator+(double a, const Matrix &b);

Matrix operator-(const Matrix &a, const Matrix &b);
Matrix operator-(const Matrix &a, double b);
Matrix operator-(double a, const Matrix &b);

Matrix operator*(const Matrix &a, const Matrix &b);
Matrix operator*(const Matrix &a, double b);
Matrix operator*(double a, const Matrix &b);

Matrix operator/(const Matrix &a, const Matrix &b);
Matrix operator/(const Matrix &a, double b);
Matrix operator/(double a, const Matrix &b);


#endif  // MATRIX_H
