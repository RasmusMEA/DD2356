#include "matrix.h"

#include <cmath>
#include <cstdio>
#include <stdexcept>

Matrix operator*(const Matrix a, const Matrix b) {
  if (a.rows() != b.rows() || a.cols() != b.cols()) {
    throw std::out_of_range("Mismatched dimensions for element_mul");
  }

  Matrix res = a;

  for (size_t y = 0; y < a.rows(); ++y) {
    for (size_t x = 0; x < a.cols(); ++x) {
      res(y, x) *= b(y, x);
    }
  }

  return res;
}

Matrix operator*(const Matrix a, double b) {
  Matrix res(a);

  for (size_t y = 0; y < a.rows(); ++y) {
    for (size_t x = 0; x < a.cols(); ++x) {
      res(y, x) *= b;
    }
  }

  return res;
}

Matrix operator*(double a, const Matrix b) { return b * a; }

Matrix operator+(const Matrix a, const Matrix b) {
  if (a.rows() != b.rows() || a.cols() != b.cols()) {
    throw std::out_of_range("Mismatched dimensions for element_add");
  }

  Matrix res = a;

  for (size_t y = 0; y < a.rows(); ++y) {
    for (size_t x = 0; x < a.cols(); ++x) {
      res(y, x) += b(y, x);
    }
  }

  return res;
}

Matrix operator+(const Matrix a, double b) {
  Matrix res(a);

  for (size_t y = 0; y < a.rows(); ++y) {
    for (size_t x = 0; x < a.cols(); ++x) {
      res(y, x) += b;
    }
  }

  return res;
}

Matrix operator+(double a, const Matrix b) { return b + a; }

Matrix operator-(const Matrix a, const Matrix b) {
  if (a.rows() != b.rows() || a.cols() != b.cols()) {
    throw std::out_of_range("Mismatched dimensions for element_sub");
  }

  Matrix res = a;

  for (size_t y = 0; y < a.rows(); ++y) {
    for (size_t x = 0; x < a.cols(); ++x) {
      res(y, x) -= b(y, x);
    }
  }

  return res;
}

Matrix operator-(const Matrix a, double b) {
  Matrix res(a);

  for (size_t y = 0; y < a.rows(); ++y) {
    for (size_t x = 0; x < a.cols(); ++x) {
      res(y, x) -= b;
    }
  }

  return res;
}

Matrix operator-(double a, const Matrix b) {
  Matrix res(b);

  for (size_t y = 0; y < b.rows(); ++y) {
    for (size_t x = 0; x < b.cols(); ++x) {
      res(y, x) = a - res(y, x);
    }
  }

  return res;
}

Matrix operator/(const Matrix a, const Matrix b) {
  if (a.rows() != b.rows() || a.cols() != b.cols()) {
    throw std::out_of_range("Mismatched dimensions for element_div");
  }

  Matrix res = a;

  for (size_t y = 0; y < a.rows(); ++y) {
    for (size_t x = 0; x < a.cols(); ++x) {
      if (b(y, x) == 0) {
        res(y, x) = 0;
      } else {
        res(y, x) /= b(y, x);
      }
    }
  }

  return res;
}

Matrix operator/(const Matrix a, double b) {
  Matrix res(a);

  bool set_all_zero = b == 0;

  for (size_t y = 0; y < a.rows(); ++y) {
    for (size_t x = 0; x < a.cols(); ++x) {
      if (set_all_zero) {
        res(y, x) = 0;
      } else {
        res(y, x) /= b;
      }
    }
  }

  return res;
}

Matrix operator/(double a, const Matrix b) {
  Matrix res(b);

  for (size_t y = 0; y < b.rows(); ++y) {
    for (size_t x = 0; x < b.cols(); ++x) {
      if (res(y, x) == 0) {
        res(y, x) = 0;
      } else {
        res(y, x) = a / res(y, x);
      }
    }
  }

  return res;
}

// Constructors
Matrix::Matrix(size_t rows, size_t cols) : m_rows(rows), m_cols(cols) {
  m_data = new double[m_rows * m_cols]{};
}

Matrix::~Matrix() { delete[] m_data; }

// copy constructor
Matrix::Matrix(const Matrix &other)
    : m_rows(other.m_rows), m_cols(other.m_cols) {
  m_data = new double[m_rows * m_cols];

  for (size_t i = 0; i < other.m_cols * other.m_rows; ++i) {
    m_data[i] = other.m_data[i];
  }
}

// move constructor
Matrix::Matrix(Matrix &&other) noexcept {
  if (this == &other) {
    return;
  }

  m_cols = std::exchange(other.m_cols, 0);
  m_rows = std::exchange(other.m_rows, 0);
  m_data = std::exchange(other.m_data, nullptr);
}

// copy assignment operator
Matrix &Matrix::operator=(const Matrix &other) {
  if (this == &other) {
    return *this;
  }

  if (m_rows * m_cols != 0) {
    delete[] m_data;
  }

  m_data = new double[m_rows * m_cols];

  for (size_t i = 0; i < other.m_cols * other.m_rows; ++i) {
    m_data[i] = other.m_data[i];
  }

  m_cols = other.m_cols;
  m_rows = other.m_rows;

  return *this;
}

// move assignment operator
Matrix &Matrix::operator=(Matrix &&other) noexcept {
  if (this == &other) {
    return *this;
  }

  delete[] m_data;

  m_cols = std::exchange(other.m_cols, 0);
  m_rows = std::exchange(other.m_rows, 0);
  m_data = std::exchange(other.m_data, nullptr);

  return *this;
}

// Accessors
size_t Matrix::rows() const { return m_rows; }

size_t Matrix::cols() const { return m_cols; }

double &Matrix::operator()(size_t row, size_t col) {
  if (row >= m_rows || col >= m_cols || row < 0 || col < 0) {
    throw std::out_of_range("Accessed element outside array");
  }

  return m_data[row * m_cols + col];
}

const double &Matrix::operator()(size_t row, size_t col) const {
  if (row >= m_rows || col >= m_cols || row < 0 || col < 0) {
    throw std::out_of_range("Accessed element outside array");
  }

  return m_data[row * m_cols + col];
}

// Misc
void Matrix::print() const {
  for (size_t y = 0; y < m_rows; ++y) {
    for (size_t x = 0; x < m_cols; ++x) {
      std::cout << m_data[y * m_cols + x] << " ";
    }
    std::cout << std::endl;
  }
}

Matrix Matrix::zeroCheck(double d) const {
  Matrix res(*this);
  for (size_t i = 0; i < m_cols * m_rows; ++i) {
    if (res.m_data[i] == 0) {
      res.m_data[i] = d;
    }
  }

  return res;
}

Matrix Matrix::rollDown() const {
  Matrix res(m_rows, m_cols);

  for (int y = 0; y < m_rows; ++y) {
    for (int x = 0; x < m_cols; ++x) {
      res.m_data[((y + 1) % m_rows) * m_cols + x] = m_data[y * m_cols + x];
    }
  }

  return res;
}

Matrix Matrix::rollUp() const {
  Matrix res(m_rows, m_cols);

  for (int y = 0; y < m_rows; ++y) {
    for (int x = 0; x < m_cols; ++x) {
      res.m_data[((y + m_rows - 1) % m_rows) * m_cols + x] =
          m_data[y * m_cols + x];
    }
  }

  return res;
}

Matrix Matrix::rollLeft() const {
  Matrix res(m_rows, m_cols);

  for (int y = 0; y < m_rows; ++y) {
    for (int x = 0; x < m_cols; ++x) {
      res.m_data[y * m_cols + ((x - 1 + m_cols) % m_cols)] =
          m_data[y * m_cols + x];
    }
  }

  return res;
}

Matrix Matrix::rollRight() const {
  Matrix res(m_rows, m_cols);

  for (int y = 0; y < m_rows; ++y) {
    for (int x = 0; x < m_cols; ++x) {
      res.m_data[y * m_cols + (x + 1) % m_cols] = m_data[y * m_cols + x];
    }
  }

  return res;
}

Matrix Matrix::maximum(const Matrix &a, const Matrix &b) {
  if (a.rows() != b.rows() || a.cols() != b.cols()) {
    throw std::out_of_range("Mismatched dimensions in maximum()");
  }

  Matrix res(a);

  for (size_t y = 0; y < a.rows(); ++y) {
    for (size_t x = 0; x < a.cols(); ++x) {
      if (b(y, x) > res(y, x)) {
        res(y, x) = b(y, x);
      }
    }
  }

  return res;
}

Matrix Matrix::maximum(Matrix a, double b) {
  for (size_t y = 0; y < a.rows(); ++y) {
    for (size_t x = 0; x < a.cols(); ++x) {
      if (b > a(y, x)) {
        a(y, x) = b;
      }
    }
  }

  return a;
}

Matrix Matrix::minimum(const Matrix &a, const Matrix &b) {
  if (a.rows() != b.rows() || a.cols() != b.cols()) {
    throw std::out_of_range("Mismatched dimensions in maximum()");
  }

  Matrix res(a);

  for (size_t y = 0; y < a.rows(); ++y) {
    for (size_t x = 0; x < a.cols(); ++x) {
      if (b(y, x) < res(y, x)) {
        res(y, x) = b(y, x);
      }
    }
  }

  return res;
}

Matrix Matrix::minimum(Matrix a, double b) {
  for (size_t y = 0; y < a.rows(); ++y) {
    for (size_t x = 0; x < a.cols(); ++x) {
      if (b < a(y, x)) {
        a(y, x) = b;
      }
    }
  }

  return a;
}

Matrix Matrix::sqrt(Matrix a) {
  for (size_t i = 0; i < a.m_rows * a.m_cols; ++i) {
    if (a.m_data[i] < 0) {
      std::clog << "Warning: sqrt with negative number" << std::endl;
    }
    a.m_data[i] = std::sqrt(a.m_data[i]);
  }

  return a;
}

Matrix Matrix::abs(Matrix a) {
  for (size_t i = 0; i < a.m_rows * a.m_cols; ++i) {
    a.m_data[i] = std::abs(a.m_data[i]);
  }

  return a;
}

Matrix Matrix::sin(Matrix a) {
  for (size_t i = 0; i < a.m_rows * a.m_cols; ++i) {
    a.m_data[i] = std::sin(a.m_data[i]);
  }
  return a;
}

Matrix Matrix::exp(Matrix a) {
  for (size_t i = 0; i < a.m_rows * a.m_cols; ++i) {
    a.m_data[i] = std::pow(M_E, a.m_data[i]);
  }
  return a;
}

Matrix Matrix::filter_lt(Matrix a, double b) {
  for (size_t i = 0; i < a.m_rows * a.m_cols; ++i) {
    a.m_data[i] = a.m_data[i] < b ? 1 : 0;
  }
  return a;
}

Matrix Matrix::ones(size_t rows, size_t cols) {
  Matrix res(rows, cols);

  for (size_t i = 0; i < rows * cols; ++i) {
    res.m_data[i] = 1;
  }

  return res;
}

double Matrix::min(const Matrix &a) {
  double minval = std::numeric_limits<double>::max();

  for (size_t i = 0; i < a.m_rows * a.m_cols; ++i) {
    if (minval > a.m_data[i]) {
      minval = a.m_data[i];
    }
  }

  return minval;
}

void Matrix::writeToFile(FILE *stream) const {
  fwrite(m_data, sizeof(double), m_rows * m_cols, stream);
}
