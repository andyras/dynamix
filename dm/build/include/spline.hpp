#ifndef __SPLINE__
#define __SPLINE__

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <iomanip>

// this class holds the data for each point along the spline
class Point {
public:
  // default constructor
  Point() : x(0.0), a(0.0), b(0.0), c(0.0), d(0.0) {}

  // copy-constructor
  Point(const Point& o) : x(o.x), a(o.a), b(o.b), c(o.c), d(o.d) {}

  // destructor
  ~Point() {}

  // swap
  friend void swap(Point& first, Point& second) {
    using std::swap;

    swap(first.x, second.x);
    swap(first.a, second.a);
    swap(first.b, second.b);
    swap(first.c, second.c);
    swap(first.d, second.d);
  }

  // assignment
  Point& operator=(Point other) {
    swap(*this, other);
    return *this;
  }

  // move constructor
  Point(Point&& other) : Point() { // initialize via default constructor, C++11 only
      swap(*this, other);
  }

  // data
  double x, a, b, c, d;
};

class Spline {
public:
  // default constructor
  Spline() : s(std::vector<Point> (0)) {}

  // copy-constructor
  Spline(const Spline& o) : s(o.s) {}

  // destructor
  ~Spline() {}

  // swap
  friend void swap(Spline& first, Spline& second) {
    using std::swap;

    swap(first.s, second.s);
  }

  // assignment
  Spline& operator=(Spline other) {
    swap(*this, other);
    return *this;
  }

  // move constructor
  Spline(Spline&& other) : Spline() { // initialize via default constructor, C++11 only
      swap(*this, other);
  }

  Spline(const char * dataFile);
  // method to get the value at a certain point
  double value(double x);
  double getFirstX();
  double getLastX();
  // method to print the contents of the spline (for debug)
  void print();
  std::vector<Point> s;
public:
};

bool comparePoints(const Point &pa, const Point &pb);

#endif
