#pragma once

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <iomanip>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/string.hpp>

// this class holds the data for each point along the spline
class Point {
public:
  // default constructor
  Point() : x(0.0), a(0.0), b(0.0), c(0.0), d(0.0) {}

  // data
  double x, a, b, c, d;

private:
  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version) {
    ar & x;
    ar & a;
    ar & b;
    ar & c;
    ar & d;
  }
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

  // method to read data from file
  void readFile(const char * dataFile);

  // method to get the value at a certain point
  double value(double x);
  double getFirstX();
  double getLastX();
  // method to print the contents of the spline (for debug)
  void print();

  // data
  std::vector<Point> s;

private:
  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version) {
    ar & s;
  }
};

bool comparePoints(const Point &pa, const Point &pb);