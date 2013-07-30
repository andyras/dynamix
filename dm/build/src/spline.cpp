#include "spline.hpp"

#define DEBUG_SPLINE

Spline::Spline(const char * dataFile) {
  if (dataFile == NULL) {
    std::cout << "ERROR: specify an input file" << std::endl;
  }
#ifdef DEBUG_SPLINE
  std::cout << "Creating spline from data file " << dataFile << "." << std::endl;
#endif

  // placeholders for x and a values
  double xTmp;
  double aTmp;

  Point p;

  // stream for input data file
  std::ifstream data(dataFile, std::ios::in);

  // for each line, read in the x and y point values
  while (data >> p.x >> p.a) {
    s.emplace_back(p);
  }

  // sort the vector of values to ensure they are ordered by x value
  std::sort(s.begin(), s.end(), comparePoints);

  data.close();

  // vectors which will be used while calculating spline parameters
  // method from:
  // http://banach.millersville.edu/~BobBuchanan/math375/CubicSpline/main.pdf
  std::vector<double> h (s.size(), 0.0);
  std::vector<double> al (s.size(), 0.0);
  std::vector<double> I (s.size(), 0.0);
  std::vector<double> mu (s.size(), 0.0);
  std::vector<double> z (s.size(), 0.0);

  //// Step 1: compute differences between x points

  for (int ii = 0; ii < (s.size() - 1); ii++) {
  std::cout << "NO BUGS HERE\n";
    h[ii] = s[ii+1].x - s[ii].x;
    // check that there are no duplicates
    if (h[ii] == 0.0) {
      std::cout << "ERROR: points number " << ii << " and " << (ii+1)
	<< " have the same x value." << std::endl;
    }
  }

  //// Step 2: compute alpha
  for (int ii = 1; ii < (s.size() - 1); ii++) {
    al[ii] = (3.0/h[ii])*(s[ii+1].a - s[ii].a) - (3.0/h[ii-1])*(s[ii].a - s[ii-1].a);
  }

  //// Step 3: set I, mu, and z for first point
  I[0] = 1.0;
  mu[0] = 0.0;
  z[0] = 0.0;

  //// Step 4: set I, mu, and z for other points
  for (int ii = 1; ii < (s.size() - 1); ii++) {
    I[ii] = 2.0*(s[ii+1].x - s[ii-1].x) - h[ii-1]*mu[ii-1];
    mu[ii] = h[ii]/I[ii];
    z[ii] = (al[ii] - h[ii-1]*z[ii-1])/I[ii];
  }

  //// Step 5: set I, c, and z for last point
  I.back() = 1;
  // second derivative at last point is zero
  s.back().c = 0;
  z.back() = 0;


  //// Step 6: set b, c, and d for other points
  for (int ii = (s.size() - 2); ii >= 0; ii--) {
    s[ii].c = z[ii] - mu[ii]*s[ii+1].c;
    s[ii].b = (s[ii+1].a - s[ii].a)/h[ii] - h[ii]*(s[ii+1].c + 2*s[ii].c)/3.0;
    s[ii].d = (s[ii+1].c - s[ii].c)/(3.0*h[ii]);
  }

  print();
}

// used with std::sort, gives vector with smallest x first
bool comparePoints(const Point &pa, const Point &pb) {
  return pa.x < pb.x;
}

double Spline::value(double x) {
  if ((x < s.front().x) || (x > s.back().x)) {
    std::cout << "ERROR: value too large." << std::endl;
  }

  // index of point/spline to use
  int index = 0;

  // s[ii].x <= x < s[ii+1].x
  for (int ii = 0; ii < (s.size() - 1); ii++) {
    if ((x >= s[ii].x) && (x < s[ii+1].x)) {
      index = ii;
      break;
    }
  }
  // the point at the end of the range is interpolated by the last spline
  // s[N-1].x <= x <= s[N].x
  if (x == s.back().x) {
    index = (s.size() - 2);
  }

  // this is how we get the cubic interpolated value
  return ((((s[index].d*x) + s[index].c)*x + s[index].b)*x + s[index].a);
}

void Spline::print() {
  std::cout << std::setw(7) << "index ";
  std::cout << std::setw(13) << "x      ";
  std::cout << std::setw(13) << "a      ";
  std::cout << std::setw(13) << "b      ";
  std::cout << std::setw(13) << "c      ";
  std::cout << std::setw(13) << "d      ";
  std::cout << std::endl;
  for (int ii = 0; ii < s.size(); ii++) {
    std::cout << std::setw(6) << std::scientific << ii << " ";
    std::cout << std::setw(13) << std::scientific << s[ii].x;
    std::cout << std::setw(13) << std::scientific << s[ii].a;
    std::cout << std::setw(13) << std::scientific << s[ii].b;
    std::cout << std::setw(13) << std::scientific << s[ii].c;
    std::cout << std::setw(13) << std::scientific << s[ii].d;
    std::cout << std::endl;
  }
}
