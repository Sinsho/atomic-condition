#include <math.h>

double foo(double x, double y, double z) {
    double m = x * y +z;
    double k = m * z;
    double val = (1-cos(x*y))/(x*x*k);
    return val;
}

