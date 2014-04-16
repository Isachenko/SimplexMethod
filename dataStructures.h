#ifndef DATDSTRUCTURES_H
#define DATDSTRUCTURES_H

#include <memory>
#include <vector>

struct Point3D {
    double x,y,z;

    Point3D() : x(0), y(0), z(0) {}

    Point3D(double a, double b, double c) :
       x(a), y(b), z(c) {}

    Point3D operator - (const Point3D &p) const {
        return Point3D(this->x - p.x, this->y - p.y, this->z - p.z);
    }

    Point3D operator + (const Point3D &p) const {
        return Point3D(this->x + p.x, this->y + p.y, this->z + p.z);
    }
};

struct Simplex2D {
    size_t p1, p2, p3;

};

#endif // DATDSTRUCTURES_H
