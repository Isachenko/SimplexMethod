#include <stdio.h>
#include <dataStructures.h>
#include <string>

using namespace std;

size_t n; // base size
size_t m; // non base size
int* baseIndexes;
int* nonBaseIndexes;
double** a;

//coefs
const double factorX = 0.5;
const double factorY = 0.5;
const double factorZ = 10;

static Point3D vectorMultVect(const Point3D &p1, const Point3D &p2) {
    return Point3D(p1.y*p2.z - p1.z*p2.y, p1.z*p2.x - p1.x*p2.z, p1.x*p2.y - p1.y*p2.x);
}

static double vectorMultScolar(const Point3D &p1, const Point3D &p2) {
    return p1.x*p2.x + p1.y*p2.y + p1.z*p2.z;
}

// create simplexMatrix
void testInition() {
    n = 3;
    m = 4;
    baseIndexes = new int[n];
    nonBaseIndexes = new int[m];
    a = new double*[n+1];
    for(size_t i = 0; i <=n; ++i) {
        a[i] = new double[m+1];
    }
    for(size_t i = 0; i < m; ++i) {
        nonBaseIndexes[i] = i;
    }
    for(size_t i = 0; i < n; ++i) {
        baseIndexes[i] = i+m;
    }
    a[0][0] = 2;
    a[0][1] = 5;
    a[0][2] = 3;
    a[0][3] = 8;
    a[0][4] = 0;

    a[1][0] = 3;
    a[1][1] = 6;
    a[1][2] = -4;
    a[1][3] = 1;
    a[1][4] = 12;

    a[2][0] = -4;
    a[2][1] = 13;
    a[2][2] = -10;
    a[2][3] = -5;
    a[2][4] = -6;

    a[3][0] = -3;
    a[3][1] = -7;
    a[3][2] = -1;
    a[3][3] = 0;
    a[3][4] = -1;

}

void parseBrs() {
    size_t pointsCount;
    m = 3;
    //moove all construction for x,y,z >=0
    Point3D center(10, 10, 10);

    //open file
    //FILE *input = fopen("./litleTest.txt", "r");
    FILE *input = fopen("./input.txt", "r");
    if (input == NULL){
        printf("file not open\n");
        throw std::string("file can't be opened");
    }

    //read points
    fscanf(input, "%d\n", &pointsCount);
    vector<Point3D> points(pointsCount);
    for(size_t i = 0; i < pointsCount; ++i){
        fscanf(input, "%lf %lf %lf", &(points[i].x), &(points[i].y), &(points[i].z));
        points[i] = points[i] + center; // move all points
    }

    fscanf(input, "%d\n", &n);

    //simplex params init
    baseIndexes = new int[n];
    nonBaseIndexes = new int[m];
    a = new double*[n+1];
    for(size_t i = 0; i <=n; ++i) {
        a[i] = new double[m+1];
    }
    for(size_t i = 0; i < m; ++i) {
        nonBaseIndexes[i] = i;
    }
    for(size_t i = 0; i < n; ++i) {
        baseIndexes[i] = i+m;
    }
    a[0][0] = -factorX;
    a[0][1] = -factorY;
    a[0][2] = -factorZ;
    a[0][3] = -(center.x*factorX + center.y*factorY + center.z*factorZ);
    //a[0][3] = 0;


    //read simplexes and make equations
    for(size_t i = 1; i <= n; ++i){
        int p1, p2, p3;
        fscanf(input, "%d %d %d", &p1, &p2, &p3);
        Point3D v1 = points[p2] - points[p1];
        Point3D v2 = points[p3] - points[p1];
        Point3D h = vectorMultVect(v1, v2); // normal
        Point3D nullV = points[p1] - center; //
        if (vectorMultScolar(h, nullV) > 0) {
            h = vectorMultVect(v2, v1);
        }
        a[i][0] = -h.x;
        a[i][1] = -h.y;
        a[i][2] = -h.z;
        a[i][3] = -vectorMultScolar(points[p1], h);
    }
    fclose(input);
}

void transformMatrix(size_t k, size_t l) {
    int tmp = baseIndexes[k-1];
    baseIndexes[k-1] = nonBaseIndexes[l];
    nonBaseIndexes[l] = tmp;
    for(size_t i = 0; i <= n; ++i) {
        if (i == k) {
            continue;
        }
        for(size_t j = 0; j <= m; ++j){
            if (j == l) {
                continue;
            }
            a[i][j] -= (a[i][l] * a[k][j]) / a[k][l];
        }
    }
    for(size_t j = 0; j <= m; ++j){
        if (j != l) {
            a[k][j] *= 1/a[k][l];
        }
    }
    for(size_t i = 0; i <= n; ++i){
        if (i != k) {
            a[i][l] *= -1/a[k][l];
        }
    }
    a[k][l] = 1/a[k][l];

}

//return 1 - if solution find, -1 - if solution doesn't exist
int simplexMethod() {
    //step 1
    while(true){
        size_t minI = 1;
        for (size_t i = 2; i <= n; ++i) {
            if (a[i][m] < a[minI][m]) {
                minI = i;
            }
        }
        if (a[minI][m] >= 0) {
            break;
        }
        size_t minJ = 0;
        for(size_t j = 1; j < m; ++j) {
            if (a[minI][j] < a[minI][minJ]) {
                minJ = j;
            }
        }
        if (a[minI][minJ] >= 0) {
            printf("solution is not exist\n");
            return -1;
        }
        transformMatrix(minI, minJ);
    }
    //step 2
    while(true) {
        size_t minJ = 0;
        for (size_t j = 1; j < m; ++j) {
            if (a[0][j] < a[0][minJ]) {
                minJ = j;
            }
        }
        if (a[0][minJ] >= 0) {
            printf("solution found\n");
            return 1;
        }
        int minI = -1;
        for(size_t i = 1; i <= n; ++i) {
            if ((a[i][minJ] > 0) && (a[i][m] > 0)) {
                if((minI == -1) || ((a[i][m]/a[i][minJ]) < (a[minI][m]/a[minI][minJ]))) {
                    minI = i;
                }
            }
        }
        if (minI == -1) {
            printf("function is not limited\n");
            return -1;
        }
        transformMatrix(minI, minJ);
    }
}

int main() {
    //testInition();
    parseBrs();
    int result = simplexMethod();
    if (result == 1) {
        printf("ans: %lf\n", a[0][m]);
    }
    return 0;
}
