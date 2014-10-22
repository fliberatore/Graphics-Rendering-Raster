#if !defined MATRIX_H
#define MATRIX_H
#include "gz.h"

typedef GzMatrix4 float[4][4];
typedef GzVector3 float[3];
typedef GzVector4 float[4];

//Matrix multiplication. Result is returned in res.
int mult(GzMatrix4 A, GzMatrix4 B, GzMatrix4 &res);

//Matrix-vector multiplication. Result is returned in res.
int mult(GzMatrix4 A, GzVector4 b, GzVector4 &res);

//Transforms a Vector4 to a Vector3. Result is returned in res.
int toVector3(GzVector4 a, GzVector3 &res);

//Transforms a Vector3 to a Vector 4. Result is returned in res.
int toVector4(GzVector3 a, GzVector4 &res);

//Dot product of two Vector4. Result is returned in res.
int dotProd(GzVector4 a, GzVector4 b, float &res);

//Cross product of two Vector4. Result is returned in res.
int crossProd(GzVector4 a, GzVector4 b, GzVector4 &res);

#endif