#include <math.h>

float dxfunc3(float x) {
	float ans;
	ans = -sinf(x + sqrt(2)) + x + sqrt(2);
	return ans;
}