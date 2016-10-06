#include <stdio.h>

#define _USE_MATH_DEFINES
#include <math.h>
#include "nrutil.h"

#define JMAX 40
#define MAXIT 30
#define JMAX_1 20
#define MAXIT_1 100

int count = 0;

float sgn(float x) {
	if (x >= 0)
		return 1.0;
	if (x < 0)
		return -1.0;
}

float rtbis(float(*func)(float), float x1, float x2, float xacc)
{
	void nrerror(char error_text[]);
	int j;
	float dx, f, fmid, xmid, rtb;

	f = (*func)(x1);
	fmid = (*func)(x2);
	if (f*fmid >= 0.0) nrerror("Root must be bracketed for bisection in rtbis");
	rtb = f < 0.0 ? (dx = x2 - x1, x1) : (dx = x1 - x2, x2);
	for (j = 1; j <= JMAX; j++) {
		count++;
		fmid = (*func)(xmid = rtb + (dx *= 0.5));
		if (fmid <= 0.0) rtb = xmid;
		if (fabs(dx) / rtb * 100  < xacc || fmid == 0.0) return rtb;
	}
	nrerror("Too many bisections in rtbis");
	return 0.0;
}

float rtflsp(float(*func)(float), float x1, float x2, float xacc)
{
	void nrerror(char error_text[]);
	int j;
	float fl, fh, xl, xh, swap, dx, del, f, rtf;

	fl = (*func)(x1);
	fh = (*func)(x2);
	if (fl*fh > 0.0) nrerror("Root must be bracketed in rtflsp");
	if (fl < 0.0) {
		xl = x1;
		xh = x2;
	}
	else {
		xl = x2;
		xh = x1;
		swap = fl;
		fl = fh;
		fh = swap;
	}
	dx = xh - xl;
	for (j = 1; j <= MAXIT; j++) {
		count++;
		rtf = xl + dx*fl / (fl - fh);
		f = (*func)(rtf);
		if (f < 0.0) {
			del = xl - rtf;
			xl = rtf;
			fl = f;
		}
		else {
			del = xh - rtf;
			xh = rtf;
			fh = f;
		}
		dx = xh - xl;
		if (fabs(del)/rtf * 100 < xacc || f == 0.0) return rtf;
	}
	nrerror("Maximum number of iterations exceeded in rtflsp");
	return 0.0;
}

float rtsec(float(*func)(float), float x1, float x2, float xacc)
{
	void nrerror(char error_text[]);
	int j;
	float fl, f, dx, swap, xl, rts;

	fl = (*func)(x1);
	f = (*func)(x2);
	if (fabs(fl) < fabs(f)) {
		rts = x1;
		xl = x2;
		swap = fl;
		fl = f;
		f = swap;
	}
	else {
		xl = x1;
		rts = x2;
	}
	for (j = 1; j <= MAXIT; j++) {
		count++;
		dx = (xl - rts)*f / (f - fl);
		xl = rts;
		fl = f;
		rts += dx;
		f = (*func)(rts);
		if (fabs(dx) / rts * 100  < xacc || f == 0.0) return rts;
	}
	nrerror("Maximum number of iterations exceeded in rtsec");
	return 0.0;
}

float rtmuller(float(*func)(float), float x1, float x2, float xacc) {

	void nrerror(char err_text[]);
	int j;
	float p_0, p_1, p_2, rt, dx;
	float f_0, f_1, f_2;
	float a, b, c;

	p_0 = x1;
	p_2 = x2;
	p_1 = (x1 + x2) / 2.0;


	for (j = 1; j <= JMAX_1; j++) {

		count++;

		f_0 = (*func)(p_0);
		f_1 = (*func)(p_1);
		f_2 = (*func)(p_2);

		c = f_2;

		b = ((p_0 - p_2)*(p_0 - p_2)*(f_1 - f_2) - (p_1 - p_2)*(p_1 - p_2)*(f_0 - f_2))
			/ ((p_0 - p_2)*(p_1 - p_2)*(p_0 - p_1));

		a = ((p_1 - p_2)*(f_0 - f_2) - (p_0 - p_2)*(f_1 - f_2))
			/ ((p_0 - p_2)*(p_1 - p_2)*(p_0 - p_1));

		p_0 = p_1;
		p_1 = p_2;

		dx = 2.0f * c / (b + sgn(b)*(float)sqrt(b*b - 4.0*a*c));

		rt = p_2 - dx;

		p_2 = rt;

		if (fabs(dx) / rt * 100  < xacc) return rt;
	}
	nrerror("Maximum number of iterations exceeded in rtnewt");
	return 0.0;

}

float rtnewt(void(*funcd)(float, float *, float *), float x1, float x2,
	float xacc)
{
	void nrerror(char error_text[]);
	int j;
	float df, dx, f, rtn;

	rtn = 0.5*(x1 + x2);
	for (j = 1; j <= JMAX_1; j++) {
		count++;
		(*funcd)(rtn, &f, &df);
		dx = f / df;
		rtn -= dx;
		if ((x1 - rtn)*(rtn - x2) < 0.0)
			nrerror("Jumped out of brackets in rtnewt");
		if (fabs(dx) / rtn * 100  < xacc) return rtn;
	}
	nrerror("Maximum number of iterations exceeded in rtnewt");
	return 0.0;
}

float rtsafe(void(*funcd)(float, float *, float *), float x1, float x2,
	float xacc)
{
	void nrerror(char error_text[]);
	int j;
	float df, dx, dxold, f, fh, fl;
	float temp, xh, xl, rts;

	(*funcd)(x1, &fl, &df);
	(*funcd)(x2, &fh, &df);
	if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))
		nrerror("Root must be bracketed in rtsafe");
	if (fl == 0.0) return x1;
	if (fh == 0.0) return x2;
	if (fl < 0.0) {
		xl = x1;
		xh = x2;
	}
	else {
		xh = x1;
		xl = x2;
	}
	rts = 0.5*(x1 + x2);
	dxold = fabs(x2 - x1);
	dx = dxold;
	(*funcd)(rts, &f, &df);
	for (j = 1; j <= MAXIT_1; j++) {
		count++;
		if ((((rts - xh)*df - f)*((rts - xl)*df - f) > 0.0)
			|| (fabs(2.0*f) > fabs(dxold*df))) {
			dxold = dx;
			dx = 0.5*(xh - xl);
			rts = xl + dx;
			if (xl == rts) return rts;
		}
		else {
			dxold = dx;
			dx = f / df;
			temp = rts;
			rts -= dx;
			if (temp == rts) return rts;
		}
		if (fabs(dx) / rts * 100  < xacc) return rts;
		(*funcd)(rts, &f, &df);
		if (f < 0.0)
			xl = rts;
		else
			xh = rts;
	}
	nrerror("Maximum number of iterations exceeded in rtsafe");
	return 0.0;
}

float bessj0(float x)
{
	float ax, z;
	double xx, y, ans, ans1, ans2;

	if ((ax = fabs(x)) < 8.0) {
		y = x*x;
		ans1 = 57568490574.0 + y*(-13362590354.0 + y*(651619640.7
			+ y*(-11214424.18 + y*(77392.33017 + y*(-184.9052456)))));
		ans2 = 57568490411.0 + y*(1029532985.0 + y*(9494680.718
			+ y*(59272.64853 + y*(267.8532712 + y*1.0))));
		ans = ans1 / ans2;
	}
	else {
		z = 8.0 / ax;
		y = z*z;
		xx = ax - 0.785398164;
		ans1 = 1.0 + y*(-0.1098628627e-2 + y*(0.2734510407e-4
			+ y*(-0.2073370639e-5 + y*0.2093887211e-6)));
		ans2 = -0.1562499995e-1 + y*(0.1430488765e-3
			+ y*(-0.6911147651e-5 + y*(0.7621095161e-6
				- y*0.934945152e-7)));
		ans = sqrt(0.636619772 / ax)*(cos(xx)*ans1 - z*sin(xx)*ans2);
	}
	return ans;
}

float bessj1(float x)
{
	float ax, z;
	double xx, y, ans, ans1, ans2;

	if ((ax = fabs(x)) < 8.0) {
		y = x*x;
		ans1 = x*(72362614232.0 + y*(-7895059235.0 + y*(242396853.1
			+ y*(-2972611.439 + y*(15704.48260 + y*(-30.16036606))))));
		ans2 = 144725228442.0 + y*(2300535178.0 + y*(18583304.74
			+ y*(99447.43394 + y*(376.9991397 + y*1.0))));
		ans = ans1 / ans2;
	}
	else {
		z = 8.0 / ax;
		y = z*z;
		xx = ax - 2.356194491;
		ans1 = 1.0 + y*(0.183105e-2 + y*(-0.3516396496e-4
			+ y*(0.2457520174e-5 + y*(-0.240337019e-6))));
		ans2 = 0.04687499995 + y*(-0.2002690873e-3
			+ y*(0.8449199096e-5 + y*(-0.88228987e-6
				+ y*0.105787412e-6)));
		ans = sqrt(0.636619772 / ax)*(cos(xx)*ans1 - z*sin(xx)*ans2);
		if (x < 0.0) ans = -ans;
	}
	return ans;
}

static void funcd(float x, float *fn, float *df)
{
	*fn = bessj0(x);
	*df = -bessj1(x);
}

void zbrak(float(*fx)(float), float x1, float x2, int n, float xb1[],
	float xb2[], int *nb)
{
	int nbb, i;
	float x, fp, fc, dx;

	nbb = 0;
	dx = (x2 - x1) / n;
	fp = (*fx)(x = x1);
	for (i = 1; i <= n; i++) {
		fc = (*fx)(x += dx);
		if (fc*fp <= 0.0) {
			xb1[++nbb] = x - dx;
			xb2[nbb] = x;
			if (*nb == nbb) return;

		}
		fp = fc;
	}
	*nb = nbb;
}

float function1(float x) {
	return (10.0*pow(M_E, -x)*sin(2 * M_PI*x) - 2);
}
float function1_d(float x) {
	return 10 * pow(M_E, -x)*((2 * M_PI)*cos(2 * M_PI*x) - sin(2 * M_PI*x));
}
float function1_cd(float x, float *fd, float *df) {
	*fd = function1(x);
	*df = function1_d(x);
}

float function2(float x) {
	return (x - pow(M_E, -x))*(x - pow(M_E, -x));
}
float function2_d(float x) {
	return 2 * pow(M_E, -2 * x)*(pow(M_E, x) + 1)*(pow(M_E, x)*x - 1);
}
float function2_cd(float x, float *fd, float *df) {
	*fd = function2(x);
	*df = function2_d(x);
}

float function3(float x) {
	return cos(x + sqrt(2.0)) + x*(x / 2 + sqrt(2));
}
float function3_d(float x) {
	return x - sin(x + sqrt(2)) + sqrt(2);
}
float function3_cd(float x, float *fn, float *df) {
	*fn = function3(x);
	*df = function3_d(x);
}

float function4(float x) {
	return x*x - 4 * x + 3;
}
float function4_d(float x) {
	return 2 * x - 4;
}
float function4_cd(float x, float *fn, float *df) {
	*fn = function4(x);
	*df = function4_d(x);
}

float function5(float x) {
	return pow(M_E, -0.005*x)*cos(sqrt(2000 - 0.01*x*x)*0.05) - 0.01;
}
float function5_d(float x) {
	return (0.0005*pow(M_E, -0.005*x)*x*sin(0.05*sqrt(2000 - 0.01*x*x))) / sqrt(2000 - 0.01*x*x)
		- 0.005*pow(M_E, -0.005*x)*cos(0.05*sqrt(2000 - 0.01*x*x));
}
float function5_cd(float x, float *fn, float *df) {
	*fn = function5(x);
	*df = function5_d(x);
}

float function32(float x) {
	float Z = 75.0;
	float R = 225.0;
	float C = 0.6*pow(10, -6);
	float L = 0.5;

	return sqrt(1.0 / (R*R) + (x*C - 1.0 / (x*L))*(x*C - 1.0 / (x*L)))
		- 1.0 / Z;
}

float function36(float x) {
	return tan(x)*35.0 - 9.81*35.0*35.0 / 
		(2.0*20.0*20.0*(0.5+0.5*cos(2.0*x))) + 1.0;
}

int main() {

	int i = 0;
	int function_type = 0;
	float interval_1 = 0.0;
	float interval_2 = 0.5*M_PI;

	float xb1[5] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
	float xb2[5] = { 0.0, 0.0, 0.0, 0.0, 0.0 };

	int nb = 2;

	float answer = 0.0;

	printf("Which method do you want to use?\n");
	printf("1.bisection\n2.linear interpolation\n3.secant\n4.newton raphson\n5.newton with bracketing\n6.muller\n");

	scanf_s("%d", &function_type);

	switch (function_type) {
	case 1:
	{
		for (i = 1; i <= nb; i++) {
			zbrak(function36, interval_1, interval_2, 10, xb1, xb2, &nb);
			answer = rtbis(function36, xb1[i], xb2[i], (float)pow(10.0, -1.0));
			printf("%f\n", answer);
		}
		printf("Number of iterations: %d\n", count);
		break;
	}
	case 2:
	{
		for (i = 1; i <= nb; i++) {
			zbrak(function36, interval_1, interval_2, 10, xb1, xb2, &nb);
			answer = rtflsp(function36, xb1[i], xb2[i], (float)pow(10.0, -1.0));
			printf("%f\n", answer);
		}
		printf("Number of iterations: %d\n", count);
		break;
	}
	case 3:
	{
		for (i = 1; i <= nb; i++) {
			zbrak(function5, interval_1, interval_2, 1, xb1, xb2, &nb);
			answer = rtsec(function5, xb1[i], xb2[i], (float)pow(10.0, -6.0));
			printf("%f\n", answer);
		}
		printf("Number of iterations: %d\n", count);
		break;
	}
	case 4:
	{
		for (i = 1; i <= nb; i++) {
			zbrak(function5, interval_1, interval_2, 1, xb1, xb2, &nb);
			answer = rtnewt(function5_cd, xb1[i], xb2[i], (float)pow(10.0, -6.0));
			printf("%f\n", answer);
		}
		printf("Number of iterations: %d\n", count);
		break;
	}
	case 5:
	{

		for (i = 1; i <= nb; i++) {
			zbrak(function5, interval_1, interval_2, 1, xb1, xb2, &nb);
			answer = rtsafe(function5_cd, xb1[i], xb2[i], (float)pow(10.0, -6.0));
			printf("%f\n", answer);
		}
		printf("Number of iterations: %d\n", count);
		break;
	}
	case 6:
	{
		for (i = 1; i <= nb; i++) {
			zbrak(function5, interval_1, interval_2, 1, xb1, xb2, &nb);
			answer = rtmuller(function5, xb1[i], xb2[i], (float)pow(10.0, -6.0));
			printf("%f\n", answer);
		}
		printf("Number of iterations: %d\n", count);
		break;
	}
	}

/*	printf("Four other functions using rtsafe.c\n");

	interval_1 = 0.1;
	interval_2 = 1.0;
	for (i = 1; i <= nb; i++) {
		zbrak(function1, interval_1, interval_2, 5, xb1, xb2, &nb);
		answer = rtsafe(function1_cd, xb1[i], xb2[i], (float)pow(10.0, -6.0));
		printf("%f\n\n\n", answer);
	}
	
	interval_1 = 0.0;
	interval_2 = 1.0;
	for (i = 1; i <= nb; i++) {
	zbrak(function2, interval_1, interval_2, 5, xb1, xb2, &nb);
	answer = rtsafe(function2_cd, xb1[i], xb2[i], (float)pow(10.0, -6.0));
	printf("%f\n", answer);
	}

	interval_1 = -2.0;
	interval_2 = -1.0;
	for (i = 1; i <= nb; i++) {
	zbrak(function3, interval_1, interval_2, 5, xb1, xb2, &nb);
	answer = rtsafe(function3_cd, xb1[i], xb2[i], (float)pow(10.0, -6.0));
	printf("%f\n", answer);
	} 

	interval_1 = 0.0;
	interval_2 = 5.0;
	nb = 2;
	for (i = 1; i <= nb; i++) {
		zbrak(function4, interval_1, interval_2, 6, xb1, xb2, &nb);
		answer = rtsafe(function4_cd, xb1[i], xb2[i], (float)pow(10.0, -6.0));
		printf("%f ", answer);
	}
	printf("\n"); */
	return 0;

}