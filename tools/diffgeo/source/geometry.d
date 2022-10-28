// Copyright 2021 KU Leuven.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
// Author: Rinaldo Wander Montalvão, PhD
//
module geometry;

import std.math, std.conv;

// Numerical Analysis 9th ed - Burden, Faires (Ch. 3 Natural Cubic Spline, Pg. 149) 
void spline(in double[] x, in double[] a, ref double[][] a2)
{
	const int n = to!int(x.length) - 1;

	double[] h = new double[](n);
	double[] A = new double[](n);

	double[] l = new double[](n + 1);
	double[] u = new double[](n + 1);
	double[] z = new double[](n + 1);

	double[] b = new double[](n);
	double[] d = new double[](n);
	double[] c = new double[](n + 1);

	a2 = new double[][](n, 3);

	// Step 1 
	foreach (i; 0 .. n)
		h[i] = x[i + 1] - x[i];

	// Step 2 
	foreach (i; 1 .. n)
		A[i] = 3.0 * (a[i + 1] - a[i]) / h[i] - 3.0 * (a[i] - a[i - 1]) / h[i - 1];

	// Step 3 
	l[0] = 1.0;
	u[0] = 0.0;
	z[0] = 0.0;

	// Step 4 
	foreach (i; 1 .. n)
	{
		l[i] = 2.0 * (x[i + 1] - x[i - 1]) - h[i - 1] * u[i - 1];
		u[i] = h[i] / l[i];
		z[i] = (A[i] - h[i - 1] * z[i - 1]) / l[i];
	}

	// Step 5 
	l[n] = 1.0;
	z[n] = 0.0;
	c[n] = 0.0;

	// Step 6 
	foreach_reverse (j; 0 .. n)
	{
		c[j] = z[j] - u[j] * c[j + 1];
		b[j] = (a[j + 1] - a[j]) / h[j] - h[j] * (c[j + 1] + 2.0 * c[j]) / 3.0;
		d[j] = (c[j + 1] - c[j]) / (3.0 * h[j]);
	}

	// Step 7 
	foreach (i; 0 .. n)
	{
		a2[i][0] = b[i];
		a2[i][1] = c[i];
		a2[i][2] = d[i];
	}
}

double splint(in double[] x, in double[] a, in double[][] a2, in double xv)
{
	const int n = to!int(x.length) - 1;

	// TODO: deal with boundary errors
	int i = n - 1;
	foreach (j; 0 .. n - 1)
	{
		if ((x[j] <= xv) && (xv <= x[j + 1]))
		{
			i = j;
			break;
		}
	}

	double y = a[i]
		+ a2[i][0] * (xv - x[i])
		+ a2[i][1] * pow(xv - x[i], 2)
		+ a2[i][2] * pow(xv - x[i], 3);

	return y;
}

// Chebyshev 
void chebfit(in double x0, in double x1, in double[] x, in double[] a, in double[][] a2, ref double[] c)
{
	const int N = to!int(c.length);

	double slope = (x1 - x0) / 2.0;
	double intercept = (x1 + x0) / 2.0;

	double[] f = new double[](N);
	foreach (k; 0 .. N)
	{
		double xk = cos(PI * (k + 0.5) / N) * slope + intercept;

		f[k] = splint(x, a, a2, xk);
	}

	foreach (j; 0 .. N)
	{
		double sum = 0.0;
		foreach (k; 0 .. N)
		{
			sum += f[k] * cos(PI * j * (k + 0.5) / N);
		}

		c[j] = 2.0 * sum / N;
	}

}

double chebeval(in double x0, in double x1, in double[] c, in double x)
{
	const int N = to!int(c.length);

	double y = (2.0 * x - x0 - x1) / (x1 - x0);

	double y2 = 2.0 * y;

	double d = 0.0, dd = 0.0;
	foreach_reverse (j; 1 .. N)
	{
		double temp = d;
		d = y2 * d - dd + c[j];
		dd = temp;
	}

	return y * d - dd + 0.5 * c[0];
}


void chebder(in double x0, in double x1, in double[] c, ref double[] cder)
{
	const int N = to!int(c.length);

	cder[N - 1] = 0.0;
	cder[N - 2] = 2.0 * (N - 1) * c[N - 1];

	foreach_reverse (j; 1 .. N - 1)
		cder[j - 1] = cder[j + 1] + 2.0 * j * c[j];

	cder[] *= 2.0 / (x1 - x0);
}


void nderiv(in double[] xx, in double[] yy, in double[][] y2, in double t, ref double[3] derivatives)
{

	double[] c0 = new double[](50);
	double[] d1 = new double[](50);
	double[] d2 = new double[](50);
	double[] d3 = new double[](50);

	chebfit(t - 1.0, t + 1.0, xx, yy, y2, c0);

	chebder(t - 1.0, t + 1.0, c0, d1);
	chebder(t - 1.0, t + 1.0, d1, d2);
	chebder(t - 1.0, t + 1.0, d2, d3);

	derivatives[0] = chebeval(t - 1.0, t + 1.0, d1, t);
	derivatives[1] = chebeval(t - 1.0, t + 1.0, d2, t);
	derivatives[2] = chebeval(t - 1.0, t + 1.0, d3, t);

}

double[3] cross(in double[3] a, in double[3] b)
{

	double[3] result;

	result[0] = a[1] * b[2] - a[2] * b[1];
	result[1] = a[2] * b[0] - a[0] * b[2];
	result[2] = a[0] * b[1] - a[1] * b[0];

	return result;
}

double dot(in double[3] a, in double[3] b)
{

	double result = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];

	return result;
}

double norm(in double[3] a)
{

	double result = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);

	return result;
}

double[3] div(in double[3] a, in double b)
{

	double[3] result;

	result[0] = a[0] / b;
	result[1] = a[1] / b;
	result[2] = a[2] / b;

	return result;
}

/*!
 * Function used to calculate a dihedral angle.
 * @param list a matrix that contains the list of vectors.
 */
double dihed(in double[4][3] list)
{
	double[3] co1;
	double[3] co2;
	double[3] co3;
	double[3] co4;

	/// Need to use a set of vectors to calculate the dihedral.

	double[3] v_12; // Vector from 1 to 2 
	double[3] v_23; // Vector from 2 to 3 
	double[3] v_34; // Vector from 3 to 4 

	double[3] n_1; // Normal vector of the plane 1-2-3.
	double[3] n_2; // normal vector of the plane 2-3-4.

	double abs_n1; // Length of the vectors.
	double abs_n2;
	double abs_n1_n2;

	double n1_n2; // Scalar product (dot product) of n_1 and n_2.
	double ss; // Scalar product of v_12 and n_2 to determine the sign of phi.

	double phi; // The dihedral angle

	foreach (j; 0 .. 3)
	{
		co1[j] = list[j][0];
		co2[j] = list[j][1];
		co3[j] = list[j][2];
		co4[j] = list[j][3];
	}

	foreach (j; 0 .. 3)
	{
		v_12[j] = co2[j] - co1[j];
		v_23[j] = co3[j] - co2[j];
		v_34[j] = co4[j] - co3[j];
	}

	n_1[0] = v_12[1] * v_23[2] - v_12[2] * v_23[1];
	n_1[1] = v_12[2] * v_23[0] - v_12[0] * v_23[2];
	n_1[2] = v_12[0] * v_23[1] - v_12[1] * v_23[0];

	abs_n1 = n_1[0] * n_1[0] + n_1[1] * n_1[1] + n_1[2] * n_1[2];

	n_2[0] = v_23[1] * v_34[2] - v_23[2] * v_34[1];
	n_2[1] = v_23[2] * v_34[0] - v_23[0] * v_34[2];
	n_2[2] = v_23[0] * v_34[1] - v_23[1] * v_34[0];

	abs_n2 = n_2[0] * n_2[0] + n_2[1] * n_2[1] + n_2[2] * n_2[2];

	abs_n1_n2 = abs_n1 * abs_n2;

	if (abs_n1_n2 < 0.0)
	{
		//writeln("warning: dihedral angle not defined.");
		phi = 0.0;
		return (phi);
	}

	n1_n2 = n_1[0] * n_2[0] + n_1[1] * n_2[1] + n_1[2] * n_2[2];

	ss = v_12[0] * n_2[0] + v_12[1] * n_2[1] + v_12[2] * n_2[2];

	if (ss >= 0.0)
	{
		phi = acos(n1_n2 / sqrt(abs_n1_n2));
	}
	else
	{
		phi = -acos(n1_n2 / sqrt(abs_n1_n2));
	}

	/*printf("the dihedral angle is (remember it is in radians) %f\n",phi);
    phi = (180.0/3.14159) * phi;
    */

	/*printf("the dihedral angle is %f\n",phi);*/

	phi = (180.0 / PI) * phi;
	return (phi);
}
