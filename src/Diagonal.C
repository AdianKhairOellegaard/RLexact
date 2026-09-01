/* Program file Diagonal.C -
* Routines to diagonalize Hermitian matrices
* Adopted from Numerical Recipes and modified to komplex numbers
* by Per Hedegaard
* Last change: SJ 08.06.16
*
============================================
*
* RLexact: The exact diagonalization package
* Christian Rischel & Kim Lefmann, 26.02.94
* Version 4.0, September 2017
*
============================================
*/

/* #include<gcc.h> */
// #include "/usr/include/sys/types.h"
#include <sys/types.h>

#include <stdlib.h>
#include <stdio.h>
#include <complex>
#include <math.h>
#include "RLexact.h"
#include <cnr.h>
#include "Functions.h"

#define SIGN(a, b) ((b) < 0 ? -fabs(a) : fabs(a))

void eigsrt(double *, komplex **, long long);

void Diagonalize(komplex **H, long long num, double *ddin, struct FLAGS *input_flags) // Naming: ddin = energies, num = Nuniq_k
{
	time_t time_single;

	long long i, j, jj;
	double *ee, *dd;
	komplex *d, *e;

	if (input_flags->TEST_DIAGONALIZE)
	{
		LogMessageCharInt("\nDiagonalization called with matrix of dimension ", num);
		LogMessageChar("Hamiltonian=\n");
		for (i = 1; i <= num; i++)
		{
			for (j = 1; j <= num; j++)
			{
				LogMessageCharDouble("( ", real(H[i][j]));
				LogMessageCharDouble(" + i ", imag(H[i][j]));
				LogMessageChar(" ) ");
			}
			LogMessageChar("\n");
		}
		LogMessageChar("\n");
		LogMessageChar("\n");
	}

	/* Allocate */
	e = kvector(1, num);
	d = kvector(1, num);
	ee = dvector(1, num + 1); // meaning ee: 1 : num+1
	dd = dvector(0, num + 1); // copy container for ddin=energies

	if (input_flags->TEST_DIAGONALIZE)
	{
		LogMessageChar("After allocation\n");
		for (int i = 1; i <= num; i++)
		{
			LogMessageCharDouble("dd[i] =", dd[i]);
			LogMessageCharDouble(",ee[i+1] =", ee[i]);
			LogMessageChar("\n");
		}
	}

	/* Move pointers to change indexation */
	// H--;
	// for (i = 1; i <= num; i++) H[i]++;
	for (jj = 0; jj <= num; jj++)
	{
		dd[jj] = ddin[jj];
	}

	dd = dd - 1;
	if (input_flags->TEST_DIAGONALIZE)
	{
		LogMessageChar("After dd=ddin\n");
		for (int i = 1; i <= num + 1; i++)
		{
			LogMessageCharInt("i =", i);
			LogMessageCharDouble(", dd[i] =", dd[i]);
			LogMessageCharDouble(",ee[i] =", ee[i]);
			LogMessageChar("\n");
		}
	}

	/* TriDiagonalize matrix */
	if (input_flags->VERBOSE_TIME_LV1)
	{
		LogMessageChar("\n");
		time_stamp(&time_single, START, "Tridiagonalizing matrix");
	}

	htred2(H, num, d, e, input_flags);

	if (input_flags->TEST_DIAGONALIZE)
	{
		LogMessageChar("After htred2\n");
		for (int i = 0; i <= num + 1; i++)
		{
			LogMessageCharInt("i =", i);
			LogMessageCharDouble(", d[i] =", real(d[i]));
			LogMessageCharDouble(" + i", imag(d[i]));
			LogMessageCharDouble(",e [i] =", real(e[i]));
			LogMessageCharDouble(" + i", imag(e[i]));
			LogMessageChar("\n");
		}
	}

	if (input_flags->VERBOSE_TIME_LV1)
		time_stamp(&time_single, STOP, "");

	ee[1] = 0;
	for (i = 1; i <= num; i++)
	{
		ee[i + 1] = abs(e[i]);
		dd[i] = real(d[i]);

		if (input_flags->TEST_DIAGONALIZE)
		{
			LogMessageCharInt("i =", i);
			LogMessageCharDouble(", dd[i] =", dd[i]);
			LogMessageCharDouble(",ee[i] =", ee[i]);
			LogMessageChar("\n");
		}
	}

	if (input_flags->TEST_DIAGONALIZE)
	{
		LogMessageChar("\n\n Tridiag hamil= \n ");
		for (i = 1; i <= num; i++)
		{
			for (j = 1; j <= num; j++)
			{
				if (j == i)
				{
					LogMessageCharDouble("( ", dd[i]);
					LogMessageChar(" ) ");
				}
				else if (j - 1 == i)
				{
					LogMessageCharDouble("( ", ee[i + 1]);
					LogMessageChar(" ) ");
				}
				else if (j == i - 1)
				{
					LogMessageCharDouble("( ", ee[j + 1]);
					LogMessageChar(" ) ");
				}
				else
				{
					LogMessageCharDouble("( ", real(zero));
					LogMessageChar(" ) ");
				}
			}
			LogMessageChar("\n");
		}

		for (i = 1; i <= num; i++)
		{
			for (j = 1; j <= num; j++)
			{
				LogMessageCharDouble("( ", real(H[i][j]));
				LogMessageCharDouble(" + i ", imag(H[i][j]));
				LogMessageChar(" ) ");
			}
			LogMessageChar("\n");
		}
	}

	// Diagonalize:

	if (input_flags->VERBOSE_TIME_LV1)
		time_stamp(&time_single, START, "Diagonalizing tridiagonal matrix");

	htqli(dd, ee, num, H);

	if (input_flags->TEST_DIAGONALIZE)
	{
		LogMessageChar("After htqli:\n");
		for (int i = 1; i <= num + 1; i++)
		{
			LogMessageCharInt("i =", i);
			LogMessageCharDouble(", dd[i] =", dd[i]);
			LogMessageCharDouble(",ee[i] =", ee[i]);
			LogMessageChar("\n");
		}
	}

	if (input_flags->VERBOSE_TIME_LV1)
		time_stamp(&time_single, STOP, "");

	if (input_flags->TEST_DIAGONALIZE)
	{
		LogMessageChar("Diagonalization done, resulting eigenvectors\n");
		for (i = 1; i <= num; i++)
		{
			for (j = 1; j <= num; j++)
			{
				LogMessageCharDouble(" ( ", real(H[i][j]));
				LogMessageCharDouble(" + i ", imag(H[i][j]));
				LogMessageChar(" ) ");
			}
			LogMessageChar("\n");
		}
	}
	/* Sort eigenvalules */
	/*  eigsrt(dd,H,num); */

	for (jj = 0; jj <= num; jj++)
		ddin[jj] = dd[jj + 1];
	// for (i = 1; i <= num; i++) H[i]--;
	// H++;

	/* Reset pointers */
	dd = dd + 1;

	if (input_flags->TEST_DIAGONALIZE)
		LogMessageChar("\nBefore deallocation\n");

	freekvector(e, 1, num);
	freekvector(d, 1, num);
	freedvector(ee, 1, num + 1);
	freedvector(dd, 0, num + 1);

	return;
}

void htred2(komplex **a, long long num, komplex *d, komplex *e, struct FLAGS *input_flags)
{
	/* Implemented by Erik - borrowed from OpenMP implementation of the Householder reduction for complex
	matrices by Andreas Honecker and Josef Sch�le.
	Code is like old Hedegaard code inspired by numerical recipes c code Householders algorithm
	for real symmetric matrices
	Original code from Andreas Honecker and Josef Sch�ler can be found at:
	http://www.theorie.physik.uni-goettingen.de/~honecker/householder/
	*/
	int l, k, j, i;
	double scale, hh, h;
	komplex g, f;
	komplex sigma;
	komplex c1, c2;
	komplex *cptr1, *cptr2, *cptr3;
	long long dim = num;
	komplex *p;

	// e[dim - 1]	//off - diagonal
	// d diagonal elements
	p = kvector(1, num);

	if (input_flags->TEST_DIAGONALIZE)
	{
		LogMessageChar("\n\nInside Diagonalization htread2. H=\n ");
		for (i = 1; i <= num; i++)
		{
			for (j = 1; j <= num; j++)
			{
				LogMessageCharDouble("( ", real(a[i][j]));
				LogMessageCharDouble(" + i ", imag(a[i][j]));
				LogMessageChar(" ) ");
			}
			LogMessageChar("\n");
		}
	}

	for (i = dim; i >= 2; i--) /* start from the lower right corner of a */
	{
		l = i - 1;
		h = scale = 0;
		if (l > 1)
		{
			for (k = 1; k <= l; k++)
			{ /* sum row of a */
				scale += abs(a[i][k]);
			}
			if (scale == 0)			  /* skip transformation if this is zero */
				e[l] = conj(a[i][l]); /* copy complex conjugate of matrix element above diagonal */
			else
			{
				for (k = 1; k <= l; k++)
				{
					a[i][k] /= scale;					 /* rescale entries of a */
					h += abs(a[i][k] * (conj(a[i][k]))); /* and accumulate norm (squared) of "u" in h (to form sigma) */
				}
				f = a[i][l]; /* store entry m[i][i-1] in f */
				/* we have to choose sigma = a[i][i-1]/a[i][i-1]^* if |a[i]|^2 not zero */
				if (abs(f * (conj(f))))
				{
					sigma = f / (conj(f));
				}
				else
					sigma = 1; /* phase doesn't matter for |a[i][i-1]|^2 = 0 */
				sigma *= h;

				/* g = +/- sqrt(sigma), where the sign is chosen such that |a[i][i-1] - g| is as large as possible */
				g = skrt(sigma);
				c1 = f + g;
				c2 = f - g;
				if (abs(c1 * (conj(c1))) > abs(c2 * (conj(c2))))
					g = -g;

				c1 = scale * g;	 /* store scale*g below diagonal */
				e[l] = conj(c1); /* and scale*g^* above */
				/* h = |a[i]|^2 - Re(g^* f) */
				h -= (real(g) * real(f) + imag(g) * imag(f));
				/* substitute a[i][i-1] by  a[i][i-1] - g = f - g (=>  a[i] = u^t) */
				a[i][l] = f - g;
				for (j = 1; j <= l; j++)
				{
					if (input_flags->find_eigenstate)
						a[j][i] = a[i][j] / h;

					p[j] = zero; /* Initialize p=0 */
				}
				for (k = 1; k <= l; k++)
				{ /* First part of A u -> p */
					for (j = 1; j <= k - 1; j++)
					{
						/* The elements above the diagonal can be related to those below it using hermiticity */
						p[j] += (conj(a[k][j])) * (conj(a[i][k])); /* p += a[k][j]^* * a[i][k]^* DET kan godt se
														  som om at der skal tages conj(a[k][j]*a[i][k])*/
					}
				}

				for (j = 1; j <= l; j++)
				{
					p[j] /= h; /* Rescale properly */
				}
				for (j = 1; j <= l; j++)
				{		   /* Second part of A u -> p */
					g = 0; /* Form an element of A*u in g */
					for (k = 1; k <= j; k++)
					{
						g += a[j][k] * (conj(a[i][k])); /* g += a[j][k]*a[i][k]^* */
					}
					p[j] += (g / h); /* Store the result in p[j] */
				}
				f = zero;
				for (j = 1; j <= l; j++)
				{
					f += a[i][j] * p[j]; /* accumulate u * p in f */
				}
				// if ((imag(f) / abs(f)) > 1e-10)
				if (imag(f) > SMALL_NUMBER)
				{
					LogMessageChar("\n\n Error to big, a(i,j)*p(j)=\n");
					for (j = 1; j <= l; j++)
					{
						LogMessageCharInt("\n a[", i);
						LogMessageCharInt("][", j);
						LogMessageCharInt("]*p[", j);
						LogMessageCharDouble("]=(", real(a[i][j]));
						LogMessageCharDouble("+ i*", imag(a[i][j]));
						LogMessageCharDouble(")*(", real(p[j]));
						LogMessageCharDouble("+ i*", imag(p[j]));
						LogMessageCharDouble(") = ", real(a[i][j] * p[j]));
						LogMessageCharDouble("+ i*", imag(a[i][j] * p[j]));
					}
					LogMessageCharDouble("\n\n imag(f)=", imag(f));
					LogMessageCharDouble("\n\n abs(f)=", abs(f));
					LogMessageChar("\n\n");
					fatalerror("Householder: K is not real, as expected error=", 1);
				}
				hh = real(f) / (h + h); /* form K */
				for (j = 1; j <= l; j++)
				{													 /* Form q and overwrite p with it */
					f = conj(a[i][j]);								 /* store entry a[i][j]^* = u[j] in f */
					g = p[j] = p[j] - hh * f;						 /* p[j] = p[j] - hh*f */
					for (k = 1; k <= j; k++)						 /* Reduce a - we can stop at j */
						a[j][k] -= (f * (conj(p[k])) + g * a[i][k]); /* a[j][k] -= (f * p[k]^* + g * a[i][k]) */
				}
			}
		}
		else
		{
			e[l] = conj(a[i][l]); /* copy complex conjugate of matrix element above diagonal */
		}
		p[i] = h; /* Now store h in p[i] - to be able to check later if it was zero */
	}

	p[1] = zero; /* In the last transformation step, there is nothing to do */

	if (input_flags->find_eigenstate)
	{
		for (i = 1; i <= dim; i++)
		{
			l = i - 1;
			if (real(p[i]))
			{ /* skip block if we did not do a transformation */
				for (j = 1; j <= l; j++)
					p[j] = 0;
				for (k = 1; k <= l; k++)
				{
					for (j = 0; j <= l; j++)	   /* Use u^t und u^t/H stored in m to form PQ */
						p[k] += a[i][k] * a[k][j]; /* g += a[i][k]*a[k][j] */
				}
				for (k = 1; k <= l; k++)
				{ /* Now we need a second loop */
					for (j = 1; j <= l; j++)
						a[k][j] -= p[k] * (conj(a[k][i])); /* a[k][j] -= g * a[k][i]^* */
				}
			}
			d[i] = a[i][i]; /* this statement must be kept in any case */
			a[i][i] = one;	/* Reset row and column of a to identity matrix for next iteration */
			for (j = 1; j <= l; j++)
				a[i][j] = a[j][i] = zero;
		}
	}
	else
	{
		for (i = 1; i <= dim; i++)
			d[i] = a[i][i]; /* this statement must be kept in any case */
	}
}

// htqli used for both MATRIX and LANCZOS!
long long htqli(double *d, double *e, long long num, komplex **z)
{ /* CRECIPES -- HEDEGARD : Tridiagonal QL with Implicit shifts*/
	long long m, l, iter, i, k;
	double s, r, p, g, f, dd, c, b;
	komplex fc;
	long long itercollect = 0;
	char errortext[100];
	iter = 0;

	for (i = 2; i <= num; i++)
		e[i - 1] = e[i];
	e[num] = 0.0;
	for (l = 1; l <= num; l++)
	{
		itercollect += iter;
		iter = 0;
		do
		{
			for (m = l; m <= num - 1; m++)
			{
				dd = fabs(d[m]) + fabs(d[m + 1]);
				if (fabs(e[m]) + dd == dd)
				{
					break;
				}
			}
			if (m != l)
			{
				if (iter++ == 500)
				{
					LogMessageCharDouble("It ends here", iter);
					fatalerror("Too many interations in TQLI", num);
				}
				g = (d[l + 1] - d[l]) / (2.0 * e[l]);
				r = sqrt((g * g) + 1.0);
				g = d[m] - d[l] + e[l] / (g + SIGN(r, g));
				s = c = 1.0;
				p = 0.0;
				for (i = m - 1; i >= l; i--)
				{
					f = s * e[i];
					b = c * e[i];
					if (fabs(f) >= fabs(g))
					{
						c = g / f;
						r = sqrt((c * c) + 1.0);
						e[i + 1] = f * r;
						c *= (s = 1.0 / r);
					}
					else
					{
						s = f / g;
						r = sqrt((s * s) + 1.0);
						e[i + 1] = g * r;
						s *= (c = 1.0 / r);
					}
					g = d[i + 1] - p;
					r = (d[i] - g) * s + 2.0 * c * b;
					p = s * r;
					d[i + 1] = g + p;
					g = c * r - b;
					/* Next loop could be omitted if eigenvectors not wanted,
					ecxept we need it for the lanczos stopping criteria */
					for (k = 1; k <= num; k++)
					{
						fc = z[k][i + 1];
						z[k][i + 1] = s * z[k][i] + c * fc;
						z[k][i] = c * z[k][i] - s * fc;
					}
				}

				d[l] = d[l] - p;
				e[l] = g;
				e[m] = 0.0;
			}
		} while (m != l);
	}
	return itercollect;
}

void eigsrt(double *d, komplex **v, long long num)
{ /* CRECIPES -- HEDEGARD */
	long long k, j, i;
	double p;
	komplex pc;

	for (i = 1; i < num; i++)
	{
		p = d[k = i];
		for (j = i + 1; j <= num; j++)
			if (d[j] <= p)
				p = d[k = j];
		if (k != i)
		{
			d[k] = d[i];
			d[i] = p;
			for (j = 1; j <= num; j++)
			{
				pc = v[j][i];
				v[j][i] = v[j][k];
				v[j][k] = pc;
			}
		}
	}
}

/* Example of call:
void Eigenvectors(double qx,double qy,komplex **cpoint){
int i,j;
komplex cp;
htred2(Hk,nn,dd,ee);

for(i=1;i<=nn;i++){
ff[i] = real(dd[i]);
gg[i] = real(ee[i]);
};

htqli(ff,gg,nn,Hk);
eigsrt(ff,Hk,nn);


}
*/
