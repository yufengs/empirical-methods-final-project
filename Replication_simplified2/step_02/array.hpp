//
// Utility classes to work with matlab arrays
//

#pragma once

#include <math.h>
#include <vector>
#include "blas.h"

#include "arraypool.hpp"

static void doAssert(const char *str)
{
	mexPrintf("Assertion Failure:");
	mexPrintf(str);
    mexErrMsgIdAndTxt("Likelihood::assert", str);
}

#define mexAssert(cond) ( (cond) || (doAssert(#cond),true) )


struct Range
{
    int lo, hi;

    Range(int i)
    {
        lo = i;
        hi = i;
    }

    Range(int _lo, int _hi)
    {
        lo = _lo;
        hi = _hi;
    }

    int size() const
    {
        return (hi-lo)+1;
    }
};

struct Array1D
{
    mwSize      ncols;
    double *    ptr;

    Array1D(const mxArray * A)
    {
        ncols = mxGetN(A);
        ptr   = mxGetPr(A);
    }

    double & operator()(int col)
    {
        mexAssert(col > 0);
        --col;
        return ptr[col];
    }

    double & operator[](int col)
    {
        mexAssert(col > 0);
        --col;
        return ptr[col];
    }
};


// Utility class to work with 2D arrays
struct Array2D
{
    mwSize      nrows;
    mwSize      ncols;
    double *    ptr;

    Array2D(int _nrows, int _ncols, double * _ptr)
    {
        nrows = _nrows;
        ncols = _ncols;
        ptr = _ptr;
    }

    Array2D(const mxArray * A)
    {
        nrows = mxGetM(A);
        ncols = mxGetN(A);
        ptr   = mxGetPr(A);
    }

    double & operator()(int col)
    {
        mexAssert(col > 0);

        return ptr[col-1];
    }

    double & operator()(int row, int col)
    {
        mexAssert(row > 0);
        mexAssert(col > 0);

        return ptr[(col-1)*nrows + (row-1)];
    }
};



// Utility class to work with arrays upto ND dimensions
// ND is not really constant, you'll need to update the code after
// increasing it.
#define ND (4)

struct IndexList;

struct ArrayND
{
    double *        ptr;
    int             ndims;
    mwSize          dims[ND];
    mwSize          offsets[ND-1];
	mwSize			numelems;

    static ArrayND create2D(int nrows, int ncols, mxArray * & out)
    {
        out = mxCreateDoubleMatrix(nrows, ncols, mxREAL);

        return ArrayND(out);
    }

    ArrayND(const mxArray * A)
    {
        ndims = mxGetNumberOfDimensions(A);
        ptr  = mxGetPr(A);

        const mwSize *mxdims = mxGetDimensions(A);

        dims[0] = ndims < 1 ? 1 : mxdims[0];
        dims[1] = ndims < 2 ? 1 : mxdims[1];
        dims[2] = ndims < 3 ? 1 : mxdims[2];
        dims[3] = ndims < 4 ? 1 : mxdims[3];

        offsets[0] = ndims < 2 ? 1 : dims[0];
        offsets[1] = ndims < 3 ? 1 : offsets[0]*dims[1];
        offsets[2] = ndims < 4 ? 1 : offsets[1]*dims[2];

		numelems = dims[0]*dims[1]*dims[2]*dims[3];
	}

    ArrayND(int _ndims, const mwSize *_dims)
    {
        ndims = _ndims;

		numelems = 1;
        for(int i = 0; i < ndims; i++)
            numelems *= _dims[i];

        ptr = (double *)arrayPool->alloc(numelems*sizeof(double));

        dims[0] = _dims[0];
        dims[1] = ndims > 1 ? _dims[1] : 1;
        dims[2] = ndims > 2 ? _dims[2] : 1;
        dims[3] = ndims > 3 ? _dims[3] : 1;

        offsets[0] = dims[0];
        offsets[1] = dims[0]*dims[1];
        offsets[2] = dims[0]*dims[1]*dims[2];
	}

    ArrayND(int m, int n)
    {
        ndims = 2;
        ptr = (double *)arrayPool->alloc(m*n*sizeof(double));

        dims[0] = m;
        dims[1] = n;
        dims[2] = 1;
        dims[3] = 1;

        offsets[0] = dims[0];
        offsets[1] = 1;
        offsets[2] = 1;

		numelems = dims[0]*dims[1]*dims[2]*dims[3];
	}

    ArrayND(int nrows, int ncols, double * _ptr)
    {
        ndims = 2;
        ptr = _ptr;

        dims[0] = nrows;
        dims[1] = ncols;
        dims[2] = 1;
        dims[3] = 1;

        offsets[0] = dims[0];
        offsets[1] = 1;
        offsets[2] = 1;

		numelems = dims[0]*dims[1]*dims[2]*dims[3];
	}

    // return number of elements in the array (sum all dimensions)
    int size() const
    {
        return numelems;
    }

    double & operator()(int x) const
    {
        mexAssert(x > 0);
        mexAssert(x <= numelems);

        return ptr[x-1];
    }

    double & operator()(int x, int y) const
    {
		int elem = (x-1) + (y-1)*offsets[0];

		mexAssert(x > 0);
        mexAssert(y > 0);
		mexAssert(elem < numelems);

        return ptr[elem];
    }

    double & operator()(int x, int y, int z) const
    {
		int elem = (x-1) + (y-1)*offsets[0] + (z-1)*offsets[1];

		mexAssert(x > 0);
        mexAssert(y > 0);
		mexAssert(elem < numelems);

        return ptr[elem];
    }

    double & operator()(int x, int y, int z, int w) const
    {
		int elem = (x-1) + (y-1)*offsets[0] + (z-1)*offsets[1] + (w-1)*offsets[2];

		mexAssert(x > 0);
        mexAssert(y > 0);
		mexAssert(elem < numelems);

        return ptr[elem];
    }
    
    ArrayND slice(Range a, Range b)
    {
        mwSize newdims[2];
        
        newdims[0] = a.size();
        newdims[1] = b.size();

        int ndims = 2;
        if(newdims[1] == 1)
            ndims -= 1;
        
        ArrayND out(ndims, newdims);
        
        int outi = 1;
        for(int i = a.lo; i <= a.hi; i++)
        {
            int outj = 1;
            for(int j = b.lo; j <= b.hi; j++)
            {
                out(outi,outj) = (*this)(i, j);
                outj++;
            }
            outi++;
        }

        return out;
    }

    ArrayND slice(Range a, Range b, Range c)
    {
        mwSize newdims[3];

        newdims[0] = a.size();
        newdims[1] = b.size();
        newdims[2] = c.size();

        int ndims = 3;
        if(newdims[2] == 1)
        {
            ndims -= 1;
            if(newdims[1] == 1)
                ndims -= 1;
        }

        ArrayND out(ndims, newdims);

        int outi = 1;
        for(int i = a.lo; i <= a.hi; i++)
        {
            int outj = 1;
            for(int j = b.lo; j <= b.hi; j++)
            {
                int outk = 1;
                for(int k = c.lo; k <= c.hi; k++)
                {
                    out(outi,outj,outk) = (*this)(i, j, k);
                    outk++;
                }
                outj++;
            }
            outi++;
        }


        return out;
    }

    void dump() const
    {
        printf("ndims:%d dims:[%d %d %d %d]\n", ndims, dims[0], dims[1], dims[2], dims[3]);
        for(int i = 0; i <= size(); i++)
            printf("   el:%d = %f\n", i, ptr[i]);
    }

    IndexList operator()(const ArrayND & a);
};

struct IndexList
{
    ArrayND &tgt;
    const ArrayND &idx;

    IndexList(ArrayND &_tgt, const ArrayND &_idx) : tgt(_tgt), idx(_idx)
    {
    }

    ArrayND operator = (double val)
    {
        ArrayND out(tgt);

        for(int i=1; i <= out.size(); i++)
        {
            if(idx(i))
                out(i) = val;
            else
                out(i) = tgt(i);
        }

        return out;
    }
};

IndexList ArrayND::operator()(const ArrayND & a)
{
    return IndexList(*this, a);
}


static double sum(const ArrayND &a)
{
    mexAssert(a.ndims == 2);

    double sum = 0.0;

    for(int col = 1; col <= a.dims[1]; col++)
    {
        for(int row = 1; row <= a.dims[0]; row++)
            sum += a(row,col);
    }

    return sum;
}


#if WITH_BLAS

static bool sameshape(const ArrayND & a, const ArrayND &b)
{
    return  a.numelems == b.numelems && a.dims[0] == b.dims[0] &&
    a.dims[1] == b.dims[1] &&
    a.dims[2] == b.dims[2] &&
    a.dims[3] == b.dims[3];
}

ArrayND sum(const ArrayND &a, int dim)
{
    mexAssert(a.ndims <= 2);

    if(dim == 2)
    {
        // sum cols
        ArrayND out(a.dims[0], 1);

        for(int i = 1; i <= a.dims[0]; i++) // rows
        {
            double sum = 0.0;

            for(int j = 1; j <= a.dims[1]; j++)     // cols
                sum += a(i,j);

            out(i,1) = sum;
        }

        return out;
    }

    // TODO: sum rows
    mexAssert(false);

    return ArrayND(1,1);
}

ArrayND operator - (const ArrayND &b)
{
    ArrayND out(b.ndims, b.dims);
    
    for(int i = 1; i <= b.size(); i++)
        out(i) = -b(i);
    
    return out;
}

ArrayND operator * (double a, const ArrayND &b)
{
    ArrayND out(b.ndims, b.dims);

    for(int i = 1; i <= b.size(); i++)
        out(i) = a * b(i);

    return out;
}

ArrayND operator + (double a, const ArrayND &b)
{
    ArrayND out(b.ndims, b.dims);

    for(int i = 1; i <= b.size(); i++)
        out(i) = a + b(i);

    return out;
}

ArrayND operator - (double a, const ArrayND &b)
{
    ArrayND out(b.ndims, b.dims);

    for(int i = 1; i <= b.size(); i++)
        out(i) = a - b(i);

    return out;
}

ArrayND operator + (const ArrayND & a, const ArrayND &b)
{
    mexAssert(sameshape(a,b));

    ArrayND out(b.ndims, b.dims);

    for(int i = 1; i <= b.size(); i++)
        out(i) = a(i) + b(i);

    return out;
}

ArrayND operator - (const ArrayND & a, const ArrayND &b)
{
    mexAssert(sameshape(a,b));

    ArrayND out(b.ndims, b.dims);

    for(int i = 1; i <= b.size(); i++)
        out(i) = a(i) - b(i);

    return out;
}


ArrayND smul(const ArrayND & a, const ArrayND &b)
{
    mexAssert(sameshape(a,b));
    
    ArrayND out(b.ndims, b.dims);
    
    for(int i = 1; i <= b.size(); i++)
        out(i) = a(i) * b(i);
    
    return out;
}

// matrix multiply
ArrayND operator * (const ArrayND & a, const ArrayND &b)
{
    /* form of op(A) & op(B) to use in matrix multiplication */
    const char *chn = "N";
    /* scalar values to use in dgemm */
    double one = 1.0, zero = 0.0;

    double * A = a.ptr;
    double * B = b.ptr;
    ptrdiff_t m = a.dims[0];
    ptrdiff_t p = a.dims[1];
    ptrdiff_t n = b.dims[1];

    if (p != b.dims[0]) {
        mexErrMsgIdAndTxt("MATLAB:matrixMultiply:matchdims",
                "Inner dimensions of matrix multiply do not match.");
        a.dump();
    }

    /* create output matrix C */
    ArrayND out(m, n);

    double *C = out.ptr;

    dgemm((char *)chn, (char*)chn, &m, &n, &p, &one, A, &m, B, &p, &zero, C, &m);

    return out;
}

// per-element divide
ArrayND operator / (const ArrayND & a, const ArrayND &b)
{
    mexAssert(sameshape(a,b));

    ArrayND out(a.ndims, a.dims);

    for(int i = 1; i <= a.size(); i++)
        out(i) = a(i) / b(i);

    return out;
}

ArrayND operator / (const ArrayND & a, const double b)
{
    ArrayND out(a.ndims, a.dims);

    for(int i = 1; i <= a.size(); i++)
        out(i) = a(i) / b;

    return out;
}


ArrayND operator < (const ArrayND & a, const double b)
{
    ArrayND out(a.ndims, a.dims);
    
    for(int i = 1; i <= a.size(); i++)
        out(i) = a(i) < b;
    
    return out;
}


ArrayND operator <= (const ArrayND & a, const double b)
{
    ArrayND out(a.ndims, a.dims);
    
    for(int i = 1; i <= a.size(); i++)
        out(i) = a(i) <= b;
    
    return out;
}


ArrayND operator == (const ArrayND & a, const double b)
{
    ArrayND out(a.ndims, a.dims);

    for(int i = 1; i <= a.size(); i++)
        out(i) = a(i) == b;

    return out;
}

ArrayND make1x1(double a)
{
    ArrayND mat(1, 1);

    mat(1) = a;

    return mat;
}

ArrayND make2x1(double a, double b)
{
    ArrayND mat(2, 1);

    mat(1) = a;
    mat(2) = b;

    return mat;
}

ArrayND zeros(int m, int n)
{
    ArrayND mat(m, n);

    for(int i = 1; i <= m*n; i++)
        mat(i) = 0.0;

    return mat;
}

ArrayND ones(int m, int n)
{
    ArrayND mat(m, n);

    for(int i = 1; i <= m*n; i++)
        mat(i) = 1.0;

    return mat;
}

ArrayND eye(int n)
{
    ArrayND mat(n, n);

    for(int j = 1; j <= n; j++)
    {
        for(int i = 1; i <= n; i++)
            mat(j,i) = (i == j ? 1.0 : 0.0);
    }

    return mat;
}

ArrayND exp(const ArrayND &a)
{
    ArrayND out(a.ndims, a.dims);

    for(int i = 1; i <= a.size(); i++)
        out(i) = ::exp( a(i) );

    return out;
}


ArrayND log(const ArrayND &a)
{
    ArrayND out(a.ndims, a.dims);

    for(int i = 1; i <= a.size(); i++)
        out(i) = ::log( a(i) );

    return out;
}

ArrayND transpose(const ArrayND &a)
{
    mexAssert(a.ndims <= 2);

    ArrayND out(a.dims[1], a.dims[0]);

    for(int j = 1; j <= a.dims[1]; j++)
    {
        for(int i = 1; i <= a.dims[0]; i++)
        {
            out(j,i) = a(i, j);
        }
    }

    return out;
}

ArrayND repmat(const ArrayND & a, int m, int n)
{
    int nrows = a.dims[0];
    int ncols = a.dims[1];

    ArrayND out(nrows*m, ncols*n);

    for(int i = 1; i <= nrows*m; i++)
    {
        for(int j = 1; j <= ncols*n; j++)
        {
            double val = a((i-1)%nrows+1,(j-1)%ncols+1);

            out(i, j) = val;
        }
    }

    return out;
}



ArrayND horzcat(std::vector<const ArrayND *> & a)
{
    mexAssert(a.size() > 0);

    int nrows = a[0]->dims[0];
    int ncols = 0;

    for(int i = 0; i < a.size(); i++)
    {
        ncols += a[i]->dims[1];
        mexAssert(a[i]->dims[0] == nrows); // all inputs have same number of rows
    }

    ArrayND out(nrows, ncols);

    for(int j = 1; j <= nrows; j++)
    {
        int col = 1;
        for(int i = 0; i < a.size(); i++)
        {
            const ArrayND & src = *a[i];

            for(int k = 1; k <= src.dims[1]; k++)
            {
                out(j, col) = src(j, k);
                col++;
            }
        }
    }

    return out;
}

ArrayND vertcat(std::vector<const ArrayND *> & a)
{
    mexAssert(a.size() > 0);

    int nrows = 0;
    int ncols = a[0]->dims[1];

    for(int i = 0; i < a.size(); i++)
    {
        nrows += a[i]->dims[0];
        mexAssert(a[i]->dims[1] == ncols); // all inputs have same number of cols
    }

    ArrayND out(nrows, ncols);

    int row = 1;
    for(int i = 0; i < a.size(); i++)
    {
        const ArrayND & src = *a[i];

        for(int k = 1; k <= src.dims[0]; k++)       //rows
        {
            for(int j = 1; j <= src.dims[1]; j++)    //cols
                out(row, j) = src(k, j);

            row++;
        }
    }

    return out;
}


ArrayND horzcat(const ArrayND &a, const ArrayND &b)
{
    std::vector<const ArrayND *> ms;

    ms.push_back(&a);
    ms.push_back(&b);

    return horzcat(ms);
}

ArrayND horzcat(const ArrayND &a, const ArrayND &b, const ArrayND &c, const ArrayND &d)
{
    std::vector<const ArrayND *> ms;

    ms.push_back(&a);
    ms.push_back(&b);
    ms.push_back(&c);
    ms.push_back(&d);

    return horzcat(ms);
}

ArrayND horzcat(const ArrayND &a, const ArrayND &b, const ArrayND &c, const ArrayND &d, const ArrayND &e)
{
    std::vector<const ArrayND *> ms;

    ms.push_back(&a);
    ms.push_back(&b);
    ms.push_back(&c);
    ms.push_back(&d);
    ms.push_back(&e);

    return horzcat(ms);
}

ArrayND vertcat(const ArrayND &a, double b)
{
    std::vector<const ArrayND *> ms;
    ArrayND _b = make1x1(b);

    ms.push_back(&a);
    ms.push_back(&_b);  

    return vertcat(ms);
}

ArrayND vertcat(const double &a, const double &b)
{
    std::vector<const ArrayND *> ms;
    ArrayND _a = make1x1(a);
    ArrayND _b = make1x1(b);

    ms.push_back(&_a);
    ms.push_back(&_b);

    return vertcat(ms);
}

ArrayND vertcat(const ArrayND &a, const ArrayND &b)
{
    std::vector<const ArrayND *> ms;

    ms.push_back(&a);
    ms.push_back(&b);

    return vertcat(ms);
}
#endif
