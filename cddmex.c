/* combined.c
 *
 * A single mexFunction dispatches based on the first input argument (a command string):
 *    "solve_lp"     : Solve an LP (primal/dual) using cddlib.
 *    "solve_lp_ds"  : Solve an LP using the Dual Simplex method.
 *    "adj_extreme"  : Compute extreme generators and adjacency from an H‐representation.
 *    "extreme"      : Compute extreme generators (V‐representation) from an H‐representation.
 *    "hull"         : Compute the H‐representation (inequalities) from a V‐representation.
 *    "version"      : Print version information.
 *
 * Compile (example):
 *   mex -I/usr/local/include/cddlib -I/usr/local/include -L/usr/local/lib -DGMPRATIONAL combined.c -lcddgmp -lgmp
 *
 * Version: 2.0
 */

#include "mex.h"
#include "setoper.h"   /* Must come before cdd.h so that set_type is defined */
#include <gmp.h>
#include "cdd.h"
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <inttypes.h>

#define CDDMEX_VERSION "2.0"
#define BUFFER_SIZE 4096

#ifndef dd_get_si2
static inline void dd_get_si2(mpq_srcptr element, int64_t *num, int64_t *den) {
    *num = mpz_get_si(mpq_numref(element));
    *den = mpz_get_si(mpq_denref(element));
}
#endif



// Debug: Log LP fields.
void LogLpFields(const dd_LPPtr lp) {
    if (!lp) {
        mexPrintf("lp is NULL.\n");
        return;
    }
    mexPrintf("=== Logging cddlib LP Fields ===\n");
    mexPrintf("lp->d = %d  (dimension including homogeneous)\n", lp->d);
    mexPrintf("lp->m = %d  (number of constraints + 1 for homogeneous)\n", lp->m);
    mexPrintf("lp->objective (objective type): %d\n", lp->objective);
    mexPrintf("lp->solver (solver type): %d\n", lp->solver);
    mexPrintf("lp->Homogeneous (is homogeneous?): %d\n", lp->Homogeneous);
    if (lp->sol) {
        mexPrintf("lp->sol (primal solution), length = %d:\n", lp->d);
        for (int i = 0; i < lp->d; i++) {
            char *valStr = mpq_get_str(NULL, 10, lp->sol[i]);
            mexPrintf("  sol[%d] = %s\n", i, valStr);
            free(valStr);
        }
    } else {
        mexPrintf("lp->sol is NULL (no primal solution).\n");
    }
    if (lp->dsol) {
        mexPrintf("lp->dsol (dual solution), length = %d:\n", lp->d);
        for (int i = 0; i < lp->d; i++) {
            char *valStr = mpq_get_str(NULL, 10, lp->dsol[i]);
            mexPrintf("  dsol[%d] = %s\n", i, valStr);
            free(valStr);
        }
    } else {
        mexPrintf("lp->dsol is NULL (no dual solution).\n");
    }
    if (lp->nbindex) {
        int nbSize = lp->d + 2;
        mexPrintf("lp->nbindex (size ~ %d):\n", nbSize);
        for (int i = 0; i < nbSize; i++) {
            mexPrintf("  nbindex[%d] = %d\n", i, lp->nbindex[i]);
        }
    } else {
        mexPrintf("lp->nbindex is NULL.\n");
    }
    {
        char *valStr = mpq_get_str(NULL, 10, lp->optvalue);
        mexPrintf("lp->optvalue = %s  (objective)\n", valStr);
        free(valStr);
    }
    if (lp->A) {
        mexPrintf("lp->A (constraint matrix):\n");
        for (dd_rowrange i = 0; i < lp->m; i++) {
            mexPrintf("  Row %d: ", i);
            for (dd_colrange j = 0; j < lp->d; j++) {
                char *valStr = mpq_get_str(NULL, 10, lp->A[i][j]);
                mexPrintf("[%d][%d] = %s  ", i, j, valStr);
                free(valStr);
            }
            mexPrintf("\n");
        }
    } else {
        mexPrintf("lp->A is NULL.\n");
    }
    if (lp->B) {
        mexPrintf("lp->B (RHS vector):\n");
        for (dd_rowrange i = 0; i < lp->m; i++) {
            char *valStr = mpq_get_str(NULL, 10, (mpq_srcptr)lp->B[i]);
            mexPrintf("  B[%d] = %s\n", i, valStr);
            free(valStr);
        }
    } else {
        mexPrintf("lp->B is NULL.\n");
    }
    mexPrintf("lp->objrow (row index for objective): %d\n", lp->objrow);
    mexPrintf("lp->rhscol (column index for RHS): %d\n", lp->rhscol);
    mexPrintf("=== End of LP Fields ===\n\n");
}


/*
 * LogSetFamily
 * ------------
 * Logs the contents of a dd_SetFamilyPtr structure.
 *
 * Input:
 *   A - a pointer to a dd_setfamily structure.
 *
 * It prints:
 *   - The total number of sets (famsize)
 *   - The overall set size (setsize)
 *   - For each set, it first counts the actual number of members and then
 *     prints each member index.
 */
void LogSetFamily(const dd_SetFamilyPtr A) {
    if (!A) {
        mexPrintf("dd_SetFamilyPtr is NULL.\n");
        return;
    }
    mexPrintf("Logging dd_SetFamilyPtr structure:\n");
    mexPrintf("  famsize: %lld\n", (long long)A->famsize);
    mexPrintf("  setsize: %lld\n", (long long)A->setsize);
    
    for (dd_bigrange i = 0; i < A->famsize; i++) {
        dd_bigrange card = set_card(A->set[i]);
        int count = 0;
        // Count the number of members in the set.
        for (dd_bigrange j = 1; j <= card; j++) {
            if (set_member(j, A->set[i])) {
                count++;
            }
        }
        mexPrintf("  Set %lld (actual count %d): ", (long long)i, count);
        // Print each member.
        for (dd_bigrange j = 1; j <= card; j++) {
            if (set_member(j, A->set[i])) {
                mexPrintf("%lld ", (long long)j);
            }
        }
        mexPrintf("\n");
    }
    mexPrintf("Finished logging dd_SetFamilyPtr structure.\n");
}


static void LogVStruct(const mxArray *vStruct) {
    if (!vStruct) {
        mexPrintf("vStruct is NULL.\n");
        return;
    }
    
    mexPrintf("Logging vStruct:\n");
    
    // List of expected field names.
    const char *fields[] = {"Vnum", "Vden", "Rnum", "Rden", "vpos", "rpos", "lin"};
    int nfields = sizeof(fields) / sizeof(fields[0]);
    
    for (int i = 0; i < nfields; i++) {
        const char *fieldName = fields[i];
        mxArray *field = mxGetField(vStruct, 0, fieldName);
        if (!field) {
            mexPrintf("Field '%s' is missing.\n", fieldName);
            continue;
        }
        mwSize m = mxGetM(field);
        mwSize n = mxGetN(field);
        mexPrintf("Field '%s': size = [%d x %d]\n", fieldName, (int)m, (int)n);
        
        // We assume the data is stored as int64_t.
        int64_t *data = (int64_t *) mxGetData(field);
        if (!data) {
            mexPrintf("  No data in field '%s'.\n", fieldName);
            continue;
        }
        
        // Print each element.
        for (mwSize r = 0; r < m; r++) {
            mexPrintf("  Row %d: ", (int)r);
            for (mwSize c = 0; c < n; c++) {
                // MATLAB stores data in column-major order: index = r + c*m
                int64_t value = data[r + c * m];
                mexPrintf("%" PRId64 " ", value);
            }
            mexPrintf("\n");
        }
    }
    
    mexPrintf("Finished logging vStruct.\n");
}


/*------------------------------------------------------------------
 * LogMatrix: Helper to display a dd_MatrixPtr in MATLAB console.
 *------------------------------------------------------------------*/
static void LogMatrix(const dd_MatrixPtr M, const char *label)
{
    if (!M) {
        mexPrintf("[LogMatrix] %s is NULL.\n\n", label);
        return;
    }

    mexPrintf("[LogMatrix] %s:\n", label);
    
    /* Basic fields */
    mexPrintf("  rowsize      = %d\n", (int)M->rowsize);
    mexPrintf("  colsize      = %d\n", (int)M->colsize);
    mexPrintf("  representation: %s\n",
        (M->representation == dd_Generator)   ? "Generator" :
        (M->representation == dd_Inequality)  ? "Inequality" :
        (M->representation == dd_Unspecified) ? "Unspecified" : "Unknown");
    mexPrintf("  numbtype       = %s\n",
        (M->numbtype == dd_Rational) ? "Rational" :
        (M->numbtype == dd_Real)     ? "Real" : "Unknown");
    mexPrintf("  objective      = %d\n", (int)M->objective);
    
    /* linset is a set of row indices corresponding to lineality or equality rows */
    mexPrintf("  linset members:");
    if (M->linset) {
        dd_rowrange maxRow = M->rowsize;
        for (dd_rowrange i = 1; i <= maxRow; i++) {
            if (set_member(i, M->linset)) {
                mexPrintf(" %d", (int)i);
            }
        }
        mexPrintf("\n");
    } else {
        mexPrintf(" (none)\n");
    }

    /* rowvec is used by cddlib to store an objective row or similar. */
    if (M->rowvec) {
        mexPrintf("  rowvec (%d entries):\n", (int)M->colsize);
        for (dd_colrange j = 0; j < M->colsize; j++) {
            char *valStr = mpq_get_str(NULL, 10, M->rowvec[j]);
            mexPrintf("    rowvec[%d] = %s\n", (int)j, valStr);
            free(valStr);
        }
    } else {
        mexPrintf("  rowvec is NULL.\n");
    }
    
    /* Print the actual matrix entries. M->matrix is of type dd_Amatrix = mpq_t** */
    mexPrintf("  matrix (rows x cols):\n");
    for (dd_rowrange i = 0; i < M->rowsize; i++) {
        mexPrintf("    Row %d:", (int)i);
        for (dd_colrange j = 0; j < M->colsize; j++) {
            char *valStr = mpq_get_str(NULL, 10, M->matrix[i][j]);
            mexPrintf(" %s", valStr);
            free(valStr);
        }
        mexPrintf("\n");
    }
    
    mexPrintf("[End of LogMatrix] %s\n\n", label);
}


/*------------------------------------------------------------------
 * LogPolyhedra: Helper to display a dd_PolyhedraPtr in MATLAB console.
 *------------------------------------------------------------------*/
static void LogPolyhedra(const dd_PolyhedraPtr P, const char *label)
{
    if (!P) {
        mexPrintf("[LogPolyhedra] %s is NULL.\n\n", label);
        return;
    }

    mexPrintf("[LogPolyhedra] %s:\n", label);

    /* Representation type */
    mexPrintf("  representation  = %s\n",
        (P->representation == dd_Inequality) ? "Inequality" :
        (P->representation == dd_Generator)  ? "Generator" :
        (P->representation == dd_Unspecified)? "Unspecified" : "Unknown");
    /* Is it homogeneous? */
    mexPrintf("  homogeneous     = %d\n", (int)P->homogeneous);
    
    /* Dimensions and matrix size */
    mexPrintf("  d (columns)     = %d\n", (int)P->d);
    mexPrintf("  m (rows)        = %d\n", (int)P->m);
    
    /* The main matrix A (H- or V-representation) */
    /* You can optionally re-use LogMatrix here, but be cautious:
       P->A is a dd_Amatrix with dimension m x d. For a quick listing: */
    mexPrintf("  A (dd_Amatrix) entries:\n");
    for (dd_rowrange i = 0; i < P->m; i++) {
        mexPrintf("    Row %d:", (int)i);
        for (dd_colrange j = 0; j < P->d; j++) {
            char *valStr = mpq_get_str(NULL, 10, P->A[i][j]);
            mexPrintf(" %s", valStr);
            free(valStr);
        }
        mexPrintf("\n");
    }

    /* Number type, e.g., dd_Rational or dd_Real */
    mexPrintf("  numbtype        = %s\n",
        (P->numbtype == dd_Rational) ? "Rational" :
        (P->numbtype == dd_Real)     ? "Real" : "Unknown");
    
    /* Possibly child data: P->child is a dd_ConePtr, used for homogenized data */
    /* We skip logging child->... fields here, but you could dive deeper if desired. */
    
    /* Memory-allocation info (not usually critical for debugging, but can be printed): */
    mexPrintf("  m_alloc         = %d\n", (int)P->m_alloc);
    mexPrintf("  d_alloc         = %d\n", (int)P->d_alloc);
    
    /* The cost vector c is stored in P->c with length d (or maybe d+1). */
    if (P->c) {
        mexPrintf("  cost vector c (length %d):\n", (int)P->d);
        for (dd_colrange k = 0; k < P->d; k++) {
            char *valStr = mpq_get_str(NULL, 10, P->c[k]);
            mexPrintf("    c[%d] = %s\n", (int)k, valStr);
            free(valStr);
        }
    } else {
        mexPrintf("  cost vector c is NULL.\n");
    }

    /* EqualityIndex is an array of length m, telling if row i is eq/strict/ineq. */
    if (P->EqualityIndex) {
        mexPrintf("  EqualityIndex (length %d):\n", (int)P->m);
        for (dd_rowrange i = 0; i < P->m; i++) {
            mexPrintf("    EqualityIndex[%d] = %d\n", (int)i, (int)P->EqualityIndex[i]);
        }
    } else {
        mexPrintf("  EqualityIndex is NULL.\n");
    }

    /* Additional flags */
    mexPrintf("  IsEmpty             = %d\n", (int)P->IsEmpty);
    mexPrintf("  NondegAssumed       = %d\n", (int)P->NondegAssumed);
    mexPrintf("  InitBasisAtBottom   = %d\n", (int)P->InitBasisAtBottom);
    mexPrintf("  RestrictedEnumeration= %d\n", (int)P->RestrictedEnumeration);
    mexPrintf("  RelaxedEnumeration  = %d\n", (int)P->RelaxedEnumeration);

    /* m1 can be m or m+1, used internally by cddlib. */
    mexPrintf("  m1                  = %d\n", (int)P->m1);
    /* AincGenerated tells if some internal data is already computed. */
    mexPrintf("  AincGenerated       = %d\n", (int)P->AincGenerated);

    /* ldim: dimension of lineality space */
    mexPrintf("  ldim                = %d\n", (int)P->ldim);
    /* n: total number of rays in the computed cone + linearity dimension */
    mexPrintf("  n (size of output)  = %lld\n", (long long)P->n);

    /* Ainc, Ared, Adom are more advanced internal structures for cdd.
       Ainc is an incidence matrix. Ared, Adom are rowsets (similar to sets).
       We show an example for printing Ared or Adom. */

    if (P->Ared) {
        mexPrintf("  Ared (redundant rows):");
        for (dd_rowrange i = 1; i <= P->m; i++) {
            if (set_member(i, P->Ared)) {
                mexPrintf(" %d", (int)i);
            }
        }
        mexPrintf("\n");
    } else {
        mexPrintf("  Ared is NULL.\n");
    }

    if (P->Adom) {
        mexPrintf("  Adom (dominant rows):");
        for (dd_rowrange i = 1; i <= P->m; i++) {
            if (set_member(i, P->Adom)) {
                mexPrintf(" %d", (int)i);
            }
        }
        mexPrintf("\n");
    } else {
        mexPrintf("  Adom is NULL.\n");
    }

    /* For Ainc (dd_Aincidence), it’s typically an array of sets. You could do more detailed logging. */

    mexPrintf("[End of LogPolyhedra] %s\n\n", label);
}


// For LP solving, convert MATLAB structure to cddlib matrix without negating A.
// (This version is used by the LP routines.)
dd_MatrixPtr FT_get_H_MatrixPtr(const mxArray *in, int negateA) {
    // Check that 'in' is a struct.
    if (!mxIsStruct(in))
        mexErrMsgTxt("Input must be a struct.");

    // Get required fields.
    mxArray *Anum = mxGetField(in, 0, "Anum");
    mxArray *Aden = mxGetField(in, 0, "Aden");
    mxArray *Bnum = mxGetField(in, 0, "Bnum");
    mxArray *Bden = mxGetField(in, 0, "Bden");

    // Validate presence and type (using int64 checks).
    if (!(Anum && Aden && Bnum && Bden))
        mexErrMsgTxt("Missing one or more required fields: Anum, Aden, Bnum, Bden.");
    if (!mxIsInt64(Anum) || !mxIsInt64(Aden) ||
        !mxIsInt64(Bnum) || !mxIsInt64(Bden))
        mexErrMsgTxt("Anum, Aden, Bnum, Bden must all be int64 arrays.");

    // Get dimensions.
    mwSize m = mxGetM(Anum);
    mwSize nA = mxGetN(Anum);
    if (mxGetM(Bnum) != m || mxGetN(Bnum) != 1)
        mexErrMsgTxt("Bnum must be m x 1 int64.");
    if (mxGetM(Bden) != m || mxGetN(Bden) != 1)
        mexErrMsgTxt("Bden must be m x 1 int64.");
    if (mxGetM(Aden) != m || mxGetN(Aden) != nA)
        mexErrMsgTxt("Dimensions of Aden must match Anum.");

    // Create the matrix with m rows and (nA+1) columns.
    dd_MatrixPtr M = dd_CreateMatrix((dd_rowrange)m, (dd_colrange)(nA + 1));
    M->representation = dd_Inequality;
    M->numbtype = dd_Rational;

    // Get raw data pointers.
    int64_t *AnumData = (int64_t*) mxGetData(Anum);
    int64_t *AdenData = (int64_t*) mxGetData(Aden);
    int64_t *BnumData = (int64_t*) mxGetData(Bnum);
    int64_t *BdenData = (int64_t*) mxGetData(Bden);

    // Fill in column 0 with B values.
    for (mwSize i = 0; i < m; i++) {
        dd_set_si2(M->matrix[i][0], BnumData[i], BdenData[i]);
    }

    // Fill in the A coefficients (with optional negation).
    for (mwSize i = 0; i < m; i++) {
        for (mwSize j = 0; j < nA; j++) {
            int64_t a_num = AnumData[i + j * m];
            int64_t a_den = AdenData[i + j * m];
            // Apply negation if required.
            dd_set_si2(M->matrix[i][j+1], negateA ? -a_num : a_num, a_den);
        }
    }

    // Process the "lin" field if present.
    mxArray *linField = mxGetField(in, 0, "lin");
    if (linField) {
        // If lin is provided as int64.
        if (mxIsInt64(linField) && mxGetM(linField) == 1) {
            mwSize nLin = mxGetN(linField);
            int64_t *linData = (int64_t*) mxGetData(linField);
            for (mwSize idx = 0; idx < nLin; idx++) {
                int64_t rowIndex = linData[idx];
                if (rowIndex < 1 || rowIndex > (int64_t)m)
                    mexWarnMsgTxt("Ignoring out-of-range lineality index.");
                else
                    set_addelem(M->linset, (long)rowIndex);
            }
        }
        // Alternatively, handle double data (as in the original FT_get_H_MatrixPtr).
        else if (mxGetNumberOfDimensions(linField) <= 2 && mxGetM(linField) == 1) {
            double *lin = mxGetPr(linField);
            mwSize nLin = mxGetN(linField);
            for (mwSize i = 0; i < nLin; i++) {
                if (lin[i] <= (double)m)
                    set_addelem(M->linset, (long)lin[i]);
                else
                    mexWarnMsgTxt("Error in the lineality vector.");
            }
        }
    }
    return M;
}


// FT_set_H_MatrixPtr: Convert a dd_MatrixPtr (H-representation) to MATLAB structure.
mxArray* FT_set_H_MatrixPtr(const dd_MatrixPtr M) {
    // Validate that M is non-NULL, has an inequality representation,
    // and uses rational numbers.
    if (!M || M->representation != dd_Inequality) {
        mexErrMsgTxt("FT_set_H_MatrixPtr: invalid or non-inequality matrix.");
    }
    
    // m is the number of constraints, n is 1 + the problem dimension.
    dd_rowrange m = M->rowsize;
    dd_colrange n = M->colsize;
    
    // Create the output MATLAB struct with fields: Anum, Aden, Bnum, Bden, and lin.
    const char *fieldNames[] = {"Anum", "Aden", "Bnum", "Bden", "lin"};
    mxArray *outStruct = mxCreateStructMatrix(1, 1, 5, fieldNames);
    
    // Allocate numeric arrays:
    // For A, dimensions: m x (n-1); for B, dimensions: m x 1.
    mxArray *AnumArr = mxCreateNumericMatrix(m, n - 1, mxINT64_CLASS, mxREAL);
    mxArray *AdenArr = mxCreateNumericMatrix(m, n - 1, mxINT64_CLASS, mxREAL);
    mxArray *BnumArr = mxCreateNumericMatrix(m, 1, mxINT64_CLASS, mxREAL);
    mxArray *BdenArr = mxCreateNumericMatrix(m, 1, mxINT64_CLASS, mxREAL);
    
    int64_t *AnumData = (int64_t*) mxGetData(AnumArr);
    int64_t *AdenData = (int64_t*) mxGetData(AdenArr);
    int64_t *BnumData = (int64_t*) mxGetData(BnumArr);
    int64_t *BdenData = (int64_t*) mxGetData(BdenArr);
    
    // Loop over each row to fill in the right-hand side (B) and the coefficients for A.
    for (dd_rowrange i = 0; i < m; i++) {
        int64_t num, den;
        // Column 0 of M holds the B coefficient.
        dd_get_si2(M->matrix[i][0], &num, &den);
        BnumData[i] = num;
        BdenData[i] = den;
        
        // For each column of A (columns 1 to n-1 in M),
        // store the NEGATIVE of the original value.
        for (dd_colrange j = 1; j < n; j++) {
            dd_get_si2(M->matrix[i][j], &num, &den);
            // The merged version stores negative of A
            AnumData[i + (j - 1) * m] = -num;
            AdenData[i + (j - 1) * m] = den;
        }
    }
    
    // Set the fields for Anum, Aden, Bnum, and Bden.
    mxSetField(outStruct, 0, "Anum", AnumArr);
    mxSetField(outStruct, 0, "Aden", AdenArr);
    mxSetField(outStruct, 0, "Bnum", BnumArr);
    mxSetField(outStruct, 0, "Bden", BdenArr);
    
    // Process the "lin" (lineality) field:
    dd_rowrange countLin = set_card(M->linset);
    if (countLin > 0) {
        // Create a 1 x countLin int64 array for the lineality indices.
        mxArray *linArr = mxCreateNumericMatrix(1, countLin, mxINT64_CLASS, mxREAL);
        int64_t *linData = (int64_t*) mxGetData(linArr);
        dd_rowrange idx = 0;
        // cddlib uses 1-based indexing.
        for (dd_rowrange r = 1; r <= m; r++) {
            if (set_member(r, M->linset)) {
                linData[idx++] = (int64_t) r;
            }
        }
        mxSetField(outStruct, 0, "lin", linArr);
    } else {
        // If no lineality exists, return an empty int64 array.
        mxArray *linEmpty = mxCreateNumericMatrix(1, 0, mxINT64_CLASS, mxREAL);
        mxSetField(outStruct, 0, "lin", linEmpty);
    }
    
    return outStruct;
}


// FT_get_V_MatrixPtr: Convert MATLAB V-structure (with Vnum/Vden and optional Rnum/Rden) to dd_MatrixPtr.
dd_MatrixPtr FT_get_V_MatrixPtr(const mxArray *in) {
    mxArray *tmpv_num = mxGetField(in, 0, "Vnum");
    mxArray *tmpv_den = mxGetField(in, 0, "Vden");
    if (!(tmpv_num && tmpv_den) ||
        (mxGetNumberOfDimensions(tmpv_num) > 2) ||
        (mxGetNumberOfDimensions(tmpv_den) > 2)) {
        return 0;
    }
    mxArray *tmpr_num = mxGetField(in, 0, "Rnum");
    mxArray *tmpr_den = mxGetField(in, 0, "Rden");
    int mr, m, n, i, j;
    int64_t *VnumData = (int64_t *) mxGetData(tmpv_num);
    int64_t *VdenData = (int64_t *) mxGetData(tmpv_den);
    int64_t *RnumData = NULL;
    int64_t *RdenData = NULL;
    if (tmpr_num && tmpr_den &&
        (mxGetNumberOfDimensions(tmpr_num) <= 2) &&
        (mxGetNumberOfDimensions(tmpr_den) <= 2)) {
        mr = mxGetM(tmpr_num);
        RnumData = (int64_t *) mxGetData(tmpr_num);
        RdenData = (int64_t *) mxGetData(tmpr_den);
    } else {
        mr = 0;
    }
    m = mxGetM(tmpv_num);
    n = mxGetN(tmpv_num) + 1;  /* +1 for homogenizing coordinate */
    dd_MatrixPtr V = dd_CreateMatrix(m + mr, n);
    for (i = 0; i < m; i++) {
        dd_set_si2(V->matrix[i][0], 1, 1);
        for (j = 0; j < n - 1; j++) {
            dd_set_si2(V->matrix[i][j+1], VnumData[i + j * m], VdenData[i + j * m]);
        }
    }
    for (i = m; i < m + mr; i++) {
        dd_set_si2(V->matrix[i][0], 0, 1);
        for (j = 0; j < n - 1; j++) {
            dd_set_si2(V->matrix[i][j+1], RnumData[(i - m) + j * mr], RdenData[(i - m) + j * mr]);
        }
    }
    V->representation = dd_Generator;
    V->numbtype = dd_Rational;
    return V;
}


// FT_set_V_MatrixPtr: Convert cddlib generator matrix (V-representation) to MATLAB structure.
mxArray *FT_set_V_MatrixPtr(const dd_MatrixPtr M) {
    mxArray *P;
    int mr, mv, m, i, j;
    int ir, iv;
    mxArray *tmpv_num, *tmpv_den;
    mxArray *tmpr_num, *tmpr_den;
    mxArray *tmpvpos, *tmprpos;
    mxArray *tmpl;
    int64_t *vnum, *vden;
    int64_t *rnum, *rden;
    int64_t *vpos, *rpos;
    int64_t *l;
    int k = 0;
    dd_rowrange ii;
    if ((M != NULL) &&
        (M->representation == dd_Generator)) 
        //(M->numbtype == dd_Rational)) 
        {
        const char *f[] = {"Vnum", "Vden", "Rnum", "Rden", "vpos", "rpos", "lin"};
        mwSize dims[] = {1};
        P = mxCreateStructArray(1, dims, 7, f);
        if (set_card(M->linset)) {
            tmpl = mxCreateNumericMatrix(1, set_card(M->linset), mxINT64_CLASS, mxREAL);
            l = (int64_t *) mxGetData(tmpl);
            for (ii = 1; ii <= M->rowsize; ii++) {
                if (set_member(ii, M->linset))
                    l[k++] = (int64_t) ii;
            }
            mxSetField(P, 0, "lin", tmpl);
        }
        mr = 0;
        mv = 0;
        m = (int) M->rowsize;
        for (i = 0; i < m; i++) {
            int64_t num0, den0;
            dd_get_si2(M->matrix[i][0], &num0, &den0);
            if (num0 == 0)
                mr++;
            else
                mv++;
        }
        tmpr_num = mxCreateNumericMatrix(mr, M->colsize - 1, mxINT64_CLASS, mxREAL);
        tmpr_den = mxCreateNumericMatrix(mr, M->colsize - 1, mxINT64_CLASS, mxREAL);
        tmpv_num = mxCreateNumericMatrix(mv, M->colsize - 1, mxINT64_CLASS, mxREAL);
        tmpv_den = mxCreateNumericMatrix(mv, M->colsize - 1, mxINT64_CLASS, mxREAL);
        tmprpos = mxCreateNumericMatrix(mr, 1, mxINT64_CLASS, mxREAL);
        tmpvpos = mxCreateNumericMatrix(mv, 1, mxINT64_CLASS, mxREAL);
        rnum = (int64_t *) mxGetData(tmpr_num);
        rden = (int64_t *) mxGetData(tmpr_den);
        vnum = (int64_t *) mxGetData(tmpv_num);
        vden = (int64_t *) mxGetData(tmpv_den);
        rpos = (int64_t *) mxGetData(tmprpos);
        vpos = (int64_t *) mxGetData(tmpvpos);
        ir = 0;
        iv = 0;
        for (i = 0; i < m; i++) {
            int64_t num0, den0;
            dd_get_si2(M->matrix[i][0], &num0, &den0);
            if (num0 == 0) {
                for (j = 0; j < (int) M->colsize - 1; j++) {
                    int64_t a_num, a_den;
                    dd_get_si2(M->matrix[i][j+1], &a_num, &a_den);
                    rnum[ir + j * mr] = a_num;
                    rden[ir + j * mr] = a_den;
                }
                rpos[ir] = (int64_t) (i + 1);
                ir++;
            } else {
                for (j = 0; j < (int) M->colsize - 1; j++) {
                    int64_t a_num, a_den;
                    dd_get_si2(M->matrix[i][j+1], &a_num, &a_den);
                    vnum[iv + j * mv] = a_num;
                    vden[iv + j * mv] = a_den;
                }
                vpos[iv] = (int64_t) (i + 1);
                iv++;
            }
        }
        mxSetField(P, 0, "Vnum", tmpv_num);
        mxSetField(P, 0, "Vden", tmpv_den);
        mxSetField(P, 0, "Rnum", tmpr_num);
        mxSetField(P, 0, "Rden", tmpr_den);
        mxSetField(P, 0, "vpos", tmpvpos);
        mxSetField(P, 0, "rpos", tmprpos);
        return P;
    }
    return 0;
}


// Convert LP solution to MATLAB rational structure.
mxArray *MB_set_LPsol_MatrixPtr(const dd_LPPtr lp) {
    if (!lp)
        return NULL;
    mxArray *P;
    const char *fields[] = {"Homogeneous", "m", "d", "optvalue", "sol", "dsol", "how"};
    mwSize dims[1] = {1};
    P = mxCreateStructArray(1, dims, 7, fields);
    mxSetField(P, 0, "Homogeneous", mxCreateLogicalScalar(lp->Homogeneous));
    mxSetField(P, 0, "m", mxCreateDoubleScalar((double)lp->m));
    mxSetField(P, 0, "d", mxCreateDoubleScalar((double)lp->d));
    mxArray *optvalueStruct = mxCreateStructMatrix(1, 1, 2, (const char *[]){"num", "den"});
    mxSetField(optvalueStruct, 0, "num", mxCreateDoubleScalar(mpz_get_d(mpq_numref(lp->optvalue))));
    mxSetField(optvalueStruct, 0, "den", mxCreateDoubleScalar(mpz_get_d(mpq_denref(lp->optvalue))));
    mxSetField(P, 0, "optvalue", optvalueStruct);
    mxArray *solNum = mxCreateDoubleMatrix(lp->d, 1, mxREAL);
    mxArray *solDen = mxCreateDoubleMatrix(lp->d, 1, mxREAL);
    double *solNumData = mxGetPr(solNum);
    double *solDenData = mxGetPr(solDen);
    for (int ii = 0; ii < lp->d; ii++) {
        solNumData[ii] = mpz_get_d(mpq_numref(lp->sol[ii]));
        solDenData[ii] = mpz_get_d(mpq_denref(lp->sol[ii]));
    }
    mxArray *solStruct = mxCreateStructMatrix(1, 1, 2, (const char *[]){"num", "den"});
    mxSetField(solStruct, 0, "num", solNum);
    mxSetField(solStruct, 0, "den", solDen);
    mxSetField(P, 0, "sol", solStruct);
    mxArray *dsolNum = mxCreateDoubleMatrix(lp->m, 1, mxREAL);
    mxArray *dsolDen = mxCreateDoubleMatrix(lp->m, 1, mxREAL);
    double *dsolNumData = mxGetPr(dsolNum);
    double *dsolDenData = mxGetPr(dsolDen);
    for (int ii = 0; ii < lp->m; ii++) {
        dsolNumData[ii] = mpz_get_d(mpq_numref(lp->dsol[ii]));
        dsolDenData[ii] = mpz_get_d(mpq_denref(lp->dsol[ii]));
    }
    mxArray *dsolStruct = mxCreateStructMatrix(1, 1, 2, (const char *[]){"num", "den"});
    mxSetField(dsolStruct, 0, "num", dsolNum);
    mxSetField(dsolStruct, 0, "den", dsolDen);
    mxSetField(P, 0, "dsol", dsolStruct);
    mxSetField(P, 0, "how", mxCreateDoubleScalar((double)lp->LPS));
    return P;
}


// Convert MATLAB structure to LP structure for LP solving.
dd_LPPtr MB_get_LP_MatrixPtr(const mxArray *in) {
    dd_MatrixPtr A;
    dd_ErrorType error = dd_NoError;
    dd_LPPtr lp;
    dd_set_global_constants();
    A = FT_get_H_MatrixPtr(in, 0);
    if (!A)
        mexErrMsgTxt("Error constructing the H matrix from the input structure.");
    mxArray *tmpObjNum = mxGetField(in, 0, "objNum");
    mxArray *tmpObjDen = mxGetField(in, 0, "objDen");
    if (!tmpObjNum || !tmpObjDen) {
        dd_FreeMatrix(A);
        mexErrMsgTxt("Missing or invalid fields 'objNum' or 'objDen' in LP structure.");
    }
    mwSize n = mxGetNumberOfElements(tmpObjNum);
    if (n != A->colsize) {
        dd_FreeMatrix(A);
        mexErrMsgTxt("Objective vector size does not match the expected dimension.");
    }
    if (!mxIsInt64(tmpObjNum) || !mxIsInt64(tmpObjDen)) {
        dd_FreeMatrix(A);
        mexErrMsgTxt("Objective numerator and denominator must be of type int64.");
    }
    int64_t *objNumData = (int64_t *)mxGetData(tmpObjNum);
    int64_t *objDenData = (int64_t *)mxGetData(tmpObjDen);
    for (mwSize j = 0; j < n; j++) {
        if (objDenData[j] == 0) {
            dd_FreeMatrix(A);
            mexErrMsgTxt("Objective denominator contains zero, which is invalid.");
        }
    }
    for (mwSize j = 0; j < n; j++) {
        dd_set_si2(A->rowvec[j], objNumData[j], objDenData[j]);
    }
    A->objective = dd_LPmax;
    lp = dd_Matrix2LP(A, &error);
    dd_FreeMatrix(A);
    if (error != dd_NoError) {
        mexErrMsgTxt("Error converting matrix to LP in dd_Matrix2LP.");
    }
    return lp;
}


mxArray * FT_set_Set(const set_type S)
{
    mxArray * tmp;
    int64_t * r;
    int i, j;
    int num;
    
    if (S) {
        num = (int)(S[0]);
        j = 0;
        for (i = 1; i <= num; i++) {
            if (set_member(i,S)) {
                j++; 
            }
        }
        tmp = mxCreateNumericMatrix(1, j, mxINT64_CLASS, mxREAL);
        r = (int64_t *) mxGetData(tmp);
        j = 0;
        for (i = 1; i <= num; i++) {
            if (set_member(i,S)) {
                r[j++] = (int64_t) i;
            }
        }	
        return tmp;
    }
    return 0;
}


// FT_set_SetFamilyPtr: Convert dd_SetFamilyPtr (adjacency info) to a MATLAB cell array.
mxArray *FT_set_SetFamilyPtr(const dd_SetFamilyPtr F) {
    if (!F) {
        return 0;
    }
    int famsize = (int)(F->famsize);
    // Create a cell array with one cell per set.
    mxArray *cellArray = mxCreateCellMatrix(famsize, 1);
    for (int i = 0; i < famsize; i++) {
        // In cddlib, F->set[i][0] holds the maximum possible index.
        int maxIndex = (int)(F->set[i][0]);
        int count = 0;
        // Count how many members are actually present.
        for (int j = 1; j <= maxIndex; j++) {
            if (set_member(j, F->set[i])) {
                count++;
            }
        }
        // Allocate a 1-by-count int64 array.
        mxArray *tmp = mxCreateNumericMatrix(1, count, mxINT64_CLASS, mxREAL);
        int64_t *data = (int64_t *) mxGetData(tmp);
        int pos = 0;
        // Fill the array with the indices of members.
        for (int j = 1; j <= maxIndex; j++) {
            if (set_member(j, F->set[i])) {
                data[pos++] = (int64_t) j;
            }
        }
        mxSetCell(cellArray, i, tmp);
    }
    return cellArray;
}


// LP solving using the default CrissCross solver.
void solve_lp(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    dd_ErrorType error = dd_NoError;
    dd_LPSolverType solver = dd_CrissCross;
    dd_LPPtr lp;
    if (nrhs != 1)
        mexErrMsgTxt("solve_lp requires exactly one input argument (the LP structure).");
    if (nlhs != 1)
        mexErrMsgTxt("solve_lp produces exactly one output argument (the LP solution).");
    dd_set_global_constants();
    lp = MB_get_LP_MatrixPtr(prhs[0]);
    if (lp == NULL)
        mexErrMsgTxt("Failed to convert MATLAB input to LP matrix. Ensure the input structure is valid.");
    // Optionally log fields:
    LogLpFields(lp);
    dd_LPSolve(lp, solver, &error);
    if (error != dd_NoError) {
        dd_WriteErrorMessages(stdout, error);
        dd_FreeLPData(lp);
        mexErrMsgTxt("Error solving LP using cddlib. Check input parameters or constraints.");
    }
    plhs[0] = MB_set_LPsol_MatrixPtr(lp);
    if (plhs[0] == NULL) {
        dd_FreeLPData(lp);
        mexErrMsgTxt("Failed to convert LP solution to MATLAB format.");
    }
    dd_FreeLPData(lp);
}


// LP solving using the Dual Simplex method.
void solve_lp_ds(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    dd_ErrorType error = dd_NoError;
    dd_LPSolverType solver = dd_DualSimplex;
    dd_LPPtr lp;
    if (nrhs != 1)
        mexErrMsgTxt("solve_lp_DS requires exactly one input argument (the LP structure).");
    if (nlhs != 1)
        mexErrMsgTxt("solve_lp_DS produces exactly one output argument (the LP solution).");
    dd_set_global_constants();
    lp = MB_get_LP_MatrixPtr(prhs[0]);
    if (lp == NULL)
        mexErrMsgTxt("Failed to convert MATLAB input to LP matrix. Ensure the input structure is valid.");
    dd_LPSolve(lp, solver, &error);
    if (error != dd_NoError) {
        dd_WriteErrorMessages(stdout, error);
        dd_FreeLPData(lp);
        mexErrMsgTxt("Error solving LP using the Dual Simplex method in cddlib. Check input constraints.");
    }
    plhs[0] = MB_set_LPsol_MatrixPtr(lp);
    if (plhs[0] == NULL) {
        dd_FreeLPData(lp);
        mexErrMsgTxt("Failed to convert LP solution to MATLAB format.");
    }
    dd_FreeLPData(lp);
}


// adj_extreme: Compute V-representation and adjacency from an H-structure.
void adj_extreme(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    dd_PolyhedraPtr P;
    dd_ErrorType err;
    dd_MatrixPtr H, V;
    dd_SetFamilyPtr A;
    
    if (nrhs == 1 && nlhs == 2 && mxIsStruct(prhs[0])) {
        dd_set_global_constants();
        
        // Get the H-matrix without any negation.
        H = FT_get_H_MatrixPtr(prhs[0], 0);
        if (!H)
            mexErrMsgTxt("Failed to construct H-matrix from input structure.");
        //LogMatrix(H, "H (original)");
        
        // Negate the inequality coefficients (columns 1 to n-1; column 0 is the RHS).
        mwSize m = H->rowsize;
        mwSize n = H->colsize;
        for (mwSize i = 0; i < m; i++) {
            for (mwSize j = 1; j < n; j++) {
                int64_t num, den;
                dd_get_si2(H->matrix[i][j], &num, &den);
                dd_set_si2(H->matrix[i][j], -num, den);
            }
        }
        
        // Convert the modified H-matrix to a polyhedron.
        P = dd_DDMatrix2Poly(H, &err);
        if (err == dd_NoError) {
            V = dd_CopyGenerators(P);
            A = dd_CopyAdjacency(P);
            //LogMatrix(V, "V (converted)");
            plhs[0] = FT_set_V_MatrixPtr(V);
            plhs[1] = FT_set_SetFamilyPtr(A);
            dd_FreeMatrix(V);
            dd_FreeSetFamily(A);
        } else {
            dd_WriteErrorMessages(stdout, err);
            mexErrMsgTxt("CDD returned an error, see above(!) for details");
        }
        dd_FreeMatrix(H);
        dd_FreePolyhedra(P);
        return;
    } else {
        mexErrMsgTxt("adj_extreme expects an H input struct and produces a V output struct and the adjacency struct");
    }
}


// extreme: Compute V-representation from an H-structure.
void extreme(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    dd_PolyhedraPtr P;
    dd_ErrorType err;
    dd_MatrixPtr H, V;

    if (nrhs == 1 && nlhs == 1 && mxIsStruct(prhs[0])) {
        dd_set_global_constants();
        
        // 1) Load H
        H = FT_get_H_MatrixPtr(prhs[0], 0);
        if (!H) {
            mexErrMsgTxt("Failed to construct H-matrix from input structure.");
        }

        /* Log H as read in. */
        //LogMatrix(H, "H (original)");

        // 2) Negate inequalities
        for (mwSize i = 0; i < H->rowsize; i++) {
            for (mwSize j = 1; j < H->colsize; j++) {
                int64_t num, den;
                dd_get_si2(H->matrix[i][j], &num, &den);
                dd_set_si2(H->matrix[i][j], -num, den);
            }
        }
        /* Log H after negation. */
        //LogMatrix(H, "H (negated)");

        // 3) Convert to polyhedron
        P = dd_DDMatrix2Poly(H, &err);
        if (err != dd_NoError) {
            dd_WriteErrorMessages(stdout, err);
            dd_FreeMatrix(H);
            dd_FreePolyhedra(P);
            mexErrMsgTxt("CDD returned an error, see above for details");
        }
        
        /* Log the polyhedra P. */
        //LogPolyhedra(P, "P after dd_DDMatrix2Poly");

        // 4) Copy generators
        V = dd_CopyGenerators(P);
        if (!V) {
            mexWarnMsgTxt("dd_CopyGenerators returned NULL, assigning empty output.");
            plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
        } else {
            /* Log V. */
            //LogMatrix(V, "V (generators)");

            mxArray *vStruct = FT_set_V_MatrixPtr(V);
            if (!vStruct) {
                mexWarnMsgTxt("FT_set_V_MatrixPtr returned NULL, assigning empty output.");
                plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
            } else {
                //LogVStruct(vStruct);
                plhs[0] = vStruct;
            }
            dd_FreeMatrix(V);
        }

        dd_FreeMatrix(H);
        dd_FreePolyhedra(P);
    } else {
        mexErrMsgTxt("extreme expects an H input struct and produces a V output struct");
    }
}



// hull: Convert V-representation (MATLAB structure) to H-representation.
void hull(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    dd_PolyhedraPtr P;
    dd_ErrorType err;
    dd_MatrixPtr H, V;

    if (nrhs == 1 && nlhs == 1 && mxIsStruct(prhs[0])) {
        V = FT_get_V_MatrixPtr(prhs[0]);
        //LogMatrix(V, "V (original)");
        dd_set_global_constants();

        P = dd_DDMatrix2Poly(V, &err);
        if (err == dd_NoError) {
            H = dd_CopyInequalities(P);
            
            // Explicitly set the representation to dd_Inequality
            H->representation = dd_Inequality;

            plhs[0] = FT_set_H_MatrixPtr(H);
            //LogMatrix(H, "H (inequality)");
            dd_FreeMatrix(H);
        } else {
            dd_WriteErrorMessages(stdout, err);
            mexErrMsgTxt("CDD returned an error, see above for details");
        }

        dd_FreeMatrix(V);
        dd_FreePolyhedra(P);
        return;
    } else {
        mexErrMsgTxt("hull expects a V input struct and produces an H output struct");
    }
}



/* 
 * adjacency: Computes the input adjacency of a polyhedron from its V-representation.
 * Input: MATLAB struct (V-structure)
 * Output: MATLAB cell array (set family)
 */
void adjacency(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
    dd_PolyhedraPtr P;
    dd_ErrorType err;
    dd_MatrixPtr V;
    dd_SetFamilyPtr A;
    
    if (mxIsStruct(prhs[0])) {
        V = FT_get_V_MatrixPtr(prhs[0]);
        //LogMatrix(V, "V (original)");
        dd_set_global_constants();
        
        P = dd_DDMatrix2Poly(V, &err);
        if (err == dd_NoError) {
            A = dd_CopyInputAdjacency(P);
            //LogSetFamily(A);
            plhs[0] = FT_set_SetFamilyPtr(A);
            dd_FreeSetFamily(A);
        } else {
            dd_WriteErrorMessages(stdout, err);
            mexErrMsgTxt("CDD returned an error, see above(!) for details");
        }
        dd_FreeMatrix(V);
        dd_FreePolyhedra(P);
    } else {
        mexErrMsgTxt("adjacency expects a V input struct");
    }
}

/* 
 * copy_h: Copies an H-representation.
 * Input: MATLAB struct (H-structure)
 * Output: MATLAB struct (copied H-structure)
 */
void copy_h(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
    dd_MatrixPtr H;
    
    if (nrhs == 1 && nlhs == 1 && mxIsStruct(prhs[0])) {
        dd_set_global_constants();
        H = FT_get_H_MatrixPtr(prhs[0], 1);
        //LogMatrix(H, "H (original)");
        plhs[0] = FT_set_H_MatrixPtr(H);
        dd_FreeMatrix(H);
        return;
    } else {
        mexErrMsgTxt("copy_h expects a H input struct and produces a H output struct");
    }
}

/* 
 * copy_v: Copies a V-representation.
 * Input: MATLAB struct (V-structure)
 * Output: MATLAB struct (copied V-structure)
 */
void copy_v(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
    dd_MatrixPtr V;
    
    if (nrhs == 1 && nlhs == 1 && mxIsStruct(prhs[0])) {
        dd_set_global_constants();
        V = FT_get_V_MatrixPtr(prhs[0]);
        plhs[0] = FT_set_V_MatrixPtr(V);
        dd_FreeMatrix(V);
        return;
    } else {
        mexErrMsgTxt("copy_v expects a V input struct and produces a V output struct");
    }
}

/* 
 * find_interior: Finds an interior point for a polyhedron in H-representation
 * using the CrissCross method.
 * Input: MATLAB struct (H-structure)
 * Output: MATLAB struct (LP solution structure)
 */
void find_interior(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
    dd_ErrorType error = dd_NoError;
    dd_LPSolverType solver = dd_CrissCross;  /* Alternatively, dd_DualSimplex */
    dd_LPPtr lp, lp1;
    dd_MatrixPtr A;
    int j;
    
    dd_set_global_constants();
    if ((A = FT_get_H_MatrixPtr(prhs[0], 1))) {
        A->objective = dd_LPmin;
        for (j = 0; j < A->colsize; j++)
            dd_set_d(A->rowvec[j], 0.0);
        lp = dd_Matrix2LP(A, &error);
        dd_FreeMatrix(A);
    } else {
        mexErrMsgTxt("Error in the setting of LP matrix.");
        return;
    }
    
    lp1 = dd_MakeLPforInteriorFinding(lp);
    dd_LPSolve(lp1, solver, &error);
    if (error != dd_NoError)
        dd_WriteErrorMessages(stdout, error);
    
    plhs[0] = MB_set_LPsol_MatrixPtr(lp1);
    
    dd_FreeLPData(lp);
    dd_FreeLPData(lp1);
}

/* 
 * find_interior_DS: Finds an interior point for a polyhedron in H-representation
 * using the Dual Simplex method.
 * Input: MATLAB struct (H-structure)
 * Output: MATLAB struct (LP solution structure)
 */
void find_interior_DS(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
    dd_ErrorType error = dd_NoError;
    dd_LPSolverType solver = dd_DualSimplex;
    dd_LPPtr lp, lp1;
    dd_MatrixPtr A;
    int j;
    
    dd_set_global_constants();
    if ((A = FT_get_H_MatrixPtr(prhs[0], 0))) {
        A->objective = dd_LPmin;
        for (j = 0; j < A->colsize; j++)
            dd_set_d(A->rowvec[j], 0.0);
        lp = dd_Matrix2LP(A, &error);
        dd_FreeMatrix(A);
    } else {
        mexErrMsgTxt("Error in the setting of LP matrix.");
        return;
    }
    
    lp1 = dd_MakeLPforInteriorFinding(lp);
    dd_LPSolve(lp1, solver, &error);
    if (error != dd_NoError)
        dd_WriteErrorMessages(stdout, error);
    
    plhs[0] = MB_set_LPsol_MatrixPtr(lp1);
    
    dd_FreeLPData(lp);
    dd_FreeLPData(lp1);
}

/* 
 * reduce_h: Removes redundant inequalities from an H-representation.
 * Input: MATLAB struct (H-structure)
 * Output: MATLAB struct (reduced H-structure) and optionally a vector of removed rows.
 */
void reduce_h(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
    dd_ErrorType err;
    dd_MatrixPtr H, H1;
    dd_rowset red;
    
    if (nrhs == 1 && nlhs >= 1 && nlhs <= 2 && mxIsStruct(prhs[0])) {
        dd_set_global_constants();
        H = FT_get_H_MatrixPtr(prhs[0], 1);
        //LogMatrix(H, "H (original)");
        
        red = dd_RedundantRows(H, &err);
        if (err == dd_NoError) {
            H1 = dd_MatrixSubmatrix(H, red);
            //LogMatrix(H1, "H1 (reduced)");
            plhs[0] = FT_set_H_MatrixPtr(H1);
            dd_FreeMatrix(H1);
            if (nlhs == 2) {
                plhs[1] = FT_set_Set(red);
            }
        } else {
            dd_WriteErrorMessages(stdout, err);
            mexErrMsgTxt("CDD returned an error, see above(!) for details");
        }
        dd_FreeMatrix(H);
        set_free(red);
        return;
    } else {
        mexErrMsgTxt("reduce_h expects an H input struct and produces an H output struct and an optional vector of removed rows");
    }
}

/* 
 * reduce_v: Removes redundant generators from a V-representation.
 * Input: MATLAB struct (V-structure)
 * Output: MATLAB struct (reduced V-structure) and optionally a vector of removed vertices.
 */
void reduce_v(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
    dd_ErrorType err;
    dd_MatrixPtr V, V1;
    dd_rowset red;
    
    if (nrhs == 1 && nlhs >= 1 && nlhs <= 2 && mxIsStruct(prhs[0])) {
        dd_set_global_constants();
        V = FT_get_V_MatrixPtr(prhs[0]);
        
        red = dd_RedundantRows(V, &err);
        if (err == dd_NoError) {
            V1 = dd_MatrixSubmatrix(V, red);
            plhs[0] = FT_set_V_MatrixPtr(V1);
            dd_FreeMatrix(V1);
            if (nlhs == 2) {
                plhs[1] = FT_set_Set(red);
            }
        } else {
            dd_WriteErrorMessages(stdout, err);
            mexErrMsgTxt("CDD returned an error, see above(!) for details");
        }
        dd_FreeMatrix(V);
        set_free(red);
        return;
    } else {
        mexErrMsgTxt("reduce_v expects a V input struct and produces a V output struct and an optional vector of removed vertices");
    }
}

/* 
 * v_hull_extreme: Computes the vertex hull extreme representation from a V-structure.
 * Input: MATLAB struct (V-structure)
 * Output: MATLAB struct (computed extreme generators in V-representation)
 */
void v_hull_extreme(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
    dd_PolyhedraPtr P, P1;
    dd_ErrorType err;
    dd_MatrixPtr H, V, V1;
    
    if (nrhs == 1 && nlhs == 1 && mxIsStruct(prhs[0])) {
        V = FT_get_V_MatrixPtr(prhs[0]);
        dd_set_global_constants();
        
        P = dd_DDMatrix2Poly(V, &err);
        if (err == dd_NoError) {
            H = dd_CopyInequalities(P);
            P1 = dd_DDMatrix2Poly(H, &err);
            if (err == dd_NoError) {
                V1 = dd_CopyGenerators(P1);
                plhs[0] = FT_set_V_MatrixPtr(V1);
                dd_FreeMatrix(V1);
            } else {
                dd_WriteErrorMessages(stdout, err);
                mexErrMsgTxt("CDD returned an error, see above(!) for details");
            }
            dd_FreePolyhedra(P1);
            dd_FreeMatrix(H);
        } else {
            dd_WriteErrorMessages(stdout, err);
            mexErrMsgTxt("CDD returned an error, see above(!) for details");
        }
        dd_FreeMatrix(V);
        dd_FreePolyhedra(P);
        return;
    } else {
        mexErrMsgTxt("v_hull_extreme expects a V input struct and produces a V output struct");
    }
}


/* 
 * The gateway function expects a command string as the first argument,
 * followed by the appropriate input structure(s).
 */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int buflen, status;
    char *input_buf;
    
    if (nrhs < 1) {
        mexErrMsgTxt("At least one input required: a command string.");
    }
    
    if (mxIsChar(prhs[0]) && (mxGetM(prhs[0]) == 1)) {
        buflen = mxGetN(prhs[0]) + 1;
        input_buf = mxCalloc(buflen, sizeof(char));
        status = mxGetString(prhs[0], input_buf, buflen);
        if (status != 0) {
            mexErrMsgTxt("Error reading command string.");
        }
        
        /* Dispatch commands */
        if (strcmp(input_buf, "solve_lp") == 0) {
            if (nrhs < 2)
                mexErrMsgTxt("solve_lp requires an LP structure as additional input.");
            solve_lp(nlhs, plhs, nrhs - 1, prhs + 1);
        }
        else if (strcmp(input_buf, "solve_lp_ds") == 0) {
            if (nrhs < 2)
                mexErrMsgTxt("solve_lp_ds requires an LP structure as additional input.");
            solve_lp_ds(nlhs, plhs, nrhs - 1, prhs + 1);
        }
        else if (strcmp(input_buf, "adj_extreme") == 0) {
            if (nrhs < 2)
                mexErrMsgTxt("adj_extreme requires an H-structure as additional input.");
            adj_extreme(nlhs, plhs, nrhs - 1, prhs + 1);
        }
        else if (strcmp(input_buf, "extreme") == 0) {
            if (nrhs < 2)
                mexErrMsgTxt("extreme requires an H-structure as additional input.");
            extreme(nlhs, plhs, nrhs - 1, prhs + 1);
        }
        else if (strcmp(input_buf, "hull") == 0) {
            if (nrhs < 2)
                mexErrMsgTxt("hull requires a V-structure as additional input.");
            hull(nlhs, plhs, nrhs - 1, prhs + 1);
        }
        else if (strcmp(input_buf, "version") == 0) {
            mexPrintf("CDDMEX version: %s\n", CDDMEX_VERSION);
            if (nlhs > 0)
                plhs[0] = mxCreateString(CDDMEX_VERSION);
        }
        else if (strcmp(input_buf, "adjacency") == 0) {
            adjacency(nlhs, plhs, nrhs - 1, prhs + 1);
        }
        else if (strcmp(input_buf, "copy_h") == 0) {
            copy_h(nlhs, plhs, nrhs - 1, prhs + 1);
        }
        else if (strcmp(input_buf, "copy_v") == 0) {
            copy_v(nlhs, plhs, nrhs - 1, prhs + 1);
        }
        else if (strcmp(input_buf, "find_interior") == 0) {
            find_interior(nlhs, plhs, nrhs - 1, prhs + 1);
        }
        else if (strcmp(input_buf, "find_interior_DS") == 0) {
            find_interior_DS(nlhs, plhs, nrhs - 1, prhs + 1);
        }
        else if (strcmp(input_buf, "reduce_h") == 0) {
            reduce_h(nlhs, plhs, nrhs - 1, prhs + 1);
        }
        else if (strcmp(input_buf, "reduce_v") == 0) {
            reduce_v(nlhs, plhs, nrhs - 1, prhs + 1);
        }
        else if (strcmp(input_buf, "v_hull_extreme") == 0) {
            v_hull_extreme(nlhs, plhs, nrhs - 1, prhs + 1);
        }
        else {
            mexErrMsgTxt("Unknown function command.");
        }
        mxFree(input_buf);
    } else {
        mexErrMsgTxt("First input must be a row vector string command.");
    }
}
