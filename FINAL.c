#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#define ITER 10000000000000000
#define MAX_LINE 100000
#define MAX_VARIABLES 5000


double logit(double z)
{   double e;
    if (z >= 0.0)
    {
        e = exp(-z);
        return 1.0 / (1.0 + e);
    }
    else
    {
        e = exp(z);
        return e / (1.0 + e);
    }
}


int Gauss_Elimination(double **A, double *b, double *x, int m)
{
    int i, j, k, pivot;
    double max,a,z;
    double f;
  // FINDING THE PIVOT

    for (i = 0; i < m; i++)
    {
        pivot = i;
        max = fabs(A[i][i]);
        for (j = i + 1; j < m; j++)
        {
            if (fabs(A[j][i]) > max)
            {
                max = fabs(A[j][i]);
                pivot = j;
            }

        }

    // CHECKING WHETHER MATRIX IS SINGULAR OR NOT
        if (max < 0.0000000001)
        {
            printf("Error: Coefficient Matrix not Compatible, Check whether the dataset is Completely Separable\n");
            return 0;
        }
        // Swapping the rows of Coefficent matrix and RHS vector
        if (pivot != i)
        {
            for (j = i; j < m; j++)
            {
                a = A[i][j];
                A[i][j] = A[pivot][j];
                A[pivot][j] = a;
            }

            z = b[i];
            b[i] = b[pivot];
            b[pivot] = z;
        }

    // Converting Matrix to Upper Triangular using Gaussian Elimination
        for (j = i + 1; j < m; j++)
        {
            f = A[j][i] / A[i][i];
            for (k = i; k < m; k++)
               { A[j][k] = A[j][k] - f * A[i][k];}
            b[j] = b[j]- f * b[i];
        }
    }
   // Performing Back Substituion to find x
    for (i = m - 1; i >= 0; i--)
    {
        x[i] = b[i];
        for (j = i + 1; j < m; j++)
            {x[i] = x[i]-A[i][j] * x[j];}
        x[i] = x[i]/A[i][i];
    }
return 1;
}

// READING THE FILE

int read_data(const char *file, int *n, int *k, double ***X, double **y)
{
    FILE *fp;
    char line[MAX_LINE];
    char *t;
    char *e;
    int observation = 0;
    int var = -1;
    int val;
    int i, j, z;
    double value;
    fp = fopen(file, "r");
    if (fp == NULL)
    {
        printf("Error opening file: %s\n", file);
        return 0;
    }
    while (fgets(line, MAX_LINE, fp) != NULL)    //  CHECKING FROM START OF THE FILE
    {
        if (strspn(line, " \t\r\n") == strlen(line))   // IF THERE IS ANY SPECIAL CHARACTER THEN GO TO THE NEXTT STEP
            {continue;}
        val = 0;
        t = strtok(line, " ,\t\r\n");  // TAKING VALUES BEFORE THE SPECIAL CHARACTERS
        while (t != NULL)
        {   errno = 0; 
            value = strtod(t, &e);    // CONVERT STRING TO DOUBLE
            if (e == t || *e != '\0' || errno == ERANGE || value != value || value == HUGE_VAL || value == -HUGE_VAL )
            {
                printf("Error: Invalid numeric value: %s\n", t);
                fclose(fp);
                return 0;
            }
            if (val == 0 && value != 0.0 && value != 1.0)
            {
                printf("Error: y must be binary (0 or 1).\n");
                fclose(fp);
                return 0;
            }
            val++;
            t = strtok(NULL, " ,\t\r\n");
        }
        if (val < 2)
        {
            printf("Error: Incomplete observation.\n");
            fclose(fp);
            return 0;
        }
        if (var == -1)
        {
            var = val - 1;

            if (var > MAX_VARIABLES)
            {
                printf("Error: Too many independent variables.\n"); // CHECKING CONDITION
                fclose(fp);
                return 0;
            }
        }
        else if (val - 1 != var)
        {
            printf("Error: Inconsistent number of columns.\n");   // CHECKING FOR MISSING GAPS IF ANY
            fclose(fp);
            return 0;
        }

        observation++;    
    }

    if (ferror(fp))
    {
        printf("Error while reading file.\n");
        fclose(fp);
        return 0;
    }

    if (observation == 0)
    {
        printf("Error: File contains no observations.\n");
        fclose(fp);
        return 0;
    }

    *n = observation;  // ALLOCATING THE NUMBER OF OBSERVATION
    *k = var;         // ALLOCATING THE NUMBER OF INDEPENDENT VAIRABLES

    *y = malloc((size_t)*n * sizeof(double));   // ALLOCATING A POINTER FOR STORING Y
    if (*y == NULL)
    {
        printf("Error allocating y.\n");
        fclose(fp);
        return 0;
    }
    *X = malloc((size_t)*n * sizeof(double *)); // ALLOCATING A POINTER FOR STORING X
    if (*X == NULL)
    {
        printf("Error allocating X.\n");
        free(*y);
        *y = NULL;
        fclose(fp);
        return 0;
    }
    for (i = 0; i < *n; i++)
    {
        (*X)[i] = malloc((size_t)*k * sizeof(double)); // ALLOCATING A POINTER FOR STORING VALUES IN ROWS OF X 
        if ((*X)[i] == NULL)
        {
            printf("Error allocating row of X.\n");

            for (j = 0; j < i; j++)
                free((*X)[j]);
            free(*X);
            free(*y);
            *X = NULL;
            *y = NULL;
            fclose(fp);
            return 0;
        }
    }
    rewind(fp);   // INITIALISING FP to use it here to read file agian from the start
    z = 0;
    while (fgets(line, MAX_LINE, fp) != NULL)  // second check the data before allocation
    {
        if (strspn(line, " \t\r\n") == strlen(line))   
            continue;
        if (z >= *n)    // checking for conssitency if blanks are also ignored
        {
            printf("Error: Inconsistent observations during second check.\n");

            for (i = 0; i < *n; i++)
                {free((*X)[i]);}
            free(*X);
            free(*y);
            *X = NULL;
            *y = NULL;
            fclose(fp);
            return 0;
        }
    t = strtok(line, " ,\t\r\n");
    (*y)[z] = strtod(t, NULL);
    for (j = 0; j < *k; j++)
    {
        t = strtok(NULL, " ,\t\r\n");
       (*X)[z][j] = strtod(t, NULL);
    }
    if (strtok(NULL, " ,\t\r\n") != NULL)
        {
            printf("Error: Too many columns in observation.\n");
            goto second_pass_error;
        }
    z++;
    }
    if (ferror(fp) || z != *n)
    {
        printf("Error: Could not read all observations.\n");
        goto second_pass_error;
    }
    fclose(fp);
    return 1;
    second_pass_error:    // LABEL
        for (i = 0; i < *n; i++)
            {free((*X)[i]);}
        free(*X);
        free(*y);
        *X = NULL;
        *y = NULL;
        fclose(fp);
        return 0;
}
//  LOGISTIC NEWTON RAPHSON METHOD
int logistic_newton_raphson(double **X, double *y, int n, int k, double *beta, double TOL)
{   int c = k + 1;
    int iter;
    int i,j,z,s,t,u;
    int q,e,m;
    double change;
    double w,residual;
    double o,v;
    double (*logit_ptr)(double) = logit;
    int (*gauss_ptr)(double **, double *, double *, int) =Gauss_Elimination;
    double *p = malloc(n * sizeof(double));
    double *score = malloc(c * sizeof(double));
    double *delta = malloc(c * sizeof(double));
    double **H = calloc((size_t)c, sizeof(double *));
    if (p == NULL || score == NULL ||
        delta == NULL || H == NULL)
    {
        printf("Error: Memory allocation failed in logistic regression.\n");
        free(score);
        free(delta);
        free(H);
        free(p);
        return 0; }
    for ( i = 0; i < c; i++)
    {    H[i] = malloc(c * sizeof(double));
                if (H[i] == NULL)
        {
            printf("Error: Memory allocation failed for Hessian row.\n");

            for (j = 0; j < i; j++)
                free(H[j]);

            free(H);
            free(p);
            free(score);
            free(delta);
            return 0;
        }
    }
    for (iter = 0; iter < ITER; iter++)    
    {   
        for (j = 0; j < n; j++)
        {
        double eta = beta[0];
        for (z = 0; z < k; z++)
        {
         eta = eta+ X[j][z] * beta[z + 1];
         }
            p[j] = logit_ptr(eta);  // CNA USE OINTER TO A FUCNTION HERE
        }
        for (s = 0; s < c; s++)
        {
            score[s] = 0.0;   // Initializing the score and Hessian
        }
        for (t = 0; t < c; t++)
        {
            for (u = 0; u < c; u++)
            {
                H[t][u] = 0.0;
            }
        }
for (q = 0; q < n; q++)
{
    w = p[q] * (1.0 - p[q]);
    residual = y[q] - p[q];

    for (e = 0; e < c; e++)
    {
        o = (e == 0) ? 1.0 : X[q][e - 1];
        score[e] = score[e] + o * residual;
        for (m = 0; m < c; m++)
        {
            if (m == 0)
            {v = 1.0;}
          else
            {v = X[q][m - 1];}
            H[e][m] =   H[e][m] + o * w * v;
        }
    }
}
    if (gauss_ptr(H, score, delta, c) == 0)
        { 
            printf("Newton-Raphson stopped at iteration %d.\n", iter + 1);
            for (i = 0; i < c; i++)
                {free(H[i]);}
            free(H);
            free(p);
            free(score);
            free(delta);
            return 0 ;
        
        }
        change = 0.0;
        for (j = 0; j < c; j++)
        {
            beta[j] = beta[j]  + delta[j];  //newton rahson method

            if (fabs(delta[j]) > change)
                change = fabs(delta[j]);
        }
        if (change < TOL)
        {
            printf("\nNewton-Raphson converged.\n");
            printf("Number of iterations: %d\n", iter + 1);
            break;
        }
    }
    if (iter == ITER)
    {
        printf("Error Maximum iterations reached.\n");
        for (j = 0; j < c; j++)
                {free(H[j]);}
            free(H);
            free(p);
            free(score);
            free(delta);
        return 0;
    }
for (i = 0; i < c; i++)
    {free(H[i]);}

free(H);
free(p);
free(score);
free(delta);
return 1;
}
int main()
{   
    char filename[5000];
    char outputfile[5000];
    int n, k;
    int i, j;
    double **X;
    double *y;
    double *beta;
    double TOL;
    FILE *out;
    printf("Warning Input file format: .txt file with numeric data in each row: y x1 x2 ... xk\n");
    printf("Warning The first column must be binary (0 or 1), followed by regression variables.\n");
    printf("Warning Use spaces, commas, or tabs to separate values; do not include  headers \n");
    printf("Enter data file name without spaces:");
    if (scanf("%4999s", filename) != 1)
    {
        printf("Error reading input filename.\n");
        return 0;
    }
     printf("Enter tolerance level (e.g., 0.000001): ");
    if (scanf("%lf", &TOL) != 1 || TOL <= 0.0)
    {
        printf("Error: Tolerance must be a positive, finite number.\n");
        return 0;
    }
    if (read_data(filename, &n, &k, &X, &y) != 1)
    {
        printf("Error reading data.\n");
        return 0;
    }
    printf("\nData successfully read.\n");
    printf("Number of observations = %d\n", n);
    printf("Number of independent variables = %d\n", k);
    beta = malloc((k+1) * sizeof(double));
    if (beta == NULL)
    {
        printf("Error: Memory allocation failed.\n");
        for (i = 0; i < n; i++)
            free(X[i]);
        free(X);
        free(y);
        return 0;
    }
    for (j = 0; j <= k; j++)
    {
        beta[j] = 0.0;
    }
if (!logistic_newton_raphson(X, y, n, k, beta, TOL))
   {
    printf("Error: Logistic regression failed to converge.\n");

    for (i = 0; i < n; i++)
    {  free(X[i]);   }
    free(X);
    free(y);
    free(beta);
    return 0;
   }
    printf("RESULTS\n");
    printf("Number of observations: %d\n", n);
    printf("Number of independent variables: %d\n", k);
    printf("\nEstimated coefficients:\n");
   printf("Intercept = %.10f\n", beta[0]);
   for (j = 1; j <= k; j++)
    {printf("Beta_%d = %.10f\n", j, beta[j]);}
    printf("\nEnter output file name: ");
    if (scanf("%4999s", outputfile) != 1)
    {
        printf("Error reading output filename.\n");
        for (i = 0; i < n; i++)
            free(X[i]);
        free(X);
        free(y);
        free(beta);
        return 0;
    }
    out = fopen(outputfile, "w");
    if (out == NULL)
    {
        printf("Error: Could not create output file.\n");
        for (i = 0; i < n; i++)
            free(X[i]);

        free(X);
        free(y);
        free(beta);

        return 0;
    }
    fprintf(out, "LOGISTIC REGRESSION RESULTS\n");
    fprintf(out, "Tolerance level = %.10g\n", TOL);
    fprintf(out, "Number of observations: %d\n", n);
    fprintf(out, "Number of independent variables: %d\n", k);
    fprintf(out, "\nEstimated coefficients:\n");
    fprintf(out,"Intercept = %.10f\n", beta[0]);
    for (j = 1; j <= k; j++)
    {fprintf(out,"Beta_%d = %.10f\n", j, beta[j]);}
    fclose(out);
    printf("\nResults successfully saved to %s\n", outputfile);
    for (i = 0; i < n; i++)
    {
        free(X[i]);
    }
    free(X);
    free(y);
    free(beta);
    return 1;
}

