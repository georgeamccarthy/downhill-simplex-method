/*****************************************************
 * downhillsimplex.c (C99)                           *
 *                                                   *
 * Candidate number: 22853                           *           
 *                                                   *
 * Finds the minimum of a particular Rosenbrock      *
 * function using the Nelder-Mead downhill simpex    *
 * method.                                           *
 *                                                   *
 * Setup at starting points:                         *
 * p0 = (0,0), p1 = (2,0), p2 = (0,2)                *
 *                                                   *
 * Compatible with other single valued functions and *
 * functions of more than 2 variables.               *
 *                                                   *
 *****************************************************/

#include <stdio.h>
// For pow & sqrt
#include <math.h>
// For boolean data types
#include <stdbool.h>
// For memcpy
#include <string.h>

// For dimensions other than two, the function and starting 
// coordinates must be modified.
#define N 2

// Contains all information relating to the N dimensional simplex.
typedef struct
{
    // Contains the coordinates of each point of the simplex.
    double p[N+1][N];
    // Stores the function evaluation at each point. 
    double y[N+1];

    // Indices of highest and lowest points.
    int h,l;

    // Position of the centroid.
    double pbar[N];

    // Position and y value of the first trial point.
    double pstar[N];
    double ystar;

    // Position and y value of the second trial point.
    double pstarstar[N];
    double ystarstar;

} simplex;

// Function declarations
double F(double p[N]);
void init_simplex(simplex *s);
void print_points(simplex *s);
void calc_highest_point(simplex *s);
void calc_lowest_point(simplex *s);
void calc_centroid(simplex *s);
void reflect(simplex *s);
void expand(simplex *s);
void contract(simplex *s);
void replace(simplex *s);
void update_y_values(simplex *s);
double standard_deviation(simplex *s);

int main(void)
{    
    simplex s;
    init_simplex(&s);
    
    // Iterations of the algorithm performed.
    int iterations = 0;
    // The algorithm will finish when it is done.
    bool done = false;
    while(!done)
    {
        // Updates the highest and lowest point.
        calc_highest_point(&s);
        calc_lowest_point(&s);
        // Calculates the centroid.
        calc_centroid(&s);
        // Performs a reflection.
        reflect(&s);

        if (s.ystar < s.y[s.l])
        {
            // Performs an expansion.
            expand(&s);

            if (s.ystarstar < s.y[s.l])
            {
                // Replaces ph with p**
                memcpy(s.p[s.h], s.pstarstar, sizeof(s.p[s.h]));
                update_y_values(&s);
            }
            else
            {
                // Replaces ph with p*
                memcpy(s.p[s.h], s.pstar, sizeof(s.p[s.h]));
                update_y_values(&s);
            }
        }
        else
        {
            // The result is true only if y* is greater than all y values
            // excluding the highest point.
            bool result = true;
            for (int i = 0; i <= N; i++)
            {
                if ((s.ystar > s.y[i] || s.h == i) && result != false)
                {
                    result = true;
                }
                else
                {
                    result = false;
                }
            }

            if (result)
            {
                if (!(s.ystar > s.y[s.h]))
                {
                    // Replaces ph with p*
                    memcpy(s.p[s.h], s.pstar, sizeof(s.p[s.h]));
                    update_y_values(&s);
                }

                // Performs a contraction.
                contract(&s);

                if (s.ystarstar > s.y[s.h])
                {
                    // Performs a replacement.
                    replace(&s);
                }
                else
                {
                    // Replaces ph with p**
                    memcpy(s.p[s.h], s.pstarstar, sizeof(s.p[s.h]));
                    update_y_values(&s);
                }
            }
            else
            {
                // Replaces ph with p*
                memcpy(s.p[s.h], s.pstar, sizeof(s.p[s.h]));
                update_y_values(&s);
            }
        }

        // Algorithm is complete if standard deviation meets the criterion.
        if (standard_deviation(&s) < pow(10, -8))
        {
            done = true;
            printf("Minimum reached.\n");
        }
        // Algorithm terminates if iterations surpass 1000.
        else if (iterations >= 1000)
        {
            done = true;
            printf("Max iterations (1000) reached.\n");
        }
        iterations = iterations + 1;
    }

    print_points(&s);
    printf("Standard deviation: %e\n", standard_deviation(&s));
    printf("Iterations taken: %d\n", iterations);

    return 0;
}

// Function to be minimised, returns function value at point.
double F(double p[N])
{
    // Requires modification to be compatible with higher than 2 dimensions
    return 100*pow(p[1] - pow(p[0], 2), 2) + pow((1 - p[0]), 2);
}

// Plots F between (x0 = -2 -> 2) with constant (x1 = 1).
/*void plot(void)
{
    int i;
    double p[N] = { 0,1 };
    double (*F_ptr)(double p[N]) = &F;

    FILE *fp;
    fp = fopen("output.txt", "w");
    
    for (i = 0; i <= 100; i++)
    {
        p[0] = -2 + (4.0/100) * (double)i;
        fprintf(fp, "%f,%f\n", p[0], (*F_ptr)(p));
    }

    fclose(fp);
}*/

// Initializes the simplex with starting points.
// Requires modification for larger than 2 dimensions.
void init_simplex(simplex *s) {
    s->p[0][0] = 0;
    s->p[0][1] = 0;
    s->y[0] = F(s->p[0]);

    s->p[1][0] = 2;
    s->p[1][1] = 0;
    s->y[1] = F(s->p[1]);

    s->p[2][0] = 0;
    s->p[2][1] = 2;
    s->y[2] = F(s->p[2]);
}

// Prints all of the points and their function evaluation.
void print_points(simplex *s)
{
    for (int i = 0; i <= N; i++)
    {
        printf("p%d F(", i);
        for (int j = 0; j <= N-2; j++)
        {
            // Allows for compatibility with N dimensions requiring
            // N coordinates for N+1 points.

            printf("%f, ", s->p[i][j]);
        }
        printf("%f", s->p[i][N-1]);
        printf(") = %e\n", F(s->p[i]));
    }
}

// Calculates the highest point.
void calc_highest_point(simplex *s)
{
    // Reset highest point.
    s->h = 0;
    for (int i = 0; i <= N; i++)
    {
        if (s->y[i] > s->y[s->h])
        {
            // Swap out highest point with next largest.
            s->h = i;
        }
    }
}

// Calculates the lowest point.
void calc_lowest_point(simplex *s)
{
    // Reset lowest point.
    s->l = 0;
    for (int i = 0; i <= N; i++)
    {
        if (s->y[i] < s->y[s->l])
        {
            // Swap out lowest point with next smallest.
            s->l = i;
        }
    }
}

// Finds the centroid excluding the highest point.
void calc_centroid(simplex *s)
{
    // Reset pbar.
    for (int i = 0; i <= N-1; i++)
    {
        s->pbar[i] = 0;
    }

    for (int i = 0; i <= N; i++)
    {
        // The highest point is excluded for calculation of the centroid
        if (i != s->h)
        {
            for (int j = 0; j <= N-1; j++)
            {
                s->pbar[j] = s->pbar[j] + (s->p[i][j] * 0.5);
            }
        }
    }
}

// Performs reflection transformation.
void reflect(simplex *s)
{
    for (int i = 0; i <= N-1; i++)
    {
        s->pstar[i] = (2 * s->pbar[i]) - s->p[s->h][i];
    }
    s->ystar = F(s->pstar);
}

// Performs expansion transformation.
void expand(simplex *s)
{
    for (int i = 0; i <= N-1; i++)
    {
        s->pstarstar[i] = (2 * s->pstar[i]) - s->pbar[i];
    }
    s->ystarstar = F(s->pstarstar);
}

// Performs contraction transformation.
void contract(simplex *s)
{
    for (int i = 0; i <= N-1; i++)
    {
        s->pstarstar[i] = (s->p[s->h][i] + s->pbar[i]) * 0.5;
    }
    s->ystarstar = F(s->pstarstar);
}

// Replaces all points with a point between it and the lowest point.
void replace(simplex *s)
{
    for (int i = 0; i <= N; i++)
    {
        if (i != s->h) // the highest point is excluded for calculation of the centroid
        {
            for (int j = 0; j <= N-1; j++)
            {
                s->p[i][j] = (s->p[i][j] + s->p[s->l][j])*0.5;
            }
        }
    }
    update_y_values(s);
}

// Updates y values when points are altered.
void update_y_values(simplex *s)
{
    for (int i = 0; i <= N; i++)
    {
        s->y[i] = F(s->p[i]);
    }
}

// Returns the standard deviation of the y-values.
double standard_deviation(simplex *s)
{
    double sd = 0;
    for (int i = 0; i <= N; i++)
    {
        // Calculate variance.
        sd = sd + pow(s->y[i] - F(s->pbar), 2)*((double)1/N);
    }
    // Return standard deviation as sqrt of variance.
    return sqrt(sd);
} 
