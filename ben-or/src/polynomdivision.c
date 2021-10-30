//////////////////////////////////////////////////////////////////////////
// Polynom Division in C
//
// Author:  Adriano De Rosa
// Created: December 23rd, 2002
// Updated: August 1st, 2020: conio.h removed for usage on linux systems. Dynamic array added
// Updated: December  13th, 2020: Time measurment added with time.h librarx      
//
//////////////////////////////////////////////////////////////////////////

#include<stdio.h>   // This is required for printf, scanf
#include<stdlib.h>  // This is required for malloc, free, realloc
#include<time.h>    // This is required for malloc, free, realloc


////////////////////////////////////////////////////////////////////
void main()
{
     
 int* poly= malloc(sizeof(int[10]));  // Dynamisches array für polynom koeffizienten
 int* q   = malloc(sizeof(int[10]));  // 

 //int poly[5], q[5];                 // so würde es statisch aussehen, z.B. für olynom mit Grad 5 .... nicht so flexibel

 int n, r, i;

 printf("\t POLYNOMDIVISION");
 printf("\n Please enter the degree of the equation: ");
 
 scanf("%d",&n);  // Einabe 

 // Arraygrößen anpassen auf Polynomgrad
 poly = realloc(poly, sizeof(int[n])); 
 q    = realloc(q, sizeof(int[n]));    

 for(i=0;i<=n;i++)
	{
 	 printf("\n Input the coefficient x[%d] = ", n-i);
 	 scanf("%d",&poly[i]);
	}

 printf("\n Enter the value of constant in (x-r) : ");
 scanf("%d",&r);

 //////////////////////////////////////////////////////////////////////////////
 clock_t begin = clock(); // Start time measurement from here on (count clocks)
 //////////////////////////////////////////////////////////////////////////////

 q[0] = poly[0];
 for(i=1;i<=n;i++)
	{
	 q[i] = (q[i-1]*r)+poly[i];
	}

 printf("\n The quotient coefficients are: \n");

 for(i=0;i<n;i++)
	{
	 printf("\t%d",q[i]);
	}

 printf("\n The remainder is: %d", q[n]);
 printf("\n");
 
 free (poly); // free the memory ...
 free (q);

 /////////////////////////////////////////////////////////////
 clock_t end = clock(); // Save the execution clocks */
 double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
 printf("Execution time: %f ms\n",time_spent*1000);
 //////////////////////////////////////////////////////////////

}


