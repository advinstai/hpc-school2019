// a simple program to compute macheps for the machine you are working on
// S.C. 
// taken from wikipedia

 #include <stdio.h>
 
 int main( int argc, char **argv )
 {
    float machEps = 1.0f;
 
    printf( "current Epsilon, 1 + current Epsilon\n" );
    do {
       printf( "%G\t%.20f\n", machEps, (1.0f + machEps) );
       machEps /= 2.0f;
       // If next epsilon yields 1, then break, because current
       // epsilon is the machine epsilon.
    }
    while ((float)(1.0 + (machEps/2.0)) != 1.0);
 
    printf( "\nCalculated Machine epsilon: %G\n", machEps );
    return 0;
 }
