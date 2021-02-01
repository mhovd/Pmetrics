# include <stdlib.h>
# include <stdio.h>
# include <time.h>
# include <math.h>
# include "c_utils.h"
/******************************************************************************/
void r8_abs ( int n, double *x )
/******************************************************************************/
/*
  Purpose:
    R8_ABS returns the absolute value of an R8.
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    14 November 2006
  Author:
    John Burkardt
  Parameters:
    Input, double X, the quantity whose absolute value is desired.
    Output, double R8_ABS, the absolute value of X.
  */
{
  int iii;
  double value;
/* Code for a pointer to double; w/ dimension ( n ) in Fortran main
  */
  for ( iii = 0; iii < n; iii++, x++ )
  {
    if ( 0.0 <= *x )
    {
      value = *x ;
    }
    else
    {
      value = -1.0 *  ( *x ) ;
    }
    *x = sqrt(value) ;
    /* fprintf(stderr, "%3.17g ", value) ; */
  }
/* Code for a pointer to double; w/ dimension ( 1 ) in Fortran main
  if ( 0.0 <= *x )
  {
    value = *x;
  }
  else
  {
    value = -1.0 *  ( *x );
  }
  *x = value;
  */
  return;
}
/******************************************************************************/
void timestamp ( void )
/******************************************************************************/
/*
  Purpose:
    TIMESTAMP prints the current YMDHMS date as a time stamp.
  Example:
    31 May 2001 09:45:54 AM
  Licensing:
    This code is distributed under the GNU LGPL license.
  Modified:
    24 September 2003
  Author:
    John Burkardt
  Parameters:
    None
  */
{
# define TIME_SIZE 40
  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;
  now = time ( NULL );
  tm = localtime ( &now );
  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );
  printf ( "%s\n", time_buffer );
  return;
# undef TIME_SIZE
}
