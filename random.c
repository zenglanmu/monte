/*
 *  Written by John Viega (viega@list.org)
 *  Free for any use.
 */

#include <stdio.h>
#include "random.h"

int secure_inited = 0;
int strong_inited = 0;
FILE *rand_file;
FILE *urand_file;

void __init_secure()
{
  rand_file  = fopen("/dev/random", "r");
  secure_inited = 1;
}

void __init_strong()
{
  urand_file = fopen("/dev/urandom", "r");
  strong_inited = 1;
}

/*
 * Return a secure random number from 0 to 2^32-1
 */
inline unsigned int secure_rand_base()
{
  unsigned int num = 0;
  int read_ret = 0;
  if(!secure_inited) __init_secure();
  while(!read_ret) read_ret = fread(&num, sizeof(unsigned int), 1, rand_file);
  return num;
}

/*
 * Return a cryptographically strong random number from 0 to 2^32-1
 */
inline unsigned int strong_rand_base()
{
  unsigned int num = 0;
  int read_ret = 0;
  if(!strong_inited) __init_strong();
  while(!read_ret) read_ret = fread(&num, sizeof(unsigned int), 1, urand_file);
  return num;
}

/*
 * Return a secure random number from 0 to 1.
 */
double secure_rand_real()
{
 double res;
   res = secure_rand_base() / (double)0xffffffff;
   return res;
}

/*
 * Return a cryptographically strong random number from 0 to 1.
 */
double strong_rand_real()
{
 double res;
  res = strong_rand_base() / (double)0xffffffff;
  return res;
}

/*
 * Return a secure random number between 0 and x-1.
 */
unsigned int secure_rand_int(unsigned int x)
{
  /* The % x is for the almost impossible situation when
   * the generated float is 1 exactly.
   */
  return ((unsigned int)(x * secure_rand_real())) % x;
}

/*
 * Return a cryptographically strong random number between 0 and x-1.
 */
unsigned int strong_rand_int(unsigned int x)
{
  /* The % x is for the almost impossible situation when
   * the generated float is 1 exactly.
   */
  return ((unsigned int) (x * strong_rand_real())) % x;
}

double rand()
{
	return strong_rand_real();
}
	
