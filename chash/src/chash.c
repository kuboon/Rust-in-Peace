#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <assert.h>


#define N 256
#define NN 64
#define NB 4  	/* 128bit 固定として規格されている(データの長さ) */
#define NBb 16
#define SEED 12345678901234567899ULL

typedef union
{
  unsigned long long int u[8];
  unsigned char d[64];
} arrayul;

typedef union
{
  unsigned long long int u[N];
  unsigned char d[N * 8];
} arrayull;


typedef union a4
{
  unsigned char ar[4];
  unsigned int n;
} array;

typedef struct a8
{
  unsigned char ar[8];
} array8;

typedef union
{
  unsigned int h[16];
  unsigned char c[64];
} array16;

typedef union aN
{
  unsigned int d[64];
  unsigned long long int u[32];
  unsigned char ar[256];
  //
} arrayn;

typedef struct pub
{
  unsigned char a[NN];
  unsigned char  b[NN];
} set;

arrayn c={0};


#define I8T char
#define U8C(v) (v##U)

#define U8V(v) ((unsigned char)(v) & U8C(0xFF))
#define ROTL8(v, n) \
  (U8V((v) << (n)) | ((v) >> (8 - (n))))


unsigned int rotate_left(unsigned int x, int n)
{
    assert(0 < n && n < 32);
    return (x << n) | (x >> (32 - n));
}

//言わずと知れたxorshift
unsigned int
xorshift (unsigned long long int u)
{
  static unsigned int y = 2463534242;
  y= u ^ y;
  y = y ^ (y << 13);
  y = y ^ (y >> 17);
  return y = y ^ (y << 15);
}

unsigned long long int
xorshift64 (unsigned long long int u)
{
  static unsigned long long int x = 88172645463325252ULL;
  x = x ^ u;
  x = x ^ (x << 13);
  x = x ^ (x >> 7);
  return x = x ^ (x << 17);
}


void rp(unsigned char* a){
 unsigned int i,j,x;


 for(i = 0; i < NN; i++)
    a[i] = i;
 
  for(i = 0; i < NN - 2; i++){
    // rand from i+1 to N-1
    j = (xorshift64(SEED) % (NN-1-i)) + i + 1;
    
    // swap a[i] and a[j]
    x = a[j];
    a[j] = a[i];
    a[i] = x;
  }
  if(a[NN-1] == NN-1){
    a[NN-1] = a[NN-2];
    a[NN-2] = NN - 1;
  }
}


//ハッシュ関数本体
arrayn
chash (unsigned char b[2048])
{
  int i, j = 0;
  arrayn n;

  unsigned char salt[NN] ={ 148, 246, 52, 251, 16, 194, 72, 150, 249, 23, 90, 107, 151, 42, 154, 124, 48, 58, 30, 24, 42, 33, 38, 10, 115, 41, 164, 16, 33, 32, 252, 143, 86, 175, 8, 132, 103, 231, 95, 190, 61, 29, 215, 75, 251, 248, 72, 48, 224, 200, 147, 93, 112, 25, 227, 223, 206, 137, 51, 88, 109, 214, 17, 172};

  unsigned char z[NN];
  unsigned char f[NN] = { 0 };
  unsigned char x0[NN]={0};  
  unsigned char inv_x[NN] = { 0 };
  unsigned char x1[NN]={0};

  
  rp(x0);
  rp(x1);
  
  for (i = 0; i < NN; i++)
      inv_x[x0[i]] = i;
    
  memset(f,0,sizeof(f));
  
  //デバッグ中なので保留
  //for (i = 0; i < NN; i++)
    //   f[i] ^= salt[i];


  /*  
   for(i=0;i<NN;i++)
     printf("%d,",f[i]);
   printf("\n\n");
  */

  //バッファを埋める回数だけ回す
  for (j = 0; j < 2048/NN; j++)
    {
      for (i = 0; i < NN; i++)
	z[i] = x0[x1[inv_x[i]]];
      
      for (i = 0; i < NN; i++)
	f[i] ^= b[j * NN + i];
      
      memcpy (x1, z, sizeof (unsigned char) * NN);
      
      for (i = 0; i < NN; i++)
	{	  
	  //mode 2(自己書き換え系)
	  f[x1[i]]+=abs(ROTL8(f[(i+1)%NN],3)-ROTL8(f[i],5));
	}
    }

  memcpy (n.ar, f, sizeof (unsigned char) * NN);
  
  
  return n;
}


//ファイル操作
array16
hash (char *filename)
{
  int i, k, n;
  array16 h = { 0 };

  unsigned char buf[2048] = { 0 };
  FILE *fp;
  arrayn a = { 0 };


  fp = fopen (filename, "rb");
  if (fp == NULL)
    {
      printf ("no file\n");
      exit (1);
    }
  
  while ((n = fread (buf, 1, 2048, fp)) > 0)
    {
      //paddaing
      if(n<2048){
	for(i=n;i<2048;i++)
	  buf[i]=0xc6;
      }
      
      a = chash (buf);
      for (k = 0; k < NN / 64; k++)
	{
	  for (i = 64/4 * k; i < 64/4 * k + 64/4; i++)
	      h.h[i - 64/4 * k] ^= a.d[i];
	}
    }
  
  
  return h;   
}


//蛇足
arrayul
crand (unsigned char u[NN])
{
  arrayn a = { 0 };
  arrayul b = { 0 };

  
  a = chash (u);
  
  memset (b.u, 0, sizeof (b.u));
  memcpy(b.d,a.ar,sizeof(unsigned char)*NN);
  
  return b;
}


int
main (int argc, char *argv[])
{
  int i;
  array16 t;
  //  time_t o;
  
  
  //  srand (clock () + time (&o));
  
  t = hash (argv[1]);
  //慎ましくここは256ビットだけ
  for (i = 0; i < 16 / 2; i++)
    printf ("%08x", t.h[i]);
  printf (" %s", argv[1]);
  printf ("\n");
  
  
  return 0;
}
