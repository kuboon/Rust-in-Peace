#include <stdio.h>
#include <time.h>
#include <stdlib.h>

#define V 3     //変数の数
#define P 10000 //最大次数
#define Q 4     //基礎体
#define N Q *Q  //定義体
#define I Q + 1 //曲線の次数
#define J 3
#define K I - 2                         //number h
#define H (K + 1) * (K + 2) / 2         //シンドローム行列の横ベクトルの長さ
#define F (J - K + 1) * (J - K + 2) / 2 //シンドローム行列の縦ベクトル
#define U 26
#define E 5
#define DEG 10 + 1


unsigned char gf[16] =
    {0, 1, 2, 4, 8, 9, 11, 15, 7, 14, 5, 10, 13, 3, 6, 12};
unsigned char fg[16] =
    {0, 1, 2, 13, 3, 10, 14, 8, 4, 5, 11, 6, 15, 12, 9, 7};

typedef struct
{

  unsigned short n[V];
  unsigned short a;

} mterm;

typedef struct
{

  mterm x[P];
  int t;

} MP;

typedef struct
{

  unsigned char z[V][15000]; //零点の数

} PO;

typedef struct
{
  unsigned short m[DEG * DEG][DEG * DEG]; // 多項式の二次元配列型
} mvx;

MP v2m(mvx x) //配列型から多項式型への変換
{
  int i, j, count = 0;
  MP o = {0};

  for (i = 0; i < DEG*DEG; i++)
  {
    for (j = 0; j < DEG*DEG; j++)
    {
      if (x.m[i][j] > 0)
      {
        o.x[count].n[0] = i;
        o.x[count].n[1] = j;
        o.x[count].a = x.m[i][j];
        count++;
      }
    }
  }
  o.t = count;

  return o;
}

mvx m2v(MP x) //多項式から配列型への変換
{
  mvx z = {0};
  int i;

  for (i = 0; i < x.t; i++)
  {
    if (x.x[i].a > 0)
      z.m[x.x[i].n[0]][x.x[i].n[1]] = x.x[i].a;
  }

  return z;
}

void msm(mvx *x)
{
  int i, j;

  for (i = 0; i < DEG; i++)
  {
    for (j = 0; j < DEG; j++)
      x->m[i][j] = rand() % 2;
  }
}

mvx madd(mvx a, mvx b)
{
  int i, j;

  for (i = 0; i < DEG; i++)
  {
    for (j = 0; j < DEG; j++)
      a.m[i][j] ^= b.m[i][j];
  }

  return a;
}

void printv(MP x)
{
  int i;

  for (i = 0; i < x.t; i++)
  {
    if (x.x[i].a > 0)
      printf("%d*x^%d*y^%d+", x.x[i].a, x.x[i].n[0], x.x[i].n[1]);
  }
}

void printm(mvx m)
{
  int i, j;

  for (i = 0; i < DEG*DEG; i++)
  {
    for (j = 0; j < DEG*DEG; j++)
      if (m.m[i][j] > 0)
        printf("%d*x^%d*y^%d+", m.m[i][j], i, j);
  }

  printf(" ordering by lex\n");
}


int mlt(int x, int y)
{

  if (x == 0 || y == 0)
    return 0;

  return ((x + y - 2) % (N - 1)) + 1;
}

int mltn(int n, int x)
{
  int i, j;

  if (n == 0)
  {
    return 1;
  }
  if (x == 0)
  {
    return 0;
  }
  else
  {
    i = x;
    for (j = 0; j < n - 1; j++)
      i = mlt(i, x);

    return i;
  }
}

mvx mmul(mvx x, mvx y)
{
  int i, j, k, l;
  mvx c = {0};

  for (l = 0; l < DEG; l++)
  {
    for (i = 0; i < DEG; i++)
    {
      for (j = 0; j < DEG; j++)
      {
        for (k = 0; k < DEG; k++)
        {
          if (x.m[j][k] > 0 && y.m[i][l] > 0)
            c.m[i + j][k + l] ^= gf[mlt(fg[x.m[j][k]], fg[y.m[i][l]])];
        }
      }
    }
  }

  return c;
}


int main(void){
time_t pp;

 srand(clock() + time(&pp));

  mvx mt = {0};
  mvx xc = {0};
  mvx cx = {0};
  MP za = {0};
  mvx tm = {0};

  msm(&mt);
  printm(mt);
  printf("\n");
  msm(&xc);
  printm(xc);
  printf("\n");
  cx = mmul(xc, mt);
  printm(cx);
  printf("\n");

  za = v2m(cx);
  printv(za);
  printf(" t=%d\n",za.t);

  tm = m2v(za);
  printm(tm);
  printf("\n");
  exit(1);

return 0;
}