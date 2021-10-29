#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include <stdbool.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <assert.h>
#include <execinfo.h>

#include "debug.c"
#include "global.h"
#include "struct.h"
#include "chash.c"


//Goppa多項式
static unsigned short g[K + 1] = {1, 0, 0, 0, 1, 0, 1};

//aに何をかけたらbになるか
unsigned short
equ(unsigned short a, unsigned short b)
{
  int i;

  for (i = 0; i < N; i++)
  {
    if (gf[mlt(fg[a], fg[i])] == b)
      break;
  }
  return i;
}

//OP型からベクトル型への変換
vec o2v(OP f)
{
  vec a = {0};
  int i;

  //#pragma omp parallel for
  for (i = 0; i < DEG; i++)
  {
    if (f.t[i].a > 0 && f.t[i].n < DEG)
      a.x[f.t[i].n] = f.t[i].a;
  }

  return a;
}

//ベクトル型からOP型への変換
OP v2o(vec a)
{
  int i, j = 0;
  OP f = {0};

  //#pragma omp parallel for
  for (i = 0; i < DEG; i++)
  {
    if (a.x[i] > 0)
    {
      f.t[j].n = i;
      f.t[j++].a = a.x[i];
    }
  }

  return f;
}

//停止コマンド
void wait(void)
{
  int a;                                     // 読み込む変数はローカルに取るべし
  printf(" (enter number and hit return) "); // 何か表示させたほうが良いだろう
  fflush(stdout);                            // just in case
  scanf("%d", &a);                           // fgets(line, LINESIZE, stdin); という手も
}

//OP型を正規化する
OP conv(OP f)
{
  vec v = {0};
  OP g = {0};

  v = o2v(f);
  g = v2o(v);

  return g;
}

void op_print_raw(const OP f)
{
  puts("op_print_raw:");
  for (int i = 0; i < DEG; i++)
  {
    if (f.t[i].a > 0)
      printf("[%d] %ux^%u\n", i, f.t[i].a, f.t[i].n);
  }
}

bool op_verify(const OP f)
{
  bool end = false;
  unsigned short n_max = 0;
  for (int i = 0; i < DEG; i++)
  {
    if (end && (f.t[i].n != 0 || f.t[i].a != 0))
    {
      op_print_raw(f);
      printf("found data after end: i=%d\n", i);
      print_trace();
      fflush(stdout);
      return false;
    }
    if (f.t[i].a == 0)
    {
      end = true;
      continue;
    }
    if (f.t[i].n + 1 <= n_max)
    {
      op_print_raw(f);
      printf("found invalid order: i=%d\n", i);
      print_trace();
      fflush(stdout);
      return false;
    }
    n_max = f.t[i].n + 1;
  }
  return true;
}



//20200816:正規化したいところだがうまく行かない
//多項式の足し算
OP oadd(OP f, OP g)
{
  f = conv(f);
  g = conv(g);
  assert(op_verify(f));
  assert(op_verify(g));

  vec a = {0}, b = {0}, c = {0};
  int i, j, k, l = 0;
  OP h = {0}, f2 = {0}, g2 = {0};

  a = o2v(f);
  b = o2v(g);

  //k=deg(o2v(f));
  //l=deg(o2v(g));

  for (i = 0; i < DEG; i++)
  {
    c.x[i] = a.x[i] ^ b.x[i];
  }
  h = v2o(c);
  //h=conv(h);
  assert(op_verify(h));
  return h;
}

//多項式の次数(degのOP型)
int odeg(OP f)
{
  int i, j = 0, k;

  if (f.t[0].a == 0)
    return 0;

  //k=terms(f);
  for (i = 0; i < DEG; i++)
  {
    if (j < f.t[i].n && f.t[i].a > 0)
      j = f.t[i].n;
  }

  return j;
}

//多項式の次数(default)
int deg(vec a)
{
  int i, n = 0, flg = 0;

  //#pragma omp parallel for
  for (i = 0; i < DEG; i++)
  {
    if (a.x[i] > 0)
    {
      n = i;
      flg = 1;
    }
  }
  if (flg == 0)
    return 0;

  return n;
}

//リーディグタームを抽出(default)
oterm LT(OP f)
{
  int i, k;
  oterm t = {0};

  //k = deg (o2v (f));
  for (i = 0; i < DEG; i++)
  {
    //printf("a=%d %d\n",f.t[i].a,f.t[i].n);
    if (f.t[i].a > 0)
    {
      t.n = f.t[i].n;
      t.a = f.t[i].a;
    }
  }

  return t;
}


//配列からベクトル表現の多項式へ変換する
vec Setvec(int n)
{
  int i;
  vec v = {0};

  for (i = 0; i < n; i++)
  {
    v.x[n - 1 - i] = c[i];
  }

  return v;
}


//配列の値を係数として多項式に設定する
OP setpol(unsigned short f[], int n)
{
  OP g;
  vec a;
  int i;

  memset(c, 0, sizeof(c));
  memcpy(c, f, 2 * n);
  a = Setvec(n);

  g = v2o(a);

  return g;
}
//ランダム多項式の生成
static void
ginit(void)
{
    int j, count = 0, k = 0;
    unsigned short gg[K + 1] = {0};

    printf("in ginit\n");

    g[K] = 1;          //xor128();
    g[0] = rand() % 2; //or N
    k = rand() % (K - 1);
    if (k > 0)
    {
        while (count < k)
        {
            printf("in whule\n");
            j = rand() % (K);
            if (j < K && j > 0 && g[j] == 0)
            {
                g[j] = rand() % 2; //or N;
                count++;
            }
        }
    }

    for (j = 0; j < K + 1; j++)
        gg[j] = g[K - j];

    memcpy(g, gg, sizeof(g));
}



//有限体の元の逆数
unsigned short
oinv(unsigned short a)
{
  int i;

  if (a == 0)
    return -1;

  for (i = 0; i < N; i++)
  {
    if (gf[mlt(fg[a], i)] == 1)
      return (unsigned short)i;
  }

  printf("no return \n");
  //  exit (1);
}

unsigned short
v2a(oterm a)
{
  int i, j;

  if (a.a == 0)
    return 0;

  //printf("aa=%d\n",a.a);
  for (j = 0; j < N; j++)
  {
    if (gf[j] == a.a && a.a > 0)
    {
      //printf("j==%d\n",j);
      return j - 1;
    }
  }
}


void printsage(vec a)
{
  int i, j, k;
  oterm b;

  printf("poly=");
  for (i = 0; i < DEG; i++)
  {
    if (a.x[i] > 0)
    {
      b.a = a.x[i];
      b.n = i;
      j = v2a(b);
      //printf("%d,==ba\n",b.a);
      //printf ("X**%d+", i); //for GF2
      printf("B('a^%d')*X**%d+", j, i); //for GF(2^m)
    }
  }
}


//多項式を項ずつ掛ける
OP oterml(OP f, oterm t)
{
  f = conv(f);
  assert(op_verify(f));
  int i, k, j;
  OP h = {0};
  vec test;
  unsigned int n;

  //f=conv(f);
  //k = deg (o2v(f));
  j = 0;
  for (i = 0; i < DEG; i++)
  {
    h.t[i].n = f.t[i].n + t.n;
    h.t[i].a = gf[mlt(fg[f.t[i].a], fg[t.a])];
  }

  h = conv(h);
  assert(op_verify(h));
  return h;
}

//多項式を単行式で割る
oterm LTdiv(OP f, oterm t)
{
  oterm tt = {0}, s = {
                      0};

  tt = LT(f);
  if (tt.n < t.n)
  {
    s.n = 0;
    s.a = 0;
  }
  else if (tt.n == t.n)
  {
    s.n = 0;
    s.a = equ(t.a, tt.a);
  }
  else if (tt.n > t.n)
  {
    s.n = tt.n - t.n;
    s.a = equ(t.a, tt.a);
    //printf("%u\n",s.a);
  }
  else if (t.n == 0 && t.a > 0)
  {
    s.a = gf[mlt(fg[tt.a], oinv(t.a))];
    s.n = tt.n;
  }

  return s;
}



//多項式の掛け算
OP omul(OP f, OP g)
{
  f = conv(f);
  g = conv(g);
  assert(op_verify(f));
  assert(op_verify(g));
  int i, count = 0, k, l;
  oterm t = {0};
  OP h = {0}, e = {0}, r = {0};
  vec c = {0};

  k = odeg(f);
  l = odeg(g);
  if (l > k)
  {
    k = l;
  }

  for (i = 0; i < k + 1; i++)
  {
    t = g.t[i];
    e = oterml(f, t);
    h = oadd(h, e);
  }
  assert(op_verify(h));
  return h;
}

//多項式の剰余を取る
OP omod(OP f, OP g)
{
  int i = 0, j, n, k;
  OP h = {0}, e = {
                  0};
  oterm a, b = {0}, c = {0};

  //n = LT(g).n;

  //  assert (("baka^\n", LT (f).n != 0));

  //  assert (("baka(A)\n", LT (g).n != 0));

  if (LT(f).n < LT(g).n)
  {
    //    exit(1);
    return f;
  }

  //printf ("in omod\n");
  //exit(1);

  k = LT(g).n;
  b = LT(g);
  //printpol(o2v(g));
  //printf("d %d\n",LT(g).n);

  //assert(("double baka\n", b.a != 0 && b.n != 0));
  while (LT(f).n > 0 && LT(g).n > 0)
  {

    c = LTdiv(f, b);
    h = oterml(g, c);
    f = oadd(f, h);
    if (odeg((f)) == 0 || odeg((g)) == 0)
    {
      //      printf("blake1\n");
      break;
    }

    if (c.n == 0 || b.n == 0)
      break;
  }

  return f;
}

//gcd
OP gcd(OP xx, OP yy)
{
  OP tt = {0}, tmp,h={0};

  h.t[0].a=1;
  h.t[0].n=0;
  if (deg(o2v(xx)) < deg(o2v(yy)))
  {
    tmp = xx;
    xx = yy;
    yy = tmp;
  }
  tt = omod(xx, yy);
  while (odeg(tt) > 0)
  {
    xx = yy;
    yy = tt;
    if(odeg(yy)>0)
    tt = omod(xx, yy);
    if(LT(tt).a==0)
    return yy;
  }
  if(LT(tt).a==0)
  {
  return yy;
  }else{
    return h;
  }
//  return yy;
}



//多項式のべき乗余
OP opowmod(OP f, OP mod, int n)
{
  int i, j = 0,l;


  //printpol(o2v(mod));
  //printf(" =mod %d\n",LT(mod).n);
  //繰り返し２乗法
  for (i = 1; i <  n+1; i++)
  {
    f = omul(f, f);
    if ((deg(o2v(f)) > deg(o2v(mod))) && odeg(mod)>0)
      f = omod(f, mod);
  }

  return f;
}



//GF(2^m) then set m in this function.
int ben_or(OP f)
{
  int i, n, flg = 0;
  OP s = {0}, u = {0}, r = {0}; 
  vec v = {0}, x = {0};
  //if GF(8192) is 2^m and m==13 or if GF(4096) and m==12 if GF(16384) is testing
  int m = E;
  // m=12 as a for GF(4096)=2^12 defined @ gloal.h or here,for example m=4 and GF(16)


  v.x[1] = 1;
  s = v2o(v);
  r = s;
  n = deg(o2v(f));

if(LT(f).n==0 && LT(f).a==1)
{
  printf("f==0\n");
  exit(1);
}
  if (n == 0)
    return -1;

  i = 0;
  
  //r(x)^{q^i} square pow mod
  while (i < K / 2 )
  {
    flg = 1;
    

    // irreducible over GH(8192) 2^13
    r = opowmod(r, f, m);

    // irreducible over GF2
    //r=omod(opow(r,2),f);

    u = oadd(r, s);
    //if (deg(o2v(u)) > 0)
      u = gcd(f, u);

    if (LT(u).n > 0)
      return -1;


      i++;
  }

  return 0;
}

//多項式の代入値
unsigned short
trace(OP f, unsigned short x)
{
  int i, d;
  unsigned short u = 0;

  d = deg(o2v(f));

  for (i = 0; i < d + 1; i++)
  {
    u ^= gf[mlt(fg[f.t[i].a], mltn(f.t[i].n, fg[x]))];
  }

  return u;
}


//多項式を表示する(default)
void printpol(vec a)
{
  int i, n;

  n = deg(a);

  //printf ("baka\n");
  assert(("baka\n", n >= 0));

  for (i = n; i > -1; i--)
  {
    if (a.x[i] > 0)
    {
      printf("%u", a.x[i]);
      //if (i > 0)
      printf("x^%d", i);
      //if (i > 0)
      printf("+");
    }
  }
  //  printf("\n");

  return;
}

OP mkpol()
{
    int i, j, k, fail, flg, l, ii = 0;
    OP w = {0};

    do
    {
        fail = 0;
        j = 0;
        k = 0;
        flg = 0;
        l = 0;
        memset(g, 0, sizeof(g));
        //memset(ta, 0, sizeof(ta));
        memset(w.t, 0, sizeof(w));
        ginit();
        ii++;
        if (ii > K/2)
        {
            printf("erro=%d\n", ii);
            exit(1);
        }

        for (i = 0; i < K; i++)
        {
            if (g[K - 1] > 0)
                flg = 1;
            if (i % 2 == 1 && g[i] > 0 && i < K)
                k++;
        }

        //偶数項だけにならないようにする
        if ((k > 0 && flg == 0) || (k > 1 && flg == 1))
        //if(k>0)
        {
            w = setpol(g, K + 1);
            j = 1;
            //if(isquad(w)==-1)
            //exit(1);
        }
        // exit(1);

    } while (j == 0);

    printpol(o2v(w));
    printf(" ==g\n");
    //exit(1);

    return w;
}




int main(void)
{
  int i, j, k, l, ii = 0;
  OP w = {0};
  unsigned short tr[N] = {0};
  unsigned short ta[N] = {0};

  j=0;
aa:

  //printf("\n");

  //既約性判定のためのBen-Orアルゴリズム。拡大体にも対応している。デフォルトでGF(8192)
  //既約多項式しか使わない。
  l = -1;
  ii = 0;


  while (l == -1)
  {
    w = mkpol();
    l = ben_or(w);
    printf("irr=%d\n", l);
    if (ii > 300)
    {
      printf("too many tryal\n");
      goto aa;
      //exit(1);
    }
    printf("ben=%d\n",ii);
    ii++;
    //
  } 
  
  //w = mkpol();

  //多項式の値が0でないことを確認
  for (i = 0; i < N; i++)
  {
    ta[i] = trace(w, i);
    if (ta[i] == 0)
    {
      printf("trace 0 @ %d\n", i);
      //fail = 1;
      exit(1);
    }
  }
  for (i = 0; i < N; i++)
  {
    tr[i] = oinv(ta[i]);
    //printf("%d,", tr[i]);
  }

  printpol(o2v(w));
  printf(" =irreducible\n");
  printsage(o2v(w));
  printf("\n");
  j++;

  if(j>100)
  exit(1);
  goto aa;
}
