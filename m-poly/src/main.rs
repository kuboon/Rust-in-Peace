use std::{i16, i32};
extern crate ndarray;
//use rand::prelude::*;
use std::mem::MaybeUninit;

pub const GF: [i8; 16] = [0, 1, 2, 4, 8, 9, 11, 15, 7, 14, 5, 10, 13, 3, 6, 12];
pub const FG: [i8; 16] = [0, 1, 2, 13, 3, 10, 14, 8, 4, 5, 11, 6, 15, 12, 9, 7];

#[derive(Copy, Clone)]
struct MTERM {
    n: [i16; 3],
    a: i32,
}

#[derive(Copy, Clone)]
struct MP {
    x: [MTERM; 1000],
    t: i32,
}

#[derive(Clone)]
struct PO {
    z: Vec<Vec<u8>>,
}

#[derive(Copy, Clone)]
struct MVX {
    m: [[i16; 100]; 100],
}

fn v2m(x: MVX) -> MP //配列型から多項式型への変換
{
    let mut count: i32;
    let mut o: MP = unsafe { MaybeUninit::uninit().assume_init() };
    count = 0;
    for it in 0..20 {
        for be in 0..20 {
            if x.m[it][be] > 0 {
                o.x[count as usize].n[0] = it as i16;
                o.x[count as usize].n[1] = be as i16;
                o.x[count as usize].a = x.m[it][be] as i32;
                count = count + 1;
            }
        }
    }
    o.t = count;

    return o;
}

fn m2v(x: MP) -> MVX //多項式から配列型への変換
{
    let mut z: MVX = unsafe { MaybeUninit::uninit().assume_init() };

    for i in 0..x.t {
        if x.x[i as usize].a > 0 {
            z.m[x.x[i as usize].n[0] as usize][x.x[i as usize].n[1] as usize] =
                x.x[i as usize].a as i16;
        }
    }

    return z;
}

fn msm(mut x: MVX) -> MVX {
    for me in 0..10 {
        for inn in 0..10 {
            x.m[me][inn] = rand::random::<i16>() % 2;
        }
    }
    return x;
}

/*
fn madd(mut a:MVX, b:MVX) -> MVX
{

  for i in 0..100
  {
    for j in 0..100
    {
      a.m[i][j] ^= b.m[i][j];
    }
  }

  return a;
}
*/

fn printv(x: MP) {
    for i in 0..x.t {
        if x.x[i as usize].a > 0 {
            print!(
                "{}*x^{}*y^{}+",
                x.x[i as usize].a, x.x[i as usize].n[0], x.x[i as usize].n[1]
            );
        }
    }
}

fn printm(m: MVX) {
    for i in 0..20 {
        for j in 0..20 {
            if m.m[i][j] > 0 {
                print!("{}*x^{}*y^{}+", m.m[i][j], i, j);
            }
        }
    }

    print!(" ordering by lex\n");
}

fn mlt(x: i32, y: i32) -> i32 {
    if x == 0 || y == 0 {
        return 0;
    }
    return ((x + y - 2) % (16 - 1)) + 1;
}

/*
fn mltn(n:i32, x:i32) -> i32
{
let mut i;
let mut j;

  if n == 0
  {
    return 1;
  }
  if x == 0
  {
    return 0;
  }
  else
  {
    i = x;
    j=0;
    while (j < (n-1)){
      i  = mlt(i, x);
      j+=1;
    }
    return i;
  }
}
*/

fn mmul(x: MVX, y: MVX) -> MVX {
    let mut c: MVX = unsafe { MaybeUninit::uninit().assume_init() };

    for l in 0..100 {
        for i in 0..100 {
            for j in 0..100 {
                for k in 0..100 {
                    if x.m[j as usize][k as usize] > 0 && y.m[i as usize][l as usize] > 0 {
                        c.m[i + j][k + l] ^= GF[mlt(
                            (FG[x.m[j as usize][k as usize] as usize]) as i32,
                            (FG[y.m[i as usize][l as usize] as usize]) as i32,
                        ) as usize] as i16;
                    }
                }
            }
        }
    }

    return c;
}

fn main() {
    let mut mt: MVX = unsafe { MaybeUninit::uninit().assume_init() };
    let mut xc: MVX = unsafe { MaybeUninit::uninit().assume_init() };
    let cx: MVX;
    let za: MP;
    let tm: MVX;
    let mut v = PO {
        z: vec![vec![1, 1, 1], vec![1, 1, 1]],
    };

    //  let mut t=Vec::new();

    println!("{:?}", v.z);
    v.z[0][1] = 2; //vec![0, 0, 0];
    //t=v.z;
    println!("{:?}", v.z);

    // return;

    mt = msm(mt);
    printm(mt);
    print!("\n");
    xc = msm(xc);
    printm(xc);
    print!("\n");
    cx = mmul(xc, mt);
    printm(cx);
    print!("\n");

    za = v2m(cx);
    printv(za);
    print!(" t={}\n", za.t);

    tm = m2v(za);
    printm(tm);
    print!("\n");
}
