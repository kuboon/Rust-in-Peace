use std::{i32};


fn main(){
 let mut gf: [i32;256]= [0;256];
 let mut fg: [i32;256]= [0;256];
 let o =256;

let normal: [i32; 14]= [
    0b1011,
    0b10011,
    0b110111,
    0b1100001,
    0b11000001,
    0b100011101, //sage
    0b1100110001,
    0b10001101111, //sage1024
    0b110000001101,
    0b1000011101011, //sage 4096
    0b10000000011011, /* Classic McEliece */
    0b110000100010001,
    0b1100000000000001,
    0b11010000000010001];

  let mut l :i32;
  let mut _m:i32;
  let mut bit:i32;
  let mut count: i32;
  let n:i32;
  let pol:i32;
 
  l=1;
  pol=normal[5];
  count=0;
  while pol > l //原始多項式の最大次数を計算する。
  {
    l = l << 1;
    count=count+1;
  }
  l = l >> 1;
  n = pol ^ l; //原始多項式の最大次数を消した多項式の残り。

  gf[0] = 0;
  bit = 1;
  for i in 1..256
  {
    gf[i]=0;
    fg[i]=0;
  
    if bit > (l - 1) //もしbitが最大次数に達したら
    {
      bit = bit - l; //bitから最大次数の項 x^n を消す。
      bit = bit ^ n; //最大次数の項以下の原始多項式を bit に xorする。
    }
    gf[i] = bit; //最初は何もしないでx^iを入れていく。
    
    bit = bit << 1; //原始多項式の次数を1上げる。
  }
println!("static const unsigned short gf[{}]={{", o);


  for j in 0..256
  {
      print!("{},", gf[j]);
      
  }

  println!("}};");
  
  
  println!("static const unsigned short fg[256]={{");
  for i in 0..256
    {
      for j in 0..256
      {
      if i == gf[j]
        {
          fg[i as usize] = j as i32;
        }
    }
  }

      for i in 0..256
    {
        print!("{},", fg[i]);
    }

    println!("}};\n");
}
