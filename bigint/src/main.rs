extern crate num_bigint;
//use std::i32;

use num_bigint::BigInt;


fn main() {
    //let a = 6;
    //let b = 1;
    //let c = 1;
    //p751=2^372*3^239-1
    let p751 = "10354717741769305252977768237866805321427389645549071170116189679054678940682478846502882896561066713624553211618840202385203911976522554393044160468771151816976706840078913334358399730952774926980235086850991501872665651576831";
    let s = "3248238549023850492380";
    let s2 = "3248238549023850492379";
    
    let n: BigInt = s.parse().unwrap();
    let n2: BigInt = s2.parse().unwrap();
    let p: BigInt = p751.parse().unwrap();

    let add = &n + &n2;
    let suv = &n - &n2;
    let mul = &n * &n2;
    let div = &n / &n2;
    let rem = &n % &n2;

    println!("add: {}", add);  // 6496477098047700984759
    println!("sub: {}", suv);  // 1
    println!("mul: {}", mul);  // 10551053671364569578520014120677744587572020
    println!("div: {}", div);  // 1
    println!("rem: {}", rem);  // 1
    println!("p751: {}", p);
}
