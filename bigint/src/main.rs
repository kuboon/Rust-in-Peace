extern crate num_bigint;
use num_bigint::BigInt;


fn main() {
    let s = "3248238549023850492380";
    let n: BigInt = s.parse().unwrap();

    println!("{}", n);    // 3248238549023850492380

/*
    let s2 = "3248238549023850492379";
    let n: BigInt = s.parse().unwrap();
    let n2: BigInt = s2.parse().unwrap();

    let add = &n + &n2;
    let suv = &n - &n2;
    let mul = &n * &n2;
    let div = &n / &n2;
    let rem = &n % &n2;

    println!("add: {}", add);  // 6496477098047700984759
    println!("sub: {}", sub);  // 1
    println!("mul: {}", mul);  // 10551053671364569578520014120677744587572020
    println!("div: {}", div);  // 1
    println!("rem: {}", rem);  // 1
*/
}

