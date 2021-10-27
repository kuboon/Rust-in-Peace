fn main() {
    let vector = vec![20, 80, 60, 40];
    let s = sum(&vector);
    println!("{:?} の総和は {}", vector, s); // vector が使える
}

fn sum(v: &Vec<i32>) -> i32 {
    let mut ret = 0;
    for &i in v {
        ret += i;
    }
    ret
}