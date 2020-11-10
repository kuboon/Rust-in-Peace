fn main() {
    let closures = [3, 7, 1, 5, 8, 9, 2].iter().map(|&i| {
        move |j| i + j
    }).collect::<Vec<_>>();
    println!("{}", closures[3](14));
}
