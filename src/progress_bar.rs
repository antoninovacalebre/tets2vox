use std::io;
use std::io::prelude::*;

pub fn pb_init(text: &str) -> u32 {
    print!("{}", text);
    io::stdout().flush().ok().expect("Could not flush stdout");
    0
}

pub fn pb_update(iter: u32, total: u32, old_nbars: u32, length: u32) -> u32 {
    let nbars = length * iter / total;
    let ten = length / 10;

    for i in old_nbars..nbars {
        let imod10 = i % ten;
        if imod10 == 0 && i != 0 {
            print!("{}", i / ten * 10);
            io::stdout().flush().ok().expect("Could not flush stdout");
        }
        print!("-");
        io::stdout().flush().ok().expect("Could not flush stdout");
    }

    if iter == total {
        println!("100");
    }

    nbars
}
