// main.rs 

mod text;
mod GF_2m_arithmetic;
mod GF_2mt_arithmetic;
mod matrix;
mod encrypt_decrypt;
mod test;
mod system;
use GF_2m_arithmetic as f1;
use GF_2mt_arithmetic as f2;
mod util;
use std::time::Instant;
use std::time::Duration;

fn main() {

    // experiment 1 attack 1
    /*let mut total_duration:Duration = Duration::new(0,0);
    for i in 0..100 {
        let start = Instant::now();
        let (public_key, private_key) = system::create_system(8, 2);
        //let mut bits = 3*64 + public_key.length*public_key.dimension*64;
        //println!("public key size of ({},{}) McEliece cryptosystem: {} bits", public_key.length, public_key.dimension, bits);
        let duration = start.elapsed();
        //println!("duration for ({},{}) McEliece cryptosystem: {:?}", public_key.length, public_key.t, duration);
        total_duration += duration;
    }
    let mut average_duration = total_duration/100;
    println!("average duration for ({},{}) McEliece cryptosystem: {:?}", 8, 2, average_duration);
    
    let mut total_duration:Duration = Duration::new(0,0);
    for i in 0..100 {
        let start = Instant::now();
        let (public_key, private_key) = system::create_system(8, 2);
        //let mut bits = 3*64 + public_key.length*public_key.dimension*64;
        //println!("public key size of ({},{}) McEliece cryptosystem: {} bits", public_key.length, public_key.dimension, bits);
        let duration = start.elapsed();
        //println!("duration for ({},{}) McEliece cryptosystem: {:?}", public_key.length, public_key.t, duration);
        total_duration += duration;
    }
    let mut average_duration = total_duration/100;
    println!("average duration for ({},{}) McEliece cryptosystem: {:?}", 8, 2, average_duration);

    let mut total_duration:Duration = Duration::new(0,0);
    for i in 0..100 {
        let start = Instant::now();
        let (public_key, private_key) = system::create_system(16, 2);
        //let mut bits = 3*64 + public_key.length*public_key.dimension*64;
        //println!("public key size of ({},{}) McEliece cryptosystem: {} bits", public_key.length, public_key.dimension, bits);
        let duration = start.elapsed();
        //println!("duration for ({},{}) McEliece cryptosystem: {:?}", public_key.length, public_key.t, duration);
        total_duration += duration;
    }
    let mut average_duration = total_duration/100;
    println!("average duration for ({},{}) McEliece cryptosystem: {:?}", 16, 2, average_duration);

    let mut total_duration:Duration = Duration::new(0,0);
    for i in 0..100 {
        let start = Instant::now();
        let (public_key, private_key) = system::create_system(32, 2);
        //let mut bits = 3*64 + public_key.length*public_key.dimension*64;
        //println!("public key size of ({},{}) McEliece cryptosystem: {} bits", public_key.length, public_key.dimension, bits);
        let duration = start.elapsed();
        //println!("duration for ({},{}) McEliece cryptosystem: {:?}", public_key.length, public_key.t, duration);
        total_duration += duration;
    }
    let mut average_duration = total_duration/100;
    println!("average duration for ({},{}) McEliece cryptosystem: {:?}", 32, 2, average_duration);

    let mut total_duration:Duration = Duration::new(0,0);
    for i in 0..100 {
        let start = Instant::now();
        let (public_key, private_key) = system::create_system(64, 2);
        //let mut bits = 3*64 + public_key.length*public_key.dimension*64;
        //println!("public key size of ({},{}) McEliece cryptosystem: {} bits", public_key.length, public_key.dimension, bits);
        let duration = start.elapsed();
        //println!("duration for ({},{}) McEliece cryptosystem: {:?}", public_key.length, public_key.t, duration);
        total_duration += duration;
    }
    let mut average_duration = total_duration/100;
    println!("average duration for ({},{}) McEliece cryptosystem: {:?}", 64, 2, average_duration);

    let mut total_duration:Duration = Duration::new(0,0);
    for i in 0..100 {
        let start = Instant::now();
        let (public_key, private_key) = system::create_system(64, 2);
        //let mut bits = 3*64 + public_key.length*public_key.dimension*64;
        //println!("public key size of ({},{}) McEliece cryptosystem: {} bits", public_key.length, public_key.dimension, bits);
        let duration = start.elapsed();
        //println!("duration for ({},{}) McEliece cryptosystem: {:?}", public_key.length, public_key.t, duration);
        total_duration += duration;
    }
    let mut average_duration = total_duration/100;
    println!("average duration for ({},{}) McEliece cryptosystem: {:?}", 64, 2, average_duration);

    let mut total_duration:Duration = Duration::new(0,0);
    for i in 0..100 {
        let start = Instant::now();
        let (public_key, private_key) = system::create_system(128, 2);
        //let mut bits = 3*64 + public_key.length*public_key.dimension*64;
        //println!("public key size of ({},{}) McEliece cryptosystem: {} bits", public_key.length, public_key.dimension, bits);
        let duration = start.elapsed();
        //println!("duration for ({},{}) McEliece cryptosystem: {:?}", public_key.length, public_key.t, duration);
        total_duration += duration;
    }
    let mut average_duration = total_duration/100;
    println!("average duration for ({},{}) McEliece cryptosystem: {:?}", 128, 2, average_duration);

    let mut total_duration:Duration = Duration::new(0,0);
    for i in 0..100 {
        let start = Instant::now();
        let (public_key, private_key) = system::create_system(256, 2);
        //let mut bits = 3*64 + public_key.length*public_key.dimension*64;
        //println!("public key size of ({},{}) McEliece cryptosystem: {} bits", public_key.length, public_key.dimension, bits);
        let duration = start.elapsed();
        //println!("duration for ({},{}) McEliece cryptosystem: {:?}", public_key.length, public_key.t, duration);
        total_duration += duration;
    }
    let mut average_duration = total_duration/100;
    println!("average duration for ({},{}) McEliece cryptosystem: {:?}", 256, 2, average_duration);

    let mut total_duration:Duration = Duration::new(0,0);
    for i in 0..100 {
        let start = Instant::now();
        let (public_key, private_key) = system::create_system(512, 2);
        //let mut bits = 3*64 + public_key.length*public_key.dimension*64;
        //println!("public key size of ({},{}) McEliece cryptosystem: {} bits", public_key.length, public_key.dimension, bits);
        let duration = start.elapsed();
        //println!("duration for ({},{}) McEliece cryptosystem: {:?}", public_key.length, public_key.t, duration);
        total_duration += duration;
    }
    let mut average_duration = total_duration/100;
    println!("average duration for ({},{}) McEliece cryptosystem: {:?}", 512, 2, average_duration);*/

    // experiment 1 attack 2 ===========================================
    let mut total_duration:Duration = Duration::new(0,0);
    for i in 0..100 {
        let start = Instant::now();
        let (public_key, private_key) = system::create_system(8, 2);
        //let mut bits = 3*64 + public_key.length*public_key.dimension*64;
        //println!("public key size of ({},{}) McEliece cryptosystem: {} bits", public_key.length, public_key.dimension, bits);
        let duration = start.elapsed();
        total_duration += duration;
    }
    let mut average_duration = total_duration/100;
    println!("average duration for ({},{}) McEliece cryptosystem: {:?}", 8, 2, average_duration);

    let mut total_duration:Duration = Duration::new(0,0);
    for i in 0..100 {
        let start = Instant::now();
        let (public_key, private_key) = system::create_system(16, 2);
        //let mut bits = 3*64 + public_key.length*public_key.dimension*64;
        //println!("public key size of ({},{}) McEliece cryptosystem: {} bits", public_key.length, public_key.dimension, bits);
        let duration = start.elapsed();
        total_duration += duration;
    }
    let mut average_duration = total_duration/100;
    println!("average duration for ({},{}) McEliece cryptosystem: {:?}", 16, 2, average_duration);

    let mut total_duration:Duration = Duration::new(0,0);
    for i in 0..100 {
        let start = Instant::now();
        let (public_key, private_key) = system::create_system(32, 2);
        //let mut bits = 3*64 + public_key.length*public_key.dimension*64;
        //println!("public key size of ({},{}) McEliece cryptosystem: {} bits", public_key.length, public_key.dimension, bits);
        let duration = start.elapsed();
        total_duration += duration;
    }
    let mut average_duration = total_duration/100;
    println!("average duration for ({},{}) McEliece cryptosystem: {:?}", 32, 2, average_duration);

    let mut total_duration:Duration = Duration::new(0,0);
    for i in 0..100 {
        let start = Instant::now();
        let (public_key, private_key) = system::create_system(64, 2);
        //let mut bits = 3*64 + public_key.length*public_key.dimension*64;
        //println!("public key size of ({},{}) McEliece cryptosystem: {} bits", public_key.length, public_key.dimension, bits);
        let duration = start.elapsed();
        total_duration += duration;
    }
    let mut average_duration = total_duration/100;
    println!("average duration for ({},{}) McEliece cryptosystem: {:?}", 64, 2, average_duration);

    let mut total_duration:Duration = Duration::new(0,0);
    for i in 0..100 {
        let start = Instant::now();
        let (public_key, private_key) = system::create_system(128, 2);
        //let mut bits = 3*64 + public_key.length*public_key.dimension*64;
        //println!("public key size of ({},{}) McEliece cryptosystem: {} bits", public_key.length, public_key.dimension, bits);
        let duration = start.elapsed();
        total_duration += duration;
    }
    let mut average_duration = total_duration/100;
    println!("average duration for ({},{}) McEliece cryptosystem: {:?}", 128, 2, average_duration);

    let mut total_duration:Duration = Duration::new(0,0);
    for i in 0..100 {
        let start = Instant::now();
        let (public_key, private_key) = system::create_system(256, 2);
        //let mut bits = 3*64 + public_key.length*public_key.dimension*64;
        //println!("public key size of ({},{}) McEliece cryptosystem: {} bits", public_key.length, public_key.dimension, bits);
        let duration = start.elapsed();
        total_duration += duration;
    }
    let mut average_duration = total_duration/100;
    println!("average duration for ({},{}) McEliece cryptosystem: {:?}", 256, 2, average_duration);

    let mut total_duration:Duration = Duration::new(0,0);
    for i in 0..10 {
        let start = Instant::now();    
        let (public_key, private_key) = system::create_system(512, 2);
        //let mut bits = 3*64 + public_key.length*public_key.dimension*64;
        //println!("public key size of ({},{}) McEliece cryptosystem: {} bits", public_key.length, public_key.dimension, bits);
        let duration = start.elapsed();
        total_duration += duration;
    }
    let mut average_duration = total_duration/10;
    println!("average duration for ({},{}) McEliece cryptosystem: {:?}", 512, 2, average_duration);

    let mut total_duration:Duration = Duration::new(0,0);
    for i in 0..10 {
        let start = Instant::now();
        let (public_key, private_key) = system::create_system(1024, 2);
        //let mut bits = 3*64 + public_key.length*public_key.dimension*64;
        //println!("public key size of ({},{}) McEliece cryptosystem: {} bits", public_key.length, public_key.dimension, bits);
        let duration = start.elapsed();
        total_duration += duration;
    }
    let mut average_duration = total_duration/10;
    println!("average duration for ({},{}) McEliece cryptosystem: {:?}", 1024, 2, average_duration);

    let mut total_duration:Duration = Duration::new(0,0);
    for i in 0..10 {
        let start = Instant::now();
        let (public_key, private_key) = system::create_system(2048, 2);
        //let mut bits = 3*64 + public_key.length*public_key.dimension*64;
        //println!("public key size of ({},{}) McEliece cryptosystem: {} bits", public_key.length, public_key.dimension, bits);
        let duration = start.elapsed();
        total_duration += duration;
    }
    let mut average_duration = total_duration/10;
    println!("average duration for ({},{}) McEliece cryptosystem: {:?}", 2048, 2, average_duration);

    //let mut message = String::from("!@#$%^&*()_+-=1234567890kjfkjhg");
    /*let mut message = (&text::text()[0..40 as usize]).to_string();
    println!("message.len() = {} characters", message.len());
    let (public_key, private_key) = system::create_system(8, 2);
    //println!("dimension = {}", public_key.dimension);
    //println!("length = {}", public_key.length);
    let iv                        = system::initialization_vector(&public_key); 
    let ciphertext                = system::encrypt(message.clone(), &public_key, &iv);
    println!("finished encrypting, starting decryption");
    let plaintext                 = system::decrypt(&ciphertext, &private_key, &iv);
    if message == plaintext {
        println!("True");
    } else {
        println!("False");
    }*/
    
   // let mut f = 99;
   // let mut is_irreducible = f1::is_irreducible(f);
   // println!("is_irreducible = {}", is_irreducible);

    /*let mut f = f1::get_irreducible_polynomial(4);
    println!("f = {}", f1::string(f));
    let mut is_irreducible = f1::is_irreducible(f);
    println!("is_irreducible = {}", is_irreducible);*/

    /*let mut f = 19 as u64;
    let mut g = vec![8,1,0,0,1];
    let mut is_irreducible = f2::is_irreducible(&g, f);
    println!("is_irreducible = {}", is_irreducible);*/

    // works
    /*let mut f = f1::get_irreducible_polynomial(10);
    println!("f = {}", f1::string(f));
    let mut is_irreducible = f1::is_irreducible(f);
    println!("is_irreducible = {}", is_irreducible);
    let mut g = f2::get_irreducible_polynomial(3, f);
    println!("g = \n{}", f2::string(&g));
    let mut is_irreducible = f2::is_irreducible(&g, f);
    println!("is_irreducible = {}", is_irreducible);*/
}


