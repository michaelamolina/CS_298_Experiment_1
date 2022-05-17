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

fn main() {

    //let mut message = String::from("!@#$%^&*()_+-=1234567890kjfkjhg");
    let mut message = (&text::text()[0..40 as usize]).to_string();
    println!("message.len() = {} characters", message.len());
    let (public_key, private_key) = system::create_system(256, 3);
    println!("dimension = {}", public_key.dimension);
    let iv                        = system::initialization_vector(&public_key); 
    let ciphertext                = system::encrypt(message.clone(), &public_key, &iv);
    println!("finished encrypting, starting decryption");
    let plaintext                 = system::decrypt(&ciphertext, &private_key, &iv);
    if message == plaintext {
        println!("True");
    } else {
        println!("False");
    }
    
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


