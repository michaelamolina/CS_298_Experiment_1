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
    let mut message = (&text::text()[0..512 as usize]).to_string();
    println!("message.len() = {} characters", message.len());
    let (public_key, private_key) = system::create_system(16, 2);
    let iv                        = system::initialization_vector(&public_key); 
    let ciphertext                = system::encrypt(message.clone(), &public_key, &iv);
    let plaintext                 = system::decrypt(&ciphertext, &private_key, &iv);
    if message == plaintext {
        println!("True");
    } else {
        println!("False");
    }
}


