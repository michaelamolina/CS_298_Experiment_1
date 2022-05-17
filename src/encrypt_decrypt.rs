// encrypt_decrypt.rs
// This file contains:
//    - methods to encrypt and decrypt a message 
//    - includes Patterson's error correcting algorithm, method to 
//      compute the syndrome, and method to compute error vector 
// Explanation: 
//    - the primitive polynomial for GF(2^mt) is represented by g
//    - the primitive polynomial for GF(2^m) is represented by f
//    - the public key is the scambled generator matrix
//    - the private key is S, G, P, L, g, and f

use crate::GF_2m_arithmetic as f1;
use crate::GF_2mt_arithmetic as f2;
use crate::matrix;
use peroxide::fuga::*;
use itertools::Itertools;

// Returns encrypted message
// generator is the public key
pub fn encrypt(message:&Matrix, generator:&Matrix, length:usize, t:u64) -> Matrix {
    //println!("message = ");
    //message.print();
    let mut mG = matrix::mod2_matrix(&(message.clone() * generator.clone()), 1 as usize, length as usize);
    //println!("mG = ");
    //mG.print();
    let mut error_vector = choose_error_vector(length as u64, t); 
    //println!("error_vector = ");
    //error_vector.print();
    let mut y = matrix::mod2_matrix(&(mG.clone() + error_vector.clone()), 1 as usize, length as usize);
    //println!("y = ");
    //y.print();
    y
}

pub fn encrypt2(message:&Matrix, generator:&Matrix, length:usize, t:u64, error_vector:&Matrix) -> Matrix {
    let mut mG = matrix::mod2_matrix(&(message.clone() * generator.clone()), 1 as usize, length as usize);
    let mut y = matrix::mod2_matrix(&(mG.clone() + error_vector.clone()), 1 as usize, length as usize);
    y
}

// Returns decrypted message
// y is the received, encrypted message
pub fn decrypt(y:&Matrix, S:&Matrix, G:&Matrix, P:&Matrix, g:&Vec<u64>, f:u64, L:&Vec<u64>, length:u64, dimension:u64) -> Matrix {
    let mut yP_inv = y.clone()*P.inv();
    //println!("yP_inv = ");
    //yP_inv.print();
    let mut error_vector = error_vector(&yP_inv, &g, f, &L);
    //println!("error_vector = ");
    //error_vector.print();
    let mut mSG = matrix::mod2_matrix(&(yP_inv + error_vector), 1 as usize, length as usize);
    //println!("mSG = ");
    //mSG.print();
    let mut find_mS = matrix::append(&G.transpose(), &mSG.transpose(), length, dimension, length, 1 as u64);
    //println!("find_mS = ");
    //find_mS.print();
    matrix::rref_for_decrypt(&mut find_mS, length, dimension+1, length, dimension);
    //println!("reduced find_mS = ");
    //find_mS.print();
    let mut mS = zeros(1 as usize, dimension as usize);
    let mut col = dimension as usize;
    for row in 0..dimension {
        mS[(0 as usize,row as usize)] = find_mS[(row as usize, col as usize)];
    }
    //println!("mS = ");
    //mS.print();
    let mut S_inv = matrix::mod2_matrix(&(S.inv()), dimension as usize, dimension as usize);
    let mut decoded = matrix::mod2_matrix(&(mS.clone() * S_inv), 1 as usize, dimension as usize);
    //println!("decrypted = ");
    //decoded.print();
    decoded
}

// Returns error vector which contains moved errors
pub fn error_vector(yP_inv:&Matrix, g:&Vec<u64>, f:u64, L:&Vec<u64>) -> Matrix {
    let mut syndrome = syndrome(yP_inv, g, f, L); 
    //println!("syndrome = \n{}", f2::string(&syndrome));
    let mut sigma = sigma(&syndrome, g, f);
    //println!("sigma = \n{}", f2::string(&sigma));
    let mut error_vector = zeros(1 as usize, L.len() as usize);
    for i in 0..L.len() {
        let mut result = f2::compute(&sigma, f, L[i]);
        if result == 0 as u64 {
            error_vector[(0,i)] = 1 as f64;
        } else {
            error_vector[(0,i)] = 0 as f64;
        }
    }
    error_vector
}

// Returns polynomial to correct errors
// Uses Patterson's error correcting algorithm
pub fn sigma(syndrome:&Vec<u64>, g:&Vec<u64>, f:u64) -> Vec<u64> {
    let syndrome_inv = f2::inverse(syndrome, g, f);
    let mut syndrome_inv_plus_x = f2::add(&syndrome_inv, &vec![0 as u64, 1 as u64]); // low order to high order
    let mut sqrt = f2::sqrt(&syndrome_inv_plus_x, &g, f);
    //println!("v(x) = \n{}", f2::string(&sqrt));
    let (mut bx, mut ax) = f2::get_bx_ax(&sqrt, &g, f);
    //println!("a(x) = {}", f2::string(&ax));
    //println!("b(x) = {}", f2::string(&bx));
    let mut bx2 = f2::square(&bx, &g, f);
    let mut ax2 = f2::square(&ax, &g, f);
    f2::shift_left(&mut bx2); // multiply b(x) by x
    let mut sigma = f2::add(&ax2, &bx2);
    sigma
}
    
// Returns syndrome which is a polynomial in GF(2^tm)
pub fn syndrome(yP_inv:&Matrix, g:&Vec<u64>, f:u64, L:&Vec<u64>) -> Vec<u64> {
    let mut sum = vec![0 as u64];
    for i in 0..L.len() {
        if yP_inv[(0,i)] == 1 as f64 {
            let mut denom = vec![L[i], 1]; // low order to high order 
            let mut inverse = f2::inverse(&denom, &g, f);
            sum = f2::add(&sum, &inverse);
        } 
    }
    sum
}
// O(

// returns vector of all Choose(width, group_size) combinations of possible error vector
pub fn all_possible_error_vectors(width:u64, group_size:u64) -> Vec<Matrix> {
    let it = (0..width).combinations(group_size as usize);
    let mut vectors = Vec::new();
    for element in it {
        let mut vec = zeros(1 as usize, width as usize);
        for index in element {
            vec[(0 as usize,index as usize)] = 1 as f64;
        }
        vectors.push(vec);
    }
    vectors
}

pub fn choose_error_vector(length:u64, t:u64) -> Matrix {
    let mut error_vectors = all_possible_error_vectors(length, t); // takes too long to do all of them
    let mut index = rand::thread_rng().gen_range(0..error_vectors.len());
    return error_vectors[index as usize].clone();
}

