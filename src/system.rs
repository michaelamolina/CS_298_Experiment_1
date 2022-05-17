// system.rs

use std::process;
use peroxide::fuga::*;
use crate::GF_2m_arithmetic  as f1;
use crate::GF_2mt_arithmetic as f2;
use crate::matrix;
use crate::encrypt_decrypt;
use crate::test;
use crate::util;


pub struct public_key {
    pub generator: Matrix,
    pub length: u64,
    pub dimension: u64, 
    pub t: u64
}

pub struct private_key {
    pub S: Matrix,
    pub G: Matrix, 
    pub P: Matrix, 
    pub L: Vec<u64>,
    pub g: Vec<u64>,
    pub f: u64,
    pub length: u64,
    pub dimension: u64
}

pub fn create_system(n:u64, t:u64) -> (public_key, private_key) {
    let m = (n as f64).log2();
    if m.fract() != 0.0 {
        println!("Error: n must be a power of 2");
        process::exit(1);
    } 
    let mut m = m as u64;
    let mut length = 2_u64.pow(m as u32); // set n = 2^m always
    test::check_precondition(length, t, m);
    let mut f = f1::get_irreducible_polynomial(m);
    let mut g = f2::get_irreducible_polynomial(t, f);
    let mut L = f1::create_L(&g, f, length as usize);
    let (mut G, mut dimension, mut length) = matrix::generator(&g, f, &L);
    let mut S = matrix::dense_nonsingular_matrix(dimension as u64);
    let mut P = matrix::permutation_matrix(length as u64); 
    let mut generator = matrix::mod2_matrix(&(S.clone() * G.clone() * P.clone()), dimension as usize, length as usize);
    let mut public_key = public_key {
        generator: generator,
        length: length,
        dimension: dimension,
        t: t
    };
    let mut private_key = private_key {
        S: S,
        G: G,
        P: P,
        L: L,
        g: g, 
        f: f,
        length: length,
        dimension: dimension
    };
    (public_key, private_key)
}

pub fn initialization_vector(public_key:&public_key) -> Matrix {
    let mut dimension = public_key.dimension;
    let mut iv = zeros(1 as usize, dimension as usize);
    let mut init_number = rand::thread_rng().gen_range(0..2_u64.pow(dimension as u32));
    for i in (0..dimension as usize).rev() { // fill matrix from right to left
        if init_number%2 == 0 {
            iv[(0,i)] = 0 as f64;
        } else {
            iv[(0,i)] = 1 as f64;
        }
        init_number = init_number >> 1;
    }
    iv
}

pub fn encrypt(message:String, public_key:&public_key, iv:&Matrix) -> Vec<u64> {
    let mut plaintext = util::string_to_vec_of_u64s(message);
    let     dimension = public_key.dimension;
    let        length = public_key.length;
    while plaintext.len() % dimension as usize != 0 { // make sure vector of 1 and 0 is divisible by dimension
        plaintext.push(0);
    }
    let mut prev_ciphertext = iv.clone(); 
    let mut      ciphertext = Vec::new(); // vector of 1 x length matrices
    let              blocks = vec_into_matrix_blocks(&plaintext, dimension); // split plaintext into 1 x dimension matrices
    for block in blocks { // iterate through matrices
        let x = xor_matrix(&block, &prev_ciphertext, dimension); // xor with previous ciphertext
        let mut curr_ciphertext = encrypt_decrypt::encrypt(&x, 
                                                            &public_key.generator,  
                                                            public_key.length as usize, 
                                                            public_key.t);
        prev_ciphertext = curr_ciphertext.clone();
        ciphertext.push(curr_ciphertext.clone()); // ciphertext is a vector of encrypted matrices
    }
    let mut ciphertext = matrix_blocks_into_vec(&ciphertext, length);
    ciphertext
}

pub fn decrypt(ciphertext:&Vec<u64>, private_key:&private_key, iv:&Matrix) -> String {
    let mut ciphertext = vec_into_matrix_blocks(&ciphertext, private_key.length); 
    let mut plaintext = Vec::new();
    let mut prev_ciphertext = iv.clone();
    for block in ciphertext {
        let mut plaintext_block = encrypt_decrypt::decrypt(&block, 
                                                             &private_key.S, 
                                                             &private_key.G, 
                                                             &private_key.P, 
                                                             &private_key.g, 
                                                             private_key.f, 
                                                             &private_key.L, 
                                                             private_key.length, 
                                                             private_key.dimension);
        let mut plaintext_block = xor_matrix(&plaintext_block, &prev_ciphertext, private_key.dimension);
        plaintext.push(plaintext_block);
        prev_ciphertext = block.clone();
    }
    let mut plaintext = matrix_blocks_into_vec(&plaintext, private_key.dimension);
    while plaintext.len() % 8 as usize != 0 {
        plaintext.pop();
    }
    let mut plaintext = util::vec_of_u64s_to_string(&plaintext);
    plaintext
}

pub fn matrix(a:&Vec<u64>, dimension:u64) -> Matrix {
    let mut m = zeros(1 as usize, dimension as usize);
    for i in 0..dimension {
        m[(0,i as usize)] = a[i as usize] as f64;
    }
    m
}

// convert vector of 1 x width matrices to one long vector 
pub fn matrix_blocks_into_vec(a:&Vec<Matrix>, width:u64) -> Vec<u64> {
    let mut v = Vec::new();
    for i in 0..a.len() {
        for j in 0..width as usize {
            v.push(a[i][(0,j)] as u64);
        }   
    }
    v
}

// takes in vector of length that is divisible by dimension
// returns vector of matrices of length equal to dimension
// precondition: a.len()%dimension == 0
pub fn vec_into_matrix_blocks(a:&Vec<u64>, dimension:u64) -> Vec<Matrix> {
    let mut blocks = Vec::new();
    let mut i = 0;
    while i < a.len() {
        let mut block = Vec::new();
        while block.len() < dimension as usize {
            block.push(a[i]);
            i += 1;
        }
        blocks.push(matrix(&block, dimension));
    }
    blocks
}

// takes in two 1 x dimension matrices
// return 1 x dimension matrix that represents their xor
pub fn xor_matrix(a:&Matrix, b:&Matrix, dimension:u64) -> Matrix {
    let mut c = zeros(1 as usize, dimension as usize);
    for i in 0..dimension as usize {
        if a[(0,i)] == 1.0 && b[(0,i)] == 1.0 {
            c[(0,i)] = 0 as f64;
        } else if a[(0,i)] == 1.0 || b[(0,i)] == 1.0 {
            c[(0,i)] = 1.0 as f64;
        } else {
            c[(0,i)] = 0.0 as f64;
        }
    }
    c
}

