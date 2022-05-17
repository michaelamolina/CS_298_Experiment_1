// test.rs

use rand::{thread_rng, Rng};
use peroxide::fuga::*;
use std::process;
use itertools::Itertools;
use rand::seq::SliceRandom;
use crate::GF_2m_arithmetic  as f1;
use crate::GF_2mt_arithmetic as f2;
use crate::matrix;
use crate::encrypt_decrypt;
use std::time::Instant;
use std::time::Duration;
use get_size::GetSize;


pub fn run_some_tests() {
    let start = Instant::now();
    // test 1
    let (mut n, mut t, mut length)  = (16, 2, 16);
    let (correct, incorrect, total) = test(n, t, length, 1);
    let duration = start.elapsed();
    //println!("duration for parameters ({},{}) = {:?}", n, t, duration);
    println!("message length for parameters ({},{}) = {}", length, t, length*8 );

    let start = Instant::now();
    // test 2
    let (mut n, mut t, mut length)  = (16, 3, 16);
    let (correct, incorrect, total) = test(n, t, length, 2);
    let duration = start.elapsed();
    //println!("duration for parameters ({},{}) = {:?}", n, t, duration);
    println!("message length for parameters ({},{}) = {}", length, t, length*8 );

    let start = Instant::now();
    // test 3
    let (mut n, mut t, mut length)  = (32, 2, 16);
    let (correct, incorrect, total) = test(n, t, length, 3);
    let duration = start.elapsed();
    //println!("duration for parameters ({},{}) = {:?}", n, t, duration);
    println!("message length for parameters ({},{}) = {}", 32, t, 32*8 );

    let start = Instant::now();
    // test 4
    let (mut n, mut t, mut length)  = (32, 3, 16);
    let (correct, incorrect, total) = test(n, t, length, 4);
    let duration = start.elapsed();
    //println!("duration for parameters ({},{}) = {:?}", n, t, duration);
    println!("message length for parameters ({},{}) = {}", 32, t, 32*8 );

    let start = Instant::now();
    // test 5
    let (mut n, mut t, mut length)  = (64, 2, 16);
    let (correct, incorrect, total) = test(n, t, length, 5); 
    let duration = start.elapsed();
    //println!("duration for parameters ({},{}) = {:?}", n, t, duration);
    println!("message length for parameters ({},{}) = {}", 64, t, 64*8 );
}




pub fn test(n:u64, t:u64, length:u64, test_number:u64) -> (u64, u64, u64) {
    let m = (n as f64).log2();
    if m.fract() != 0.0 {
        println!("Error: n must be a power of 2");
        process::exit(1);
    } 
    let mut m = m as u64;
    check_precondition(length, t, m);
    let mut f = f1::get_irreducible_polynomial(m);
    let mut g = f2::get_irreducible_polynomial(t, f);
    //println!("test {} polynomials:", test_number);
    //println!("{}", f1::string(f));
    //println!("{}", f2::string(&g));
    let mut L = f1::create_L(&g, f, length as usize);
    let (mut G, mut dimension, mut length) = matrix::generator(&g, f, &L);
    let mut S = matrix::dense_nonsingular_matrix(dimension as u64);
    let mut P = matrix::permutation_matrix(length as u64); 
    let mut generator = matrix::mod2_matrix(&(S.clone() * G.clone() * P.clone()), dimension as usize, length as usize);
    // create public and private keys here 
    let (correct, incorrect, total) = run_unit_tests(dimension, length, t, &generator, &S, &G, &P, &g, f, &L); 
    //println!("test {} results:\n  {}/{} correct, {}/{} incorrect\n", test_number, correct, total, incorrect, total);
    return (correct, incorrect, total);
}


pub fn run_unit_tests(dimension:u64, length:u64, t:u64, generator:&Matrix, S:&Matrix, G:&Matrix, P:&Matrix, g:&Vec<u64>, f:u64, L:&Vec<u64>) -> (u64, u64, u64) {
    let mut messages      = all_possible_messages(dimension);
    let mut error_vectors = all_possible_error_vectors(length, t); // takes too long to do all of them
    let mut error_indices: Vec<u32> = (0..error_vectors.len() as u32).collect();
    let mut number_of_error_vectors = 10;
    let mut chosen_error_indices:Vec<u32> = error_indices.choose_multiple(&mut rand::thread_rng(), number_of_error_vectors).cloned().collect();
    let mut correct_count = 0;
    let mut incorrect_count = 0;

    /*let public_key = public_key {
        generator: generator.clone(),
        length: length,
        dimension: dimension, 
        t: t
    };

    let private_key = private_key {
        S: S.clone(),
        G: G.clone(), 
        P: P.clone(), 
        L: L.clone(),
        g: g.clone(),
        f: f,
        length: length,
        dimension: length
    };*/

    let mut k = dimension;
    let mut n = length;
    let mut S_size = k*k;
    let mut P_size = n*n;
    let mut G_size = k*n;
    let mut L_size = n;
    let mut generator_size = k*n;

    let mut private_key_size = generator_size;
    let mut public_key_size = S_size + P_size + G_size + L_size;
    //println!("\npublic_key_size for parameters ({},{}) = {}", length, t, public_key_size*8);
    //println!("private_key_size for parameters ({},{}) = {}", length, t,  private_key_size*8);
    let mut count = 0;
    for message in &messages {
        for error_index in &chosen_error_indices {
            let mut y = encrypt_decrypt::encrypt2(message, &generator, length as usize, t, &error_vectors[*error_index as usize]);
            let mut decrypted = encrypt_decrypt::decrypt(&y, &S, &G, &P, &g, f, &L, length, dimension);
            if is_equal(&message, &decrypted, dimension as u64) {
                correct_count += 1;
            } else {
                incorrect_count += 1;
            }
            count += 1;
            if count > 0 {
                break;
            }
        }
        if count > 0 {
            break;
        }
    } 
    (correct_count, incorrect_count, messages.len() as u64 * number_of_error_vectors as u64)// * number_of_error_vectors as u64)
}


pub fn check_precondition(length:u64, t:u64, m:u64) {
    if length <= t*m {
        println!("Error: length must be > t*m");
        process::exit(1);
    }
}

// check if two vectors are equal
pub fn is_equal(a:&Matrix, b:&Matrix, width:u64) -> bool {
    for i in 0..width as usize {
        if a[(0,i)] != b[(0,i)] {
            return false;
        } 
    }
    return true;
}

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

// returns vector of all possible binary vectors of width width
pub fn all_possible_messages(width:u64) -> Vec<Matrix> {
    let mut vectors = Vec::new();
    for i in 0..2_u64.pow(width as u32) {
        let mut vec = f1::vec(i, 2_u64.pow(width as u32));
        let mut matrix = zeros(1 as usize, width as usize);
        for j in 0..vec.len() {
            matrix[(0 as usize, j as usize)] = vec[j as usize];
        } 
        vectors.push(matrix);
    }
    vectors
}




