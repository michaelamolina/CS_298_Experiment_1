// GF_2m_arithmetic.rs 
// Used as "f1" for "field 1"
// This file contains:
//    - arithmetic operations add, multiply, reduce, inverse, power, and gcd for elements of GF(2^m)
//    - method to generate a primitive polynomial of degree m over GF(2)
// Explanation: 
//    - elements of GF(2^m) are polynomials of degree m-1 over GF(2).
//    - since these polynomials can be represented by binary strings, 
//      they can be represented by by binary u64s, with higher order 
//      coefficients corresponding to higher order bits.
//    - elements of GF(2^m) are represented by the names a, b, or u, 
//      and the primitive polynomial is represented by f

use crate::GF_2mt_arithmetic as f2;
use std::process;
use rand::{thread_rng, Rng};
use peroxide::fuga::*;
use itertools::Itertools;

// Returns gcd of two polynomials a and b 
pub fn gcd(a:u64, b:u64) -> u64 {
    let mut u = a.clone();
    let mut v = b.clone();
    let mut g1 = 1 as u64;
    let mut g2 = 0 as u64;
    let mut h1 = 0 as u64;
    let mut h2 = 1 as u64;
    while u != 0 {
        let mut j:i64 = deg(u) as i64 - deg(v) as i64;
        if j < 0 {
            swap(&mut u, &mut v);
            swap(&mut g1, &mut g2);
            swap(&mut h1, &mut h2);
            j *= -1;
        }
        let mut p1 = v.clone();
        let mut p2 = g2.clone();
        let mut p3 = h2.clone();
        shift_left_i(&mut p1, j as u64);
        shift_left_i(&mut p2, j as u64);
        shift_left_i(&mut p3, j as u64);
        u = add(u, p1);
        g1 = add(g1, p2);
        h1 = add(h1, p3);
    }
    let mut d = v.clone();
    d
}
// O(degree of bigger)

// Returns the inverse of a modulo f
pub fn inverse(a:u64, f:u64) -> u64 {
    if a == 0 as u64 {
        println!("Error in f1::inverse - a was 0");
        process::exit(1);
    }
    let mut u = a.clone();
    let mut v = f.clone();
    let mut g1 = 1;
    let mut g2 = 0;
    while u != 1 {
        if u == 0 {
            println!("Error in f1: polynomial was not irreducible");
            println!("a = {:#066b}, f = {:#066b}", a, f);
            process::exit(1);
        } 
        let mut j:i64 = deg(u) as i64 - deg(v) as i64;
        if j < 0 as i64 {
            swap(&mut u, &mut v);            
            swap(&mut g1, &mut g2);            
            j *= -1 as i64
        }
        let mut p1 = v.clone();
        let mut p2 = g2.clone();
        shift_left_i(&mut p1, j as u64);
        shift_left_i(&mut p2, j as u64);
        u = add(u, p1);
        g1 = add(g1, p2);
    }
    return g1;
}
// O(m)

// Right-to-left shift-and-add field multiplication 
// Multiplies two polynomials a, b while keeping result less than degree of f
// Precondition: a, b in GF(2^m)
pub fn multiply(a:u64, b:u64, f:u64) -> u64 {
    let mut m = deg(f);
    let mut b = b.clone();
    let mut c = 0;
    if a%2 == 1 {
        c = b;
    }
    for i in 1..m {
        shift_left(&mut b);
        if deg(b) == m {
            b = add(b, rz(f));
            b = add(b, 2_u64.pow(m as u32) as u64);
        }
        if bit(a, i) == 1 {
            c = add(c, b);
        }
    }
    c
}
// O(m)

// Reduces polynomial u modulo f
pub fn reduce(u:&mut u64, f:u64) {
    let m = deg(f);
    let n = deg(*u);
    let mut rz = rz(f);
    shift_left_i(&mut rz, 63-m);
    for i in (m..64).rev() {
        if 2_u64.pow(i as u32) & (*u) != 0 {
            *u = add(*u, rz);
            *u = add(*u, 2_u64.pow(i as u32));
        } 
        shift_right(&mut rz);
    }
}
// O(n)

// Return u^i modulo f for u a polynomial in GF(2^m)
pub fn pow(u:u64, i:u64, f:u64) -> u64 {
    let mut p = 1;
    for i in 0..i {
        p = multiply(p,u,f);
    }
    p
}
// O(m*i)


// Adds two polynomials a and b over GF(2) - same as bitwise XOR
pub fn add(a:u64, b:u64) -> u64 {
    a ^ b
}
// ========= helper methods ===============================================
// Returns the vector representation of polynomial u
// The low order bits go into high order vector index
pub fn vec(u:u64, f:u64) -> Vec<f64> {
    let mut vec = vec![0 as f64; deg(f) as usize];
    let mut u = u.clone();
    for i in (0..vec.len()).rev() {
        if u%2 == 0 {
            vec[i] = 0 as f64;
        } else {
            vec[i] = 1 as f64;
        }
        shift_right(&mut u);
    }
    vec
}
// O(m)


// Swaps two polynomial references
fn swap(a:&mut u64, b:&mut u64) {
    let mut temp = *a;
    *a = *b;
    *b = temp;
}



// Returns true if coefficient at position bit is 1 in polynomial a, false otherwise
fn bit(a:u64, bit:u64) -> u64 {
    if a & (1 << bit) != 0 {
        return 1;
    } else {
        return 0;
    }
}

// Returns poynomial f, except with highest order coefficient removed
// Used with the following facts in multiply, reduce:
// f(z) = z^m + r(z) 
// z^m mod f(z) = r(z)
pub fn rz(f:u64) -> u64 {
    let mut x = 1;
    shift_left_i(&mut x, deg(f));
    add(f, x)   
}


// Returns degree of polynomial a
pub fn deg(a:u64) -> u64 {
    if a == 0 {
        return 0;
    }
    return ((a as f64).log2().floor()+1 as f64) as u64 - 1;
}

// Returns polynomial a as binary string
pub fn string(a:u64) -> String {
    return format!("{:#066b}", a);
}

// Returns polynomial a as a text string: {0,1}x^m-1 + ... + {0,1}x + {0,1}
pub fn poly_string(a:u64) -> String {
    let mut result = String::from("");
    let mut deg = deg(a) as i64;
    while deg >= 0 {
        if a & 2_u64.pow(deg as u32) != 0 {
            if deg != 0 {
                result.push_str(&format!("1x^{} + ", deg));
            } else {
                result.push_str(&format!("1x^{}", deg));
            }
        } else {
            if deg != 0 {
                result.push_str(&format!("0x^{} + ", deg));
            } else {
                result.push_str(&format!("0x^{}", deg));
            }
        }
        deg -= 1;
    }
    result
}

// Divide polynomial a by x
pub fn shift_right(a:&mut u64) -> u64 {
    *a = *a >> 1;
    a.clone()
}

// Multiply polynomial a by x
pub fn shift_left(a:&mut u64) -> u64 {
    *a = *a << 1;
    a.clone()
}

// Divide polynomial a by x^i
pub fn shift_right_i(a:&mut u64, i:u64) -> u64 {
    for i in 0..i {
        shift_right(a); 
    }
    a.clone()
}

// Multiply polynomial a by x^i 
pub fn shift_left_i(a:&mut u64, i:u64) -> u64 {
    for i in 0..i {
        shift_left(a); 
    }
    a.clone()
}

// ========== polynomial methods ==================================

// Returns true if f is irreducible over GF(2), false otherwise
// max degree of f is m = 63
pub fn is_irreducible(f:u64) -> bool {
    let mut m = deg(f);
    let mut px = 2 as u64;
    for i in 0..m {
        px = multiply(px, px, f);
    }
    if px == 2 as u64 {
        return true;
    } else {
        return false;
    }
}


// Returns a primitive (irreducible) polynomial of degree m
pub fn get_irreducible_polynomial(m:u64) -> u64 {
    let mut f = 2_u64.pow(m as u32); // degree m 
    f = add(f, 1 as u64); // constant is 1
    let mut index = rand::thread_rng().gen_range(1..m); // choose 3rd index
    f = add(f, 2_u64.pow(index as u32));
    let mut count = 0;
    while !is_irreducible(f) { 
        f = add(f, 2_u64.pow(index as u32)); // undo third index 
        let mut index = rand::thread_rng().gen_range(1..m); // choose new 3rd index
        f = add(f, 2_u64.pow(index as u32)); // add 3rd index
        count += 1;
        if count > 100 {
            println!("Could not find primitive polynomial over GF(2) of degree {}", m);
            process::exit(1);
        }
    }   
    return f
}

// Returns length random elements of GF(2^m) in a shuffled vector
// If one of them is a root, polynomial g was not irreducible
pub fn create_L(g:&Vec<u64>, f:u64, length:usize) -> Vec<u64> {
    let mut L = Vec::new();
    for i in 0..2_u64.pow(deg(f) as u32) {
        if f2::compute(g, f, i) == 0 {
            println!("Error in choosing L: g(x) is not irreducible over GF(2^{})", deg(f));
            process::exit(1);
        }
        L.push(i);
    }
    L.shuffle(&mut thread_rng());   
    return (&L[0..length]).to_vec();
}
// O(mt^2) * 2^m + 2^m

