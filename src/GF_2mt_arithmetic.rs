// GF_2mt_arithmetic.rs
// Used as "f2" for "field 2"
// This file contains:
//    - arithmetic operations add, multiply, reduce, inverse, square root, and gcd for elements of GF(2^mt)
//    - method to generate a primitive polynomial of degree t over GF(2^m)
// Explanation: 
//    - elements of GF(2^mt) are polynomials of degree t-1 over GF(2^m).
//    - they are represented as vectors of u64s, with the u64 being the coefficient in GF(2^m),
//      and low order indices corresponding to low order coefficients
//    - elements of GF(2^mt) are represented by the names a or b, with the primitive polynomial 
//      represented by g, and the primitive polynomial for the field GF(2^m) represented by f
//    - for visualization, the higher order coefficients are on the left so that a multiply by x
//      is a shift left 

use crate::GF_2m_arithmetic as f1;
use std::process;
use std::mem;
use rand::{thread_rng, Rng};
use peroxide::fuga::*;
use itertools::Itertools;

// Returns gcd of a and b
pub fn gcd(a:&Vec<u64>, b:&Vec<u64>, f:u64) -> Vec<u64> {
    if is_zero(&a) || is_zero(&b) {
        return vec![0];
    }
    let mut u = a.clone();
    let mut v = b.clone();
    let mut g1 = vec![1 as u64];
    let mut g2 = vec![0 as u64];
    let mut h1 = vec![0 as u64];
    let mut h2 = vec![1 as u64];
    let mut t = deg(&b);
    let mut count = 0;
    while !is_zero(&u) {
        /*if count%100 == 0 {
            println!("f2::gcd iteration {}", count);
        }*/
        if count > 5000 {
            println!("Error: max iterations in f2 gcd");
            process::exit(1);
        }
        if is_zero(&v) {
            println!("Error in f2 gcd: g = {}", string(&b));
            process::exit(1);
            //return vec![0];
        }            
        let mut j:i64 = deg(&u) as i64 - deg(&v) as i64;
        if j<0 as i64 {
            mem::swap(&mut u, &mut v);
            mem::swap(&mut g1, &mut g2);
            mem::swap(&mut h1, &mut h2);
            j *= -1 as i64;
        }
        let mut p1 = v.clone();
        let mut p2 = g2.clone();
        let mut p3 = h2.clone();
        shift_left_i(&mut p1, j as usize);
        shift_left_i(&mut p2, j as usize);
        shift_left_i(&mut p3, j as usize);
        let mut u_coeff = u[u.len()-1];
        let mut v_coeff = v[v.len()-1];
        let mut coeff = f1::multiply(f1::inverse(v_coeff, f), u_coeff, f); 
        p1 = multiply_by_constant(&p1, coeff, f);
        p2 = multiply_by_constant(&p2, coeff, f);
        p3 = multiply_by_constant(&p3, coeff, f);
        let mut len1 = u.len();
        u  = add(&u, &p1);
        let mut len2 = u.len();
        g1 = add(&g1, &p2);
        h1 = add(&h1, &p3);
        count += 1;
    }
    let mut d = v;
    d    
}
// O(t*tm*m)

// Returns (b(x), a(x)) such that v(x)*b(x) = a(x) mod g(x)
// Pass in (v(x), g(x)) and use extended euclidean algorithm
// Used in decode 
pub fn get_bx_ax(vx:&Vec<u64>, gx:&Vec<u64>, f:u64) -> (Vec<u64>, Vec<u64>) {
    let mut u = vx.clone();
    let mut v = gx.clone();
    let mut g1 = vec![1 as u64];
    let mut g2 = vec![0 as u64];
    let mut h1 = vec![0 as u64];
    let mut h2 = vec![1 as u64];
    let mut t = deg(&gx);
    let mut u_bound = ((t as f64)/(2 as f64)).floor() as u64;
    let mut g1_bound = ((t as f64 - 1 as f64)/(2 as f64)).floor() as u64;
    let mut count = 0;
    while !is_zero(&u) {
        if count > 100 {
            println!("Error: max iterations in f2 get_bx_ax");
            process::exit(1);
        }
        if is_zero(&v) {
            println!("Error in f2 get_bx_ax: g = {}", string(&gx));
            process::exit(1);
        }            
        let mut j:i64 = deg(&u) as i64 - deg(&v) as i64;
        if j<0 as i64 {
            mem::swap(&mut u, &mut v);
            mem::swap(&mut g1, &mut g2);
            mem::swap(&mut h1, &mut h2);
            j *= -1 as i64;
        }
        let mut p1 = v.clone();
        let mut p2 = g2.clone();
        let mut p3 = h2.clone();
        shift_left_i(&mut p1, j as usize);
        shift_left_i(&mut p2, j as usize);
        shift_left_i(&mut p3, j as usize);
        let mut u_coeff = u[u.len()-1];
        let mut v_coeff = v[v.len()-1];
        let mut coeff = f1::multiply(f1::inverse(v_coeff, f), u_coeff, f); 
        p1 = multiply_by_constant(&p1, coeff, f);
        p2 = multiply_by_constant(&p2, coeff, f);
        p3 = multiply_by_constant(&p3, coeff, f);
        u  = add(&u, &p1);
        g1 = add(&g1, &p2);
        h1 = add(&h1, &p3);
        if deg(&u) > u_bound || deg(&g1) > g1_bound {
            if deg(&v) <= u_bound && deg(&g2) <= g1_bound { // check prev 
                break;
            }
        }
        count += 1;
    }
    return (g2, v);
}
// O(t*tm*m)

// Returns the inverse of a modulo g
// Coefficient arithmetic done in GF(2^m) using f
pub fn inverse(a:&Vec<u64>, g:&Vec<u64>, f:u64) -> Vec<u64> {
    let mut u = a.clone();
    let mut v = g.clone();
    let mut g1 = vec![1 as u64];
    let mut g2 = vec![0 as u64];
    let mut count = 0;
    while !is_constant(&u) {
        if count > 100 {
            println!("Error: max iterations in f2 inverse");
            process::exit(1);
        }
        if is_zero(&u) || is_zero(&v) {
            println!("Error in f2: g = {}", string(&g));
            process::exit(1);
        }            
        let mut j:i64 = deg(&u) as i64 - deg(&v) as i64;
        if j<0 as i64 {
            mem::swap(&mut u, &mut v);
            mem::swap(&mut g1, &mut g2);
            j *= -1 as i64;
        }
        let mut p1 = v.clone();
        let mut p2 = g2.clone();
        shift_left_i(&mut p1, j as usize);
        shift_left_i(&mut p2, j as usize);
        let mut u_coeff = u[u.len()-1];
        let mut v_coeff = v[v.len()-1];
        let mut p1_coeff = f1::multiply(f1::inverse(v_coeff, f), u_coeff, f); 
        p1 = multiply_by_constant(&p1, p1_coeff, f);
        p2 = multiply_by_constant(&p2, p1_coeff, f);
        u = add(&u, &p1);
        g1 = add(&g1, &p2);
        count += 1;
    }
    let mut inv = f1::inverse(u[0], f);
    g1 = multiply_by_constant(&g1, inv, f);
    return g1;
}
// O(t*tm*m)


// Reduces polynomial a modulo polynomial g, coefficients reduced modulo f 
pub fn reduce(a:&mut Vec<u64>, g:&Vec<u64>, f:u64) { 
    //println!("in reduce");
    let mut deg_a = deg(&a) as usize;
    let mut deg_g = deg(&g) as usize;
    //println!("deg_a = {}", deg_a);
    //println!("deg_g = {}", deg_g);
    if deg_a < deg_g {
        return;
    }
    let mut rz = g.clone();
    rz.pop();
    //println!("rz = {}", string(&rz));
    shift_left_i(&mut rz, deg_a - deg_g);
    for i in (deg_g..deg_a+1).rev() {
        //println!("in loop: i = {}", i);
        if a[i] != 0 {
            //println!("a[i] != 0");
            let mut to_add = multiply_by_constant(&rz, a[i], f);
            //println!("to_add = {}", string(&to_add));
            //println!("old a = {}", string(&a));
            *a = add(&a, &to_add);
            //println!("new a = {}", string(&a));
            a.pop();
            //println!("after pop: a.len() = {}", a.len());
            //println!("new a = {}", string(&a));
        } else {
            a.pop();
        }
        shift_right(&mut rz);
    }
}
// O(t*tm)

// Returns a^2 modulo g - do not reduce because 
// this operation is used in decoding and you don't 
// reduce it there 
pub fn square(a:&Vec<u64>, g:&Vec<u64>, f:u64) -> Vec<u64> {
    let mut b = Vec::new();
    let mut i = 0;
    while i < a.len() {
        b.push(f1::multiply(a[i],a[i],f));
        b.push(0 as u64);
        i += 1;
    }
    b.pop(); 
    //let mut test = b.clone();
    //reduce(&mut test, &g, f);
    //println!("test = {}", string(&test));
    //reduce(&mut b, &g, f);
    b
}
// O(t*m)

// Returns square root of a modulo g, coefficients reduced modulo f
pub fn sqrt(a:&Vec<u64>, g:&Vec<u64>, f:u64) -> Vec<u64> {  
    //println!("in sqrt");
    let mut t = deg(&g);
    let mut m = f1::deg(f);
    let mut n = t*m-1;
    let mut result = a.clone();
    for i in 0..n {
        result = square(&result, &g, f);
        reduce(&mut result, &g, f);
        //println!("actual = {}", string(&result));
    }
    //println!("end sqrt");
    result
}
// O(tmn)

// Adds two polynomials a and b in GF(2^tm)        
pub fn add(a:&Vec<u64>, b:&Vec<u64>) -> Vec<u64> {
    let mut c = Vec::new();
    let mut i = 0;
    while i < a.len() && i < b.len() { 
        c.push(f1::add(a[i], b[i]));
        i += 1;
    }
    while i < a.len() {
        c.push(a[i]);
        i += 1;
    }
    while i < b.len() {
        c.push(b[i]);
        i += 1;
    }
    while c[c.len()-1] == 0 && c.len() > 1 {
        c.pop();
    }
    c
}
// O(t)


// Plugs x, an element of GF(2^m), into g(x)
// Coefficients are reduced modulo f
pub fn compute(g:&Vec<u64>, f:u64, x:u64) -> u64 {
    let mut sum = 0;
    let mut power = deg(g) as usize;
    while power > 0 {
        //let mut variable = f1::pow(x, power as u64, f);
        //let mut coeff = g[power];
        //let mut prod = f1::multiply(coeff, variable, f);
        sum = f1::add(f1::multiply(g[power], f1::pow(x, power as u64, f), f), sum);
        power -= 1;
    }
    sum = f1::add(sum, g[power]); // add constant term
    sum
}
// O(mt^2 + t)

// ======= helper methods =======================================
// Multiplies polynomial a by constant c and reduces coefficients modulo f
pub fn multiply_by_constant(a:&Vec<u64>, c:u64, f:u64) -> Vec<u64> {
    let mut b = a.clone();
    for i in 0..a.len() {
        b[i] = f1::multiply(b[i], c, f);
    }
    b
}
// O(t*m)

// Returns true if polynomial is the constant 1, false otherwise
pub fn is_one(a:&Vec<u64>) -> bool {
    return a.len() == 1 && a[a.len()-1] == 1;
}

// Returns true is polynomial a is only a constant value
pub fn is_constant(a:&Vec<u64>) -> bool {
    if a.len() == 1 && a[0] != 0 {
        return true;
    }
    return false;
}

// Returns true is polynomial a is zero, false otherwise
pub fn is_zero(a:&Vec<u64>) -> bool {
    for i in 0..a.len() {
        if a[i] != 0 {
            return false;
        } 
    }
    return true;
}

// Multiplies polynomial a by x
pub fn shift_left(a:&mut Vec<u64>) {
    a.insert(0,0 as u64);
}      

// Divides polynomial a by x
pub fn shift_right(a:&mut Vec<u64>) {
    a.remove(0);
}

// Multiplies polynomial a by x^i
pub fn shift_left_i(a:&mut Vec<u64>, i:usize) {
    let mut i = i;
    while i > 0 {
        shift_left(a);
        i -= 1;
    }
}


// Returns degree of polynomial g
pub fn deg(g:&Vec<u64>) -> u64 {
    return g.len() as u64-1;
}



// Returns polynomial a as a string: {GF(2^m)}x^{t-1} + ... + {GF(2^m)}x + {GF(2^m)}
pub fn string(a:&Vec<u64>) -> String {
    let mut result = String::from("");
    let mut deg = deg(&a) as i64;
    for u in a.iter().rev() {
        result.push_str(&format!("x^{}: {:#066b}", deg, u));
        if deg != 0 {
            result.push_str("\n");
        }
        deg -= 1;
    }
    result
}
// ============ polynomial methods ======================================

pub fn get_random_polynomial(t:u64, f:u64) -> Vec<u64> {
    let mut m = f1::deg(f);
    let mut g = vec![1 as u64];
    shift_left_i(&mut g, t as usize);
    for index in 0..t {
        let mut constant = rand::thread_rng().gen_range(1..2_u64.pow(m as u32));
        g[index as usize] = constant;
    }
    g
}

pub fn get_irreducible_polynomial(t:u64, f:u64) -> Vec<u64> {
    let mut g = get_random_polynomial(t, f);
    let mut count = 0;
    while !is_irreducible(&g, f) {
        g = get_random_polynomial(t, f);
        count += 1;
    }
    println!("{} tries to find irreducible g", count);
    g
}

// Returns a primitive (irreducible) polynomial of degree t over GF(2^m)
/*pub fn get_irreducible_polynomial(t:u64, f:u64) -> Vec<u64> {
    let mut m = f1::deg(f);
    let mut g = vec![1 as u64];
    shift_left_i(&mut g, t as usize);
    let mut constant = rand::thread_rng().gen_range(1..2_u64.pow(m as u32));
    g[0] = constant;
    let mut index = rand::thread_rng().gen_range(1..t) as usize; // choose third index
    g[index] = rand::thread_rng().gen_range(1..2_u64.pow(m as u32));
    let mut count = 0;
    while !is_irreducible(&g, f) { // || !is_primitive(&g, f) { // takes too long for deg(g) > 2 
        //println!("was not irreducible: \n{}", string(&g));
        g[0] = 0; // clear constant
        g[index] = 0; // clear index
        constant = rand::thread_rng().gen_range(1..2_u64.pow(m as u32));
        g[0] = constant;
        index = rand::thread_rng().gen_range(1..t) as usize; // choose third index
        g[index] = rand::thread_rng().gen_range(1..2_u64.pow(m as u32));
        count += 1;
        if count > 500 {
            println!("Could not find primitive polynomial over GF(2^{}) of degree {}", m, t);
            process::exit(1);
        }   
    }
    println!("success after {} tries", count);
    g
}*/


// Returns true if g is irreducible over GF(2^m), false otherwise
/*pub fn is_irreducible(g:&Vec<u64>, f:u64) -> bool {
    //return true; 
    /*for i in 0..2_u64.pow(f1::deg(f) as u32) {
        if compute(g, f, i) == 0 {
            return false;
        }
    }*/
    let mut n = deg(&g);
    let mut q = 2_u64.pow(f1::deg(f) as u32); // q = 2^m
    for i in 1..(n as f64/2 as f64).floor() as u64 + 1 {
        let mut exp = q.pow(i as u32);
        let mut rhs = vec![1 as u64]; 
        shift_left_i(&mut rhs, exp as usize); // x^{q^i}
        let mut rhs = add(&mut rhs, &vec![0 as u64, 1 as u64]); // x^{q^i}-x
        let mut gcd = gcd(&g, &rhs, f);
        if !is_one(&gcd) {
            return false;
        } 
    }
    return true;
}*/

/*fn multiply_by_x_mod_gx(b: &Vec<u64>, g: &Vec<u64>, f:u64) -> Vec<u64> {
    let mut b = b.clone();
    let t = deg(g);
    let mut rz = g.clone();
    rz.pop();
    while rz[rz.len()-1] == 0 {
        rz.pop(); // did this wrong in reduce? 
    }
    shift_left(&mut b);   // left shift of b
    if deg(&b) == t { // followed by conditional addition of rz to b
        let mut temp_rz = rz.clone();
        temp_rz = multiply_by_constant(&temp_rz, b[b.len()-1], f);
        b = add(&b, &temp_rz);
        b.pop(); // cancel out high order coeff 
    }
    b
}*/

// Returns true if g is irreducible over GF(2^m), false otherwise
/*pub fn is_irreducible(g:&Vec<u64>, f:u64) -> bool {
    let mut t = deg(&g);
    let mut n = 2_u64.pow(f1::deg(f) as u32); // n = 2^m
    let mut m = f1::deg(f);
    let mut q = 2; 
    let mut rhs = vec![1 as u64]; 
    for i in 1..(t as f64/2 as f64).floor() as u64 + 1 { // for 1 to t/2
        for number in 0..(q*m) {
            rhs = multiply_by_x_mod_gx(&rhs, &g, f);
        } 
        let mut rhs_plus_x = add(&mut rhs, &vec![0 as u64, 1 as u64]);
        // rhs_plus_x is reduced
        let mut gcd = gcd(&g, &rhs_plus_x, f);
        if !is_one(&gcd) {
            return false;
        } 
    }
    return true;
}*/

pub fn is_irreducible(g:&Vec<u64>, f:u64) -> bool {
    //println!("\nis_irreducible: g = \n{}", string(&g));
    let mut t = deg(&g);
    let mut m = f1::deg(f);
    //let mut q = 2;
    let mut exp = 1;
    let mut rhs = vec![0,1]; // rhs = x
    //let limit = (2_u64.pow(m as u32)).pow((t as f64/2 as f64).floor() as u64);
    let t_over_2 = (t as f64/2 as f64).floor() as u64;
    //println!("t = {}", t);
    //println!("t_over_2 = {}", t_over_2);
    let limit = 2_u64.pow((m*t_over_2) as u32);
    //println!("limit = {}", limit);
    let mut j = 0;
    for i in 0..m { // square rhs m times
        rhs = square(&rhs, &g, f);
        reduce(&mut rhs, &g, f);
        exp *= 2;
    }
    j += 1;
    while j <= t_over_2 {
        //println!("j = {}", j);
        let mut rhs_plus_x = add(&mut rhs, &vec![0 as u64, 1 as u64]);
        let mut gcd = gcd(&g, &rhs_plus_x, f);
        if !is_one(&gcd) {
            return false;
        } 
        for i in 0..m { // square rhs m times
            rhs = square(&rhs, &g, f);
            reduce(&mut rhs, &g, f);
            exp *= 2;
        }
        j += 1;
    }
    return true;
}


