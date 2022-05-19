// matrix.rs
// This file contains:
//    - method to create generator matrix for a binary irreducible Goppa code
//      using an irreducible polynomial g 
//    - method to generate the matrices S and P for the private key, 
//      used to scramble the generator matrix
//    - method used to convert matrix of elements of GF(2^m) to binary matrix
//    - method to multiply two matrices using field multiplication from GF(2^m)
// Explanation: 
//    - the generator matrix for the code is found by taking the transpose of the nullspace 
//      of the parity check matrix for the code
//    - the parity check matrix is H = xyz where x,y,and z are matrices formed 
//      using g, f, and L

use crate::GF_2m_arithmetic as f1;
use crate::GF_2mt_arithmetic as f2;
use peroxide::fuga::*;
use std::process;
use rand_chacha::ChaCha20Rng;
use rand_chacha::ChaCha8Rng;

// Returns the generator matrix for the binary irreducible Goppa code corresponding
// to the irreducible polynomial g 
// Returns (generator, height, width)
pub fn generator(g:&Vec<u64>, f:u64, L:&Vec<u64>) -> (Matrix, u64, u64) {
    let t  = f2::deg(g);
    let n  = L.len();
    let m  = f1::deg(f); 
    let x  = x(g, t);
    //println!("x = ");
    //x.print();
    let y  = y(L, t, n as u64, f);
    //println!("y = ");
    //y.print();
    let z  = z(&g, &L, n as u64, f);
    //println!("z = ");
    //z.print();
    let xy = multiply(&x, &y, t as usize, t as usize, t as usize, n, f);
    let h  = multiply(&xy, &z, t as usize, n, n, n, f);
    //println!("H = xyz = ");
    //h.print();
    let mut h  = binary_matrix(&h, t as usize, n, m as usize, f);
    //println!("binary H = ");
    //h.print();
    let (generator, dimension, length) = nullspace_t(&mut h, t*m as u64, n as u64); 
    return (generator, dimension, length);
}   

// Append two matrices together
pub fn append(a:&Matrix, b:&Matrix, h1:u64, w1:u64, h2:u64, w2:u64) -> Matrix {
    let mut m = zeros(h1 as usize, (w1+w2) as usize);
    let mut r1 = 0 as usize;
    let mut c1 = 0 as usize;
    while r1 < h1 as usize {
        c1 = 0;
        while c1 < w1 as usize {
            m[(r1,c1)] = a[(r1,c1)];
            c1 += 1;
        }
        r1 += 1;
    }
    let mut r2 = 0 as usize;
    while r2 < h2 as usize {
        let mut c2 = 0 as usize;
        while c2 < w2 as usize {
            m[(r2, c1+c2)] = b[(r2,c2)];
            c2 += 1;
        }
        r2 += 1;
    }
    m
}


// Returns dense nonsingular matrix S for use in private key
// Make sure determinant is 1 or -1 so inverse contains whole numbers
// Also, if matrix is too large with high density, library cannot find the matrix
// Also Rust library sometimes crashes while finding determinant
/*pub fn dense_nonsingular_matrix(k:u64) -> Matrix {
    let mut rng = ChaCha20Rng::from_entropy();
    //let mut rng = ChaCha8Rng::seed_from_u64(1);
    let size:u64 = k*k;
    let min_density = (size as f64/2 as f64).ceil() as u64 + 1; // 1/2
    let density = min_density; // crashes for higher density or lower density
    let mut det:f64 = 0.0;
    let mut m = zeros(k as usize, k as usize);
    let mut iteration = 0;
    loop {
        if iteration > 100 {
            println!("Error in matrix.rs:: Rust libarary did not find nonsingular matrix");
            process::exit(1);
        } 
        //println!("dense_nonsingular: iteration = {}", iteration);
        m = zeros(k as usize, k as usize);
        let mut count = 0;
        let mut inner_iteration = 0;
        while count < density {
            if inner_iteration > 100 {
                break;
            }
            //println!("count < density");   // gets stuck here ?
            //let c = rand::thread_rng().gen_range(0..k);
            //let r = rand::thread_rng().gen_range(0..k);
            let c = rng.gen_range(0..k);
            let r = rng.gen_range(0..k);
            if m[(c as usize, r as usize)] != 1 as f64 {
                m[(c as usize, r as usize)] = 1 as f64;
                count += 1;
            }
            inner_iteration += 1;
        }
        //if m.det() != 0.0 {
        //    break;
       // }
        if m.det() == 1 as f64 || m.det() == -1 as f64 { // Rust library sometimes crashes here
                break;
        }
        iteration += 1;
    }
    //println!("end dense nonsingular");
    m
}*/

pub fn dense_nonsingular_matrix(k:u64) -> Matrix {
    let mut m = zeros(k as usize, k as usize);
    for i in 0..k {
        for j in 0..k {
            if i == j {
                m[(i as usize, j as usize)] = 1 as f64;
            }
        }
    }
    //m.print();
    let mut density = k;
    let mut min_density = ((k as f64*k as f64)/2 as f64).floor() as u64;
    let mut no_density_increase = 0;
    while density < min_density {
        for row in 0..k {
            let mut first_density = density;
            let mut random_row = rand::thread_rng().gen_range(0..k);
            if row != random_row {
                let mut possible = true;
                for i in 0..k {
                    if m[(row as usize, i as usize)] == 1.0 && m[(random_row as usize, i as usize)] == 1.0 {
                        possible = false;
                        break;
                    }
                }
                if possible {
                    for i in 0..k {
                        let mut before = m[(row as usize, i as usize)];
                        m[(row as usize, i as usize)] = m[(row as usize, i as usize)] + m[(random_row as usize, i as usize)]; 
                        let mut after = m[(row as usize, i as usize)];
                        if after > before {
                            density += 1;
                        }
                    }
                }
            }
            let mut last_density = density;
            if first_density == last_density {
                no_density_increase += 1;
            } else {
                no_density_increase = 0;
            }
        }
        if no_density_increase > k*3 {
            //println!("breaking: no_density_increase = {}", no_density_increase);
            break;
        }
    }
    //m.print();
    //println!("density = {}/{}", density, k*k);
    m
}

// Returns random permutation matrix P for use in private key
pub fn permutation_matrix(n:u64) -> Matrix  {
    let mut ones: Vec<u64> = (0..n).collect();
    ones.shuffle(&mut thread_rng());
    let mut m = zeros(n as usize,n as usize);
    let mut i:usize = 0;
    for i in 0..n {
        m[(i as usize,ones[i as usize] as usize)] = 1 as f64;
    }
    m
}

// Convert a matrix of elements of GF(2^m) to a binary matrix
// by converting each element of GF(2^m) to a vertical binary string where
// the top number is the high order bit
pub fn binary_matrix(h:&Matrix, t:usize, d:usize, m:usize, f:u64) -> Matrix {
    let mut b = zeros(t*m, d);
    for c in 0..d {
        for r in 0..t {
            let mut vec = f1::vec(h[(r,c)] as u64, f); // high order to low order
            for i in 0..m {
                b[(r*m+i, c)] = vec[i];
            }
        }
    }
    b
}


// Multiplies two matrices together
// Since the elements of the matrix are elements of GF(2^m), regular
// multiplication cannot be used and field multiplication in GF(2^m) must be used 
pub fn multiply(x:&Matrix, y:&Matrix, xrows:usize, xcols:usize, yrows:usize, ycols:usize, f:u64) -> Matrix {
    if xcols != yrows {
        println!("matrix multiplication error");
        process::exit(1);
    }
    let mut z = zeros(xrows, ycols);
    for r in 0..xrows {
        for c in 0..ycols {
            for t in 0..xcols {
                z[(r,c)] = f1::add(z[(r,c)] as u64, f1::multiply(x[(r,t)] as u64, y[(t,c)] as u64, f)) as f64;
            }
        }
    }
    z
}             
// O(xrows * ycols * xcols * m)


// Return z in (xyz = parity check matrix)
// It is a diagonal matrix where the elements are 
// the inverse of g(l) modulo g(x) for l in L
pub fn z(g:&Vec<u64>, L:&Vec<u64>, d:u64, f:u64) -> Matrix {
    let mut z = zeros(d as usize, d as usize);
    for r in 0..d {
        for c in 0..d {
            if r == c {
                z[(r as usize, c as usize)] = f1::inverse(f2::compute(&g, f, L[r as usize]), f) as f64;
            }
        }
    }
    z
}
// O(mt^2 * length)

// Returns y in (xyz = parity check matrix)
// It contains the elements of L raised to powers and reduced modulo f
pub fn y(L:&Vec<u64>, t:u64, n:u64, f:u64) -> Matrix {
    let mut L = L.clone();
    let mut y = zeros(t as usize, n as usize);
    for c in 0..n {
        y[(0,c as usize)] = 1 as f64;
    }
    for r in 1..t {
        for c in 0..n {
            y[(r as usize,c as usize)] = f1::pow(L[c as usize],r,f) as f64;
        }
    }
    y
}
// O(mt * length)


// Returns x in (xyz = parity check matrix)
// It is a lower triangular matrix made 
// using the coefficients of irreducible polynomial g(x)
pub fn x(g:&Vec<u64>, t:u64) -> Matrix {
    let mut g = g.clone();
    g.remove(0);
    let mut x = zeros(t as usize, t as usize);
    for r in (0..t).rev() {
        for c in 0..t {
            x[(r as usize, c as usize)] = g[c as usize] as f64;
        }
        g.remove(0);
        g.push(0);
    }
    x
}
// O(t)

// Modified rref, used in last decryption step
pub fn rref_for_decrypt(M:&mut Matrix, m:u64, n:u64, rows_to_include:u64, columns_to_include:u64) {
    let mut lead:usize = 0;
    let mut rowCount:usize = rows_to_include as usize;
    let mut columnCount:usize = columns_to_include as usize;
    for r in 0..rowCount {
        if columnCount <= lead {
            return 
        }
        let mut i = r;
        while M[(i, lead)] == 0 as f64{
            i = i+1;
            if rowCount == i {
                i = r;
                lead = lead + 1;
                if columnCount == lead {
                    return 
                }
            }
        }
        if i != r {
            unsafe {
                M.swap(i, r, Row);
            }
        }
        for i in 0..rowCount {
            if i != r && M[(i, lead)] != 0 as f64 {
                for c in 0..(n as  usize) {
                    M[(i,c)] = f1::add(M[(i,c)] as u64, M[(r,c)] as u64) as f64;    
                }
            }
        }
        lead = lead + 1;
    } 
}

// Puts a matrix into reduced row echelon form modulo 2
// Modified from wikipedia rref
pub fn rref_mod_2(M:&mut Matrix, m:u64, n:u64) {
    let mut lead:usize = 0;
    let mut rowCount:usize = m as usize;
    let mut columnCount:usize = n as usize;
    for r in 0..rowCount {
        if columnCount <= lead {
            return 
        }
        let mut i = r;
        while M[(i, lead)] == 0 as f64{
            i = i+1;
            if rowCount == i {
                i = r;
                lead = lead + 1;
                if columnCount == lead {
                    return 
                }
            }
        }
        if i != r {
            unsafe {
                M.swap(i, r, Row);
            }
        }
        for i in 0..rowCount {
            if i != r && M[(i, lead)] != 0 as f64 {
                // add row r to row i
                for c in 0..columnCount {
                    M[(i,c)] = f1::add(M[(i,c)] as u64, M[(r,c)] as u64) as f64;    
                }
            }
        }
        lead = lead + 1;
    } 
}
// O(tm by length)


// Reduces an element of a matrix modulo 2
pub fn mod2(a:f64) -> f64 {
    let result = (a.round() as i64).rem_euclid(2) as f64;
    result
}

// Reduces all the elements of a matrix modulo 2
pub fn mod2_matrix(m:&Matrix, rows:usize, cols:usize) -> Matrix {
    let mut result = m.clone();
    for row in 0..rows {
        for col in 0..cols {
            result[(row, col)] = mod2(result[(row,col)]);
        }
    }
    result
}

// Returns true if row is all zeros, false otherwise
pub fn row_is_zero(M:&Matrix, i:u64) -> bool {
    let mut row = M.row(i as usize);
    for e in row {
        if e != 0 as f64 {
            return false;
        }
    }
    return true;
}
// Check if a column is all zeros except a 1 at position col
pub fn is_col(vec:&Vec<f64>, col:usize) -> bool {
    for i in 0..vec.len() {
        if i == col {
            if vec[i] != 1  as f64 {
                return false;
            }
        } else {
            if vec[i] != 0 as f64 {
                return false;
            }
        }
    }
    return true;
}

// Find column to swap with in swap_columns()
pub fn find_col(M:&Matrix, rows:usize, cols:usize, i:usize) -> (bool, usize) {
    for c in 0..cols {
        let mut vec = M.col(c);
        if is_col(&vec, i) {
            return (true,c );
        }
    }
    return (false,0);
}

// Function to swap columns in rref 
/*pub fn swap_columns(M:&mut Matrix, rows:u64, cols:u64) {
    for row in 0..rows {
        for col in 0..cols {
            if row == col && M[(row as usize, col as usize)] != 1 as f64 {
                let (mut is_col, index) = find_col(&M, rows as usize, cols as usize, col as usize);
                if is_col {
                    unsafe {
                        M.swap(row as usize, index as usize, Col);
                    }
                } else {
                }
            }
        }
    }
}*/

// Returns the transpose of the nullspace of matrix M along with its height, width
// Arithmetic done modulo 2
pub fn nullspace_t(M:&mut Matrix, rows:u64, cols:u64) -> (Matrix, u64, u64) {  
    let mut swap = Vec::new();
    rref_mod_2(M, rows, cols);
    // swap columns
    for row in 0..rows {
        for col in 0..cols {
            if row == col && M[(row as usize, col as usize)] != 1 as f64 {
                let (mut is_col, index) = find_col(&M, rows as usize, cols as usize, col as usize);
                if is_col {
                    unsafe {
                        M.swap(row as usize, index as usize, Col);
                        swap.push((row as usize, index as usize));
                    }
                } else {
                }
            }
        }
    }
    // end swap columns
    //println!("H = ");
    //M.print();
    // find where identity matrix ends
    let mut last_nonzero_row = rows-1;
    for row in (0..rows).rev() {
        if !row_is_zero(&M, row) {
            last_nonzero_row = row;
            break;
        }
    }
    let mut first_column = last_nonzero_row+1;
    // if rref form is an identity matrix, there is no nullspace
    if (cols - first_column) < 1 {
        println!("Error: Nullspace is zero or too small");
        process::exit(1);
    }   
    // nullspace has one row for each column and one column for each column after the identity matrix
    let mut nullspace = zeros(cols as usize, cols as usize-first_column as usize);
    for row in 0..rows {
        for col in 0..(cols as usize - first_column as usize) {
            nullspace[(row as usize, col as usize)] = M[(row as usize, col as usize +first_column as usize)];
        }
    }
    // fill in identity matrix
    let mut row = first_column as usize;
    let mut col = 0 as usize;
    while col < (cols as usize - first_column as usize) {
        nullspace[(row, col)] = 1 as f64;
        col += 1;
        row += 1;
    }
    // swap rows back in reverse order
    if swap.len() != 0 {
        for &(r1, r2) in swap.iter().rev() {
            unsafe {
                nullspace.swap(r1,r2,Row);
            }
        }
    }
    // end swap rows back
    //println!("nullspace(H) =");
    //nullspace.print();
    let mut nullspace = nullspace.transpose();
    return (nullspace, cols as u64 - first_column as u64, cols as u64);
}
