// util.rs

// convert string to vector of 0 and 1 of type u64
pub fn string_to_vec_of_u64s(s:String) -> Vec<u64> {
    let mut string_of_1_and_0 = String::from("");
    for character in s.clone().into_bytes() { // character is number in ascii
        let mut character_as_binary_string = format!("{:b}", character);
        while character_as_binary_string.len() < 8 {
            let mut new_str = String::from("0");
            new_str.push_str(&character_as_binary_string);
            character_as_binary_string = new_str;
        }
        string_of_1_and_0.push_str(&character_as_binary_string);
    }
    let mut vec_of_u64s = Vec::new();
    for ch in string_of_1_and_0.chars() {
        if ch == '1' {
            vec_of_u64s.push(1);
        } else if ch == '0' {
            vec_of_u64s.push(0);
        }
    }
    if s.len()*8 != vec_of_u64s.len() {
        use std::process;
        println!("Error in util.rs: translating string to vector did not work");
        process::exit(1);
    }
    vec_of_u64s
}

pub fn vec_of_u64s_to_string(v:&Vec<u64>) -> String {
    let mut vec_of_chars = Vec::new();
    let mut binary_string = String::from("");
    for i in 0..v.len() {
        if v[i] == 0 {
            binary_string += "0";
        } else if v[i] == 1 {
            binary_string += "1";
        }
        if binary_string.len() == 8 {
            vec_of_chars.push(binary_string.clone());
            binary_string = String::from("");
        }
    }
    let mut characters = Vec::new();
    for i in 0..vec_of_chars.len() {
        let bin = vec_of_chars[i].clone();
        characters.push(isize::from_str_radix(&bin, 2).unwrap() as u8);
    }
    while characters[characters.len()-1] == 0 { // get rid of nulls 
        characters.pop();
    }
    let mut characters = &characters[..];
    let st = std::str::from_utf8(&characters).unwrap();
    st.to_string()
}