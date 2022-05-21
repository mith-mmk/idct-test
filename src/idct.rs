
use core::f32::consts::PI;
use core::f64::consts::PI as PI_f64;

pub fn idct_f64(f :&[i32]) -> Vec<u8> {
    let vals :Vec<u8> = (0..64).map(|i| {
        let (x,y) = ((i%8) as f64,(i/8) as f64);
        // IDCT from CCITT Rec. T.81 (1992 E) p.27 A3.3
        let mut val: f64=0.0;
        for u in 0..8 {
            let cu = if u == 0 {1.0 / 2.0_f64.sqrt()} else {1.0};
            for v in 0..8 {
                let cv = if v == 0 {1.0_f64/ 2.0_64_f64.sqrt()} else {1.0};
                val += cu * cv * (f[v*8 + u] as f64)
                    * ((2.0 * x + 1.0_64) * u as f64 * PI_f64 / 16.0_f64).cos()
                    * ((2.0 * y + 1.0_64) * v as f64 * PI_f64 / 16.0_f64).cos();
            }
        }
        val = val / 4.0;

        // level shift from CCITT Rec. T.81 (1992 E) p.26 A3.1
        let v = val.round() as i32 + 128;
        v.clamp(0,255) as u8
    }).collect();
    vals
}

pub fn idct(f :&[i32]) -> Vec<u8> {
    let vals :Vec<u8> = (0..64).map(|i| {
        let (x,y) = ((i%8) as f32,(i/8) as f32);
        // IDCT from CCITT Rec. T.81 (1992 E) p.27 A3.3
        let mut val: f32=0.0;
        for u in 0..8 {
            let cu = if u == 0 {1.0 / 2.0_f32.sqrt()} else {1.0};
            for v in 0..8 {
                let cv = if v == 0 {1.0_f32 / 2.0_f32.sqrt()} else {1.0};
                val += cu * cv * (f[v*8 + u] as f32)
                    * ((2.0 * x + 1.0) * u as f32 * PI / 16.0_f32).cos()
                    * ((2.0 * y + 1.0) * v as f32 * PI / 16.0_f32).cos();
            }
        }
        val = val / 4.0;

        // level shift from CCITT Rec. T.81 (1992 E) p.26 A3.1
        let v = val.round() as i32 + 128;
        v.clamp(0,255) as u8
    }).collect();
    vals
}

// LL&M 2D IDCT is 1D IDCT * 1D IDCT 
pub fn llm_idct(f: &[i32]) -> Vec<u8> {
    let m1 = 0.5411961;   // α √2cos(3π/8)
    let m2 = 1.306562965; // β √2cos(3π/8)
    let m3 = 1.414213562; // γ v2
    let m4 = 0.831469612; // η cos(3π/16)
    let m5 = 0.555570233; // θ sin(3π/16)
    let m6 = 0.98078528;  // δ cos(π/16)
    let m7 = 0.195090322; // ε sin(π/16)
    let m0 = 0.125; // √2/4 * √2/4

    let mut ff = [0_f32;64];
    for j in 0..8 {
        let i = j * 8;
        let f0 = f[0 + i] as f32;  // X0
        let f1 = f[1 + i] as f32;  // X1
        let f2 = f[2 + i] as f32;  // X2
        let f3 = f[3 + i] as f32;  // X3
        let f4 = f[4 + i] as f32;  // X4
        let f6 = f[6 + i] as f32;  // X5
        let f5 = f[5 + i] as f32;  // X6
        let f7 = f[7 + i] as f32;  // X7

        // implement batterfly mutilply

        // even part
        // part 2

        let y0 = f0 + f4;
        let y1 = f0 - f4;
        let y2 = m1 * f2 - m2 * f6;
        let y3 = m1 * f6 + m2 * f2;

        // part3

        let x0 = y0 + y3;
        let x1 = y1 + y2;
        let x2 = y1 - y2;
        let x3 = y0 - y3;

        // odds part

        // part 1
        let z4 = f1 - f7;
        let z5 = f3 * m3;
        let z6 = f5 * m3;
        let z7 = f1 + f7;

        // part 2
        let y4 = z4 + z6;
        let y5 = z7 - z5;
        let y6 = z4 - z6;
        let y7 = z7 + z5;

        // part 3
        let x4 = y4 * m4 - y7 * m5;
        let x5 = y5 * m6 - y6 * m7;
        let x6 = y6 * m6 + y5 * m7;
        let x7 = y7 * m4 + y4 * m5;

        // last part  multiply √2 / 4 after parts

        ff[0 + i] = x0 + x7;   // x0
        ff[7 + i] = x0 - x7;   // x1
        ff[1 + i] = x1 + x6;   // x2
        ff[6 + i] = x1 - x6;   // x3
        ff[2 + i] = x2 + x5;   // x4
        ff[5 + i] = x2 - x5;   // x5
        ff[3 + i] = x3 + x4;   // x6
        ff[4 + i] = x3 - x4;   // x7
    }
    for i in 0..8 {
        let f0 = ff[0 * 8 + i];
        let f1 = ff[1 * 8 + i];
        let f2 = ff[2 * 8 + i];
        let f3 = ff[3 * 8 + i];
        let f4 = ff[4 * 8 + i];
        let f5 = ff[5 * 8 + i];
        let f6 = ff[6 * 8 + i];
        let f7 = ff[7 * 8 + i];

        // odds part
        // part 1

        // even part
        // part 2

        let y0 = f0 + f4;
        let y1 = f0 - f4;
        let y2 = m1 * f2 - m2 * f6;
        let y3 = m1 * f6 + m2 * f2;

        // part3

        let x0 = y0 + y3;
        let x1 = y1 + y2;
        let x2 = y1 - y2;
        let x3 = y0 - y3;

        // odds part
        // part 1
        let z4 = f1 - f7;
        let z5 = f3 * m3;
        let z6 = f5 * m3;
        let z7 = f1 + f7;

        // part 2
        let y4 = z4 + z6;
        let y5 = z7 - z5;
        let y6 = z4 - z6;
        let y7 = z7 + z5;


        let x4 = y4 * m4 - y7 * m5;
        let x5 = y5 * m6 - y6 * m7;
        let x6 = y6 * m6 + y5 * m7;
        let x7 = y7 * m4 + y4 * m5;

        ff[0 * 8 + i] = (x0 + x7) * m0;
        ff[7 * 8 + i] = (x0 - x7) * m0;
        ff[1 * 8 + i] = (x1 + x6) * m0; 
        ff[6 * 8 + i] = (x1 - x6) * m0;
        ff[2 * 8 + i] = (x2 + x5) * m0;
        ff[5 * 8 + i] = (x2 - x5) * m0;
        ff[3 * 8 + i] = (x3 + x4) * m0; 
        ff[4 * 8 + i] = (x3 - x4) * m0;  
    }
    
    let val = ff.iter().map(|i| ((*i + 128.5) as i32).clamp(0,255) as u8).collect();
    val
}
// AAN
pub fn fast_idct(f: &[i32]) -> Vec<u8> {
    let mut _f  = [0_f32;64];
    let mut vals = [0_u8;64];
    let m0 = 1.847759;
    let m1 = 1.4142135;
    let m3 = 1.4142135;
    let m5 = 0.76536685;
    let m2 = m0 - m5;
    let m4 = m0 + m5;

    let s0 = 0.35355338;
    let s1 = 0.49039263;
    let s2 = 0.46193975;
    let s3 = 0.4157348;
    let s4 = 0.35355338;
    let s5 = 0.2777851;
    let s6 = 0.19134171;
    let s7 = 0.09754512;
    
    for i in 0..8 {
        let g0 = f[0*8 + i] as f32 * s0;
        let g1 = f[4*8 + i] as f32 * s4;
        let g2 = f[2*8 + i] as f32 * s2;
        let g3 = f[6*8 + i] as f32 * s6;
        let g4 = f[5*8 + i] as f32 * s5;
        let g5 = f[1*8 + i] as f32 * s1;
        let g6 = f[7*8 + i] as f32 * s7;
        let g7 = f[3*8 + i] as f32 * s3;
    
        let f0 = g0;
        let f1 = g1;
        let f2 = g2;
        let f3 = g3;
        let f4 = g4 - g7;
        let f5 = g5 + g6;
        let f6 = g5 - g6;
        let f7 = g4 + g7;
    
        let e0 = f0;
        let e1 = f1;
        let e2 = f2 - f3;
        let e3 = f2 + f3;
        let e4 = f4;
        let e5 = f5 - f7;
        let e6 = f6;
        let e7 = f5 + f7;
        let e8 = f4 + f6;
    
        let d0 = e0;
        let d1 = e1;
        let d2 = e2 * m1;
        let d3 = e3;
        let d4 = e4 * m2;
        let d5 = e5 * m3;
        let d6 = e6 * m4;
        let d7 = e7;
        let d8 = e8 * m5;
    
        let c0 = d0 + d1;
        let c1 = d0 - d1;
        let c2 = d2 - d3;
        let c3 = d3;
        let c4 = d4 + d8;
        let c5 = d5 + d7;
        let c6 = d6 - d8;
        let c7 = d7;
        let c8 = c5 - c6;
    
        let b0 = c0 + c3;
        let b1 = c1 + c2;
        let b2 = c1 - c2;
        let b3 = c0 - c3;
        let b4 = c4 - c8;
        let b5 = c8;
        let b6 = c6 - c7;
        let b7 = c7;
        
        _f[0 * 8 + i] = b0 + b7;
        _f[1 * 8 + i] = b1 + b6;
        _f[2 * 8 + i] = b2 + b5;
        _f[3 * 8 + i] = b3 + b4;
        _f[4 * 8 + i] = b3 - b4;
        _f[5 * 8 + i] = b2 - b5;
        _f[6 * 8 + i] = b1 - b6;
        _f[7 * 8 + i] = b0 - b7; 
    }
    
    for i in 0..8 {
        let g0 = _f[i*8 + 0] as f32 * s0;
        let g1 = _f[i*8 + 4] as f32 * s4;
        let g2 = _f[i*8 + 2] as f32 * s2;
        let g3 = _f[i*8 + 6] as f32 * s6;
        let g4 = _f[i*8 + 5] as f32 * s5;
        let g5 = _f[i*8 + 1] as f32 * s1;
        let g6 = _f[i*8 + 7] as f32 * s7;
        let g7 = _f[i*8 + 3] as f32 * s3;
    
        let f0 = g0;
        let f1 = g1;
        let f2 = g2;
        let f3 = g3;
        let f4 = g4 - g7;
        let f5 = g5 + g6;
        let f6 = g5 - g6;
        let f7 = g4 + g7;
    
        let e0 = f0;
        let e1 = f1;
        let e2 = f2 - f3;
        let e3 = f2 + f3;
        let e4 = f4;
        let e5 = f5 - f7;
        let e6 = f6;
        let e7 = f5 + f7;
        let e8 = f4 + f6;
    
        let d0 = e0;
        let d1 = e1;
        let d2 = e2 * m1;
        let d3 = e3;
        let d4 = e4 * m2;
        let d5 = e5 * m3;
        let d6 = e6 * m4;
        let d7 = e7;
        let d8 = e8 * m5;
    
        let c0 = d0 + d1;
        let c1 = d0 - d1;
        let c2 = d2 - d3;
        let c3 = d3;
        let c4 = d4 + d8;
        let c5 = d5 + d7;
        let c6 = d6 - d8;
        let c7 = d7;
        let c8 = c5 - c6;
    
        let b0 = c0 + c3;
        let b1 = c1 + c2;
        let b2 = c1 - c2;
        let b3 = c0 - c3;
        let b4 = c4 - c8;
        let b5 = c8;
        let b6 = c6 - c7;
        let b7 = c7;
        
        vals[i * 8 + 0] = ((b0 + b7 + 128.5) as i32).clamp(0,255) as u8;
        vals[i * 8 + 1] = ((b1 + b6 + 128.5) as i32).clamp(0,255) as u8;
        vals[i * 8 + 2] = ((b2 + b5 + 128.5) as i32).clamp(0,255) as u8;
        vals[i * 8 + 3] = ((b3 + b4 + 128.5) as i32).clamp(0,255) as u8;
        vals[i * 8 + 4] = ((b3 - b4 + 128.5) as i32).clamp(0,255) as u8;
        vals[i * 8 + 5] = ((b2 - b5 + 128.5) as i32).clamp(0,255) as u8;
        vals[i * 8 + 6] = ((b1 - b6 + 128.5) as i32).clamp(0,255) as u8;
        vals[i * 8 + 7] = ((b0 - b7 + 128.5) as i32).clamp(0,255) as u8;
    }
    vals.to_vec()
}

// AAN
pub fn fast_idct_f64(f: &[i32]) -> Vec<u8> {
    let mut _f  = [0_f64;64];
    let mut vals = [0_u8;64];
    let m0 = 1.847759_f64;
    let m1 = 1.4142135_f64;
    let m3 = 1.4142135_f64;
    let m5 = 0.76536685_f64;
    let m2 = m0 - m5;
    let m4 = m0 + m5;

    let s0 = 0.35355338_f64;
    let s1 = 0.49039263_f64;
    let s2 = 0.46193975_f64;
    let s3 = 0.4157348_f64;
    let s4 = 0.35355338_f64;
    let s5 = 0.2777851_f64;
    let s6 = 0.19134171_f64;
    let s7 = 0.09754512_f64;
    
    for i in 0..8 {
        let g0 = f[0*8 + i] as f64 * s0;
        let g1 = f[4*8 + i] as f64 * s4;
        let g2 = f[2*8 + i] as f64 * s2;
        let g3 = f[6*8 + i] as f64 * s6;
        let g4 = f[5*8 + i] as f64 * s5;
        let g5 = f[1*8 + i] as f64 * s1;
        let g6 = f[7*8 + i] as f64 * s7;
        let g7 = f[3*8 + i] as f64 * s3;
    
        let f0 = g0;
        let f1 = g1;
        let f2 = g2;
        let f3 = g3;
        let f4 = g4 - g7;
        let f5 = g5 + g6;
        let f6 = g5 - g6;
        let f7 = g4 + g7;
    
        let e0 = f0;
        let e1 = f1;
        let e2 = f2 - f3;
        let e3 = f2 + f3;
        let e4 = f4;
        let e5 = f5 - f7;
        let e6 = f6;
        let e7 = f5 + f7;
        let e8 = f4 + f6;
    
        let d0 = e0;
        let d1 = e1;
        let d2 = e2 * m1;
        let d3 = e3;
        let d4 = e4 * m2;
        let d5 = e5 * m3;
        let d6 = e6 * m4;
        let d7 = e7;
        let d8 = e8 * m5;
    
        let c0 = d0 + d1;
        let c1 = d0 - d1;
        let c2 = d2 - d3;
        let c3 = d3;
        let c4 = d4 + d8;
        let c5 = d5 + d7;
        let c6 = d6 - d8;
        let c7 = d7;
        let c8 = c5 - c6;
    
        let b0 = c0 + c3;
        let b1 = c1 + c2;
        let b2 = c1 - c2;
        let b3 = c0 - c3;
        let b4 = c4 - c8;
        let b5 = c8;
        let b6 = c6 - c7;
        let b7 = c7;
        
        _f[0 * 8 + i] = b0 + b7;
        _f[1 * 8 + i] = b1 + b6;
        _f[2 * 8 + i] = b2 + b5;
        _f[3 * 8 + i] = b3 + b4;
        _f[4 * 8 + i] = b3 - b4;
        _f[5 * 8 + i] = b2 - b5;
        _f[6 * 8 + i] = b1 - b6;
        _f[7 * 8 + i] = b0 - b7; 
    }
    
    for i in 0..8 {
        let g0 = _f[i*8 + 0] as f64 * s0;
        let g1 = _f[i*8 + 4] as f64 * s4;
        let g2 = _f[i*8 + 2] as f64 * s2;
        let g3 = _f[i*8 + 6] as f64 * s6;
        let g4 = _f[i*8 + 5] as f64 * s5;
        let g5 = _f[i*8 + 1] as f64 * s1;
        let g6 = _f[i*8 + 7] as f64 * s7;
        let g7 = _f[i*8 + 3] as f64 * s3;
    
        let f0 = g0;
        let f1 = g1;
        let f2 = g2;
        let f3 = g3;
        let f4 = g4 - g7;
        let f5 = g5 + g6;
        let f6 = g5 - g6;
        let f7 = g4 + g7;
    
        let e0 = f0;
        let e1 = f1;
        let e2 = f2 - f3;
        let e3 = f2 + f3;
        let e4 = f4;
        let e5 = f5 - f7;
        let e6 = f6;
        let e7 = f5 + f7;
        let e8 = f4 + f6;
    
        let d0 = e0;
        let d1 = e1;
        let d2 = e2 * m1;
        let d3 = e3;
        let d4 = e4 * m2;
        let d5 = e5 * m3;
        let d6 = e6 * m4;
        let d7 = e7;
        let d8 = e8 * m5;
    
        let c0 = d0 + d1;
        let c1 = d0 - d1;
        let c2 = d2 - d3;
        let c3 = d3;
        let c4 = d4 + d8;
        let c5 = d5 + d7;
        let c6 = d6 - d8;
        let c7 = d7;
        let c8 = c5 - c6;
    
        let b0 = c0 + c3;
        let b1 = c1 + c2;
        let b2 = c1 - c2;
        let b3 = c0 - c3;
        let b4 = c4 - c8;
        let b5 = c8;
        let b6 = c6 - c7;
        let b7 = c7;
        
        vals[i * 8 + 0] = ((b0 + b7 + 128.5) as i32).clamp(0,255) as u8;
        vals[i * 8 + 1] = ((b1 + b6 + 128.5) as i32).clamp(0,255) as u8;
        vals[i * 8 + 2] = ((b2 + b5 + 128.5) as i32).clamp(0,255) as u8;
        vals[i * 8 + 3] = ((b3 + b4 + 128.5) as i32).clamp(0,255) as u8;
        vals[i * 8 + 4] = ((b3 - b4 + 128.5) as i32).clamp(0,255) as u8;
        vals[i * 8 + 5] = ((b2 - b5 + 128.5) as i32).clamp(0,255) as u8;
        vals[i * 8 + 6] = ((b1 - b6 + 128.5) as i32).clamp(0,255) as u8;
        vals[i * 8 + 7] = ((b0 - b7 + 128.5) as i32).clamp(0,255) as u8;
    }
    vals.to_vec()
}


pub fn print_vec(f:&[u8]) -> String {
    let mut str = "".to_string();
    for i in 0..8{
        str += &format!("{:3} {:3} {:3} {:3} {:3} {:3} {:3} {:3}\n",
            f[i*8],f[i*8+1],f[i*8+2],f[i*8+3],f[i*8+4],f[i*8+5],f[i*8+6],f[i*8+7]
        );
    }
    str
}

// method1 https://note.com/mith_mmk/n/n6f57f007453b
fn idct1(f :&[i32]) -> Vec<u8> {
    let vals :Vec<u8> = (0..64).map(|i| {
        let (x,y) = ((i%8) as f32,(i/8) as f32);
        // IDCT from CCITT Rec. T.81 (1992 E) p.27 A3.3
        let mut val: f32=0.0;
        for u in 0..8 {
            let cu = if u == 0 {1.0 / 2.0_f32.sqrt()} else {1.0};
            for v in 0..8 {
                let cv = if v == 0 {1.0_f32 / 2.0_f32.sqrt()} else {1.0};
                val += cu * cv * (f[v*8 + u] as f32)
                    * ((2.0 * x + 1.0) * u as f32 * PI / 16.0_f32).cos()
                    * ((2.0 * y + 1.0) * v as f32 * PI / 16.0_f32).cos();
            }
        }
        val = val / 4.0;

        // level shift from CCITT Rec. T.81 (1992 E) p.26 A3.1
        let v = val.round() as i32 + 128;
        if v < 0 {0} else if v > 255 {255} else {v as u8}
    }).collect();
    vals
}

// method2 https://note.com/mith_mmk/n/n6f57f007453b
// use cos tables equal idct
/*
fn idct2 (f :&[i32]) -> Vec<u8> {
    let cos_table :[[f32;8];8] = 
       [[ 1.0,  0.98078528,  0.92387953,  0.83146961,  0.70710678, 0.55557023,  0.38268343,  0.19509032],
        [ 1.0,  0.83146961,  0.38268343, -0.19509032, -0.70710678,  -0.98078528, -0.92387953, -0.5555702],
        [ 1.0,  0.55557023, -0.38268343, -0.98078528, -0.70710678, 0.19509032,  0.92387953,  0.83146961],
        [ 1.0,  0.19509032, -0.92387953, -0.55557023,  0.70710678, 0.83146961, -0.38268343, -0.98078528],
        [ 1.0, -0.19509032, -0.92387953,  0.55557023,  0.70710678, -0.83146961, -0.38268343,  0.98078528],
        [ 1.0, -0.55557023, -0.38268343,  0.98078528, -0.70710678,  -0.19509032,  0.92387953, -0.83146961],
        [ 1.0, -0.83146961,  0.38268343,  0.19509032, -0.70710678,  0.98078528, -0.92387953,  0.55557023],
        [ 1.0, -0.98078528,  0.92387953, -0.83146961,  0.70710678, -0.55557023,  0.38268343, -0.19509032]];
    let vals :Vec<u8> = (0..64).map(|i| {
    let (x,y) = ((i%8) as usize,(i/8) as usize);
        // IDCT from CCITT Rec. T.81 (1992 E) p.27 A3.3
        let mut val = 0.0;
        for u in 0..8 {
            for v in 0..8 {
                let cucv :f32 = if u == 0 && v ==0 {0.5} 
                        else if  u==0 || v==0 {1.0 / 2.0_f32.sqrt()}
                        else {1.0};
                val += cucv * f[v*8 + u] as f32 * cos_table[x][u] * cos_table[y][v];
            }
        }
        val = val / 4.0;

        // level shift from CCITT Rec. T.81 (1992 E) p.26 A3.1
        let v = val.round() as isize + 128 ;
        if v < 0 {0} else if v > 255 {255} else {v as u8}
    }).collect();
    vals
}
*/

// method3 https://note.com/mith_mmk/n/n6f57f007453b
// normalize matrix multiply
fn idct3 (f :&[i32]) -> Vec<u8> {
    let c_table :[[f32;8];8] = 
    [[ 0.70710678,  0.98078528,  0.92387953,  0.83146961,  0.70710678, 0.55557023,  0.38268343,  0.19509032],
    [ 0.70710678,  0.83146961,  0.38268343, -0.19509032, -0.70710678, -0.98078528, -0.92387953, -0.55557023],
    [ 0.70710678,  0.55557023, -0.38268343, -0.98078528, -0.70710678, 0.19509032,  0.92387953,  0.83146961],
    [ 0.70710678,  0.19509032, -0.92387953, -0.55557023,  0.70710678, 0.83146961, -0.38268343, -0.98078528],
    [ 0.70710678, -0.19509032, -0.92387953,  0.55557023,  0.70710678, -0.83146961, -0.38268343,  0.98078528],
    [ 0.70710678, -0.55557023, -0.38268343,  0.98078528, -0.70710678, -0.19509032,  0.92387953, -0.83146961],
    [ 0.70710678, -0.83146961,  0.38268343,  0.19509032, -0.70710678, 0.98078528, -0.92387953,  0.55557023],
    [ 0.70710678, -0.98078528,  0.92387953, -0.83146961,  0.70710678, -0.55557023,  0.38268343, -0.19509032]];
    let vals :Vec<u8> = (0..64).map(|i| {
    let (x,y) = ((i%8) as usize,(i/8) as usize);
        // IDCT from CCITT Rec. T.81 (1992 E) p.27 A3.3
        let mut val = 0.0;
        for u in 0..8 {
            for v in 0..8 {
                val += f[v*8 + u] as f32 * c_table[x][u] * c_table[y][v];
            }
        }
        val = val / 4.0;

        // level shift from CCITT Rec. T.81 (1992 E) p.26 A3.1
        let v = val.round() as isize + 128 ;
        if v < 0 {0} else if v > 255 {255} else {v as u8}
    }).collect();
    vals
}

// method4 https://note.com/mith_mmk/n/n6f57f007453b
fn idct4 (f :&[i32]) -> Vec<u8> {
    let c_table :[[f32;8];8] = 
    [[ 0.70710678,  0.98078528,  0.92387953,  0.83146961,  0.70710678, 0.55557023,  0.38268343,  0.19509032],
    [ 0.70710678,  0.83146961,  0.38268343, -0.19509032, -0.70710678, -0.98078528, -0.92387953, -0.55557023],
    [ 0.70710678,  0.55557023, -0.38268343, -0.98078528, -0.70710678, 0.19509032,  0.92387953,  0.83146961],
    [ 0.70710678,  0.19509032, -0.92387953, -0.55557023,  0.70710678, 0.83146961, -0.38268343, -0.98078528],
    [ 0.70710678, -0.19509032, -0.92387953,  0.55557023,  0.70710678, -0.83146961, -0.38268343,  0.98078528],
    [ 0.70710678, -0.55557023, -0.38268343,  0.98078528, -0.70710678, -0.19509032,  0.92387953, -0.83146961],
    [ 0.70710678, -0.83146961,  0.38268343,  0.19509032, -0.70710678, 0.98078528, -0.92387953,  0.55557023],
    [ 0.70710678, -0.98078528,  0.92387953, -0.83146961,  0.70710678, -0.55557023,  0.38268343, -0.19509032]];
    let vals :Vec<u8> = (0..64).map(|i| {
    let (x,y) = ((i%8) as usize,(i/8) as usize);
        // IDCT from CCITT Rec. T.81 (1992 E) p.27 A3.3
        let mut val = 0.0;
        for u in 0..8 {
            let mut uval = 0.0;
            for v in 0..8 {
                uval += f[v*8 + u] as f32 * c_table[y][v];
            }
            val += uval * c_table[x][u];
        }
        val = val / 4.0;

        // level shift from CCITT Rec. T.81 (1992 E) p.26 A3.1
        let v = val.round() as isize + 128 ;
        if v < 0 {0} else if v > 255 {255} else {v as u8}
    }).collect();
    vals
}

// method5 https://note.com/mith_mmk/n/n6f57f007453b
// fixed number
fn idct5 (f :&[i32]) -> Vec<u8> {
    let c_table :[[i32;8];8] = // * 256 << 8 
       [[ 181,  251,  237,  213,  181,  142,   98,   50],
        [ 181,  213,   98,  -50, -181, -251, -237, -142],
        [ 181,  142,  -98, -251, -181,   50,  237,  213],
        [ 181,   50, -237, -142,  181,  213,  -98, -251],
        [ 181,  -50, -237,  142,  181, -213,  -98,  251],
        [ 181, -142,  -98,  251, -181,  -50,  237, -213],
        [ 181, -213,   98,   50, -181,  251, -237,  142],
        [ 181, -251,  237, -213,  181, -142,   98,  -50]];
    let vals :Vec<u8> = (0..64).map(|i| {
    let (x,y) = ((i%8) as usize,(i/8) as usize);
        // IDCT from CCITT Rec. T.81 (1992 E) p.27 A3.3
        let mut val = 0;
        for u in 0..8 {
            let mut uval = 0;
            for v in 0..8 {
                uval += f[v*8 + u]  * c_table[y][v];
            }
            val += uval * c_table[x][u];
        }
        val = val / (4 * 256);

        // level shift from CCITT Rec. T.81 (1992 E) p.26 A3.1
        let v = val as isize + 128 ;
        if v < 0 {0} else if v > 255 {255} else {v as u8}
    }).collect();
    vals
}

// method6 https://note.com/mith_mmk/n/n6f57f007453b
// extend loop
fn idct6 (f :&[i32]) -> Vec<u8> {
    let c_table :[[f32;8];8] = 
    [[ 0.70710678,  0.98078528,  0.92387953,  0.83146961,  0.70710678, 0.55557023,  0.38268343,  0.19509032],
    [ 0.70710678,  0.83146961,  0.38268343, -0.19509032, -0.70710678, -0.98078528, -0.92387953, -0.55557023],
    [ 0.70710678,  0.55557023, -0.38268343, -0.98078528, -0.70710678, 0.19509032,  0.92387953,  0.83146961],
    [ 0.70710678,  0.19509032, -0.92387953, -0.55557023,  0.70710678, 0.83146961, -0.38268343, -0.98078528],
    [ 0.70710678, -0.19509032, -0.92387953,  0.55557023,  0.70710678, -0.83146961, -0.38268343,  0.98078528],
    [ 0.70710678, -0.55557023, -0.38268343,  0.98078528, -0.70710678, -0.19509032,  0.92387953, -0.83146961],
    [ 0.70710678, -0.83146961,  0.38268343,  0.19509032, -0.70710678, 0.98078528, -0.92387953,  0.55557023],
    [ 0.70710678, -0.98078528,  0.92387953, -0.83146961,  0.70710678, -0.55557023,  0.38268343, -0.19509032]];
    let vals :Vec<u8> = (0..64).map(|i| {
    let (x,y) = ((i%8) as usize,(i/8) as usize);
        // IDCT from CCITT Rec. T.81 (1992 E) p.27 A3.3
        let mut val = 0.0;
        for u in 0..8 {
            let mut uval = 0.0;
            uval += f[0*8 + u] as f32 * c_table[y][0];
            uval += f[1*8 + u] as f32 * c_table[y][1];
            uval += f[2*8 + u] as f32 * c_table[y][2];
            uval += f[3*8 + u] as f32 * c_table[y][3];
            uval += f[4*8 + u] as f32 * c_table[y][4];
            uval += f[5*8 + u] as f32 * c_table[y][5];
            uval += f[6*8 + u] as f32 * c_table[y][6];
            uval += f[7*8 + u] as f32 * c_table[y][7];
            val += uval * c_table[x][u];
        }
        val = val / 4.0;

        // level shift from CCITT Rec. T.81 (1992 E) p.26 A3.1
        let v = val.round() as isize + 128 ;
        if v < 0 {0} else if v > 255 {255} else {v as u8}
    }).collect();
    vals
}

// method7 https://note.com/mith_mmk/n/n6f57f007453b
// use matrix symmetry X
pub fn idct7 (f :&[i32]) -> Vec<u8> {
    let c_table :[[f32;8];8] = 
    [[ 0.70710678,  0.98078528,  0.92387953,  0.83146961,  0.70710678, 0.55557023,  0.38268343,  0.19509032],
    [ 0.70710678,  0.83146961,  0.38268343, -0.19509032, -0.70710678, -0.98078528, -0.92387953, -0.55557023],
    [ 0.70710678,  0.55557023, -0.38268343, -0.98078528, -0.70710678, 0.19509032,  0.92387953,  0.83146961],
    [ 0.70710678,  0.19509032, -0.92387953, -0.55557023,  0.70710678, 0.83146961, -0.38268343, -0.98078528],
    [ 0.70710678, -0.19509032, -0.92387953,  0.55557023,  0.70710678, -0.83146961, -0.38268343,  0.98078528],
    [ 0.70710678, -0.55557023, -0.38268343,  0.98078528, -0.70710678, -0.19509032,  0.92387953, -0.83146961],
    [ 0.70710678, -0.83146961,  0.38268343,  0.19509032, -0.70710678, 0.98078528, -0.92387953,  0.55557023],
    [ 0.70710678, -0.98078528,  0.92387953, -0.83146961,  0.70710678, -0.55557023,  0.38268343, -0.19509032]];
    let mut vals :Vec<u8> = (0..64).map(|_| 0).collect();
    for i in 0..32 {
        let (x,y) = ((i%4) as usize,(i/4) as usize);
        // IDCT from CCITT Rec. T.81 (1992 E) p.27 A3.3
        let mut val1 = 0.0;
        let mut val2 = 0.0;
        let mut plus_minus = 1.0;
        for u in 0..8 {
            let mut uval1 = 0.0;
            uval1 += f[0*8 + u] as f32 * c_table[y][0];
            uval1 += f[1*8 + u] as f32 * c_table[y][1];
            uval1 += f[2*8 + u] as f32 * c_table[y][2];
            uval1 += f[3*8 + u] as f32 * c_table[y][3];
            uval1 += f[4*8 + u] as f32 * c_table[y][4];
            uval1 += f[5*8 + u] as f32 * c_table[y][5];
            uval1 += f[6*8 + u] as f32 * c_table[y][6];
            uval1 += f[7*8 + u] as f32 * c_table[y][7];

            val1 += uval1 * c_table[x][u];
            val2 += uval1 * c_table[x][u] * plus_minus;
            plus_minus *= -1.0;
        }
        val1 = val1 / 4.0;
        val2 = val2 / 4.0;

        // level shift from CCITT Rec. T.81 (1992 E) p.26 A3.1
        let v = val1.round() as isize + 128 ;
        vals[y *8 + x] = if v < 0 {0} else if v > 255 {255} else {v as u8};
        let v = val2.round() as isize + 128 ;
        vals[y *8 + 7-x] = if v < 0 {0} else if v > 255 {255} else {v as u8};
    }
    vals
}

// method8 https://note.com/mith_mmk/n/n6f57f007453b
// use matrix symmetry X,Y
pub fn idct8 (f :&[i32]) -> Vec<u8> {
    let c_table :[[f32;8];8] = 
    [[ 0.70710678,  0.98078528,  0.92387953,  0.83146961,  0.70710678, 0.55557023,  0.38268343,  0.19509032],
    [ 0.70710678,  0.83146961,  0.38268343, -0.19509032, -0.70710678, -0.98078528, -0.92387953, -0.55557023],
    [ 0.70710678,  0.55557023, -0.38268343, -0.98078528, -0.70710678, 0.19509032,  0.92387953,  0.83146961],
    [ 0.70710678,  0.19509032, -0.92387953, -0.55557023,  0.70710678, 0.83146961, -0.38268343, -0.98078528],
    [ 0.70710678, -0.19509032, -0.92387953,  0.55557023,  0.70710678, -0.83146961, -0.38268343,  0.98078528],
    [ 0.70710678, -0.55557023, -0.38268343,  0.98078528, -0.70710678, -0.19509032,  0.92387953, -0.83146961],
    [ 0.70710678, -0.83146961,  0.38268343,  0.19509032, -0.70710678, 0.98078528, -0.92387953,  0.55557023],
    [ 0.70710678, -0.98078528,  0.92387953, -0.83146961,  0.70710678, -0.55557023,  0.38268343, -0.19509032]];
    let mut vals :Vec<u8> = (0..64).map(|_| 0).collect();
    for i in 0..16 {
        let (x,y) = ((i%4) as usize,(i/4) as usize);
        // IDCT from CCITT Rec. T.81 (1992 E) p.27 A3.3
        let mut val11 = 0.0;
        let mut val12 = 0.0;
        let mut val21 = 0.0;
        let mut val22 = 0.0;
        let mut plus_minus = 1.0;
        for u in 0..8 {
            let mut uval1 = 0.0;
            let mut uval2 = 0.0;
            uval1 += f[0*8 + u] as f32 * c_table[y][0];
            uval2 += f[0*8 + u] as f32 * c_table[y][0];

            uval1 += f[1*8 + u] as f32 * c_table[y][1];
            uval2 += f[1*8 + u] as f32 *-c_table[y][1];

            uval1 += f[2*8 + u] as f32 * c_table[y][2];
            uval2 += f[2*8 + u] as f32 * c_table[y][2];

            uval1 += f[3*8 + u] as f32 * c_table[y][3];
            uval2 += f[3*8 + u] as f32 *-c_table[y][3];

            uval1 += f[4*8 + u] as f32 * c_table[y][4];
            uval2 += f[4*8 + u] as f32 * c_table[y][4];

            uval1 += f[5*8 + u] as f32 * c_table[y][5];
            uval2 += f[5*8 + u] as f32 *-c_table[y][5];

            uval1 += f[6*8 + u] as f32 * c_table[y][6];
            uval2 += f[6*8 + u] as f32 * c_table[y][6];

            uval1 += f[7*8 + u] as f32 * c_table[y][7];
            uval2 += f[7*8 + u] as f32 *-c_table[y][7];

            val11 += uval1 * c_table[x][u];
            val12 += uval1 * c_table[x][u] * plus_minus;
            val21 += uval2 * c_table[x][u];
            val22 += uval2 * c_table[x][u] * plus_minus;
            plus_minus *= -1.0;
        }
        val11 /= 4.0;
        val12 /= 4.0;
        val21 /= 4.0;
        val22 /= 4.0;

        // level shift from CCITT Rec. T.81 (1992 E) p.26 A3.1
        let v = val11.round() as isize + 128 ;
        vals[y *8 + x] = if v < 0 {0} else if v > 255 {255} else {v as u8};
        let v = val12.round() as isize + 128 ;
        vals[y *8 + 7-x] = if v < 0 {0} else if v > 255 {255} else {v as u8};
        let v = val21.round() as isize + 128 ;
        vals[(7 - y) *8 + x] = if v < 0 {0} else if v > 255 {255} else {v as u8};
        let v = val22.round() as isize + 128 ;
        vals[(7 - y) *8 + 7-x] = if v < 0 {0} else if v > 255 {255} else {v as u8};
    }
    vals
}

// method9 https://note.com/mith_mmk/n/n6f57f007453b
// use calculate same time
pub fn idct9 (f :&[i32]) -> Vec<u8> {
    let c_table :[[f32;8];8] = 
    [[ 0.70710678,  0.98078528,  0.92387953,  0.83146961,  0.70710678, 0.55557023,  0.38268343,  0.19509032],
    [ 0.70710678,  0.83146961,  0.38268343, -0.19509032, -0.70710678, -0.98078528, -0.92387953, -0.55557023],
    [ 0.70710678,  0.55557023, -0.38268343, -0.98078528, -0.70710678, 0.19509032,  0.92387953,  0.83146961],
    [ 0.70710678,  0.19509032, -0.92387953, -0.55557023,  0.70710678, 0.83146961, -0.38268343, -0.98078528],
    [ 0.70710678, -0.19509032, -0.92387953,  0.55557023,  0.70710678, -0.83146961, -0.38268343,  0.98078528],
    [ 0.70710678, -0.55557023, -0.38268343,  0.98078528, -0.70710678, -0.19509032,  0.92387953, -0.83146961],
    [ 0.70710678, -0.83146961,  0.38268343,  0.19509032, -0.70710678, 0.98078528, -0.92387953,  0.55557023],
    [ 0.70710678, -0.98078528,  0.92387953, -0.83146961,  0.70710678, -0.55557023,  0.38268343, -0.19509032]];
    let mut vals :Vec<u8> = (0..64).map(|_| 0).collect();
    for i in 0..16 {
        let (x,y) = ((i%4) as usize,(i/4) as usize);
        // IDCT from CCITT Rec. T.81 (1992 E) p.27 A3.3
        let mut val11 = 0.0;
        let mut val12 = 0.0;
        let mut val21 = 0.0;
        let mut val22 = 0.0;
        let mut plus_minus = 1.0;
        for u in 0..8 {
            let temp1 = f[0*8 + u] as f32 * c_table[y][0] + f[2*8 + u] as f32 * c_table[y][2]
                      + f[4*8 + u] as f32 * c_table[y][4] + f[6*8 + u] as f32 * c_table[y][6];

            let temp2 = f[1*8 + u] as f32 * c_table[y][0] + f[3*8 + u] as f32 * c_table[y][2]
                      + f[5*8 + u] as f32 * c_table[y][4] + f[7*8 + u] as f32 * c_table[y][6];

            let uval1 = temp1 + temp2;
            let uval2 = temp1 - temp2;
          
            val11 += uval1 * c_table[x][u];
            val12 += uval1 * c_table[x][u] * plus_minus;
            val21 += uval2 * c_table[x][u];
            val22 += uval2 * c_table[x][u] * plus_minus;
            plus_minus *= -1.0;
        }
        val11 /= 4.0;
        val12 /= 4.0;
        val21 /= 4.0;
        val22 /= 4.0;

        // level shift from CCITT Rec. T.81 (1992 E) p.26 A3.1
        let v = val11.round() as isize + 128 ;
        vals[y *8 + x] = if v < 0 {0} else if v > 255 {255} else {v as u8};
        let v = val12.round() as isize + 128 ;
        vals[y *8 + 7-x] = if v < 0 {0} else if v > 255 {255} else {v as u8};
        let v = val21.round() as isize + 128 ;
        vals[(7 - y) *8 + x] = if v < 0 {0} else if v > 255 {255} else {v as u8};
        let v = val22.round() as isize + 128 ;
        vals[(7 - y) *8 + 7-x] = if v < 0 {0} else if v > 255 {255} else {v as u8};
    }
    vals
}

// AP-922 method10 https://note.com/mith_mmk/n/n6f57f007453b
pub fn ap922_idct(f :&[i32]) -> Vec<u8> {
    let g4 = 0.707106781186548 as f32;
    let g:[[f32;7];4]  = [
    /* row 0, 4 */
      [
        0.1733799806652680, /* g1 * 0.25 * g4 */
        0.1633203706095470, /* g2 * 0.25 * g4 */
        0.1469844503024200, /* g3 * 0.25 * g4 */
        0.1250000000000000, /* g4 * 0.25 * g4 */
        0.0982118697983878, /* g5 * 0.25 * g4 */
        0.0676495125182746, /* g6 * 0.25 * g4 */
        0.0344874224103679, /* g7 * 0.25 * g4 */
      ],
      /* row 1, 7 */
      [
        0.2404849415639110, /* g1 * 0.25 * g1 */
        0.2265318615882220, /* g2 * 0.25 * g1 */
        0.2038732892122290, /* g3 * 0.25 * g1 */
        0.1733799806652680, /* g4 * 0.25 * g1 */
        0.1362237766939550, /* g5 * 0.25 * g1 */
        0.0938325693794663, /* g6 * 0.25 * g1 */
        0.0478354290456362, /* g7 * 0.25 * g1 */
      ],
      /* row 2, 6 */
      [
        0.2265318615882220, /* g1 * 0.25 * g2 */
        0.2133883476483180, /* g2 * 0.25 * g2 */
        0.1920444391778540, /* g3 * 0.25 * g2 */
        0.1633203706095470, /* g4 * 0.25 * g2 */
        0.1283199917898340, /* g5 * 0.25 * g2 */
        0.0883883476483185, /* g6 * 0.25 * g2 */
        0.0450599888754343, /* g7 * 0.25 * g2 */
      ],
      /* row 3, 5 */
      [
        0.2038732892122290, /* g1 * 0.25 * g3 */
        0.1920444391778540, /* g2 * 0.25 * g3 */
        0.1728354290456360, /* g3 * 0.25 * g3 */
        0.1469844503024200, /* g4 * 0.25 * g3 */
        0.1154849415639110, /* g5 * 0.25 * g3 */
        0.0795474112858021, /* g6 * 0.25 * g3 */
        0.0405529186026822, /* g7 * 0.25 * g3 */
      ]];
    let t:[f32;3] = [
        0.414213562373095  /* t1 = g6/g2 */,
        0.198912367379658 /* t2 = g7/g1 */,
        0.668178637919299 /* t3 = g5/g3 */,
    ];
    
    let row2idx = [0,1,2,3,0,3,2,1];
    let mut _f = [[0_f32;8];8];
    let mut vals :Vec<u8> = (0..64).map(|_| 0).collect();

    for i in 0..8 {
        let idx = row2idx[i];
        /* P */
        let p = [
                f[0 +i] as f32,  /* 1 0 0 0 0 0 0 0 */
                f[2*8+i] as f32,  /* 0 0 1 0 0 0 0 0 */
                f[4*8+i] as f32,  /* 0 0 0 0 1 0 0 0 */
                f[6*8+i] as f32,  /* 0 0 0 0 0 0 1 0 */
                f[1*8+i] as f32,  /* 0 1 0 0 0 0 0 0 */
                f[3*8+i] as f32,  /* 0 0 0 1 0 0 0 0 */
                f[5*8+i] as f32,  /* 0 0 0 0 0 1 0 0 */
                f[7*8+i] as f32];  /* 0 0 0 0 0 0 0 1 */
        let tmp = [
            p[0] * g[idx][3],
            p[1] * g[idx][1],
            p[1] * g[idx][5],
            p[2] * g[idx][3],
            p[3] * g[idx][5],
            p[3] * g[idx][1]];
        /* M */
        let m = [
            tmp[0] + tmp[1] + tmp[3] + tmp[4],                                         /*  g4  g2  g4  g6  0  0  0  0 */
            tmp[0] + tmp[2] - tmp[3] - tmp[5],                                         /*  g4  g6 -g4 -g2  0  0  0  0 */
            tmp[0] - tmp[2] - tmp[3] + tmp[5],                                         /*  g4 -g6 -g4  g2  0  0  0  0 */
            tmp[0] - tmp[1] + tmp[3] - tmp[4],                                         /*  g4 -g2  g4 -g6  0  0  0  0 */
            p[4] * g[idx][0] + p[5] * g[idx][2] + p[6] * g[idx][4] + p[7] * g[idx][6], /*  0  0  0  0  g1  g3  g5  g7 */
            p[4] * g[idx][2] - p[5] * g[idx][6] - p[6] * g[idx][0] - p[7] * g[idx][4], /*  0  0  0  0  g3 -g7 -g1 -g5 */
            p[4] * g[idx][4] - p[5] * g[idx][0] + p[6] * g[idx][6] + p[7] * g[idx][2], /*  0  0  0  0  g5 -g1  g7  g3 */
            p[4] * g[idx][6] - p[5] * g[idx][4] + p[6] * g[idx][2] - p[7] * g[idx][0], /*  0  0  0  0  g7 -g5  g3 -g1 */
        ];
        /* A */
        _f[0][i]  = m[0] + m[4];  /*  1  0  0  0  1  0  0  0 */
        _f[1][i]  = m[1] + m[5];  /*  0  1  0  0  0  1  0  0 */
        _f[2][i]  = m[2] + m[6];  /*  0  0  1  0  0  0  1  0 */
        _f[3][i]  = m[3] + m[7];  /*  0  0  0  1  0  0  0  1 */
        _f[4][i]  = m[3] - m[7];  /*  0  0  0  1  0  0  0 -1 */
        _f[5][i]  = m[2] - m[6];  /*  0  0  1  0  0  0 -1  0 */
        _f[6][i]  = m[1] - m[5];  /*  0  1  0  0  0 -1  0  0 */
        _f[7][i]  = m[0] - m[4];  /*  1  0  0  0 -1  0  0  0 */
      }
    
      // column
      // add 8*26 = 208
      // mul 8*8  = 64
      // C = A F E B D P
      for i in 0..8 {
        /* P */
        let p = [
            _f[i][0],  /* 1 0 0 0 0 0 0 0 */
            _f[i][2],  /* 0 0 1 0 0 0 0 0 */
            _f[i][4],  /* 0 0 0 0 1 0 0 0 */
            _f[i][6],  /* 0 0 0 0 0 0 1 0 */
            _f[i][1],  /* 0 1 0 0 0 0 0 0 */
            _f[i][3],  /* 0 0 0 1 0 0 0 0 */
            _f[i][5],  /* 0 0 0 0 0 1 0 0 */
            _f[i][7],  /* 0 0 0 0 0 0 0 1 */
        ];
        /* D */
        /* g4  0  0  0  0  0  0  0 */
        /*  0  0 g4  0  0  0  0  0 */
        /*  0 g2  0  0  0  0  0  0 */
        /*  0  0  0 g2  0  0  0  0 */
        /*  0  0  0  0 g1  0  0  0 */
        /*  0  0  0  0  0  0  0 g1 */
        /*  0  0  0  0  0 g3  0  0 */
        /*  0  0  0  0  0  0 g3  0 */
        let d = [p[0],p[2],p[1],p[3], p[4],p[7],p[5],p[6]];
    
        /* B t1=g6/g2, t2=g7/g1, t3=g5/g3 */
        /*  1  1  0  0  0  0  0  0 */
        /*  1 -1  0  0  0  0  0  0 */
        /*  0  0  1 t1  0  0  0  0 */
        /*  0  0 t1 -1  0  0  0  0 */
        /*  0  0  0  0  1 t2  0  0 */
        /*  0  0  0  0 t2 -1  0  0 */
        /*  0  0  0  0  0  0  1 t3 */
        /*  0  0  0  0  0  0 t3 -1 */
        let b = [
                      d[0] +        d[1],
                      d[0] -        d[1],
                      d[2] + t[0] * d[3],
               t[0] * d[2] -        d[3],
                      d[4] + t[1] * d[5],
               t[1] * d[4] -        d[5],
                      d[6] + t[2] * d[7],
               t[2] * d[6] -        d[7],
        ];
    
        /* E */
        let e = [
            b[0] + b[2], /* 1  0  1  0  0  0  0  0 */
            b[1] + b[3], /* 0  1  0  1  0  0  0  0 */
            b[1] - b[3], /* 0  1  0 -1  0  0  0  0 */
            b[0] - b[2], /* 1  0 -1  0  0  0  0  0 */
            b[4] + b[6], /* 0  0  0  0  1  0  1  0 */
            b[4] - b[6], /* 0  0  0  0  1  0 -1  0 */
            b[5] + b[7], /* 0  0  0  0  0  1  0  1 */
            b[5] - b[7], /* 0  0  0  0  0  1  0 -1 */
        ];
        /* F g=g4*/
        let _f = [
            e[0],               /* 1  0  0  0  0  0  0  0 */
            e[1],               /* 0  1  0  0  0  0  0  0 */
            e[2],               /* 0  0  1  0  0  0  0  0 */
            e[3],               /* 0  0  0  1  0  0  0  0 */
            e[4],               /* 0  0  0  0  1  0  0  0 */
            g4 * (e[5] + e[6]), /* 0  0  0  0  0  g  g  0 */
            g4 * (e[5] - e[6]), /* 0  0  0  0  0  g -g  0 */
            e[7],               /* 0  0  0  0  0  0  0  1 */
        ];
        /* A */
        let v = (_f[0] + _f[4]).round()  as isize + 128;
        vals[i*8+0]  = v.clamp(0,255) as u8;
        let v = (_f[1] + _f[5]).round()  as isize + 128;
        vals[i*8+1]  = v.clamp(0,255) as u8;
        let v = (_f[2] + _f[6]).round()  as isize + 128;
        vals[i*8+2]  = v.clamp(0,255) as u8;
        let v = (_f[3] + _f[7]).round()  as isize + 128;
        vals[i*8+3]  = v.clamp(0,255) as u8;
        let v = (_f[3] - _f[7]).round()  as isize + 128;
        vals[i*8+4]  = v.clamp(0,255) as u8;
        let v = (_f[2] - _f[6]).round()  as isize + 128;
        vals[i*8+5]  = v.clamp(0,255) as u8;
        let v = (_f[1] - _f[5]).round()  as isize + 128;
        vals[i*8+6]  = v.clamp(0,255) as u8;
        let v = (_f[0] - _f[4]).round()  as isize + 128;
        vals[i*8+7]  = v.clamp(0,255) as u8;
    }

    vals
}