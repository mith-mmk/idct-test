//use core::f32::consts::PI;

pub fn fdct(f :&[u8]) -> Vec<f32> {

    let cos = [
        1.0 , 0.9807_8528 , 0.9238_7953 , 0.8314_6961 , 0.7071_0678 , 0.5555_7023 , 0.3826_8343 , 0.1950_9032 ,
        0.0 ,-0.1950_9032 ,-0.3826_8343 ,-0.5555_7023 ,-0.7071_0678 ,-0.8314_6961 ,-0.9238_7953 ,-0.9807_8528 ,
       -1.0 ,-0.9807_8528 ,-0.9238_7953 ,-0.8314_6961 ,-0.7071_0678 ,-0.5555_7023 ,-0.3826_8343 ,-0.1950_9032 ,
       -0.0 , 0.1950_9032 , 0.3826_8343 , 0.5555_7023 , 0.7071_0678 , 0.8314_6961 , 0.9238_7953 , 0.9807_8528 ,
    ];

    let vals :Vec<f32> = (0..64).map(|i| {
        let (u,v) = ((i%8) ,(i/8) );
        // DCT from CCITT mec. T.81 (1992 E) p.27 A3.3
        let mut val: f32=0.0;
        for y in 0..8 {
            for x in 0..8 {
                val +=  (f[y * 8 + x] as f32 - 128.0)
                    * cos[((2*x+1)*u) % 32]
                    * cos[((2*y+1)*v) % 32]
            }
        }
        let cu = if u == 0 {1.0 / 2.0_f32.sqrt()} else {1.0};
        let cv = if v == 0 {1.0_f32 / 2.0_f32.sqrt()} else {1.0};
        val = cu * cv * val / 4.0;

        val
    }).collect();
    vals
}

pub fn llm_fdct(f:&[u8]) -> Vec<f32> {
    let m0 = 0.7071067811865475;
    let m1 = 1.3870398453221475;
    let m2 = 1.3065629648763766;
    let m3 = 1.1758756024193588;
    let m5 = 0.7856949583871023;
    let m6 = 0.5411961001461971;
    let m7 = 0.2758993792829431;
    let mut zz = [0_f32;64];

    for j in 0..8 {
        let i = j * 8;
        let f0 = f[i + 0] as f32 - 128.0;
        let f1 = f[i + 1] as f32 - 128.0;
        let f2 = f[i + 2] as f32 - 128.0;
        let f3 = f[i + 3] as f32 - 128.0;
        let f4 = f[i + 4] as f32 - 128.0;
        let f5 = f[i + 5] as f32 - 128.0;
        let f6 = f[i + 6] as f32 - 128.0;
        let f7 = f[i + 7] as f32 - 128.0;

        let a0 = f0 + f7;
        let a7 = f0 - f7;
        let a1 = f1 + f6;
        let a6 = f1 - f6;
        let a2 = f2 + f5;
        let a5 = f2 - f5;
        let a3 = f3 + f4;
        let a4 = f4 - f4;


        let c0 = a0 + a3;
        let c3 = a0 - a3;
        let c1 = a1 + a2;
        let c2 = a1 - a2;

        zz[i + 0] = c0 + c1;
        zz[i + 4] = c0 - c1;
        zz[i + 2] = c2 * m6 + c3 * m2;
        zz[i + 6] = c3 * m6 - c2 * m2;

        let c3 = a4 * m3 + a7 * m5;
        let c0 = a7 * m3 - a4 * m5;
        let c2 = a5 * m1 + a6 * m7;
        let c1 = a6 * m1 - a5 * m7;

        zz[i + 5] = c3 - c1;
        zz[i + 3] = c0 - c2;

        let d0 = (c0 + c2) * m0;
        let d3 = (c3 + c1) * m0;

        zz[i + 1] = d0 + d3;
        zz[i + 7] = d0 - d3;

//        for j in 0..8 {
//            zz[i + j] *= m0 * 0.5;
//        }
    }

    for i in 0..8 {
        let f0 = zz[i + 0 * 8];
        let f1 = zz[i + 1 * 8];
        let f2 = zz[i + 2 * 8];
        let f3 = zz[i + 3 * 8];
        let f4 = zz[i + 4 * 8];
        let f5 = zz[i + 5 * 8];
        let f6 = zz[i + 6 * 8];
        let f7 = zz[i + 7 * 8];

        let a0 = f0 + f7;
        let a7 = f0 - f7;
        let a1 = f1 + f6;
        let a6 = f1 - f6;
        let a2 = f2 + f5;
        let a5 = f2 - f5;
        let a3 = f3 + f4;
        let a4 = f4 - f4;

        let c0 = a0 + a3;
        let c3 = a0 - a3;
        let c1 = a1 + a2;
        let c2 = a1 - a2;

        zz[i + 0 * 8] = c0 + c1;
        zz[i + 4 * 8] = c0 - c1;
        zz[i + 2 * 8] = c2 * m6 + c3 * m2;
        zz[i + 6 * 8] = c3 * m6 - c2 * m2;

        let c3 = a4 * m3 + a7 * m5;
        let c0 = a7 * m3 - a4 * m5;
        let c2 = a5 * m1 + a6 * m7;
        let c1 = a6 * m1 - a5 * m7;

        zz[i + 5 * 8] = c3 - c1;
        zz[i + 3 * 8] = c0 - c2;

        let d0 = (c0 + c2) * m0;
        let d3 = (c3 + c1) * m0;

        zz[i + 1 * 8] = d0 + d3;
        zz[i + 7 * 8] = d0 - d3;

        // 1/(2√2) mutliply with quantize 

        for j in 0..8 {
//            zz[i + j * 8] *= m0 * 0.5 * m0 * 0.5;
            zz[i + j * 8] *= 0.125;
        }
    }

    zz.to_vec()
}

pub fn print_vec_f32(f:&[f32]) -> String {
    let mut str = "".to_string();
    for i in 0..8{
        str += &format!("{}, {}, {}, {}, {}, {}, {}, {},\n",
            f[i*8],f[i*8+1],f[i*8+2],f[i*8+3],f[i*8+4],f[i*8+5],f[i*8+6],f[i*8+7]
        );
    }
    str
}