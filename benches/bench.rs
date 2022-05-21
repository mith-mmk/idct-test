
use criterion::{black_box, criterion_group, criterion_main, Criterion};

use idct_test::idct;
use idct_test::fdct;

const zz :[i32;64] = [
      568,     0,     0,    -4,    -4,     0,     4,     0 ,
      -27,     9,    -4,    -4,     0,    -5,     5,    -5 ,
      -49,    -4,     4,     4,     0,     0,     0,     0 ,
      -12,    -4,     0,     0,     5,     0,     0,     0 ,
      -14,    -5,     0,     0,     0,     0,     0,     0 ,
       -5,     0,     0,     0,     0,     0,     0,     0 ,
       -5,     0,     0,     0,     0,     0,     0,     0 ,
        0,     0,     0,     0,     0,     0,     0,     1 ];

const z :[u8;64] = [
            128,   128,   128,   128,   128,   128,   128,   128 ,
            128,   255,   255,   255,   255,   255,   255,   255 ,
            128,   255,   255,   255,   255,   255,   255,   255 ,
            128,   255,   255,   255,   255,   255,   255,   255 ,
            128,   255,   128,   128,   128,   128,   128,   128 ,
            128,   255,   128,     0,     0,     0,     0,     0 ,
            128,   255,   128,     0,     0,     0,     0,     0 ,
            128,   255,   128,     0,     0,     0,     0,     0 
        ];

fn idct(c: &mut Criterion) {
    c.bench_function(
        "standard IDCT",
        |b| b.iter(|| idct::idct(black_box(&zz)))
    );
}

fn ap922_idct(c: &mut Criterion) {
    c.bench_function(
        "AP922 IDCT",
        |b| b.iter(|| idct::ap922_idct(black_box(&zz)))
    );
}

fn aan_idct(c: &mut Criterion) {
    c.bench_function(
        "AAN IDCT",
        |b| b.iter(|| idct::fast_idct(black_box(&zz)))
    );
}

fn aan_idct_f64(c: &mut Criterion) {
    c.bench_function(
        "AAN f64 IDCT",
        |b| b.iter(|| idct::fast_idct_f64(black_box(&zz)))
    );
}

fn llm_idct(c: &mut Criterion) {
    c.bench_function(
        "LLM IDCT",
        |b| b.iter(|| idct::llm_idct(black_box(&zz)))
    );
}

fn llm_fdct(c: &mut Criterion) {
    c.bench_function(
        "LLM FDCT",
        |b| b.iter(|| fdct::llm_fdct(black_box(&z)))
    );
}

fn std_fdct(c: &mut Criterion) {
    c.bench_function(
        "Standard FDCT",
        |b| b.iter(|| fdct::fdct(black_box(&z)))
    );
}

criterion_group!(benches,idct,ap922_idct, aan_idct, aan_idct_f64, llm_idct,std_fdct,llm_fdct);
criterion_main!(benches);