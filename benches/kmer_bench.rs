use criterion::{black_box, criterion_group, criterion_main, Criterion};
use skesa_rs::kmer::Kmer;
use skesa_rs::large_int::LargeInt;
use skesa_rs::read_holder::ReadHolder;

fn bench_kmer_from_str(c: &mut Criterion) {
    c.bench_function("kmer_from_str_21", |b| {
        b.iter(|| Kmer::from_kmer_str(black_box("ACGTACGTACGTACGTACGTA")))
    });
}

fn bench_kmer_revcomp(c: &mut Criterion) {
    let kmer = Kmer::from_kmer_str("ACGTACGTACGTACGTACGTA");
    c.bench_function("kmer_revcomp_21", |b| {
        b.iter(|| black_box(&kmer).revcomp(21))
    });
}

fn bench_kmer_oahash(c: &mut Criterion) {
    let kmer = Kmer::from_kmer_str("ACGTACGTACGTACGTACGTA");
    c.bench_function("kmer_oahash_21", |b| {
        b.iter(|| black_box(&kmer).oahash())
    });
}

fn bench_large_int_revcomp(c: &mut Criterion) {
    let x = LargeInt::<1>::from_kmer("ACGTACGTACGTACGTACGTA");
    c.bench_function("largeint1_revcomp_21", |b| {
        b.iter(|| black_box(x).revcomp(21))
    });
}

fn bench_read_holder_push(c: &mut Criterion) {
    c.bench_function("readholder_push_100bp", |b| {
        b.iter(|| {
            let mut rh = ReadHolder::new(false);
            rh.push_back_str(black_box(
                "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
            ));
            rh
        })
    });
}

fn bench_kmer_extraction(c: &mut Criterion) {
    let mut rh = ReadHolder::new(false);
    rh.push_back_str("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");

    c.bench_function("kmer_extraction_21_from_100bp", |b| {
        b.iter(|| {
            let mut ki = rh.kmer_iter(21);
            let mut count = 0;
            while !ki.at_end() {
                black_box(ki.get());
                ki.advance();
                count += 1;
            }
            count
        })
    });
}

fn bench_sorted_counter(c: &mut Criterion) {
    let data_dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data");
    let fasta = data_dir.join("small_test.fasta");
    let rg = skesa_rs::reads_getter::ReadsGetter::new(
        &[fasta.to_str().unwrap().to_string()],
        false,
    )
    .unwrap();
    let reads = rg.reads().to_vec();

    c.bench_function("sorted_counter_200reads_k21", |b| {
        b.iter(|| {
            skesa_rs::sorted_counter::count_kmers_sorted(
                black_box(&reads),
                21,
                2,
                true,
                32,
            )
        })
    });
}

criterion_group!(
    benches,
    bench_kmer_from_str,
    bench_kmer_revcomp,
    bench_kmer_oahash,
    bench_large_int_revcomp,
    bench_read_holder_push,
    bench_kmer_extraction,
    bench_sorted_counter,
);
criterion_main!(benches);
