fn main() {
    // Only compile C++ FFI when the "ffi" feature is enabled
    if std::env::var("CARGO_FEATURE_FFI").is_err() {
        return;
    }

    #[cfg(feature = "ffi")]
    compile_ffi();
}

#[cfg(feature = "ffi")]
fn compile_ffi() {
    let boost_include = std::env::var("BOOST_INCLUDE")
        .unwrap_or_else(|_| "/usr/include".to_string());
    let boost_lib = std::env::var("BOOST_LIB")
        .unwrap_or_else(|_| "/usr/lib/x86_64-linux-gnu".to_string());

    cc::Build::new()
        .cpp(true)
        .flag("-std=c++11")
        .file("ffi/skesa_wrapper.cpp")
        .file("ffi/hash_test_wrapper.cpp")
        .file("SKESA/glb_align.cpp")
        .include("SKESA")
        .include(&boost_include)
        .define("NO_NGS", None)
        .flag("-O3")
        .flag("-msse4.2")
        .flag("-pthread")
        .flag("-fPIC")
        .flag("-Wno-format-y2k")
        .flag("-finline-functions")
        .flag("-fstrict-aliasing")
        .flag("-fomit-frame-pointer")
        .compile("skesa_ffi");

    println!("cargo:rustc-link-search=native={}", boost_lib);
    println!("cargo:rustc-link-lib=static=boost_program_options");
    println!("cargo:rustc-link-lib=static=boost_iostreams");
    println!("cargo:rustc-link-lib=static=boost_regex");
    println!("cargo:rustc-link-lib=static=boost_timer");
    println!("cargo:rustc-link-lib=static=boost_chrono");
    println!("cargo:rustc-link-lib=static=boost_system");
    println!("cargo:rustc-link-lib=dylib=z");
    println!("cargo:rustc-link-lib=dylib=pthread");
    println!("cargo:rustc-link-lib=dylib=rt");
    println!("cargo:rustc-link-lib=dylib=dl");
    println!("cargo:rustc-link-lib=dylib=stdc++");

    println!("cargo:rerun-if-changed=ffi/skesa_wrapper.cpp");
    println!("cargo:rerun-if-changed=ffi/skesa_wrapper.h");
    println!("cargo:rerun-if-changed=ffi/hash_test_wrapper.cpp");
    println!("cargo:rerun-if-changed=ffi/hash_test_wrapper.h");
    println!("cargo:rerun-if-changed=SKESA/glb_align.cpp");
}
