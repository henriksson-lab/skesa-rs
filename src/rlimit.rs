//! Optional address-space cap for development and testing.
//!
//! When `SKESA_RS_RLIMIT_GB` is set in the environment, applies
//! `setrlimit(RLIMIT_AS, gb * 1e9)` at process startup. Past the cap,
//! allocations fail with `ENOMEM` and the process panics — which is dramatic
//! but better than the kernel OOM-killing the whole host while iterating on
//! memory bugs.
//!
//! Production users should not set this env var. The user-facing `--memory`
//! flag drives the chunked sorted-counter (see `sorted_counter.rs`); the
//! rlimit here is only a hard backstop.
//!
//! Linux only. No-op on every other platform.

#[cfg(target_os = "linux")]
const RLIMIT_AS: i32 = 9;

#[cfg(target_os = "linux")]
#[repr(C)]
struct Rlimit {
    rlim_cur: u64,
    rlim_max: u64,
}

#[cfg(target_os = "linux")]
unsafe extern "C" {
    fn setrlimit(resource: i32, rlim: *const Rlimit) -> i32;
}

/// Honor `SKESA_RS_RLIMIT_GB` if set. Logs to stderr on success or failure.
/// Silently ignored on non-Linux platforms.
pub fn apply_from_env() {
    #[cfg(target_os = "linux")]
    {
        let Ok(gb_str) = std::env::var("SKESA_RS_RLIMIT_GB") else {
            return;
        };
        let gb: u64 = match gb_str.parse() {
            Ok(g) if g > 0 => g,
            _ => {
                eprintln!(
                    "SKESA_RS_RLIMIT_GB={gb_str:?}: not a positive integer, ignoring"
                );
                return;
            }
        };
        let bytes = gb.saturating_mul(1_000_000_000);
        let rlim = Rlimit {
            rlim_cur: bytes,
            rlim_max: bytes,
        };
        let rc = unsafe { setrlimit(RLIMIT_AS, &rlim) };
        if rc == 0 {
            // Silent on success — the env var is dev/testing only and emitting
            // to stderr breaks tests that snapshot exact CLI error output.
            // Set SKESA_RS_RLIMIT_VERBOSE=1 to confirm the cap was applied.
            if std::env::var_os("SKESA_RS_RLIMIT_VERBOSE").is_some() {
                eprintln!("SKESA_RS_RLIMIT_GB: address space capped at {gb} GB");
            }
        } else {
            eprintln!("SKESA_RS_RLIMIT_GB: setrlimit failed (rc={rc})");
        }
    }
}
