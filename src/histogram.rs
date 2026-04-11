//! Histogram analysis utilities for k-mer count distributions.
//! Port of FindValleyAndPeak, HistogramRange, GetAverageCount, CalculateGenomeSize
//! from common_util.hpp.

/// A histogram bin: (count_value, frequency)
pub type Bins = Vec<(i32, usize)>;

/// Find a valley and peak in a histogram using a simple heuristic.
/// Returns the valley index, or -1 if not found.
pub fn find_valley_and_peak(bins: &Bins, rlimit: i32) -> i32 {
    let slope_len: i32 = 5;
    let mut peak = rlimit.min(bins.len() as i32 - slope_len - 1);

    while peak >= slope_len {
        let mut is_max = true;
        for i in 1..=slope_len {
            if !is_max {
                break;
            }
            is_max = bins[(peak + i) as usize].1 < bins[peak as usize].1;
        }
        for i in 1..=slope_len {
            if !is_max {
                break;
            }
            is_max = bins[(peak - i) as usize].1 < bins[peak as usize].1;
        }
        if is_max {
            break;
        }
        peak -= 1;
    }

    if peak < slope_len {
        return -1;
    }

    let mut valley = 0i32;
    for i in 1..=peak {
        if bins[i as usize].1 < bins[valley as usize].1 {
            valley = i;
        }
    }
    if valley == peak {
        return -1;
    }

    for i in valley..bins.len() as i32 {
        if bins[i as usize].1 > bins[peak as usize].1 {
            peak = i;
        }
    }

    if (bins[valley as usize].1 as f64) < 0.7 * bins[peak as usize].1 as f64 {
        valley
    } else {
        -1
    }
}

/// Find the main range in a histogram. Returns (valley, rlimit).
/// valley == -1 if not found.
pub fn histogram_range(bins: &Bins) -> (i32, i32) {
    let min_num = 100usize;
    let mut gsize: usize = 0;
    for bin in bins {
        if bin.1 >= min_num {
            gsize += bin.0 as usize * bin.1;
        }
    }

    let mut rl = 0i32;
    let mut gs: usize = 0;
    for bin in bins {
        gs += bin.0 as usize * bin.1;
        if rl < bins.len() as i32 - 1 {
            rl += 1;
        }
        if gs as f64 > 0.8 * gsize as f64 {
            break;
        }
    }

    let mut valley = -1i32;
    let mut rlimit = rl;
    let mut genome: usize = 0;
    let mut genome_vol: usize = 0;

    loop {
        let v = find_valley_and_peak(bins, rl);

        let mut g: usize = 0;
        let mut g_vol: usize = 0;
        let start = if v >= 0 { v } else { 0 };
        for i in start..=rl {
            g_vol += bins[i as usize].0 as usize * bins[i as usize].1;
            g += bins[i as usize].1;
        }

        if (v >= 0 && g > genome) || g_vol > genome_vol {
            valley = v;
            rlimit = rl;
            genome = g;
            genome_vol = g_vol;
        }

        if v < 0 {
            break;
        }
        rl = v;
    }

    (valley, rlimit)
}

/// Get average k-mer count from histogram
pub fn get_average_count(bins: &Bins) -> f64 {
    let (mut valley, rlimit) = histogram_range(bins);
    if valley < 0 {
        valley = 0;
    }

    let mut genome: usize = 0;
    let mut kmers: usize = 0;
    for i in valley..=rlimit {
        genome += bins[i as usize].1;
        kmers += bins[i as usize].0 as usize * bins[i as usize].1;
    }

    if genome > 0 {
        kmers as f64 / genome as f64
    } else {
        0.0
    }
}

/// Calculate genome size from histogram
pub fn calculate_genome_size(bins: &Bins) -> usize {
    let (mut valley, rlimit) = histogram_range(bins);
    if valley < 0 {
        valley = 0;
    }
    let mut genome: usize = 0;
    for i in valley..=rlimit {
        genome += bins[i as usize].1;
    }
    genome
}
