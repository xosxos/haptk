//
// https://github.com/fwcd/hpc-smith-waterman
//
pub const G_INIT: i16 = 2;
pub const G_EXT: i16 = 2;
pub const WEIGHT_IF_EQ: i16 = 3;

// fn _run() {
//     let naive_engine = NaiveSmithWaterman;

//     let db = b"ACATAAGATAGATA";
//     let query = b"TAAG";

//     let (start, stop) = naive_engine.align(db, query);

//     println!("primer part: {:?}", &db[start..=stop]);
// }

/// A facility that computes the alignment of two sequences.
pub trait Engine {
    /// The engine's name.
    fn name(&self) -> String;

    /// Aligns the given two sequences.
    fn align<'a>(&self, database: &'a [u8], query: &'a [u8]) -> (usize, usize);
}

/// An engine that computes alignments using the Smith-Waterman-Algorithm (naively) on the CPU.
pub struct NaiveSmithWaterman;

impl NaiveSmithWaterman {
    fn weight(d: u8, q: u8) -> i16 {
        if d == q {
            WEIGHT_IF_EQ
        } else {
            -WEIGHT_IF_EQ
        }
    }
}

impl Engine for NaiveSmithWaterman {
    fn name(&self) -> String {
        "Naive (CPU)".to_owned()
    }

    fn align<'a>(&self, database: &'a [u8], query: &'a [u8]) -> (usize, usize) {
        let db_len = database.len();
        let query_len = query.len();
        let height = db_len + 1;
        let width = query_len + 1;
        let size = height * width;

        // Create scoring matrix h, helper matrix f and a
        // helper matrix p that tracks the previous index

        let mut h = vec![0; size];
        let mut f = vec![0; size];
        let mut p = vec![0; size];

        // Perform scoring stage (dynamic programming-style)

        for i in 1..=db_len {
            // We don't need to store e as a matrix since we iterate
            // from left to right (thus we only need the last value)
            let mut e_here: i16 = 0;

            for j in 1..=query_len {
                // Compute indices for the neighboring cells
                let here = i * width + j;
                let above = (i - 1) * width + j;
                let left = i * width + j - 1;
                let above_left = (i - 1) * width + j - 1;

                // Compute helper values
                e_here = (e_here - G_EXT).max(h[left] - G_INIT);
                f[here] = (f[above] - G_EXT).max(h[above] - G_INIT);

                // Compute value and the remember the index the maximum came from
                // (we need this later for the traceback phase)
                let (max_origin, max_value) = [
                    (0, 0),
                    (
                        above_left,
                        h[above_left] + Self::weight(database[i - 1], query[j - 1]),
                    ),
                    (left, e_here),
                    (above, f[here]),
                ]
                .into_iter()
                .max_by_key(|&(_, x)| x)
                .unwrap();

                h[here] = max_value;
                p[here] = max_origin;
            }
        }

        // Perform traceback stage (using the previously computed scoring matrix h)

        let mut i = (0..size).max_by_key(|&i| h[i]).unwrap();
        let mut database_indices = Vec::new();

        while i > 0 && h[i] > 0 {
            database_indices.push((i / width) - 1);
            i = p[i];
        }

        let first = database_indices[database_indices.len() - 1];
        let last = database_indices[0];

        (first, last)
    }
}
