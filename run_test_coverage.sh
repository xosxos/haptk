#!/bin/bash
 CARGO_INCREMENTAL=0 RUSTFLAGS='-Cinstrument-coverage' LLVM_PROFILE_FILE='cargo-test-%p-%m.profraw' cargo test --tests --features="clap" --target-dir ./target/coverage
 mkdir -p ./target/coverage
 grcov . --binary-path ./target/coverage/ \
   -s . \
   -t html \
   --threads 60 \
   --branch \
   --ignore "*target*" \
   --ignore-not-existing \
   -o ./target/coverage/

rm -r test_coverage
mkdir test_coverage
cd target/coverage
mv html/* ../../test_coverage
