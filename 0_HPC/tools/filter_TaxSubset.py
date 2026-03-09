#!/usr/bin/env python3
import argparse
import pandas as pd
import re
import os
import csv

def detect_hit_col(cols):
    if "TargetName" in cols:
        return "TargetName"
    if "targetName" in cols:
        return "targetName"
    raise SystemExit("ERROR: HMMER file must contain 'TargetName' or 'targetName'.")


def main():
    ap = argparse.ArgumentParser(description="Subset a taxonomy table down to HMMER-positive contigs.")
    ap.add_argument("--mode", required=True, choices=["PA", "NS"], help="Dataset type")
    ap.add_argument("--hmmer", required=True, help="Path to combined HMMER hits CSV")
    ap.add_argument("--taxonomy", required=True, help="Path to taxonomy CSV")
    ap.add_argument("--output", required=True, help="Path for output filtered taxonomy CSV")
    ap.add_argument("--chunksize", type=int, default=100_000)
    ap.add_argument("--hmmer_contig_col", default=None)
    ap.add_argument("--tax_contig_col", default="nt_id")
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)

 # Load HMMER hits and build contig set
    hits_df = pd.read_csv(args.hmmer, dtype=str)
    hit_col = args.hmmer_contig_col or detect_hit_col(hits_df.columns)

    hit_series = hits_df[hit_col].dropna().astype(str)

    hit_set = set(hit_series.unique())
    if not hit_set:
        raise SystemExit("ERROR: No contigs found in HMMER hits after processing.")

 # Stream taxonomy and filter
    first = True
    kept_rows = 0

    reader = pd.read_csv(
        args.taxonomy,
        chunksize=args.chunksize,
        dtype=str,
        sep=",",
        on_bad_lines="skip"
     )

    for chunk in reader:
        if args.tax_contig_col not in chunk.columns:
            raise SystemExit(
                f"ERROR: taxonomy missing '{args.tax_contig_col}'. Columns: {list(chunk.columns)}"
            )

        chunk[args.tax_contig_col] = chunk[args.tax_contig_col].astype(str)
        sub = chunk[chunk[args.tax_contig_col].isin(hit_set)]
        if sub.empty:
            continue

        kept_rows += len(sub)
        sub.to_csv(args.output, index=False, mode="w" if first else "a", header=first)
        first = False

    if first:
        cols = pd.read_csv(args.taxonomy, nrows=0, engine="python").columns
        pd.DataFrame(columns=cols).to_csv(args.output, index=False)

    print(f"✅ Finished. HMMER contigs: {len(hit_set)} | Taxonomy rows kept: {kept_rows} | Output: {args.output}")


if __name__ == "__main__":
    main()
