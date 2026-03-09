import pandas as pd
import argparse
import re

def main():
    p = argparse.ArgumentParser(description="Filter counts by HMMER hits (PA or NS)")
    p.add_argument("--mode", required=True, choices=["PA", "NS"], help="Dataset type")
    p.add_argument("--hmmer", required=True, help="Path to HMMER hits CSV")
    p.add_argument("--counts", required=True, help="Path to counts CSV")
    p.add_argument("--output", required=True, help="Path for output CSV")
    p.add_argument("--chunksize", type=int, default=100_000)
    args = p.parse_args()

    # Load HMMER hits (TargetName or targetName)
    hits_df = pd.read_csv(args.hmmer)
    name_col = "TargetName" if "TargetName" in hits_df.columns else (
               "targetName" if "targetName" in hits_df.columns else None)
    if name_col is None:
        raise SystemExit("ERROR: HMMER file must contain 'TargetName' or 'targetName'.")

    if args.mode == "PA":
        # clean suffix _<digits>, match counts by target_id
        hit_list = hits_df[name_col].astype(str).apply(
            lambda x: re.sub(r"_[0-9]+$", "", x)
        ).unique()
        id_col = "target_id"
    else:  # NS
        # no cleaning, match counts by nt_id
        hit_list = hits_df[name_col].astype(str).unique()
        id_col = "nt_id"

    hit_set = set(hit_list)

    # Stream-filter counts
    first = True
    for chunk in pd.read_csv(args.counts, chunksize=args.chunksize):
        # Check for column-of-interest
        if id_col not in chunk.columns:
            raise SystemExit(f"ERROR: counts missing '{id_col}'")

        # Subset rows by hit IDs
        chunk[id_col] = chunk[id_col].astype(str)
        sub = chunk[chunk[id_col].isin(hit_set)]
        if sub.empty:
            continue

        # --- NEW: trim columns for NS mode ---
        if args.mode == "NS":
            # core columns you want to keep
            ns_keep = [
                "nt_id_sample",
                "nt_id",
                "sample_name",
                "raw_counts_A",
                "raw_counts_B",
                "raw_counts_C",
                "transcripts_L_A",
                "transcripts_L_B",
                "transcripts_L_C",
            ]
            # include any unnamed index column if present
            unnamed_cols = [c for c in sub.columns if c.startswith("Unnamed")]
            keep_cols = unnamed_cols + [c for c in ns_keep if c in sub.columns]
            sub = sub[keep_cols]

        sub.to_csv(
            args.output,
            index=False,
            mode="w" if first else "a",
            header=first
        )
        first = False

    print(f"✅ Finished writing filtered counts to {args.output}")

if __name__ == "__main__":
    main()

