# Single Jupyter cell: CSV -> .dna (FASTA) converter (DNA-only)
# Paste this entire cell into one Jupyter notebook code cell and run.
# Configure the variables below (INPUT_CSV, OUTPUT_DIR) and run the cell.
# Run as ipynb
# Requirements: pandas (`pip install pandas`)

import os
import re
from typing import Tuple
import pandas as pd

# ---------------------------
# User-configurable variables
# ---------------------------
INPUT_CSV = "example.csv"   # set your CSV path here
OUTPUT_DIR = "pDNA_outputs"            # set desired output directory here
NAME_COL = "name"                     # CSV column for names
BACKBONE_COL = "backbone"             # CSV column for backbone sequences
INSERT_COL = "insert"                 # CSV column for insert sequences
FILE_EXT = ".dna"                     # output extension (keeps .dna by default)
FASTA_MODE = True                     # True -> write FASTA; False -> write raw sequence
LINE_WIDTH = 80                       # FASTA line width

# ---------------------------
# Helper functions
# ---------------------------

DNA_CHARS = set("ACGTacgt")

def sanitize_filename(name: str) -> str:
    """Make a filesystem-safe filename from name."""
    name = str(name).strip()
    name = re.sub(r'[\/\\\:\*\?\"<>\|]', "_", name)
    name = re.sub(r'\s+', '_', name)
    return name or "untitled"

def is_valid_dna_sequence(seq: str) -> bool:
    """Return True only if seq is non-empty and contains only A/T/G/C (case-insensitive)."""
    if not isinstance(seq, str):
        return False
    s = seq.strip()
    if s == "":
        return False
    return all(ch in DNA_CHARS for ch in s)

def write_fasta(path: str, header: str, sequence: str, line_width: int = 80) -> None:
    """Write sequence to path in FASTA format with given header."""
    seq = sequence.upper().replace(" ", "").replace("\n", "")
    with open(path, "w", encoding="utf-8") as f:
        f.write(f">{header}\n")
        for i in range(0, len(seq), line_width):
            f.write(seq[i:i+line_width] + "\n")

def unique_output_path(directory: str, filename: str) -> str:
    """Return a non-conflicting path by appending numeric suffix if needed."""
    base, ext = os.path.splitext(filename)
    candidate = os.path.join(directory, filename)
    if not os.path.exists(candidate):
        return candidate
    i = 1
    while True:
        new_name = f"{base}_{i}{ext}"
        candidate = os.path.join(directory, new_name)
        if not os.path.exists(candidate):
            return candidate
        i += 1

# ---------------------------
# Core processing function
# ---------------------------

def combine_and_save_dna_sequences(
    input_csv_path: str,
    output_dir: str,
    name_col: str = "name",
    backbone_col: str = "backbone",
    insert_col: str = "insert",
    file_ext: str = ".dna",
    fasta_mode: bool = True,
    line_width: int = 80
) -> Tuple[int, int]:
    """
    Read CSV and produce one FASTA (.dna) file per valid row.
    Returns (created_count, skipped_count).
    """
    if not os.path.isfile(input_csv_path):
        raise FileNotFoundError(f"Input CSV not found: {input_csv_path}")

    os.makedirs(output_dir, exist_ok=True)

    df = pd.read_csv(input_csv_path, dtype=str).fillna("")

    required = {name_col, backbone_col, insert_col}
    if not required.issubset(set(df.columns)):
        missing = required - set(df.columns)
        raise ValueError(f"CSV missing required columns: {', '.join(sorted(missing))}")

    created = 0
    skipped = 0

    for idx, row in df.iterrows():
        raw_name = (row.get(name_col) or "").strip()
        backbone = (row.get(backbone_col) or "").strip()
        insert = (row.get(insert_col) or "").strip()

        if raw_name == "":
            print(f"[WARN] row {idx}: empty name -> skipped")
            skipped += 1
            continue

        if not is_valid_dna_sequence(backbone):
            print(f"[WARN] row {idx} '{raw_name}': invalid backbone (only A/T/G/C allowed) -> skipped")
            skipped += 1
            continue

        if not is_valid_dna_sequence(insert):
            print(f"[WARN] row {idx} '{raw_name}': invalid insert (only A/T/G/C allowed) -> skipped")
            skipped += 1
            continue

        combined = (backbone + insert).upper()
        safe_name = sanitize_filename(raw_name)
        out_fname = safe_name + (file_ext if file_ext.startswith(".") else f".{file_ext}")
        out_path = unique_output_path(output_dir, out_fname)

        try:
            if fasta_mode:
                header = f"{raw_name} | row={idx}"
                write_fasta(out_path, header, combined, line_width=line_width)
            else:
                with open(out_path, "w", encoding="utf-8") as f:
                    f.write(combined)
            print(f"[OK] row {idx} '{raw_name}' -> {out_path}")
            created += 1
        except OSError as e:
            print(f"[ERROR] row {idx} '{raw_name}': write failed: {e}")
            skipped += 1

    print(f"Done: created {created}, skipped {skipped}")
    return created, skipped

# ---------------------------
# Optional: create example CSV for quick testing (comment/uncomment as needed)
# ---------------------------

if not os.path.exists(INPUT_CSV):
    df_example = pd.DataFrame([
        {"name": "construct1", "backbone": "ATGCGTAC", "insert": "GGTTA"},
        {"name": "construct 2", "backbone": "atgc", "insert": "ttaa"},
        {"name": "bad_seq", "backbone": "ATGXB", "insert": "ATGC"},
        {"name": "", "backbone": "ATGC", "insert": "ATGC"},
    ])
    df_example.to_csv(INPUT_CSV, index=False)
    print(f"Example CSV created at: {INPUT_CSV}")

# ---------------------------
# Run conversion
# ---------------------------

created_count, skipped_count = combine_and_save_dna_sequences(
    input_csv_path=INPUT_CSV,
    output_dir=OUTPUT_DIR,
    name_col=NAME_COL,
    backbone_col=BACKBONE_COL,
    insert_col=INSERT_COL,
    file_ext=FILE_EXT,
    fasta_mode=FASTA_MODE,
    line_width=LINE_WIDTH
)

print(f"Summary: {created_count} created, {skipped_count} skipped. Output dir: {os.path.abspath(OUTPUT_DIR)}")
