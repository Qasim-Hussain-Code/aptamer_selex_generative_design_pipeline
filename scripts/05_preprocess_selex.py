"""Stream strict FASTQ records into SQLite counts and preserve every rejection count."""
import argparse
from collections import Counter
import gzip
from pathlib import Path
import sqlite3
import os
from common import ROOT, dataset, rows, table, atomic_text, sha256
from sequences import normalize, reverse_complement


def trim(read, d):
    read = normalize(read)
    if set(read) - set("ACGT"):
        return None, "ambiguous_base"
    fwd, rev = d["forward_primer"], d["reverse_primer"]
    orientation = "forward"
    if fwd not in read or rev not in read:
        opposite = reverse_complement(read)
        if fwd in opposite and rev in opposite:
            read, orientation = opposite, "reverse"
        else:
            return None, "primer_mismatch"
    start = read.find(fwd) + len(fwd)
    end = read.find(rev, start)
    if end < start or read.count(fwd) != 1:
        return None, "primer_order_or_multiple"
    seq = read[start:end]
    if not int(d["minimum_length"]) <= len(seq) <= int(d["maximum_length"]):
        return None, "length_outside_published_range"
    if any(primer in seq for primer in [d["fold_forward"], d["fold_reverse"]]):
        return None, "constant_region_remnant"
    return seq, orientation


def fastq(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="ascii") as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq, plus, qual = f.readline().rstrip(), f.readline(), f.readline().rstrip()
            if not header.startswith("@") or not plus.startswith("+") or not seq or len(seq) != len(qual):
                raise ValueError("Malformed FASTQ record; input is not partially accepted")
            yield seq


def preprocess(path, out, d, rnd, qc_path):
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    db = out.with_name(out.name + ".sqlite.tmp")
    db.unlink(missing_ok=True)
    conn = sqlite3.connect(db)
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("CREATE TABLE counts (sequence TEXT PRIMARY KEY, n INTEGER NOT NULL) WITHOUT ROWID")
    reasons, batch = Counter(), Counter()
    raw = 0
    try:
        for read in fastq(path):
            raw += 1
            seq, reason = trim(read, d)
            reasons[reason] += 1
            if seq is not None:
                batch[seq] += 1
            if raw % 10000 == 0:
                conn.executemany("INSERT INTO counts VALUES (?,?) ON CONFLICT(sequence) DO UPDATE SET n=n+excluded.n", batch.items())
                conn.commit()
                batch.clear()
        conn.executemany("INSERT INTO counts VALUES (?,?) ON CONFLICT(sequence) DO UPDATE SET n=n+excluded.n", batch.items())
        conn.commit()
        retained = reasons["forward"] + reasons["reverse"]
        unique = conn.execute("SELECT count(*) FROM counts").fetchone()[0]
        if not retained:
            raise ValueError("No reads pass primary-source primer filters; orientation must be resolved")
        def records():
            for seq, count in conn.execute("SELECT sequence,n FROM counts ORDER BY sequence"):
                yield dict(sequence=seq, read_count=count, total_reads=retained, relative_frequency=count/retained, round=rnd, target=d["target"])
        table(out, records(), ["sequence", "read_count", "total_reads", "relative_frequency", "round", "target"])
        verify_counts(out, retained, d)
        record = dict(target=d["target"], round=rnd, raw_reads=raw, retained_reads=retained,
              unique_sequences=unique, duplicate_fraction=1-unique/retained,
              ambiguous_base_count=reasons["ambiguous_base"], discarded_reads=raw-retained,
              primer_mismatch=reasons["primer_mismatch"], primer_order_or_multiple=reasons["primer_order_or_multiple"],
              length_outside_published_range=reasons["length_outside_published_range"], constant_region_remnant=reasons["constant_region_remnant"],
              reverse_orientation_reads=reasons["reverse"], compressed_raw_bytes=Path(path).stat().st_size,
              compact_bytes=out.stat().st_size, sqlite_peak_bytes=db.stat().st_size,
              counts_sha256=sha256(out))
        if sum(reasons.values()) != raw:
            raise AssertionError("Unaccounted reads")
        table(qc_path, [record])
        return record
    except BaseException:
        out.unlink(missing_ok=True)
        raise
    finally:
        conn.close()
        db.unlink(missing_ok=True)


def verify_counts(path, expected, d):
    total = 0
    previous = None
    for r in rows(path):
        seq = r["sequence"]
        if previous is not None and seq <= previous:
            raise ValueError("Counts must be unique and sorted")
        if set(seq) - set("ACGT") or not int(d["minimum_length"]) <= len(seq) <= int(d["maximum_length"]):
            raise ValueError("Invalid processed sequence")
        n = int(r["read_count"])
        if n <= 0 or int(r["total_reads"]) != expected or abs(float(r["relative_frequency"]) - n/expected) > 1e-12:
            raise ValueError("Count/frequency mismatch")
        total += n
        previous = seq
    if total != expected:
        raise ValueError("Retained count mismatch")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True, type=Path)
    p.add_argument("--target", required=True, choices=["tg2", "integrin"])
    p.add_argument("--round", required=True, type=int)
    a = p.parse_args()
    preprocess(a.input, ROOT / f"data/processed/{a.target}/{a.round}.counts.tsv.gz", dataset(a.target), a.round,
               ROOT / f"results/qc/{a.target}_{a.round}.tsv")


if __name__ == "__main__":
    main()
