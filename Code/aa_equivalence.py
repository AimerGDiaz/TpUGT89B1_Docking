#!/usr/bin/env python3

import argparse
import re
import sys
from pathlib import Path


def parse_cigar(cigar: str):
    tokens = re.findall(r"(\d+)([MID])", cigar)
    if not tokens:
        raise ValueError(f"No valid CIGAR operations found in: {cigar}")
    return [(int(length), op) for length, op in tokens]


def parse_positions(pos_string: str):
    try:
        return [int(x.strip()) for x in pos_string.split(",") if x.strip()]
    except ValueError:
        raise ValueError("Residue positions must be a comma-separated list of integers, for example 26,27,30,144")


def parse_alignment_line(line: str, sep: str = "\t"):
    fields = line.rstrip("\n").split(sep)
    if len(fields) < 13:
        raise ValueError(
            f"Expected at least 13 tab-separated columns, got {len(fields)}.\n"
            f"Line was:\n{line}"
        )

    return {
        "query": fields[0],
        "target": fields[1],
        "fident": fields[2],
        "alnlen": fields[3],
        "qstart": int(fields[6]),
        "qend": int(fields[7]),
        "tstart": int(fields[8]),
        "tend": int(fields[9]),
        "cigar": fields[12],
    }


def build_full_map(qstart: int, tstart: int, cigar: str):
    """
    Assumed convention:
      M = aligned residue in both query and target
      I = residue only in query (gap in target)
      D = residue only in target (gap in query)
    """
    ops = parse_cigar(cigar)
    qpos = qstart
    tpos = tstart
    rows = []

    for length, op in ops:
        if op == "M":
            for _ in range(length):
                rows.append({
                    "query_pos": qpos,
                    "target_pos": tpos,
                    "state": "match"
                })
                qpos += 1
                tpos += 1

        elif op == "I":
            for _ in range(length):
                rows.append({
                    "query_pos": qpos,
                    "target_pos": None,
                    "state": "query_only"
                })
                qpos += 1

        elif op == "D":
            for _ in range(length):
                rows.append({
                    "query_pos": None,
                    "target_pos": tpos,
                    "state": "target_only"
                })
                tpos += 1

        else:
            raise ValueError(f"Unsupported CIGAR operation: {op}")

    return rows


def residue_equivalency(line: str, direction: str, positions):
    aln = parse_alignment_line(line)
    full_map = build_full_map(aln["qstart"], aln["tstart"], aln["cigar"])

    if direction not in {"QtoT", "TtoQ"}:
        raise ValueError("Direction must be either QtoT or TtoQ")

    output = []
    for row in full_map:
        if direction == "QtoT":
            source_pos = row["query_pos"]
            equivalent_pos = row["target_pos"]
        else:
            source_pos = row["target_pos"]
            equivalent_pos = row["query_pos"]

        output.append({
            "query": aln["query"],
            "target": aln["target"],
            "direction": direction,
            "source_pos": source_pos,
            "equivalent_pos": equivalent_pos,
            "state": row["state"],
            "query_pos": row["query_pos"],
            "target_pos": row["target_pos"],
        })

    result = []
    for p in positions:
        matches = [r for r in output if r["source_pos"] == p]
        if matches:
            result.append(matches[0])
        else:
            result.append({
                "query": aln["query"],
                "target": aln["target"],
                "direction": direction,
                "source_pos": p,
                "equivalent_pos": None,
                "state": "outside_alignment_or_not_present",
                "query_pos": None,
                "target_pos": None,
            })

    return result


def read_first_valid_line(filepath: str):
    valid_lines = []
    with open(filepath, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            valid_lines.append(line)

    if not valid_lines:
        raise ValueError(f"No valid alignment lines found in {filepath}")

    if len(valid_lines) > 1:
        print(
            f"Warning: {filepath} contains {len(valid_lines)} alignment lines. "
            f"Using the first non-empty non-comment line only.",
            file=sys.stderr
        )

    return valid_lines[0]


def main():
    parser = argparse.ArgumentParser(
        description="Map residue equivalencies from a Foldseek/MMseqs tabular alignment line."
    )
    parser.add_argument("m8_file", help="Path to tabular alignment file, for example pdb.m8")
    parser.add_argument("positions", help="Comma-separated residue positions, for example 26,27,30,144")
    parser.add_argument("direction", choices=["QtoT", "TtoQ"], help="Mapping direction")
    args = parser.parse_args()

    positions = parse_positions(args.positions)
    line = read_first_valid_line(args.m8_file)
    result = residue_equivalency(line, args.direction, positions)

    print("\t".join([
        "query",
        "target",
        "direction",
        "source_pos",
        "equivalent_pos",
        "state",
        "query_pos",
        "target_pos"
    ]))

    for r in result:
        print("\t".join([
            str(r["query"]),
            str(r["target"]),
            str(r["direction"]),
            str(r["source_pos"]),
            "NA" if r["equivalent_pos"] is None else str(r["equivalent_pos"]),
            str(r["state"]),
            "NA" if r["query_pos"] is None else str(r["query_pos"]),
            "NA" if r["target_pos"] is None else str(r["target_pos"]),
        ]))


if __name__ == "__main__":
    main()
