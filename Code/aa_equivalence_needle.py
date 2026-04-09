#!/usr/bin/env python3

import argparse
import re
import sys


def parse_positions(pos_string):
    try:
        return [int(x.strip()) for x in pos_string.split(",") if x.strip()]
    except ValueError:
        raise ValueError(
            "Residue positions must be a comma-separated list of integers, for example 26,27,30,144"
        )


def normalize_direction(direction):
    d = direction.strip()
    aliases = {
        "QtoT": "QtoT",
        "TtoQ": "TtoQ",
        "1to2": "QtoT",
        "2to1": "TtoQ",
        "seq1to2": "QtoT",
        "seq2to1": "TtoQ",
        "querytotarget": "QtoT",
        "targettoquery": "TtoQ",
    }
    if d not in aliases:
        raise ValueError(
            "Direction must be one of: QtoT, TtoQ, 1to2, 2to1, seq1to2, seq2to1"
        )
    return aliases[d]


def parse_needle_pair_file(filepath):
    """
    Parse EMBOSS Needle output in pair format.

    Returns:
      seq1_name, seq2_name, seq1_aln, seq2_aln
    """

    seq1_name = None
    seq2_name = None
    matched_lines = []

    line_re = re.compile(r"^(\S+)\s+(\d+)\s+([A-Za-z\-]+)\s+(\d+)\s*$")

    with open(filepath, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")

            # Header names
            if line.startswith("# 1:"):
                seq1_name = line.split(":", 1)[1].strip()
            elif line.startswith("# 2:"):
                seq2_name = line.split(":", 1)[1].strip()

            m = line_re.match(line)
            if m:
                name, start, frag, end = m.groups()
                matched_lines.append(
                    {
                        "name": name,
                        "start": int(start),
                        "frag": frag,
                        "end": int(end),
                    }
                )

    if len(matched_lines) == 0:
        raise ValueError("No aligned sequence lines found. Is this an EMBOSS Needle pair-format file?")

    if len(matched_lines) % 2 != 0:
        raise ValueError(
            f"Unexpected odd number of sequence lines ({len(matched_lines)}). "
            "The Needle output may be truncated or not in pair format."
        )

    seq1_frags = []
    seq2_frags = []

    # In Needle pair format, matched lines alternate:
    # seq1 line, seq2 line, seq1 line, seq2 line, ...
    for i in range(0, len(matched_lines), 2):
        a = matched_lines[i]
        b = matched_lines[i + 1]

        seq1_frags.append(a["frag"])
        seq2_frags.append(b["frag"])

    seq1_aln = "".join(seq1_frags)
    seq2_aln = "".join(seq2_frags)

    if len(seq1_aln) != len(seq2_aln):
        raise ValueError("Parsed aligned strings have different lengths, which should not happen.")

    if seq1_name is None:
        seq1_name = "seq1"
    if seq2_name is None:
        seq2_name = "seq2"

    return seq1_name, seq2_name, seq1_aln, seq2_aln


def build_full_map(seq1_aln, seq2_aln):
    """
    Build residue-level equivalency map from two aligned strings.

    seq1 is treated as Q
    seq2 is treated as T
    """

    qpos = 0
    tpos = 0
    rows = []

    for a, b in zip(seq1_aln, seq2_aln):
        q_here = None
        t_here = None

        if a != "-":
            qpos += 1
            q_here = qpos

        if b != "-":
            tpos += 1
            t_here = tpos

        if a != "-" and b != "-":
            state = "match"
        elif a != "-" and b == "-":
            state = "query_only"
        elif a == "-" and b != "-":
            state = "target_only"
        else:
            # This should never occur in a valid pairwise alignment
            state = "double_gap"

        rows.append(
            {
                "query_pos": q_here,
                "target_pos": t_here,
                "query_aa": a,
                "target_aa": b,
                "state": state,
            }
        )

    return rows


def residue_equivalency(filepath, positions, direction):
    direction = normalize_direction(direction)

    seq1_name, seq2_name, seq1_aln, seq2_aln = parse_needle_pair_file(filepath)
    full_map = build_full_map(seq1_aln, seq2_aln)

    output = []
    for row in full_map:
        if direction == "QtoT":
            source_pos = row["query_pos"]
            equivalent_pos = row["target_pos"]
            source_aa = row["query_aa"] if row["query_pos"] is not None else "-"
            equivalent_aa = row["target_aa"] if row["target_pos"] is not None else "-"
        else:
            source_pos = row["target_pos"]
            equivalent_pos = row["query_pos"]
            source_aa = row["target_aa"] if row["target_pos"] is not None else "-"
            equivalent_aa = row["query_aa"] if row["query_pos"] is not None else "-"

        output.append(
            {
                "seq1_name": seq1_name,
                "seq2_name": seq2_name,
                "direction": direction,
                "source_pos": source_pos,
                "equivalent_pos": equivalent_pos,
                "source_aa": source_aa,
                "equivalent_aa": equivalent_aa,
                "query_pos": row["query_pos"],
                "target_pos": row["target_pos"],
                "query_aa": row["query_aa"],
                "target_aa": row["target_aa"],
                "state": row["state"],
            }
        )

    result = []
    for p in positions:
        matches = [r for r in output if r["source_pos"] == p]
        if matches:
            result.append(matches[0])
        else:
            result.append(
                {
                    "seq1_name": seq1_name,
                    "seq2_name": seq2_name,
                    "direction": direction,
                    "source_pos": p,
                    "equivalent_pos": None,
                    "source_aa": None,
                    "equivalent_aa": None,
                    "query_pos": None,
                    "target_pos": None,
                    "query_aa": None,
                    "target_aa": None,
                    "state": "outside_alignment_or_not_present",
                }
            )

    return result


def main():
    parser = argparse.ArgumentParser(
        description="Map residue equivalencies from EMBOSS Needle pair output."
    )
    parser.add_argument(
        "needle_file",
        help="EMBOSS Needle output file in pair format"
    )
    parser.add_argument(
        "positions",
        help="Comma-separated residue positions, for example 26,27,30,144"
    )
    parser.add_argument(
        "direction",
        help="QtoT or TtoQ. Also accepts 1to2 and 2to1"
    )

    args = parser.parse_args()

    positions = parse_positions(args.positions)
    result = residue_equivalency(args.needle_file, positions, args.direction)

    print("\t".join([
        "seq1_name",
        "seq2_name",
        "direction",
        "source_pos",
        "equivalent_pos",
        "source_aa",
        "equivalent_aa",
        "state",
        "query_pos",
        "target_pos",
        "query_aa",
        "target_aa"
    ]))

    for r in result:
        print("\t".join([
            str(r["seq1_name"]),
            str(r["seq2_name"]),
            str(r["direction"]),
            str(r["source_pos"]),
            "NA" if r["equivalent_pos"] is None else str(r["equivalent_pos"]),
            "NA" if r["source_aa"] is None else str(r["source_aa"]),
            "NA" if r["equivalent_aa"] is None else str(r["equivalent_aa"]),
            str(r["state"]),
            "NA" if r["query_pos"] is None else str(r["query_pos"]),
            "NA" if r["target_pos"] is None else str(r["target_pos"]),
            "NA" if r["query_aa"] is None else str(r["query_aa"]),
            "NA" if r["target_aa"] is None else str(r["target_aa"]),
        ]))


if __name__ == "__main__":
    main()
