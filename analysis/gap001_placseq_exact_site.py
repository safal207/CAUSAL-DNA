#!/usr/bin/env python3
"""Query processed adult hypothalamus PLAC-seq contacts at exact rs1421085 site.

The analysis is intentionally bounded: it asks whether the 10-kb contact bin
containing the exact mouse rs1421085 ortholog participates in a published
significant H3K27ac+ PLAC-seq interaction whose opposite anchor overlaps the
Irx3 promoter-side bin. It does not infer allele dependence or causality.
"""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

TARGET_CHROM = "chr8"
TARGET_SITE0 = 91374371
TARGET_TILE = (91374000, 91374500)
TARGET_BIN = (91370000, 91380000)
# NCBI Gene 16373, GRCm38.p6: 91798511..91801654, complement (1-based).
IRX3_GENE = (91798510, 91801654)  # 0-based half-open
IRX3_TSS0 = 91801653
IRX3_PROMOTER_BIN = (91800000, 91810000)


def overlap(a0: int, a1: int, b0: int, b1: int) -> bool:
    return max(a0, b0) < min(a1, b1)


def distance(a0: int, a1: int, b0: int, b1: int) -> int:
    if overlap(a0, a1, b0, b1):
        return 0
    if a1 <= b0:
        return b0 - a1
    return a0 - b1


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def parse_bedpe(path: Path):
    rows=[]
    header=None
    preview=[]
    with open_text(path) as fh:
        for line in fh:
            line=line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields=line.split("\t")
            if len(preview)<5:
                preview.append(fields)
            if len(fields)<6:
                continue
            try:
                s1,e1,s2,e2=map(int,[fields[1],fields[2],fields[4],fields[5]])
            except ValueError:
                if header is None:
                    header=fields
                continue
            rows.append({
                "chrom1":fields[0],"start1":s1,"end1":e1,
                "chrom2":fields[3],"start2":s2,"end2":e2,
                "extra":fields[6:],
            })
    if not rows:
        raise RuntimeError(f"no BEDPE-like rows parsed from {path}")
    return header,preview,rows


def anchor_overlap(row, which: int, chrom: str, interval: tuple[int,int]) -> bool:
    return row[f"chrom{which}"]==chrom and overlap(row[f"start{which}"],row[f"end{which}"],*interval)


def serial(row):
    return {
        "anchor1":[row['chrom1'],row['start1'],row['end1']],
        "anchor2":[row['chrom2'],row['start2'],row['end2']],
        "extra":row['extra'],
    }


def query_contacts(rows):
    target_rows=[]
    exact_irx3=[]
    promoter_rows=[]
    for r in rows:
        t1=anchor_overlap(r,1,TARGET_CHROM,TARGET_BIN)
        t2=anchor_overlap(r,2,TARGET_CHROM,TARGET_BIN)
        p1=anchor_overlap(r,1,TARGET_CHROM,IRX3_PROMOTER_BIN)
        p2=anchor_overlap(r,2,TARGET_CHROM,IRX3_PROMOTER_BIN)
        if t1 or t2:
            target_rows.append(serial(r))
        if p1 or p2:
            promoter_rows.append(serial(r))
        if (t1 and p2) or (t2 and p1):
            exact_irx3.append(serial(r))
    return target_rows, exact_irx3, promoter_rows


def parse_peak_probe(path: Path):
    preview=[]
    overlaps=[]
    numeric_rows=0
    with open_text(path) as fh:
        for line in fh:
            line=line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            f=line.split("\t")
            if len(preview)<8:
                preview.append(f)
            if len(f)<3:
                continue
            try:
                start,end=int(f[1]),int(f[2])
            except ValueError:
                continue
            numeric_rows+=1
            if f[0]==TARGET_CHROM and overlap(start,end,*TARGET_BIN):
                overlaps.append(f)
    return {"preview":preview,"numeric_rows":numeric_rows,"rows_overlapping_target_10kb_bin":overlaps}


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--contacts",type=Path,required=True)
    ap.add_argument("--peaks",type=Path,required=True)
    ap.add_argument("--fasted-loops",type=Path,required=False)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()

    header,preview,rows=parse_bedpe(args.contacts)
    target_rows, exact_irx3, promoter_rows=query_contacts(rows)
    peak_probe=parse_peak_probe(args.peaks)

    fasted=None
    if args.fasted_loops and args.fasted_loops.exists():
        # This file is an EdgeR differential-loop table, not the primary
        # significant-contact authority. Preserve only format/target-overlap probe.
        fp=parse_peak_probe(args.fasted_loops)
        fasted=fp

    result={
        "analysis":"GAP-001 exact rs1421085-bin adult hypothalamus PLAC-seq query",
        "claim_scope":"wild-type processed contact-map evidence only; no allele-specific or causal promotion",
        "authority":{
            "series":"GSE295850",
            "file":args.contacts.name,
            "assay":"H3K27ac+ PLAC-seq",
            "published_threshold":"FDR <1%, total read count >=12, observed/expected >2; n=10",
            "assembly":"GRCm38/mm10",
        },
        "coordinates":{
            "exact_site_0_based":[TARGET_CHROM,TARGET_SITE0,TARGET_SITE0+1],
            "target_tile_0_based":[TARGET_CHROM,*TARGET_TILE],
            "target_contact_bin_10kb":[TARGET_CHROM,*TARGET_BIN],
            "irx3_gene_0_based":[TARGET_CHROM,*IRX3_GENE],
            "irx3_tss_0_based":IRX3_TSS0,
            "irx3_promoter_contact_bin_10kb":[TARGET_CHROM,*IRX3_PROMOTER_BIN],
            "target_to_irx3_tss_bp":IRX3_TSS0-TARGET_SITE0,
        },
        "contact_file":{
            "header":header,
            "preview":preview,
            "parsed_contacts":len(rows),
            "contacts_with_target_bin_anchor":len(target_rows),
            "contacts_with_irx3_promoter_bin_anchor":len(promoter_rows),
            "target_bin_to_irx3_promoter_bin_contacts":len(exact_irx3),
            "target_bin_contacts":target_rows,
            "target_to_irx3_contacts":exact_irx3,
        },
        "h3k27ac_anchor_probe":peak_probe,
        "fasted_differential_loop_probe":fasted,
        "interpretation":{
            "positive_rule":"At least one published significant GSE295850 contact must contain the exact-site 10-kb bin on one anchor and the Irx3 promoter 10-kb bin on the other.",
            "negative_rule":"Zero matching significant contacts is evidence of no called contact at this resolution/threshold, not proof that physical contact never occurs.",
            "causal_boundary":"Even a positive wild-type contact cannot show T>C allele dependence or mediation."
        }
    }
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print('PARSED_CONTACTS',len(rows))
    print('TARGET_BIN_CONTACTS',len(target_rows))
    print('IRX3_PROMOTER_CONTACTS',len(promoter_rows))
    print('TARGET_TO_IRX3',len(exact_irx3))
    print('H3K27AC_TARGET_BIN_ROWS',len(peak_probe['rows_overlapping_target_10kb_bin']))
    if exact_irx3:
        print('EXACT_SITE_3D_CONTACT_SUPPORT: POSITIVE_AT_10KB_BIN_RESOLUTION')
        for r in exact_irx3: print(json.dumps(r))
    else:
        print('EXACT_SITE_3D_CONTACT_SUPPORT: NOT_CALLED_AT_PUBLISHED_THRESHOLD')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
