#!/usr/bin/env python3
"""Select the official Nature Codebook motif archive by the full-atlas contract.

Nature's article page may be reduced for automated clients and Springer exposes
multiple ZIP media objects for the exact DOI.  Filename/member heuristics are
not sufficient provenance.  This selector therefore evaluates every exact-DOI
motif-like ZIP with the *same motif parser* used by the downstream normalizer
and accepts a source only when exactly one candidate yields the declared number
of unique TFs with no conflicting duplicate representatives.

For CDNA-001 the declared contract is 1,421 unique human TF motifs.  Passing
this gate identifies the source archive; it is not biological evidence and does
not establish TF occupancy or causality.
"""
from __future__ import annotations

import argparse
import io
import json
from pathlib import Path
import zipfile

import gap001_fetch_codebook_atlas as fetcher
import gap001_prepare_codebook_atlas as prepare


def inspect_archive(data: bytes, expected_tfs: int) -> dict[str, object]:
    """Apply the downstream one-TF/one-PWM parsing contract without writing files."""
    parsed: dict[str, tuple[list[list[float]], str]] = {}
    duplicate_count = 0
    conflicting_duplicate_count = 0
    source_member_count = 0
    unparsed_member_count = 0

    with zipfile.ZipFile(io.BytesIO(data)) as zf:
        for info in zf.infolist():
            if info.is_dir():
                continue
            source_member_count += 1
            motifs = prepare.parse_member(info.filename, zf.read(info))
            if not motifs:
                unparsed_member_count += 1
                continue
            for tf, rows in motifs:
                if len(rows) < 2 or any(len(row) != 4 for row in rows):
                    continue
                if tf in parsed:
                    duplicate_count += 1
                    if prepare.matrix_signature(parsed[tf][0]) != prepare.matrix_signature(rows):
                        conflicting_duplicate_count += 1
                else:
                    parsed[tf] = (rows, info.filename)

    parsed_unique_tfs = len(parsed)
    return {
        "expected_unique_tfs": expected_tfs,
        "parsed_unique_tfs": parsed_unique_tfs,
        "source_member_count": source_member_count,
        "unparsed_member_count": unparsed_member_count,
        "duplicate_count": duplicate_count,
        "conflicting_duplicate_count": conflicting_duplicate_count,
        "full_atlas_contract_pass": (
            parsed_unique_tfs == expected_tfs and conflicting_duplicate_count == 0
        ),
    }


def candidate_summary(candidate: dict[str, object]) -> dict[str, object]:
    return {k: v for k, v in candidate.items() if k not in {"data", "members"}}


def select_by_contract(candidates: list[dict[str, object]], expected_tfs: int) -> dict[str, object]:
    """Select exactly one motif-like candidate that satisfies the full atlas contract."""
    evaluated: list[dict[str, object]] = []
    for candidate in candidates:
        if not candidate.get("looks_motif_like", False):
            candidate["atlas_contract"] = {
                "expected_unique_tfs": expected_tfs,
                "full_atlas_contract_pass": False,
                "not_evaluated_reason": "member inventory is not motif-like",
            }
            evaluated.append(candidate)
            continue
        contract = inspect_archive(bytes(candidate["data"]), expected_tfs)
        candidate["atlas_contract"] = contract
        evaluated.append(candidate)

    passing = [
        candidate
        for candidate in evaluated
        if bool(candidate.get("atlas_contract", {}).get("full_atlas_contract_pass"))
    ]
    if len(passing) != 1:
        inventory = [candidate_summary(x) for x in evaluated]
        raise RuntimeError(
            "full-atlas source selection requires exactly one exact-DOI candidate "
            f"to satisfy the {expected_tfs}-TF parser contract; passing={len(passing)} "
            f"inventory={json.dumps(inventory, sort_keys=True)}"
        )
    return passing[0]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--article-url", default=fetcher.DEFAULT_ARTICLE)
    parser.add_argument("--doi", default=fetcher.DEFAULT_DOI)
    parser.add_argument("--expected-tfs", type=int, default=1421)
    parser.add_argument("--springer-probe-max", type=int, default=40)
    parser.add_argument("--output-zip", type=Path, required=True)
    parser.add_argument("--provenance-json", type=Path, required=True)
    args = parser.parse_args()

    page_bytes, article_final = fetcher.fetch_bytes(args.article_url)
    page_text = page_bytes.decode("utf-8", errors="replace")
    candidates: list[dict[str, object]] = []
    article_anchor_error: str | None = None

    # Prefer an explicitly labelled Supplementary Data 1 archive when Nature
    # exposes it, but still require the same 1,421-TF parser contract.
    try:
        explicit_url = fetcher.resolve_supplementary_data_1(article_final, page_text)
        explicit_data, explicit_final = fetcher.fetch_bytes(explicit_url)
        explicit_members = fetcher.zip_members_from_bytes(explicit_data)
        motif_named, motif_ext, looks_like = fetcher.motif_likeness(explicit_members)
        candidates.append({
            "index": "article_explicit_anchor",
            "origin": "article_explicit_anchor",
            "url": explicit_url,
            "final_url": explicit_final,
            "data": explicit_data,
            "sha256": fetcher.sha256_bytes(explicit_data),
            "size_bytes": len(explicit_data),
            "member_count": len(explicit_members),
            "motif_named_members": motif_named,
            "motif_extension_members": motif_ext,
            "looks_motif_like": looks_like,
            "members": explicit_members,
        })
    except RuntimeError as exc:
        article_anchor_error = str(exc)

    # Exact DOI fallback: enumerate Springer's ZIP media objects, then let the
    # parser contract—not filename guessing—decide which archive is the atlas.
    for candidate in fetcher.discover_springer_zip_candidates(args.doi, args.springer_probe_max):
        candidate["origin"] = "springer_exact_doi_media_object"
        candidates.append(candidate)

    # Deduplicate identical bytes that may be reachable from both routes.
    unique: dict[str, dict[str, object]] = {}
    for candidate in candidates:
        unique.setdefault(str(candidate["sha256"]), candidate)
    candidates = list(unique.values())

    chosen = select_by_contract(candidates, args.expected_tfs)
    archive = bytes(chosen["data"])
    args.output_zip.parent.mkdir(parents=True, exist_ok=True)
    args.output_zip.write_bytes(archive)
    persisted_members = fetcher.validate_zip(args.output_zip)

    inventory = [candidate_summary(x) for x in candidates]
    provenance = {
        "article_requested_url": args.article_url,
        "article_final_url": article_final,
        "doi": args.doi,
        "supplement_label": "Supplementary Data 1 / full Codebook motif atlas",
        "selection_method": "exact DOI candidates + downstream 1421-TF parser contract",
        "article_anchor_error": article_anchor_error,
        "expected_unique_tfs": args.expected_tfs,
        "candidate_count": len(candidates),
        "candidate_inventory": inventory,
        "selected_origin": chosen.get("origin"),
        "selected_index": chosen.get("index"),
        "selected_url": chosen.get("url"),
        "selected_final_url": chosen.get("final_url"),
        "selected_sha256": fetcher.sha256_bytes(archive),
        "selected_size_bytes": len(archive),
        "selected_member_count": len(persisted_members),
        "selected_atlas_contract": chosen.get("atlas_contract"),
        "downstream_guard": "normalizer independently re-runs the same 1,421 unique-TF contract",
        "claim_scope": "source selection/provenance only; not TF occupancy or causality",
    }
    args.provenance_json.parent.mkdir(parents=True, exist_ok=True)
    args.provenance_json.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(provenance, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
