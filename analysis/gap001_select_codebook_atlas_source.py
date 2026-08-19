#!/usr/bin/env python3
"""Select the official Nature Codebook motif archive by the full-atlas contract.

Nature's article page may be reduced for automated clients and Springer exposes
multiple ZIP media objects for the exact DOI. Filename/member heuristics are
not sufficient provenance. This selector evaluates exact-DOI motif-like ZIPs
with the same parser used by the downstream normalizer and accepts a source
only when exactly one candidate covers the declared number of *unique TFs*.

Supplementary Data 1 can contain multiple distinct published motif
representatives for one TF. Those are not conflicts: all representatives are
preserved downstream and scanned independently. The atlas contract is therefore
1,421 unique TFs plus >=1 parsed motif record per TF, not one-TF/one-PWM.

Passing this gate identifies a source archive; it is not biological evidence
and does not establish TF occupancy or causality.
"""
from __future__ import annotations

import argparse
import io
import json
import zipfile

import gap001_fetch_codebook_atlas as fetcher
import gap001_prepare_codebook_atlas as prepare


def _text_head(data: bytes, limit: int = 500) -> str:
    return data[:limit].decode("utf-8", errors="replace").replace("\x00", "\\0")


def inspect_archive(data: bytes, expected_tfs: int) -> dict[str, object]:
    """Apply downstream parsing semantics without writing normalized files."""
    by_tf: dict[str, int] = {}
    motif_record_count = 0
    source_member_count = 0
    unparsed_member_count = 0
    unparsed_examples: list[dict[str, object]] = []

    with zipfile.ZipFile(io.BytesIO(data)) as zf:
        for info in zf.infolist():
            if info.is_dir():
                continue
            source_member_count += 1
            member_data = zf.read(info)
            motifs = prepare.parse_member(info.filename, member_data)
            if not motifs:
                unparsed_member_count += 1
                if len(unparsed_examples) < 5:
                    unparsed_examples.append({
                        "name": info.filename,
                        "size": info.file_size,
                        "text_head": _text_head(member_data),
                    })
                continue
            for tf, rows in motifs:
                if len(rows) < 2 or any(len(row) != 4 for row in rows):
                    continue
                motif_record_count += 1
                by_tf[tf] = by_tf.get(tf, 0) + 1

    parsed_unique_tfs = len(by_tf)
    multiple = {tf: count for tf, count in by_tf.items() if count > 1}
    return {
        "expected_unique_tfs": expected_tfs,
        "parsed_unique_tfs": parsed_unique_tfs,
        "parsed_motif_records": motif_record_count,
        "source_member_count": source_member_count,
        "unparsed_member_count": unparsed_member_count,
        "unparsed_examples": unparsed_examples,
        "multiple_representative_tf_count": len(multiple),
        "max_representatives_per_tf": max(by_tf.values(), default=0),
        "full_atlas_contract_pass": (
            parsed_unique_tfs == expected_tfs and motif_record_count >= parsed_unique_tfs
        ),
        "contract_semantics": "unique TF coverage; multiple published motif representatives are preserved",
    }


def candidate_summary(candidate: dict[str, object]) -> dict[str, object]:
    return {k: v for k, v in candidate.items() if k not in {"data", "members"}}


def select_by_contract(candidates: list[dict[str, object]], expected_tfs: int) -> dict[str, object]:
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
        candidate["atlas_contract"] = inspect_archive(bytes(candidate["data"]), expected_tfs)
        evaluated.append(candidate)

    passing = [
        candidate for candidate in evaluated
        if bool(candidate.get("atlas_contract", {}).get("full_atlas_contract_pass"))
    ]
    if len(passing) != 1:
        inventory = [candidate_summary(x) for x in evaluated]
        raise RuntimeError(
            "full-atlas source selection requires exactly one exact-DOI candidate "
            f"to satisfy the {expected_tfs}-unique-TF parser contract; passing={len(passing)} "
            f"inventory={json.dumps(inventory, sort_keys=True)}"
        )
    return passing[0]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--article-url", default=fetcher.DEFAULT_ARTICLE)
    parser.add_argument("--doi", default=fetcher.DEFAULT_DOI)
    parser.add_argument("--expected-tfs", type=int, default=1421)
    parser.add_argument("--springer-probe-max", type=int, default=40)
    parser.add_argument("--output-zip", required=True)
    parser.add_argument("--provenance-json", required=True)
    args = parser.parse_args()

    from pathlib import Path
    output_zip = Path(args.output_zip)
    provenance_json = Path(args.provenance_json)

    page_bytes, article_final = fetcher.fetch_bytes(args.article_url)
    page_text = page_bytes.decode("utf-8", errors="replace")
    candidates: list[dict[str, object]] = []
    article_anchor_error: str | None = None

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

    for candidate in fetcher.discover_springer_zip_candidates(args.doi, args.springer_probe_max):
        candidate["origin"] = "springer_exact_doi_media_object"
        candidates.append(candidate)

    unique: dict[str, dict[str, object]] = {}
    for candidate in candidates:
        unique.setdefault(str(candidate["sha256"]), candidate)
    candidates = list(unique.values())

    chosen = select_by_contract(candidates, args.expected_tfs)
    archive = bytes(chosen["data"])
    output_zip.parent.mkdir(parents=True, exist_ok=True)
    output_zip.write_bytes(archive)
    persisted_members = fetcher.validate_zip(output_zip)

    provenance = {
        "article_requested_url": args.article_url,
        "article_final_url": article_final,
        "doi": args.doi,
        "supplement_label": "Supplementary Data 1 / full Codebook motif atlas",
        "selection_method": "exact DOI candidates + downstream unique-TF parser contract",
        "article_anchor_error": article_anchor_error,
        "expected_unique_tfs": args.expected_tfs,
        "candidate_count": len(candidates),
        "candidate_inventory": [candidate_summary(x) for x in candidates],
        "selected_origin": chosen.get("origin"),
        "selected_index": chosen.get("index"),
        "selected_url": chosen.get("url"),
        "selected_final_url": chosen.get("final_url"),
        "selected_sha256": fetcher.sha256_bytes(archive),
        "selected_size_bytes": len(archive),
        "selected_member_count": len(persisted_members),
        "selected_atlas_contract": chosen.get("atlas_contract"),
        "downstream_guard": "normalizer independently re-runs the same 1,421 unique-TF contract and preserves all representatives",
        "claim_scope": "source selection/provenance only; not TF occupancy or causality",
    }
    provenance_json.parent.mkdir(parents=True, exist_ok=True)
    provenance_json.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(provenance, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
