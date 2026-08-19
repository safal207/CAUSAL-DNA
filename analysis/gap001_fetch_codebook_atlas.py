#!/usr/bin/env python3
"""Fetch Nature Codebook Supplementary Data 1 with provenance.

The Nature paper states that Supplementary Data 1 contains representative PWMs
for all 1,421 human TFs with characterized sequence specificity. The preferred
route resolves the explicitly labelled link from article HTML. Some automated
clients receive a reduced Nature page without supplementary anchors; in that
case this utility probes only Springer media-object ZIPs belonging to the exact
DOI and selects a candidate only when its member inventory is motif-like.

The downstream atlas normalizer independently requires exactly 1,421 unique TFs,
so media-object discovery cannot silently turn a different ZIP into the atlas.
This utility is provenance/download plumbing only; it makes no biological claim.
"""
from __future__ import annotations

import argparse
import hashlib
from html.parser import HTMLParser
import html
import io
import json
from pathlib import Path
import re
from urllib.error import HTTPError, URLError
from urllib.parse import urljoin
from urllib.request import Request, urlopen
import zipfile

DEFAULT_ARTICLE = "https://www.nature.com/articles/s41586-026-10798-9"
DEFAULT_DOI = "10.1038/s41586-026-10798-9"
USER_AGENT = "Mozilla/5.0 (compatible; CAUSAL-DNA/0.1; +https://github.com/safal207/CAUSAL-DNA)"
MOTIF_EXTENSIONS = {".pwm", ".pfm", ".meme", ".jaspar", ".transfac", ".txt"}


class AnchorCollector(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.current_href: str | None = None
        self.current_text: list[str] = []
        self.anchors: list[tuple[str, str]] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag.lower() == "a":
            data = dict(attrs)
            self.current_href = data.get("href")
            self.current_text = []

    def handle_data(self, data: str) -> None:
        if self.current_href is not None:
            self.current_text.append(data)

    def handle_endtag(self, tag: str) -> None:
        if tag.lower() == "a" and self.current_href is not None:
            self.anchors.append((self.current_href, " ".join(self.current_text).strip()))
            self.current_href = None
            self.current_text = []


def resolve_supplementary_data_1(article_url: str, page_html: str) -> str:
    parser = AnchorCollector()
    parser.feed(page_html)
    candidates: list[tuple[int, str, str]] = []
    for href, text in parser.anchors:
        href = html.unescape(href)
        norm_text = " ".join(text.lower().split())
        norm_href = href.lower()
        score = 0
        if "supplementary data 1" in norm_text:
            score += 100
        if "supplementary" in norm_text and "data" in norm_text:
            score += 20
        if norm_href.endswith(".zip") or ".zip?" in norm_href:
            score += 10
        if "mediaobject" in norm_href or "springernature" in norm_href or "static-content.springer" in norm_href:
            score += 5
        if score:
            candidates.append((score, urljoin(article_url, href), text))
    if not candidates:
        raise RuntimeError("could not resolve Supplementary Data 1 link from article HTML")
    candidates.sort(key=lambda x: x[0], reverse=True)
    best = candidates[0]
    if best[0] < 100:
        raise RuntimeError(f"no anchor explicitly labelled Supplementary Data 1; best={best}")
    return best[1]


def fetch_bytes(url: str, timeout: int = 90) -> tuple[bytes, str]:
    req = Request(url, headers={"User-Agent": USER_AGENT, "Accept": "*/*"})
    with urlopen(req, timeout=timeout) as response:
        return response.read(), response.geturl()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def zip_members_from_bytes(data: bytes) -> list[dict[str, object]]:
    if not zipfile.is_zipfile(io.BytesIO(data)):
        raise RuntimeError("candidate payload is not a ZIP archive")
    members: list[dict[str, object]] = []
    with zipfile.ZipFile(io.BytesIO(data)) as zf:
        bad = zf.testzip()
        if bad is not None:
            raise RuntimeError(f"corrupt ZIP member: {bad}")
        for info in zf.infolist():
            if not info.is_dir():
                members.append({"name": info.filename, "size": info.file_size})
    if not members:
        raise RuntimeError("ZIP is empty")
    return members


def validate_zip(path: Path) -> list[dict[str, object]]:
    return zip_members_from_bytes(path.read_bytes())


def springer_media_prefix(doi: str) -> tuple[str, str]:
    """Derive Springer media-object prefix from the exact Nature DOI.

    Modern Nature article DOIs encode years as three digits after the journal
    code (for example ``026`` for 2026), while older/synthetic forms may use a
    four-digit year.  The media-object filename uses the four-digit year.
    """
    match = re.fullmatch(r"10\.1038/s(\d+)-(\d{3,4})-(\d+)-([\w-]+)", doi)
    if not match:
        raise ValueError(f"unsupported Nature DOI shape for media-object discovery: {doi}")
    journal, year_token, article_no, _ = match.groups()
    if len(year_token) == 3:
        year = str(2000 + int(year_token))
    elif len(year_token) == 4:
        year = year_token
    else:  # guarded by regex; retained as an explicit provenance invariant
        raise ValueError(f"unsupported Nature DOI year token: {year_token}")
    stem = f"{journal}_{year}_{int(article_no)}"
    base = f"https://static-content.springer.com/esm/art%3A{doi.replace('/', '%2F')}/MediaObjects"
    return base, stem


def motif_likeness(members: list[dict[str, object]]) -> tuple[int, int, bool]:
    names = [str(x["name"]).lower() for x in members]
    motif_named = sum("motif" in name or "pwm" in name for name in names)
    motif_ext = sum(Path(name).suffix.lower() in MOTIF_EXTENSIONS for name in names)
    # Supplementary Data 1 can be either many per-TF files or one/few aggregate
    # motif files. Require explicit motif/PWM naming for small archives; large
    # motif-like text collections also qualify.
    looks_like = motif_named > 0 or motif_ext >= 100
    return motif_named, motif_ext, looks_like


def discover_springer_zip_candidates(doi: str, max_index: int = 40) -> list[dict[str, object]]:
    base, stem = springer_media_prefix(doi)
    found: list[dict[str, object]] = []
    for index in range(1, max_index + 1):
        url = f"{base}/{stem}_MOESM{index}_ESM.zip"
        try:
            data, final_url = fetch_bytes(url, timeout=30)
        except HTTPError as exc:
            if exc.code in {403, 404, 410}:
                continue
            raise
        except URLError:
            continue
        if not zipfile.is_zipfile(io.BytesIO(data)):
            continue
        members = zip_members_from_bytes(data)
        motif_named, motif_ext, looks_like = motif_likeness(members)
        found.append({
            "index": index,
            "url": url,
            "final_url": final_url,
            "data": data,
            "sha256": sha256_bytes(data),
            "size_bytes": len(data),
            "member_count": len(members),
            "motif_named_members": motif_named,
            "motif_extension_members": motif_ext,
            "looks_motif_like": looks_like,
            "members": members,
        })
    return found


def select_motif_zip(candidates: list[dict[str, object]]) -> dict[str, object]:
    motif_like = [x for x in candidates if x["looks_motif_like"]]
    if len(motif_like) != 1:
        inventory = [
            {k: v for k, v in x.items() if k not in {"data", "members"}}
            for x in candidates
        ]
        raise RuntimeError(
            "Springer media-object discovery did not yield exactly one motif-like ZIP; "
            f"motif_like={len(motif_like)} inventory={json.dumps(inventory, sort_keys=True)}"
        )
    return motif_like[0]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--article-url", default=DEFAULT_ARTICLE)
    parser.add_argument("--doi", default=DEFAULT_DOI)
    parser.add_argument("--output-zip", type=Path, required=True)
    parser.add_argument("--provenance-json", type=Path, required=True)
    parser.add_argument("--springer-probe-max", type=int, default=40)
    args = parser.parse_args()

    page_bytes, article_final = fetch_bytes(args.article_url)
    page_text = page_bytes.decode("utf-8", errors="replace")
    resolution_method = "article_explicit_anchor"
    discovery_inventory: list[dict[str, object]] = []

    try:
        supplement_url = resolve_supplementary_data_1(article_final, page_text)
        archive, supplement_final = fetch_bytes(supplement_url)
        members = zip_members_from_bytes(archive)
    except RuntimeError as anchor_error:
        candidates = discover_springer_zip_candidates(args.doi, args.springer_probe_max)
        discovery_inventory = [
            {k: v for k, v in x.items() if k not in {"data", "members"}}
            for x in candidates
        ]
        chosen = select_motif_zip(candidates)
        archive = bytes(chosen["data"])
        supplement_url = str(chosen["url"])
        supplement_final = str(chosen["final_url"])
        members = list(chosen["members"])
        resolution_method = "springer_exact_doi_media_object_discovery"
        anchor_error_text = str(anchor_error)
    else:
        anchor_error_text = None

    args.output_zip.parent.mkdir(parents=True, exist_ok=True)
    args.output_zip.write_bytes(archive)
    # Re-read persisted bytes so provenance covers exactly what downstream sees.
    members = validate_zip(args.output_zip)

    provenance = {
        "article_requested_url": args.article_url,
        "article_final_url": article_final,
        "doi": args.doi,
        "supplement_label": "Supplementary Data 1 / motif ZIP candidate",
        "resolution_method": resolution_method,
        "article_anchor_error": anchor_error_text,
        "supplement_requested_url": supplement_url,
        "supplement_final_url": supplement_final,
        "sha256": sha256_bytes(archive),
        "size_bytes": len(archive),
        "member_count": len(members),
        "members": members,
        "springer_zip_inventory": discovery_inventory,
        "expected_scope_from_article": "representative motifs for all 1,421 human TFs",
        "downstream_guard": "atlas normalizer must independently parse exactly 1,421 unique TFs",
        "claim_scope": "download/provenance only; not biological evidence",
    }
    args.provenance_json.parent.mkdir(parents=True, exist_ok=True)
    args.provenance_json.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in provenance.items() if k not in {"members"}}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
