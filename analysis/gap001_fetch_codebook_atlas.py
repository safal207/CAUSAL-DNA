#!/usr/bin/env python3
"""Fetch Nature Codebook Supplementary Data 1 with provenance.

The Nature paper states that Supplementary Data 1 contains representative PWMs
for all 1,421 human TFs with characterized sequence specificity. This utility
resolves the Supplementary Data 1 link from the article HTML, downloads the
archive, validates it as ZIP, and records a SHA256/member inventory.

This is a provenance/download utility only. It makes no biological claim.
"""
from __future__ import annotations

import argparse
import hashlib
from html.parser import HTMLParser
import html
import json
from pathlib import Path
from urllib.parse import urljoin
from urllib.request import Request, urlopen
import zipfile

DEFAULT_ARTICLE = "https://www.nature.com/articles/s41586-026-10798-9"
USER_AGENT = "Mozilla/5.0 (compatible; CAUSAL-DNA/0.1; +https://github.com/safal207/CAUSAL-DNA)"


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


def fetch_bytes(url: str) -> tuple[bytes, str]:
    req = Request(url, headers={"User-Agent": USER_AGENT, "Accept": "*/*"})
    with urlopen(req, timeout=90) as response:
        return response.read(), response.geturl()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def validate_zip(path: Path) -> list[dict[str, object]]:
    if not zipfile.is_zipfile(path):
        raise RuntimeError(f"downloaded Supplementary Data 1 is not a ZIP archive: {path}")
    members: list[dict[str, object]] = []
    with zipfile.ZipFile(path) as zf:
        bad = zf.testzip()
        if bad is not None:
            raise RuntimeError(f"corrupt ZIP member: {bad}")
        for info in zf.infolist():
            if not info.is_dir():
                members.append({"name": info.filename, "size": info.file_size})
    if not members:
        raise RuntimeError("Supplementary Data 1 ZIP is empty")
    return members


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--article-url", default=DEFAULT_ARTICLE)
    parser.add_argument("--output-zip", type=Path, required=True)
    parser.add_argument("--provenance-json", type=Path, required=True)
    args = parser.parse_args()

    page_bytes, article_final = fetch_bytes(args.article_url)
    page_text = page_bytes.decode("utf-8", errors="replace")
    supplement_url = resolve_supplementary_data_1(article_final, page_text)
    archive, supplement_final = fetch_bytes(supplement_url)

    args.output_zip.parent.mkdir(parents=True, exist_ok=True)
    args.output_zip.write_bytes(archive)
    members = validate_zip(args.output_zip)

    provenance = {
        "article_requested_url": args.article_url,
        "article_final_url": article_final,
        "supplement_label": "Supplementary Data 1",
        "supplement_requested_url": supplement_url,
        "supplement_final_url": supplement_final,
        "sha256": sha256_bytes(archive),
        "size_bytes": len(archive),
        "member_count": len(members),
        "members": members,
        "expected_scope_from_article": "representative motifs for all 1,421 human TFs",
        "claim_scope": "download/provenance only; not biological evidence",
    }
    args.provenance_json.parent.mkdir(parents=True, exist_ok=True)
    args.provenance_json.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(provenance, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
