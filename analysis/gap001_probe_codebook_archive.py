#!/usr/bin/env python3
"""Inspect an official Codebook supplementary ZIP before extending parsers.

Outputs member-extension counts, a deterministic member-name sample, and short
text previews for small candidate motif files. This is format discovery only;
it never interprets a motif or changes the 1,421-TF full-atlas contract.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
import zipfile

TEXT_EXTS = {".txt", ".pwm", ".pfm", ".meme", ".jaspar", ".transfac", ".tsv", ".csv"}


def preview_text(data: bytes, limit_lines: int = 20) -> list[str]:
    for encoding in ("utf-8-sig", "latin-1"):
        try:
            text = data.decode(encoding)
            break
        except UnicodeDecodeError:
            text = ""
    if not text:
        return ["<binary/unreadable>"]
    return [line[:500] for line in text.splitlines()[:limit_lines]]


def inspect_archive(path: Path, sample_count: int = 30, preview_count: int = 8) -> dict[str, object]:
    if not zipfile.is_zipfile(path):
        raise ValueError(f"not a ZIP: {path}")
    with zipfile.ZipFile(path) as zf:
        infos = [x for x in zf.infolist() if not x.is_dir()]
        ext_counts = Counter(Path(x.filename).suffix.lower() or "<none>" for x in infos)
        names = sorted(x.filename for x in infos)
        previews = []
        for name in names:
            if len(previews) >= preview_count:
                break
            info = zf.getinfo(name)
            if Path(name).suffix.lower() not in TEXT_EXTS or info.file_size > 2_000_000:
                continue
            previews.append({
                "name": name,
                "size": info.file_size,
                "lines": preview_text(zf.read(name)),
            })
        return {
            "archive": str(path),
            "member_count": len(infos),
            "extension_counts": dict(sorted(ext_counts.items())),
            "member_name_sample": names[:sample_count],
            "text_previews": previews,
            "claim_scope": "archive format inspection only; no motif or biological inference",
        }


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--zip", type=Path, required=True)
    p.add_argument("--output-json", type=Path, required=True)
    p.add_argument("--sample-count", type=int, default=30)
    p.add_argument("--preview-count", type=int, default=8)
    args = p.parse_args()
    result = inspect_archive(args.zip, args.sample_count, args.preview_count)
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
