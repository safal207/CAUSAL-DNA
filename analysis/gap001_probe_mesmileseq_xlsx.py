#!/usr/bin/env python3
"""Inspect meSMiLE-seq Supplementary Data 5 XLSX without spreadsheet deps.

The goal is deliberately narrow: expose workbook sheet names plus a small,
faithful row sample so a downstream methylation-class parser can be written
against the real published table structure. XLSX is an OpenXML ZIP; this tool
uses only zipfile + ElementTree and never infers biological classes.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
from xml.etree import ElementTree as ET
import zipfile

NS = {
    "main": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
    "rel": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
    "pkgrel": "http://schemas.openxmlformats.org/package/2006/relationships",
}
CELL_RE = re.compile(r"([A-Z]+)(\d+)$")


def column_index(cell_ref: str) -> int:
    match = CELL_RE.fullmatch(cell_ref)
    if not match:
        raise ValueError(f"invalid cell reference: {cell_ref}")
    letters = match.group(1)
    value = 0
    for ch in letters:
        value = value * 26 + (ord(ch) - ord("A") + 1)
    return value - 1


def read_shared_strings(zf: zipfile.ZipFile) -> list[str]:
    name = "xl/sharedStrings.xml"
    if name not in zf.namelist():
        return []
    root = ET.fromstring(zf.read(name))
    values: list[str] = []
    for item in root.findall("main:si", NS):
        parts = [node.text or "" for node in item.findall(".//main:t", NS)]
        values.append("".join(parts))
    return values


def workbook_sheet_targets(zf: zipfile.ZipFile) -> list[tuple[str, str]]:
    workbook = ET.fromstring(zf.read("xl/workbook.xml"))
    rels = ET.fromstring(zf.read("xl/_rels/workbook.xml.rels"))
    targets = {
        rel.attrib["Id"]: rel.attrib["Target"]
        for rel in rels.findall("pkgrel:Relationship", NS)
    }
    out: list[tuple[str, str]] = []
    for sheet in workbook.findall("main:sheets/main:sheet", NS):
        title = sheet.attrib["name"]
        rel_id = sheet.attrib[f"{{{NS['rel']}}}id"]
        target = targets[rel_id].lstrip("/")
        if not target.startswith("xl/"):
            target = "xl/" + target
        out.append((title, target))
    return out


def cell_value(cell: ET.Element, shared: list[str]) -> object:
    cell_type = cell.attrib.get("t")
    if cell_type == "inlineStr":
        return "".join(node.text or "" for node in cell.findall(".//main:t", NS))
    value_node = cell.find("main:v", NS)
    if value_node is None or value_node.text is None:
        return ""
    raw = value_node.text
    if cell_type == "s":
        idx = int(raw)
        return shared[idx]
    if cell_type == "b":
        return raw == "1"
    if cell_type in {"str", "e"}:
        return raw
    try:
        number = float(raw)
        return int(number) if number.is_integer() else number
    except ValueError:
        return raw


def read_rows(zf: zipfile.ZipFile, sheet_path: str, shared: list[str], max_rows: int) -> list[list[object]]:
    root = ET.fromstring(zf.read(sheet_path))
    rows: list[list[object]] = []
    for row in root.findall("main:sheetData/main:row", NS):
        values: dict[int, object] = {}
        max_col = -1
        for cell in row.findall("main:c", NS):
            ref = cell.attrib.get("r")
            if not ref:
                continue
            idx = column_index(ref)
            values[idx] = cell_value(cell, shared)
            max_col = max(max_col, idx)
        rows.append([values.get(i, "") for i in range(max_col + 1)] if max_col >= 0 else [])
        if len(rows) >= max_rows:
            break
    return rows


def inspect_xlsx(path: Path, max_rows: int = 25) -> dict[str, object]:
    if not zipfile.is_zipfile(path):
        raise ValueError(f"not an XLSX/OpenXML ZIP: {path}")
    with zipfile.ZipFile(path) as zf:
        required = {"xl/workbook.xml", "xl/_rels/workbook.xml.rels"}
        missing = required.difference(zf.namelist())
        if missing:
            raise ValueError(f"missing XLSX workbook parts: {sorted(missing)}")
        shared = read_shared_strings(zf)
        sheets = []
        for title, target in workbook_sheet_targets(zf):
            sheets.append({
                "title": title,
                "path": target,
                "sample_rows": read_rows(zf, target, shared, max_rows),
            })
        return {
            "workbook": str(path),
            "shared_string_count": len(shared),
            "sheet_count": len(sheets),
            "sheets": sheets,
            "claim_scope": "workbook structure inspection only; no methylation class inferred",
        }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--xlsx", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--max-rows", type=int, default=25)
    args = parser.parse_args()
    result = inspect_xlsx(args.xlsx, args.max_rows)
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
