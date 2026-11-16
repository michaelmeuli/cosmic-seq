"""Plain text of the lab handbook files (`data/Arbeitshandbücher`), .docx and legacy .doc.

No Word converter is installed on the cluster, so both formats are read directly:
.docx is a zip whose `word/document.xml` holds the text; .doc is an OLE2 file in which Word stores
the text as runs of cp1252 bytes or of UTF-16LE code units, which are collected as they stand.
The result is good enough to find primer names and sequences, not to reproduce the layout.
"""
from __future__ import annotations

import re
import zipfile
from pathlib import Path

_PARAGRAPH_END = re.compile(r"</w:p>|</w:tr>|<w:br[^>]*>")
_CELL_END = re.compile(r"</w:tc>|<w:tab[^>]*>")
_TAG = re.compile(r"<[^>]+>")
_ENTITIES = {"&amp;": "&", "&lt;": "<", "&gt;": ">", "&quot;": '"', "&apos;": "'"}

# bytes that can be part of a text run in a .doc: printable ASCII, Latin-1 letters, and the cp1252
# typographic quotes and dashes (0x91-0x97) that Word puts around 5'/3'
_DOC_BYTE = rb"[\x20-\x7e\x91-\x97\xa0-\xff]"
# Word's own separators inside a run: tab, paragraph mark (\r), cell mark (\x07), line break (\x0b)
_DOC_SEP = rb"[\t\r\x07\x0b]"
_DOC_RUN_8BIT = re.compile(rb"(?:" + _DOC_BYTE + rb"|" + _DOC_SEP + rb"){6,}")
_DOC_RUN_16BIT = re.compile(rb"(?:(?:" + _DOC_BYTE + rb"|" + _DOC_SEP + rb")\x00|[\x18-\x1d\x13\x14\x22\x26]\x20){6,}")
MIN_RUN_LETTERS = 4
MIN_TEXT_FRACTION = 0.8  # share of ordinary characters a .doc line needs to count as text


def docx_text(path: Path) -> str:
    """Text of a .docx: one line per paragraph or table row, cells separated by tabs."""
    with zipfile.ZipFile(path) as z:
        parts = [z.read(name).decode("utf-8", "ignore") for name in z.namelist()
                 if re.fullmatch(r"word/(document|footnotes|endnotes)\.xml", name)]
    xml = "\n".join(parts)
    xml = _CELL_END.sub("\t", xml)
    xml = _PARAGRAPH_END.sub("\n", xml)
    text = _TAG.sub("", xml)
    for entity, char in _ENTITIES.items():
        text = text.replace(entity, char)
    return _tidy(text)


def doc_text(path: Path) -> str:
    """Text of a legacy .doc: every 8-bit and UTF-16LE text run of the file, in file order."""
    raw = path.read_bytes()
    runs: list[tuple[int, str]] = []
    covered: list[tuple[int, int]] = []
    for m in _DOC_RUN_16BIT.finditer(raw):
        runs.append((m.start(), m.group().decode("utf-16-le", "ignore")))
        covered.append((m.start(), m.end()))
    for m in _DOC_RUN_8BIT.finditer(raw):
        if any(start <= m.start() < end for start, end in covered):
            continue
        runs.append((m.start(), m.group().decode("cp1252", "ignore")))
    text = "\n".join(text for _, text in sorted(runs))
    text = text.replace("\r", "\n").replace("\x07", "\t").replace("\x0b", "\n")
    return _tidy("\n".join(line for line in text.split("\n") if _is_text(line)))


def _is_text(line: str) -> bool:
    """False for the binary noise between the text runs of a .doc (few letters, or mostly non-ASCII characters)."""
    stripped = line.strip()
    if sum(c.isalpha() for c in stripped) < MIN_RUN_LETTERS:
        return False
    return sum(c.isascii() or c in "äöüÄÖÜéèàç°µ’‘`´–—" for c in stripped) >= MIN_TEXT_FRACTION * len(stripped)


def _tidy(text: str) -> str:
    lines = [re.sub(r"[ \xa0]+", " ", line).strip(" ") for line in text.split("\n")]
    return "\n".join(line for line in lines if line.strip())


def handbook_text(path: Path) -> str:
    """Text of a handbook file; .docx and .doc only."""
    suffix = path.suffix.lower()
    if suffix == ".docx":
        return docx_text(path)
    if suffix == ".doc":
        return doc_text(path)
    raise ValueError(f"not a Word file: {path}")


def handbook_files(root: Path) -> list[Path]:
    """Every .docx/.doc below `root`, without Word and LibreOffice lock files (`~$...`, `.~lock...`)."""
    return sorted(p for p in root.rglob("*") if p.is_file() and p.suffix.lower() in {".docx", ".doc"}
                  and not p.name.startswith(("~$", ".~lock")))
