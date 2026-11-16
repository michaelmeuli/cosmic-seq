#!/usr/bin/env python
"""Check the curated primer table against the lab handbooks and write docs/PRIMERS.md.

`primers.tsv` and `amplicons.tsv` were written by reading the mycobacterial Sanger protocols in
`data/Arbeitshandbücher` (notation differs too much between protocols to parse them blindly). This script
keeps the tables honest:

1. every primer sequence must occur in the text of the protocol file it names (as written or as reverse
   complement; case, spaces and line breaks ignored) - exit code 1 otherwise;
2. every oligo written with 5'/3' ends in a mycobacterial protocol that is not in `primers.tsv` is listed
   (probes, real-time PCR and fragment-length assays, or primers still to be curated);
3. each amplicon is run in silico on M. tuberculosis H37Rv, to compare with the product length of the protocol;
4. docs/PRIMERS.md is written from the tables and these three results.

Output directory: handbook_oligos.tsv (every oligo found, with its status), primer_check.tsv.

Usage: 02_handbook_primers.py [--handbook DIR] [--h37rv DIR] [--no-h37rv] [--out DIR] [--md FILE]
"""
from __future__ import annotations

import argparse
import csv
import datetime
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # repo root: the cosmic_seq package

from cosmic_seq import AMPLICONS_TSV, H37RV_ACCESSION, OUTPUT_PRIMERS, PRIMERS_MD, PRIMERS_TSV, REPO  # noqa: E402
from cosmic_seq.handbook import handbook_files, handbook_text  # noqa: E402
from cosmic_seq.insilico import IUPAC, MAX_MISMATCHES, MAX_PRODUCT, Product, revcomp  # noqa: E402
from cosmic_seq.kansasii_paths import HANDBOOK_DIR, MANIFEST_SCRIPT, NCBI_DATASET  # noqa: E402
from cosmic_seq.tables import Amplicon, Primer, genome_fasta, genome_products, load_amplicons, load_primers  # noqa: E402

# Protocols searched for oligos that are not in primers.tsv: the mycobacterial chapters of the routine
# handbook (by words of the file name) and the folder of mycobacterial research protocols.
ROUTINE_DIR = "32.43000-Routine"
MYCOBACTERIAL_NAME = re.compile(r"Mtb|MTB|Mykobakt|NTM|16S rRNA|16S rDNA|M\. ?chimaera|M\. ?leprae")
RESEARCH_DIR = "NZM Mykobakteriologie"
# an oligo as the protocols write it: 5' marker, optional label (FAM-, LC610-), bases in groups, 3' marker
_BASES = "ACGTRYSWKMBDHVNacgtryswkmbdhvn"
# The 5 must not be the end of a position number of a sequence map ("... 755 - att cga ...").
_QUOTE = "’'`´‘′"
OLIGO = re.compile(r"(?<![0-9])5\s?(?:[" + _QUOTE + r"]\s*[-–—]?|[-–—])\s*((?:[A-Za-z0-9]{2,6}\s?[-–—]\s?)?)((?:[" + _BASES + r"]\s?){15,60})")
MIN_OLIGO, MAX_OLIGO = 15, 45
SEQUENCING = re.compile(r"CycleSequencing|Sequenzierung")
# primer names of the reagent table of a research protocol ("Primer f: EC16S", "Primer J10F")
PRIMER_NAME = re.compile(r"Primer(?: [fr])?:? ([A-Za-z0-9][^\t\n/(]{1,30}?)\s*(?:\t|$)", re.MULTILINE)


def squeeze(text: str) -> str:
    """Upper case without whitespace: the form in which a primer is searched in a protocol."""
    return re.sub(r"\s+", "", text).upper()


def in_scope(path: Path) -> bool:
    if RESEARCH_DIR in path.parent.name:
        return True
    return ROUTINE_DIR in path.parts and bool(MYCOBACTERIAL_NAME.search(path.name))


def find_oligos(text: str) -> list[tuple[int, str]]:
    """(line number, sequence) of every oligo written with a 5' end. `.doc` files may break the line after the 5."""
    lines = text.split("\n")
    found: list[tuple[int, str]] = []
    for number, line in enumerate(lines, start=1):
        joined = line
        if re.search(r"(?<![0-9])5\s?[" + _QUOTE + r"]?\s*$", line) and number < len(lines):
            joined = line + lines[number]
        for match in OLIGO.finditer(joined):
            sequence = squeeze(match.group(2))
            if MIN_OLIGO <= len(sequence) <= MAX_OLIGO and len(set(sequence) & set("ACGT")) >= 3:
                found.append((number, sequence))
    return found


def status_of(sequence: str, primers: dict[str, Primer]) -> str:
    """Name of the curated primer an oligo of a protocol belongs to, or ''."""
    for primer in primers.values():
        for known in (primer.sequence, revcomp(primer.sequence)):
            if sequence == known or (len(sequence) >= MIN_OLIGO and (sequence in known or known in sequence)):
                return primer.name
    return ""


def describe(products: list[Product]) -> str:
    """Product length(s) on the reference genome, with the primer mismatches when there are any."""
    if not products:
        return "no product"
    lengths = sorted({p.length for p in products})
    text = " / ".join(f"{n} bp" for n in lengths)
    worst = max(products, key=lambda p: p.mismatches_forward + p.mismatches_reverse)
    if worst.mismatches_forward or worst.mismatches_reverse:
        text += f" ({worst.mismatches_forward} + {worst.mismatches_reverse} mismatches)"
    return text


def md_escape(text: str) -> str:
    return text.replace("|", "\\|")


def write_markdown(path: Path, primers: dict[str, Primer], amplicons: dict[str, Amplicon],
                   h37rv: dict[str, list[Product]] | None, others: list[tuple[str, str, int, str]],
                   scope_files: list[Path], unsequenced: list[tuple[str, list[str]]]) -> None:
    out: list[str] = []
    today = datetime.date.today().isoformat()
    protocols = sorted({p.source_file for p in primers.values()})
    out += ["# Primers of the mycobacterial Sanger assays", "",
            f"Written by `scripts/02_handbook_primers.py` on {today} from `primers.tsv` and `amplicons.tsv`. Do not edit by hand.",
            "", f"Source: {len(protocols)} protocols in `data/Arbeitshandbücher` "
            f"({sum('ausserhalb AHB' not in f for f in protocols)} chapters of the molecular handbook, "
            f"{sum('ausserhalb AHB' in f for f in protocols)} research protocols kept outside it, protocol `NZM <n>`). "
            "Every sequence below was found in the text of its protocol by the script.", "",
            f"**{len(primers)} primers, {len(amplicons)} primer pairs.** Only assays that end in Sanger sequencing are listed; "
            "probes and real-time or fragment-length assays are in the last section.", ""]

    out += ["## Primers", "", "| Primer | Target | Scope | Direction | Sequence 5'→3' | Role | Protocol |", "|---|---|---|---|---|---|---|"]
    for p in primers.values():
        out.append(f"| {md_escape(p.name)} | {p.target} | {p.scope} | {p.direction} | `{p.sequence}` | {p.role} | {p.protocol} |")
    out.append("")

    out += ["## Amplicons", "",
            "\"Protocol\" is the product length written in the protocol; \"H37Rv\" is the in-silico product on "
            f"*M. tuberculosis* H37Rv ({H37RV_ACCESSION}), primers included, at most {MAX_MISMATCHES} mismatches per primer. "
            "Pairs for other species give no product there.", "",
            "| Amplicon | Target | Forward + reverse | Protocol | H37Rv | Use | Sequenced with | Protocol no. |", "|---|---|---|---|---|---|---|---|"]
    differing: list[str] = []
    for a in amplicons.values():
        documented = f"{a.documented_bp} bp" if a.documented_bp else "–"
        observed = "not run" if h37rv is None else describe(h37rv[a.name])
        if h37rv is not None and a.documented_bp and h37rv[a.name] and all(p.length != a.documented_bp for p in h37rv[a.name]):
            differing.append(f"{a.name} ({a.documented_bp} bp in the protocol, {describe(h37rv[a.name])} on H37Rv)")
        out.append(f"| {a.name} | {a.target} | {md_escape(a.forward)} + {md_escape(a.reverse)} | {documented} | {observed} | {a.use} | {a.sequenced_with} | {a.protocol} |")
    out.append("")
    if differing:
        out += ["Product length differs from the protocol for: " + "; ".join(differing) + ". "
                "Some protocols give the length for another species (16S-23S ITS: *M. fortuitum*; 23S MCLR: *M. abscessus*).", ""]

    out += ["## Which protocol says what", "", "| Protocol | File | Primers |", "|---|---|---|"]
    by_file: dict[str, list[Primer]] = {}
    for p in primers.values():
        by_file.setdefault(p.source_file, []).append(p)
    for source, group in by_file.items():
        out.append(f"| {group[0].protocol} | {md_escape(Path(source).name)} | {', '.join(md_escape(p.name) for p in group)} |")
    out.append("")

    out += ["## Other oligos in the mycobacterial protocols", "",
            f"{len(scope_files)} protocols were searched (routine chapters with a mycobacterial title and the research folder). "
            "Oligos written there with 5'/3' ends that are not in `primers.tsv`: probes, real-time PCR and "
            "fragment-length assays, variants of a listed primer.", ""]
    if others:
        out += ["| Protocol file | Line | Sequence |", "|---|---|---|"]
        for source, _, line, sequence in others:
            out.append(f"| {md_escape(Path(source).name)} | {line} | `{sequence}` |")
    else:
        out.append("None.")
    out.append("")

    out += ["## Research protocols without a primer sequence", "",
            "Protocols of the research folder in which no primer sequence was found (none written, or written without 5'/3' ends), "
            "so no reference is cut for them. "
            "\"Primers named\" are the names in the reagent table of the protocol.", ""]
    if unsequenced:
        out += ["| Protocol file | Primers named |", "|---|---|"]
        for source, names in unsequenced:
            out.append(f"| {md_escape(Path(source).name)} | {', '.join(md_escape(n) for n in names) or '–'} |")
    else:
        out.append("None.")
    out.append("")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(out), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--handbook", type=Path, default=HANDBOOK_DIR, help="handbook directory")
    parser.add_argument("--primers", type=Path, default=PRIMERS_TSV)
    parser.add_argument("--amplicons", type=Path, default=AMPLICONS_TSV)
    parser.add_argument("--h37rv", type=Path, default=NCBI_DATASET / H37RV_ACCESSION, help="H37Rv genome directory")
    parser.add_argument("--no-h37rv", action="store_true", help="skip the in-silico PCR on H37Rv")
    parser.add_argument("--out", type=Path, default=OUTPUT_PRIMERS)
    parser.add_argument("--md", type=Path, default=PRIMERS_MD)
    args = parser.parse_args()

    primers = load_primers(args.primers)
    amplicons = load_amplicons(args.amplicons, primers)
    args.out.mkdir(parents=True, exist_ok=True)

    # 1. every primer in its protocol
    texts: dict[str, str] = {}
    missing: list[Primer] = []
    with (args.out / "primer_check.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["primer", "sequence", "found_as", "source_file"])
        for primer in primers.values():
            source = args.handbook / primer.source_file
            if primer.source_file not in texts:
                texts[primer.source_file] = squeeze(handbook_text(source)) if source.is_file() else ""
            text = texts[primer.source_file]
            found_as = "written" if primer.sequence in text else "reverse complement" if revcomp(primer.sequence) in text else ""
            if not found_as and set(primer.sequence) - set("ACGT"):
                # degenerate bases may be spelled out in the protocol ([C,T] for Y)
                pattern = "".join(b if b in "ACGT" else f"(?:{b}|\\[{','.join(IUPAC[b])}\\]|\\[{','.join(reversed(IUPAC[b]))}\\])"
                                  for b in primer.sequence)
                found_as = "written" if re.search(pattern, text) else ""
            if not found_as:
                missing.append(primer)
            writer.writerow([primer.name, primer.sequence, found_as or "NOT FOUND", primer.source_file])

    # 2. oligos of the mycobacterial protocols that are not curated
    scope_files = [f for f in handbook_files(args.handbook) if in_scope(f)]
    others: list[tuple[str, str, int, str]] = []
    unsequenced: list[tuple[str, list[str]]] = []
    curated_names = {squeeze(p.name) for p in primers.values()}
    with (args.out / "handbook_oligos.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["source_file", "line", "sequence", "primer"])
        for source in scope_files:
            relative = str(source.relative_to(args.handbook))
            seen: set[str] = set()
            text = handbook_text(source)
            oligos = find_oligos(text)
            if RESEARCH_DIR in source.parent.name and not oligos and SEQUENCING.search(text):
                names = sorted({m.group(1).strip() for m in PRIMER_NAME.finditer(text)})
                # a protocol that only names primers whose sequence is curated from another file is covered
                if not names or any(squeeze(n) not in curated_names for n in names):
                    unsequenced.append((relative, names))
            for line, sequence in oligos:
                name = status_of(sequence, primers)
                writer.writerow([relative, line, sequence, name or "not in primers.tsv"])
                if not name and sequence not in seen:
                    seen.add(sequence)
                    others.append((relative, name, line, sequence))

    # 3. product of every pair on H37Rv
    h37rv: dict[str, list[Product]] | None = None
    if not args.no_h37rv:
        fasta = genome_fasta(args.h37rv)
        if fasta is None:
            print(f"WARNING: no genome in {args.h37rv} (run scripts/01 first); H37Rv column left out", file=sys.stderr)
        else:
            h37rv = genome_products(fasta, primers, amplicons, MAX_MISMATCHES, MAX_PRODUCT)

    write_markdown(args.md, primers, amplicons, h37rv, others, scope_files, unsequenced)
    subprocess.run(["bash", str(MANIFEST_SCRIPT), str(args.out), str(REPO), "scripts/02_handbook_primers.py",
                    str(args.handbook), str(args.primers), str(args.amplicons)], check=True)

    print(f"{len(primers)} primers, {len(amplicons)} amplicons, {len(texts)} protocol files")
    print(f"{len(scope_files)} mycobacterial protocols searched, {len(others)} oligo(s) not in primers.tsv, "
          f"{len(unsequenced)} research protocol(s) without a primer sequence")
    print(f"wrote {args.md}")
    if missing:
        print(f"ERROR: {len(missing)} primer sequence(s) not found in their protocol:", file=sys.stderr)
        for primer in missing:
            print(f"  {primer.name}\t{primer.sequence}\t{primer.source_file}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
