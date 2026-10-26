"""Readers for the EIGENSTRAT files AADR (the Allen Ancient DNA Resource) publishes.

Ported from just-prs `ancestry/eigenstrat.py` (docs/PROVENANCE.md). One release is four files:

- `.ind`: one line per individual, `genetic_id sex group`, in the order of the `.geno` records.
- `.snp`: one line per site, `snp_id chrom genetic_position pos allele1 allele2`. The genotypes count
  `allele1`. AADR's README says it is the hg19 reference allele; the panel builder checks that against
  a GRCh37 reference rather than trusting it (`just_vcf.ancestry.aadr`).
- `.anno`: tab-separated, with a header row, one row per individual. Columns are found by their header
  name, never by position, so a release that inserts a column is read right or refused.
- `.geno`: the genotypes, in one of three layouts, told apart by the file's first bytes:
  - `tgeno` (convertf's `transpose_packed`, what AADR v66 ships): a 48-byte ASCII header
    `TGENO <n_ind> <n_snp> <hash> <hash>`, then one record of `ceil(n_snp / 4)` bytes per individual.
  - `geno` (convertf's `packedancestrymap`): a header `GENO <n_ind> <n_snp> <hash> <hash>` padded to
    one record, then one record per site; a record is `max(48, ceil(n_ind / 4))` bytes.
  - `text` (plain EIGENSTRAT): one line per site, one digit per individual.

  A packed byte holds four genotypes, two bits each, the first in the high bits; a genotype is the
  count of `allele1` (0, 1, 2), 3 is missing. Text writes 9 for missing. The `tgeno` layout was
  checked on AADR v66.p1's files (their sizes are exactly the header plus the records); the `geno`
  layout follows convertf's description and was tested here on files written to it, not on a
  published one. The two hashes in a packed header are not checked.
"""

import math
from pathlib import Path

import numpy as np
import polars as pl

from just_vcf.ancestry._projection import DEFAULT_MISSING

MISSING = int(DEFAULT_MISSING)
"""The value a missing genotype gets in what `read_genotypes` returns."""

GENO_LAYOUTS: frozenset[str] = frozenset({"tgeno", "geno", "text"})

TGENO_HEADER_BYTES = 48

ANNO_COLUMNS: dict[str, str] = {
    "genetic_id": "Genetic ID",
    "individual_id": "Individual ID",
    "date_bp": "Date mean in BP",
    "group_id": "Group ID",
    "political_entity": "Political Entity",
    "latitude": "Latitude",
    "longitude": "Longitude",
    "data_type": "Data type",
    "assessment": "ASSESSMENT",
}
"""What `read_anno` returns, and the header each comes from: the header itself, or its first words
followed by a space (AADR's headers carry long descriptions). An exact match wins, so `ASSESSMENT` is
not confused with the `ASSESSMENT WARNINGS: …` column after it."""


def read_ind(path: Path) -> pl.DataFrame:
    """A `.ind` file as `genetic_id, sex, group`, in the `.geno` record order."""
    rows = [line.split() for line in path.read_text().splitlines() if line.strip()]
    bad = [i for i, row in enumerate(rows) if len(row) != 3]
    if bad:
        raise ValueError(f"{path.name}: line {bad[0] + 1} has {len(rows[bad[0]])} fields, not 3")
    return pl.DataFrame(rows, schema={"genetic_id": pl.String, "sex": pl.String, "group": pl.String}, orient="row")


def read_snp(path: Path) -> pl.DataFrame:
    """A `.snp` file as `snp_id, chrom, pos, allele1, allele2`, in file order. `allele1` is counted."""
    frame = pl.read_csv(
        path,
        has_header=False,
        separator="\x1f",
        schema={"line": pl.String},
        quote_char=None,
    )
    parts = frame.select(pl.col("line").str.strip_chars().str.replace_all(r"\s+", " ").str.split(" ").alias("f"))
    widths = parts.select(pl.col("f").list.len()).to_series()
    if parts.height and (widths != 6).any():
        first = int((widths != 6).arg_max())
        raise ValueError(f"{path.name}: line {first + 1} has {widths[first]} fields, not 6")
    return parts.select(
        pl.col("f").list.get(0).alias("snp_id"),
        pl.col("f").list.get(1).alias("chrom"),
        pl.col("f").list.get(3).cast(pl.Int64).alias("pos"),
        pl.col("f").list.get(4).alias("allele1"),
        pl.col("f").list.get(5).alias("allele2"),
    )


def read_anno(path: Path) -> pl.DataFrame:
    """An `.anno` file's columns `ANNO_COLUMNS` names. `date_bp`, `latitude` and `longitude` are
    floats, null where the file writes `..` or anything else that is not a number."""
    frame = pl.read_csv(path, separator="\t", has_header=True, infer_schema=False, quote_char=None, encoding="utf8-lossy")
    headers = [name.strip().strip('"') for name in frame.columns]
    chosen = {}
    for key, wanted in ANNO_COLUMNS.items():
        exact = [i for i, h in enumerate(headers) if h == wanted]
        found = exact or [i for i, h in enumerate(headers) if h.startswith(wanted + " ")]
        if len(found) != 1:
            raise ValueError(f"{path.name}: {len(found)} columns headed {wanted!r}; need exactly one")
        chosen[key] = frame.columns[found[0]]
    numeric = {"date_bp", "latitude", "longitude"}
    return frame.select(
        pl.col(column).cast(pl.Float64, strict=False).alias(key) if key in numeric else pl.col(column).alias(key)
        for key, column in chosen.items()
    )


def geno_layout(path: Path) -> str:
    """Which of `GENO_LAYOUTS` a `.geno` file is, from its first bytes. Anything else is refused."""
    with path.open("rb") as handle:
        head = handle.read(64)
    if head.startswith(b"TGENO"):
        return "tgeno"
    if head.startswith(b"GENO"):
        return "geno"
    first_line = head.split(b"\n", 1)[0]
    if first_line and all(c in b"0129" for c in first_line.rstrip(b"\r")):
        return "text"
    raise ValueError(f"{path.name}: not a .geno layout this reads (TGENO, GENO or text); starts {head[:8]!r}")


def read_genotypes(path: Path, n_ind: int, n_snp: int, individuals: np.ndarray) -> np.ndarray:
    """The genotypes of `individuals` (row numbers of the `.ind` file) at every site.

    Returns an `(n_snp, len(individuals))` int8 matrix of counts of `allele1`, `MISSING` where the
    file has no genotype. `n_ind` and `n_snp` are the `.ind` and `.snp` line counts; a packed header
    or a file size that disagrees with them is refused.
    """
    individuals = np.asarray(individuals, dtype=np.int64)
    layout = geno_layout(path)
    if layout == "text":
        return _read_text(path, n_ind, n_snp, individuals)
    _check_header(path, layout, n_ind, n_snp)
    if layout == "tgeno":
        return _read_tgeno(path, n_ind, n_snp, individuals)
    return _read_packed(path, n_ind, n_snp, individuals)


def _check_header(path: Path, layout: str, n_ind: int, n_snp: int) -> None:
    with path.open("rb") as handle:
        fields = handle.read(TGENO_HEADER_BYTES).split(b"\0", 1)[0].split()
    stated_ind, stated_snp = int(fields[1]), int(fields[2])
    if stated_ind != n_ind:
        raise ValueError(f"{path.name}: its header says {stated_ind} individuals, the .ind has {n_ind}")
    if stated_snp != n_snp:
        raise ValueError(f"{path.name}: its header says {stated_snp} sites, the .snp has {n_snp}")
    if layout == "tgeno":
        expected = TGENO_HEADER_BYTES + n_ind * math.ceil(n_snp / 4)
    else:
        expected = (n_snp + 1) * max(48, math.ceil(n_ind / 4))
    actual = path.stat().st_size
    if actual != expected:
        raise ValueError(f"{path.name}: {actual} bytes, but a {layout} file of {n_ind} individuals and {n_snp} sites has {expected}")


def _unpack(raw: np.ndarray) -> np.ndarray:
    """Packed bytes (any shape, last axis the bytes) to two-bit values, four per byte, high bits first."""
    out = np.stack([(raw >> 6) & 3, (raw >> 4) & 3, (raw >> 2) & 3, raw & 3], axis=-1)
    return out.reshape(*raw.shape[:-1], raw.shape[-1] * 4)


def _read_tgeno(path: Path, n_ind: int, n_snp: int, individuals: np.ndarray) -> np.ndarray:
    """Individual-major: only the chosen individuals' records are read."""
    rlen = math.ceil(n_snp / 4)
    out = np.empty((len(individuals), n_snp), dtype=np.int8)
    with path.open("rb") as handle:
        for row, i in enumerate(individuals):
            handle.seek(TGENO_HEADER_BYTES + int(i) * rlen)
            raw = np.frombuffer(handle.read(rlen), dtype=np.uint8)
            out[row] = _unpack(raw)[:n_snp]
    out[out == 3] = MISSING
    return np.ascontiguousarray(out.T)


def _read_packed(path: Path, n_ind: int, n_snp: int, individuals: np.ndarray) -> np.ndarray:
    """Site-major: every record is read, the chosen individuals' two bits taken from each."""
    rlen = max(48, math.ceil(n_ind / 4))
    records = np.fromfile(path, dtype=np.uint8)[rlen:].reshape(n_snp, rlen)
    shifts = (6 - 2 * (individuals % 4)).astype(np.uint8)
    out = ((records[:, individuals // 4] >> shifts) & 3).astype(np.int8)
    out[out == 3] = MISSING
    return out


def _read_text(path: Path, n_ind: int, n_snp: int, individuals: np.ndarray) -> np.ndarray:
    out = np.empty((n_snp, len(individuals)), dtype=np.int8)
    with path.open("rb") as handle:
        count = 0
        for count, line in enumerate(handle, start=1):
            digits = np.frombuffer(line.rstrip(b"\r\n"), dtype=np.uint8)
            if digits.size != n_ind:
                raise ValueError(f"{path.name}: line {count} has {digits.size} genotypes, the .ind has {n_ind} individuals")
            if count > n_snp:
                break
            out[count - 1] = digits[individuals] - ord("0")
    if count != n_snp:
        raise ValueError(f"{path.name}: {count} lines, the .snp has {n_snp} sites")
    bad = (out != 9) & ((out < 0) | (out > 2))
    if bad.any():
        raise ValueError(f"{path.name}: a genotype that is not 0, 1, 2 or 9")
    out[out == 9] = MISSING
    return out
