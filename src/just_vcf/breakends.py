"""Breakend ALTs (VCF 4.2-4.5 §5.4): how one is read, and how it turns when an end changes strand.

A breakend ALT is one of four forms, `t[p[`, `t]p]`, `]p]t` and `[p[t` (`t` the bases replacing REF,
`p` the mate's `contig:position`), or a single breakend, `t.` or `.t`. Each form states two ends of a
novel adjacency, each a position and a side: the side of the position whose sequence takes part in the
join. VCF 4.3's four sentences give them: `t[p[` "piece extending to the right of p is joined after t",
so the record keeps the sequence left of POS and the mate the sequence right of p; `t]p]` "reverse comp
piece extending left of p is joined after t" (left, left); `]p]t` "piece extending to the left of p is
joined before t" (right, left); `[p[t` "reverse comp piece extending right of p is joined before t"
(right, right). A single breakend has the record's side only: `t.` left, `.t` right.

`concordance` reads the form and the mate to meet breakends (truvari's `bnd_match`); `liftover` turns
the form when an end lands on a minus-strand chain block (`turned`).
"""

import re
from dataclasses import dataclass

import polars as pl

from just_vcf.schema import CONTIG_ALIASES

BREAKEND_FORMS: tuple[str, ...] = ("t[p[", "t]p]", "]p]t", "[p[t")
"""VCF 4.2-4.5's four breakend ALT forms (`t` the bases at POS, `p` the mate's `contig:position`). They
are truvari's four `(direction, strand)` pairs (`bnd_direction_strand`): `t[p[` left/direct, `[p[t`
left/complement, `t]p]` right/direct, `]p]t` right/complement."""

SIDES: dict[str, tuple[str, str]] = {"t[p[": ("L", "R"), "t]p]": ("L", "L"), "]p]t": ("R", "L"), "[p[t": ("R", "R")}
"""Each form's (record's side, mate's side): where the sequence that takes part in the join lies,
left (`L`) or right (`R`) of POS and of p."""

FORM_OF: dict[tuple[str, str], str] = {sides: form for form, sides in SIDES.items()}

SINGLE_BREAKEND = r"^\.[A-Za-z]+$|^[A-Za-z]+\.$"
"""A single breakend (`.A`, `A.`): a join to sequence the file does not place."""

_BASES_FIRST = r"^([A-Za-z]+)([\[\]])([^\[\]]+):(\d+)([\[\]])$"
_BASES_LAST = r"^([\[\]])([^\[\]]+):(\d+)([\[\]])([A-Za-z]+)$"
"""A breakend ALT with its bases before (`t[p[`, `t]p]`) or after (`]p]t`, `[p[t`) the mate. The
contig runs to the last `:` before the digits, so a name holding `:` (`HLA-A*01:01:01:01`) reads;
truvari's `[^\\[\\]:]+` would refuse it."""

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def breakend_parts(alt: pl.Expr) -> list[pl.Expr]:
    """`_form`, `_mate_chrom` (canonical) and `_mate_pos` of a breakend ALT, all null where the ALT is
    none of the four forms: brackets that do not pair, or no `contig:position` between them."""
    first, last = alt.str.extract_groups(_BASES_FIRST), alt.str.extract_groups(_BASES_LAST)
    f_open, f_close = first.struct.field("2"), first.struct.field("5")
    l_open, l_close = last.struct.field("1"), last.struct.field("4")
    position = pl.coalesce(first.struct.field("4"), last.struct.field("3")).cast(pl.Int64, strict=False)
    form = (
        pl.when(position.is_null())
        .then(None)
        .when(f_open.is_not_null() & (f_open == f_close))
        .then(pl.when(f_open == "[").then(pl.lit("t[p[")).otherwise(pl.lit("t]p]")))
        .when(l_open.is_not_null() & (l_open == l_close))
        .then(pl.when(l_open == "[").then(pl.lit("[p[t")).otherwise(pl.lit("]p]t")))
    )
    contig = mate_contig(alt).replace(CONTIG_ALIASES)
    return [form.alias("_form"), pl.when(form.is_not_null()).then(contig).alias("_mate_chrom"), pl.when(form.is_not_null()).then(position).alias("_mate_pos")]


def mate_contig(alt: pl.Expr) -> pl.Expr:
    """The mate's contig as a breakend ALT spells it, null on any other allele."""
    return pl.coalesce(alt.str.extract(_BASES_FIRST, 3), alt.str.extract(_BASES_LAST, 2))


@dataclass(frozen=True)
class Breakend:
    """One breakend ALT, read."""

    form: str
    """One of `BREAKEND_FORMS`, or `t.` / `.t` for a single breakend."""
    bases: str
    """`t`: REF's base, with any inserted bases (after it in `t[p[`/`t]p]`/`t.`, before it otherwise)."""
    mate_contig: str | None
    """The mate's contig as written (null on a single breakend)."""
    mate_pos: int | None

    @property
    def side(self) -> str:
        """Which side of POS the record keeps: `L` or `R`."""
        return "L" if self.form in ("t[p[", "t]p]", "t.") else "R"


def parse_breakend(alt: str) -> Breakend | None:
    """A breakend ALT read, or None where it is none of the four forms or a single breakend."""
    if re.fullmatch(SINGLE_BREAKEND, alt):
        return Breakend("t." if alt.endswith(".") else ".t", alt.strip("."), None, None)
    first = re.fullmatch(_BASES_FIRST, alt)
    if first and first.group(2) == first.group(5):
        return Breakend("t[p[" if first.group(2) == "[" else "t]p]", first.group(1), first.group(3), int(first.group(4)))
    last = re.fullmatch(_BASES_LAST, alt)
    if last and last.group(1) == last.group(4):
        return Breakend("[p[t" if last.group(1) == "[" else "]p]t", last.group(5), last.group(2), int(last.group(3)))
    return None


def turned(breakend: Breakend, own_minus: bool, mate_minus: bool, mate_contig: str | None, mate_pos: int | None) -> str:
    """The ALT a breakend becomes when its own position lands on a minus-strand block (`own_minus`)
    and its mate on one (`mate_minus`), with the mate at `mate_contig:mate_pos`.

    A minus-strand block reverses the sequence around a position, so the side an end keeps changes:
    the record's side flips with `own_minus`, the mate's with `mate_minus`, and the form is the one
    stating the new pair of sides (`FORM_OF`). `t` is read along the record's side away from the join;
    when that side flips, it is read the other way, so it is reverse-complemented: REF's base becomes
    its complement (the target's base at the new POS) and lands on the other side of the inserted
    bases, which are reverse-complemented too. A mate's flip alone leaves `t` as it is.
    """
    bases = breakend.bases.upper()
    if own_minus:
        bases = bases[::-1].translate(_COMPLEMENT)
    own = breakend.side if not own_minus else ("R" if breakend.side == "L" else "L")
    if breakend.mate_contig is None:
        return f"{bases}." if own == "L" else f".{bases}"
    mate = SIDES[breakend.form][1]
    if mate_minus:
        mate = "R" if mate == "L" else "L"
    form = FORM_OF[(own, mate)]
    p = f"{mate_contig}:{mate_pos}"
    return {"t[p[": f"{bases}[{p}[", "t]p]": f"{bases}]{p}]", "]p]t": f"]{p}]{bases}", "[p[t": f"[{p}[{bases}"}[form]
