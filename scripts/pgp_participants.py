"""Fetch PGP-HMS genotype and variant data grouped by participant, for cross-source test setups.

One person with several independent facets (a 23andMe export, a PGP-sequenced VCF, a Complete Genomics
`var` file) is the strongest test just-vcf can get. Importing each facet correctly means they converge
on the same genotypes where they overlap. This script gathers, per participant:

- **Arvados datasets.** PGP's own sequencing lives in per-participant projects on the public Arvados
  cluster `ac2it`. We take the VCF collections (freeBayes + VarScan2, hg19) and the Complete Genomics
  deliveries, keeping only `var`, `masterVarBeta`, `vcfBeta` and `summary` from those. BAMs, reports and
  bulky CGI evidence/coverage files are left.
- **Participant uploads** of genotype/variant data types from the public listing (see `pgp_leech.py`).

Everything lands under `<out>/participants/<hu>/<source>-<id>/<file>`. A file already on disk from the
older by-data-type layout is moved into place, never downloaded twice. Within a participant, two items
with the same file name and byte size are one file, fetched once and recorded as a duplicate of it. An
upload whose listed size equals an Arvados file of the same person is also treated as a duplicate.

The Arvados API is read with the anonymous token the cluster itself publishes in its public config.

    uv run python scripts/pgp_participants.py --dry-run
    uv run python scripts/pgp_participants.py --only huF1DC30 --only hu3F7E26
"""

import json
import logging
import os
import re
import shutil
from collections import Counter
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path
from urllib.parse import quote

import httpx
import typer
from dotenv import load_dotenv

import pgp_leech as pl

log = logging.getLogger("pgp_participants")

ARVADOS_API = "https://ac2it.arvadosapi.com/arvados/v1"
KEEP_WEB = "https://{host}.collections.ac2it.arvadosapi.com/_/{path}"
PARTICIPANT = re.compile(r"\bhu[0-9A-F]{6}\b")
# From a Complete Genomics delivery, the small files that carry calls: the native variant list, the
# per-locus master list, CGI's own VCF export, and the run summary.
CGI_KEEP = re.compile(r"^(?:(?:var|masterVarBeta)-GS\d+-ASM\.tsv|vcfBeta-GS\d+-ASM\.vcf)\.bz2$|^summary-GS\d+-ASM\.tsv$")
DEFAULT_MAX_BYTES = 2_500_000_000


@dataclass(frozen=True)
class Item:
    """One file we want for one participant."""

    participant: str
    source: str  # "arvados" or "upload"
    item_id: str  # Arvados collection uuid, or listing download id
    item_name: str
    kind: str  # "vcf", "cgi", or the upload's self-reported data type
    file_name: str | None  # unknown for an upload until it is fetched
    size: int
    url: str | None  # keep-web URL for Arvados files; uploads resolve through pgp_leech.fetch_row

    @property
    def rel_dir(self) -> Path:
        return Path(self.participant) / f"{self.source}-{pl._safe(self.item_id)}"


@dataclass(frozen=True)
class ItemResult:
    item: Item
    outcome: str
    path: str | None = None
    bytes_on_disk: int | None = None
    reason: str | None = None
    duplicate_of: str | None = None
    fetched_at: str | None = None

    def __bool__(self) -> bool:
        return self.outcome in ("downloaded", "already_present", "moved")


def parse_manifest(text: str) -> list[tuple[str, int]]:
    """(path, size) of every file in an Arvados manifest. Stream names prefix paths; `\\040` is a space."""
    files: list[tuple[str, int]] = []
    for line in text.splitlines():
        parts = line.split()
        if not parts:
            continue
        stream = parts[0].removeprefix(".").lstrip("/")
        for tok in parts[1:]:
            m = re.fullmatch(r"(\d+):(\d+):(.+)", tok)
            if m is None:
                continue
            name = m.group(3).replace("\\040", " ")
            files.append((f"{stream}/{name}" if stream else name, int(m.group(2))))
    return files


def classify_collection(name: str) -> str | None:
    """What a dataset collection holds, from its name: "vcf", "cgi", or None to skip it."""
    if name.endswith(".vcf") or name.endswith(" vcf"):
        return "vcf"
    if "-DID" in name or "Raw sequencing data" in name or re.search(r"GS\d+-ASM", name):
        return "cgi"
    return None


class Arvados:
    def __init__(self, fetcher: pl.Fetcher) -> None:
        self.fetcher = fetcher
        config = fetcher.get(f"{ARVADOS_API}/config")
        config.raise_for_status()
        self.headers = {"Authorization": f"Bearer {config.json()['Users']['AnonymousUserToken']}"}

    def list(self, resource: str, select: list[str], filters: list | None = None, page: int = 200) -> list[dict]:
        items: list[dict] = []
        while True:
            params = {
                "select": json.dumps(select),
                "limit": str(page),
                "offset": str(len(items)),
                "order": json.dumps(["uuid asc"]),
                "count": "none",
            }
            if filters:
                params["filters"] = json.dumps(filters)
            # POST with `_method=GET`: an owner filter over ~350 projects is too long for a query string (414).
            self.fetcher._wait(ARVADOS_API)
            resp = self.fetcher.client.post(f"{ARVADOS_API}/{resource}", data={**params, "_method": "GET"}, headers=self.headers)
            resp.raise_for_status()
            batch = resp.json()["items"]
            items.extend(batch)
            if len(batch) < page:
                return items


def arvados_items(arv: Arvados) -> tuple[list[Item], Counter[str]]:
    """Every wanted file in the participants' dataset projects, and a count of what was skipped, by reason."""
    groups = arv.list("groups", ["uuid", "name"])
    projects = {g["uuid"]: PARTICIPANT.search(g["name"]).group(0)
                for g in groups if g["name"].endswith(("datasets", "dataset")) and PARTICIPANT.search(g["name"])}
    colls = arv.list("collections", ["uuid", "name", "owner_uuid", "portable_data_hash", "manifest_text"],
                     [["owner_uuid", "in", sorted(projects)]], page=50)
    items: list[Item] = []
    skipped: Counter[str] = Counter()
    for c in colls:
        kind = classify_collection(c["name"] or "")
        if kind is None:
            skipped[f"collection not vcf/cgi ({(c['name'] or '').rsplit('.', 1)[-1][:12]})"] += 1
            continue
        host = c["portable_data_hash"].replace("+", "-")
        for path, size in parse_manifest(c["manifest_text"] or ""):
            base = path.rsplit("/", 1)[-1]
            if kind == "cgi" and not CGI_KEEP.match(base):
                skipped["cgi file not in keep list"] += 1
                continue
            items.append(Item(projects[c["owner_uuid"]], "arvados", c["uuid"], c["name"], kind, base, size,
                              KEEP_WEB.format(host=host, path=quote(path))))
    return items, skipped


# Reads and alignments uploaded under a variant-sounding data type: Veritas uploads, for one, are a
# "<id> - VCF" row plus one "<id>.chrN.bam - BAM" row per chromosome, all typed "Veritas Genetics".
READS_BY_NAME = re.compile(r"\.bam\b|- BAM\b|\.bai\b|\.cram\b|\.crai\b|fastq|\.fq\b", re.IGNORECASE)


def upload_items(rows: list[pl.ListingRow]) -> list[Item]:
    return [Item(r.participant, "upload", str(r.download_id), r.name, r.data_type, None, r.listed_size, None)
            for r in rows if not pl.NON_VARIANT_TYPES.search(r.data_type) and not READS_BY_NAME.search(r.name)]


def plan(items: list[Item]) -> tuple[list[Item], dict[Item, Item]]:
    """Deduplicate within each participant. Returns the items to fetch and {duplicate: kept}.

    Order decides which copy is kept: an original collection before one named "copy of …", then by
    collection id. Arvados files are keyed by (name, size); an upload, whose file name is unknown until
    fetched, matches an Arvados file of the same person by size alone.
    """
    def order(i: Item) -> tuple:
        return (i.participant, i.source != "arvados", (i.item_name or "").startswith("copy of"), i.item_id, i.file_name or "")

    kept: list[Item] = []
    dup: dict[Item, Item] = {}
    by_name: dict[tuple[str, str, int], Item] = {}
    by_size: dict[tuple[str, int], Item] = {}
    for i in sorted(items, key=order):
        if i.source == "arvados":
            key = (i.participant, i.file_name or "", i.size)
            if key in by_name:
                dup[i] = by_name[key]
                continue
            by_name[key] = i
            by_size.setdefault((i.participant, i.size), i)
        elif (i.participant, i.size) in by_size:
            dup[i] = by_size[(i.participant, i.size)]
            continue
        kept.append(i)
    return kept, dup


def select_participants(items: list[Item], only: list[str], min_facets: int) -> list[str]:
    """Participants in priority order: with an Arvados VCF, then with CGI, then by number of items."""
    per: dict[str, Counter[str]] = {}
    for i in items:
        per.setdefault(i.participant, Counter())[i.kind if i.source == "arvados" else "upload"] += 1
    if only:
        return [p for p in only if p in per]
    chosen = [p for p, c in per.items() if c["vcf"] or c["cgi"] or c["upload"] >= min_facets]
    return sorted(chosen, key=lambda p: (not per[p]["vcf"], not per[p]["cgi"], -sum(per[p].values()), p))


def migrate_by_type_layout(out: Path, participants_root: Path) -> int:
    """Move files from `<out>/<data type>/<hu>/<download id>/` into `participants/<hu>/upload-<id>/`."""
    moved = 0
    for f in sorted(out.glob("*/hu*/*/*")):
        if not f.is_file() or f.name.endswith(".part") or participants_root in f.parents:
            continue
        hu, download_id = f.parent.parent.name, f.parent.name
        if not (PARTICIPANT.fullmatch(hu) and download_id.isdigit()):
            continue
        dest = participants_root / hu / f"upload-{download_id}" / f.name
        dest.parent.mkdir(parents=True, exist_ok=True)
        if dest.exists():
            if dest.stat().st_size == f.stat().st_size:
                f.unlink()
                continue
            raise SystemExit(f"refusing to overwrite {dest} with a different {f}")
        shutil.move(f, dest)
        moved += 1
    for d in sorted((p for p in out.glob("*/hu*/*") if p.is_dir()), reverse=True):
        if participants_root not in d.parents and not any(d.iterdir()):
            d.rmdir()
    for d in sorted((p for p in out.glob("*/hu*") if p.is_dir()), reverse=True):
        if participants_root not in d.parents and not any(d.iterdir()):
            d.rmdir()
    for d in sorted(p for p in out.iterdir() if p.is_dir() and p != participants_root):
        if not any(d.iterdir()):
            d.rmdir()
    legacy = out / "manifest.jsonl"
    if legacy.exists():
        legacy.replace(out / "manifest.legacy-by-type.jsonl")
    return moved


def fetch_arvados(fetcher: pl.Fetcher, item: Item, root: Path) -> ItemResult:
    dest = root / item.rel_dir / pl._safe(item.file_name or "")
    now = datetime.now(UTC).isoformat()
    if dest.exists() and dest.stat().st_size == item.size:
        return ItemResult(item, "already_present", str(dest), item.size)
    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        with fetcher.open(item.url or "", follow=True) as resp:
            if resp.status_code in (401, 403):
                return ItemResult(item, "not_public", reason=f"HTTP {resp.status_code}")
            resp.raise_for_status()
            written = pl._write_body(resp, dest)
    except httpx.HTTPError as err:
        return ItemResult(item, "http_error", reason=f"{type(err).__name__}: {err}")
    part = dest.with_name(dest.name + ".part")
    if written != item.size:
        return ItemResult(item, "size_mismatch", str(part), written, reason=f"expected {item.size}, got {written}")
    part.replace(dest)
    return ItemResult(item, "downloaded", str(dest), written, fetched_at=now)


def fetch_upload(fetcher: pl.Fetcher, item: Item, row: pl.ListingRow, root: Path, max_bytes: int) -> list[ItemResult]:
    results = pl.fetch_row(fetcher, row, root, max_bytes=max_bytes, dry_run=False, target_dir=root / item.rel_dir)
    return [ItemResult(Item(**{**asdict(item), "file_name": r.file_name}), r.outcome.value, r.path, r.bytes_on_disk,
                       r.reason, fetched_at=r.fetched_at) for r in results]


def write_manifest(path: Path, results: list[ItemResult]) -> None:
    """Merge results into the manifest; newest record per (participant, source, id, file) wins; sorted."""
    def key(rec: dict) -> tuple:
        it = rec["item"]
        return (it["participant"], it["source"], it["item_id"], it["file_name"] or "")

    merged: dict[tuple, dict] = {}
    if path.exists():
        for line in path.read_text().splitlines():
            if line.strip():
                rec = json.loads(line)
                merged[key(rec)] = rec
    for r in results:
        rec = asdict(r)
        merged[key(rec)] = rec
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text("".join(json.dumps(merged[k], sort_keys=True) + "\n" for k in sorted(merged)))
    tmp.replace(path)


app = typer.Typer(add_completion=False, help=__doc__.split("\n\n")[0])


@app.command()
def main(
    out: Path | None = typer.Option(None, help="Root. Default: $JUST_VCF_SAMPLES_DIR/pgp-hms, else data/input/pgp-hms."),
    only: list[str] = typer.Option([], "--only", help="Only these participants (hu ids). Repeatable."),
    limit: int | None = typer.Option(None, help="At most this many participants, in priority order."),
    min_facets: int = typer.Option(2, help="Also take participants without PGP sequencing who uploaded this many variant files."),
    max_bytes: int = typer.Option(DEFAULT_MAX_BYTES, help="Skip single files larger than this."),
    no_uploads: bool = typer.Option(False, "--no-uploads", help="Arvados datasets only."),
    delay: float = typer.Option(pl.DEFAULT_DELAY_S, help="Seconds between requests to one host."),
    dry_run: bool = typer.Option(False, "--dry-run", help="Plan and report; download and move nothing."),
) -> None:
    load_dotenv()
    if out is None:
        samples = os.environ.get("JUST_VCF_SAMPLES_DIR", "")
        out = Path(samples) / "pgp-hms" if samples else pl.FALLBACK_OUT
    root = out / "participants"
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    logging.getLogger("httpx").setLevel(logging.WARNING)
    fetcher = pl.Fetcher(delay)

    arv_items, skipped = arvados_items(Arvados(fetcher))
    for reason, n in sorted(skipped.items()):
        log.info("arvados: skipped %d: %s", n, reason)
    rows: dict[str, pl.ListingRow] = {}
    if not no_uploads:
        resp = fetcher.get(f"{pl.PGP_BASE}{pl.LISTING_PATH}")
        resp.raise_for_status()
        listed, _ = pl.parse_listing(resp.text)
        rows = {str(r.download_id): r for r in listed}
    items = arv_items + upload_items(list(rows.values()))

    people = select_participants(items, only, min_facets)
    if limit is not None:
        people = people[:limit]
    chosen = set(people)
    kept, dup = plan([i for i in items if i.participant in chosen])
    too_big = [i for i in kept if i.size > max_bytes]
    todo = [i for i in kept if i.size <= max_bytes]
    by_kind = Counter(i.kind if i.source == "arvados" else "upload" for i in todo)
    log.info("plan: %d participant(s); %d file(s), %.1f GB (%s); %d duplicate(s) folded; %d over the size cap",
             len(people), len(todo), sum(i.size for i in todo) / 1e9,
             ", ".join(f"{k} {v}" for k, v in sorted(by_kind.items())), len(dup), len(too_big))
    if dry_run:
        for p in people[:15]:
            mine = [i for i in todo if i.participant == p]
            log.info("  %s: %s", p, ", ".join(sorted(Counter(i.kind for i in mine).elements())))
        return

    root.mkdir(parents=True, exist_ok=True)
    log.info("migrated %d file(s) from the by-data-type layout", migrate_by_type_layout(out, root))
    manifest = root / "manifest.jsonl"
    results: list[ItemResult] = [ItemResult(d, "duplicate", reason="same person, same file", duplicate_of=str(root / k.rel_dir))
                                 for d, k in dup.items()]
    results += [ItemResult(i, "skipped_too_large", reason=f"{i.size} bytes > cap {max_bytes}") for i in too_big]
    order = {p: n for n, p in enumerate(people)}
    for n, item in enumerate(sorted(todo, key=lambda i: (order[i.participant], i.source, i.item_id, i.file_name or "")), 1):
        got = [fetch_arvados(fetcher, item, root)] if item.source == "arvados" else fetch_upload(fetcher, item, rows[item.item_id], root, max_bytes)
        for r in got:
            if not r:
                log.warning("%s %s-%s: %s %s", item.participant, item.source, item.item_id, r.outcome, r.reason or "")
        results += got
        if n % 10 == 0:
            write_manifest(manifest, results)
    write_manifest(manifest, results)
    typer.echo(" · ".join(f"{k}: {v}" for k, v in sorted(Counter(r.outcome for r in results).items())))
    typer.echo(f"manifest: {manifest}")


if __name__ == "__main__":
    app()
