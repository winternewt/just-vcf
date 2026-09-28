"""Download consented-public genotype files from the Harvard Personal Genome Project (PGP-HMS).

The public listing at https://my.pgp-hms.org/public_genetic_data has one row per participant file. Each
row's download link redirects to an Arvados collection page, which lists the actual files. This script
walks that chain politely (robots.txt asks for a 5 s crawl delay), writes each file under
`<out>/<data type>/<participant>/<download id>/<file name>`, and keeps `<out>/manifest.jsonl` sorted by
download id so reruns resume instead of starting over.

Licence: the PGP Researcher FAQ puts all public genomic and health data under CC0
(https://github.com/PGPHarvard/FAQ/wiki/Researcher-FAQ). Data types are self-reported by participants.

It is a dev tool, not library code: httpx is a dev dependency and nothing in `src/` imports this.

    uv run python scripts/pgp_leech.py --dry-run --limit 5
    uv run python scripts/pgp_leech.py --data-type 23andMe --max-bytes 50000000
"""

import enum
import json
import logging
import os
import re
import time
from collections import Counter
from contextlib import AbstractContextManager
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import quote, unquote, urljoin

import httpx
import typer
from dotenv import load_dotenv

log = logging.getLogger("pgp_leech")

PGP_BASE = "https://my.pgp-hms.org"
LISTING_PATH = "/public_genetic_data"
USER_AGENT = "just-vcf-pgp-leech (+https://github.com/winternewt)"
# robots.txt on my.pgp-hms.org: "Crawl-delay: 5".
DEFAULT_DELAY_S = 5.0
# Downloads land in `$JUST_VCF_SAMPLES_DIR/pgp-hms`; without that variable, in the git-ignored data/ tree.
FALLBACK_OUT = Path(__file__).resolve().parent.parent / "data" / "input" / "pgp-hms"
# Array exports are single-digit MB. The listing also carries a handful of multi-GB files under
# consumer data types, so the default cap keeps a first run to the files this project is after.
DEFAULT_MAX_BYTES = 200_000_000
CHUNK = 1 << 20


@dataclass(frozen=True)
class ListingRow:
    """One row of the public listing, as the page states it."""

    download_id: int
    participant: str
    # Early participants are also listed by their original enrolment number, e.g. `huF06AD0, PGP113`.
    pgp_alias: str | None
    participant_numeric_id: int
    published: str
    data_type: str
    source: str
    name: str
    locator: str
    listed_size: int


@dataclass(frozen=True)
class CollectionFile:
    """One file inside the Arvados collection a download link redirects to."""

    name: str
    size: int
    url: str


class Outcome(enum.StrEnum):
    DOWNLOADED = "downloaded"
    ALREADY_PRESENT = "already_present"
    SKIPPED_TOO_LARGE = "skipped_too_large"
    DRY_RUN = "dry_run"
    SIZE_MISMATCH = "size_mismatch"
    HTTP_ERROR = "http_error"
    # 401/403 from the file host: the listing shows the row but its collection is not shared publicly.
    NOT_PUBLIC = "not_public"
    UNEXPECTED_RESPONSE = "unexpected_response"
    EMPTY_COLLECTION = "empty_collection"


@dataclass(frozen=True)
class FileResult:
    """What happened to one listing row (or one file of it). Truthy only when the bytes are on disk."""

    row: ListingRow
    outcome: Outcome
    file_name: str | None = None
    path: str | None = None
    bytes_on_disk: int | None = None
    reason: str | None = None
    fetched_at: str | None = None

    def __bool__(self) -> bool:
        return self.outcome in (Outcome.DOWNLOADED, Outcome.ALREADY_PRESENT)


class _ListingParser(HTMLParser):
    """Collects the `<tr data-file-row …>` rows of the listing table.

    Cells are identified by their `data-summarize-as` attribute rather than by column position, so a
    column added to the page does not shift every field by one.
    """

    def __init__(self) -> None:
        super().__init__()
        self.rows: list[dict[str, str]] = []
        self._row: dict[str, str] | None = None
        self._cell: str | None = None
        self._text: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        a = {k: v or "" for k, v in attrs}
        if tag == "tr" and "data-file-row" in a:
            self._row = {
                "participant_numeric_id": a.get("data-participant-id", ""),
                "locator": a.get("data-locator", ""),
                "listed_size": a.get("data-file-size", ""),
            }
        elif self._row is not None and tag == "td":
            self._cell = a.get("data-summarize-as", "")
            self._text = []
        elif self._row is not None and tag == "a" and self._cell == "size":
            self._row["download_path"] = a.get("href", "")

    def handle_endtag(self, tag: str) -> None:
        if self._row is None:
            return
        if tag == "td" and self._cell is not None:
            text = " ".join("".join(self._text).split())
            key = {
                "participant": "participant",
                "file-source": "source",
                "name": "name",
            }.get(self._cell)
            if key is not None:
                self._row[key] = text
            elif self._cell == "list-distinct":
                # Two cells share this marker: the published date comes first, the data type second.
                self._row["data_type" if "published" in self._row else "published"] = text
            self._cell = None
        elif tag == "tr":
            self.rows.append(self._row)
            self._row = None

    def handle_data(self, data: str) -> None:
        if self._cell is not None:
            self._text.append(data)


_DOWNLOAD_ID = re.compile(r"^/user_file/download/(\d+)$")
_PARTICIPANT = re.compile(r"^(hu[0-9A-F]+)(?:,\s*(PGP\d+))?$")


def parse_listing(html: str) -> tuple[list[ListingRow], Counter[str]]:
    """Parse the listing page. Returns the rows and a count of rows dropped, by reason."""
    parser = _ListingParser()
    parser.feed(html)
    rows: list[ListingRow] = []
    dropped: Counter[str] = Counter()
    for raw in parser.rows:
        m = _DOWNLOAD_ID.match(raw.get("download_path", ""))
        if m is None:
            dropped["no download link"] += 1
            continue
        if not raw.get("listed_size", "").isdigit():
            dropped["no file size"] += 1
            continue
        who = _PARTICIPANT.match(raw.get("participant", ""))
        if who is None:
            dropped["participant cell is not `huXXXXXX[, PGPnn]`"] += 1
            continue
        rows.append(
            ListingRow(
                download_id=int(m.group(1)),
                participant=who.group(1),
                pgp_alias=who.group(2),
                participant_numeric_id=int(raw["participant_numeric_id"]) if raw.get("participant_numeric_id", "").isdigit() else -1,
                published=raw.get("published", ""),
                data_type=raw.get("data_type", ""),
                source=raw.get("source", ""),
                name=raw.get("name", ""),
                locator=raw.get("locator", ""),
                listed_size=int(raw["listed_size"]),
            )
        )
    rows.sort(key=lambda r: r.download_id)
    return rows, dropped


class _CollectionParser(HTMLParser):
    """Reads the file listing of an Arvados collection page: `<li>… <size> … <a class="item" href>`."""

    def __init__(self) -> None:
        super().__init__()
        self.items: list[tuple[str, str]] = []  # (text before the link, href)
        self._in_li = False
        self._li_text: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        a = {k: v or "" for k, v in attrs}
        if tag == "li":
            self._in_li = True
            self._li_text = []
        elif tag == "a" and self._in_li and "item" in a.get("class", "").split():
            self.items.append(("".join(self._li_text), a.get("href", "")))

    def handle_endtag(self, tag: str) -> None:
        if tag == "li":
            self._in_li = False

    def handle_data(self, data: str) -> None:
        if self._in_li:
            self._li_text.append(data)


def parse_collection(html: str, base_url: str) -> list[CollectionFile]:
    """Files listed on an Arvados collection page, with absolute URLs, sorted by name."""
    parser = _CollectionParser()
    parser.feed(html)
    files: list[CollectionFile] = []
    for text, href in parser.items:
        digits = re.findall(r"\d+", text.replace("\xa0", " "))
        if not digits or href.endswith("/"):
            continue  # a subdirectory entry carries no size; nested collections are not expected here
        files.append(CollectionFile(name=unquote(href.removeprefix("./")), size=int(digits[-1]), url=urljoin(base_url, href)))
    files.sort(key=lambda f: f.name)
    return files


def _safe(component: str) -> str:
    """One path component, with separators and leading dots neutralised."""
    cleaned = component.replace("/", "_").replace("\\", "_").strip()
    return cleaned.lstrip(".") or "_"


class Fetcher:
    """Sequential HTTP client that waits `delay_s` between requests to the same host."""

    def __init__(self, delay_s: float) -> None:
        self.delay_s = delay_s
        self._last: dict[str, float] = {}
        self.client = httpx.Client(headers={"User-Agent": USER_AGENT}, timeout=httpx.Timeout(60.0, read=300.0))

    def _wait(self, url: str) -> None:
        host = httpx.URL(url).host
        elapsed = time.monotonic() - self._last.get(host, 0.0)
        if elapsed < self.delay_s:
            time.sleep(self.delay_s - elapsed)
        self._last[host] = time.monotonic()

    def get(self, url: str, *, follow: bool = True) -> httpx.Response:
        self._wait(url)
        return self.client.get(url, follow_redirects=follow)

    def open(self, url: str, *, follow: bool) -> AbstractContextManager[httpx.Response]:
        """A streamed GET; the caller decides from the headers whether to read the body."""
        self._wait(url)
        return self.client.stream("GET", url, follow_redirects=follow)


def _write_body(resp: httpx.Response, dest: Path) -> int:
    """Write a response body to `dest.part`. Returns bytes written; the caller renames on success."""
    part = dest.with_name(dest.name + ".part")
    written = 0
    with part.open("wb") as fh:
        for chunk in resp.iter_bytes(CHUNK):
            fh.write(chunk)
            written += len(chunk)
    return written


def _attachment_name(resp: httpx.Response) -> str | None:
    m = re.search(r'filename\*?=(?:UTF-8\'\')?"?([^";]+)"?', resp.headers.get("content-disposition", ""))
    return unquote(m.group(1)) if m else None


def _finish(row: ListingRow, name: str, dest: Path, written: int, expected: int) -> FileResult:
    now = datetime.now(UTC).isoformat()
    part = dest.with_name(dest.name + ".part")
    if written != expected:
        return FileResult(row, Outcome.SIZE_MISMATCH, name, str(part), written,
                          reason=f"expected {expected} bytes, received {written}; left as .part", fetched_at=now)
    part.replace(dest)
    return FileResult(row, Outcome.DOWNLOADED, name, str(dest), written, fetched_at=now)


def _already_present(row: ListingRow, target_dir: Path) -> list[FileResult] | None:
    """Offline resume check: the row's directory already holds complete files totalling the listed size."""
    if not target_dir.is_dir():
        return None
    files = sorted(f for f in target_dir.iterdir() if f.is_file() and not f.name.endswith(".part"))
    if not files or sum(f.stat().st_size for f in files) != row.listed_size:
        return None
    return [FileResult(row, Outcome.ALREADY_PRESENT, f.name, str(f), f.stat().st_size) for f in files]


def fetch_row(fetcher: Fetcher, row: ListingRow, out: Path, *, max_bytes: int, dry_run: bool) -> list[FileResult]:
    """Download every file behind one listing row.

    The download link answers in one of three ways, all seen live: older uploads are served directly
    (200 with a Content-Disposition attachment); newer ones redirect to an Arvados collection, landing
    either on its directory page (`…/_/`, which lists the files) or on the file itself.
    """
    if row.listed_size > max_bytes:
        return [FileResult(row, Outcome.SKIPPED_TOO_LARGE, reason=f"listed {row.listed_size} bytes > cap {max_bytes}")]
    target_dir = out / _safe(row.data_type) / _safe(row.participant) / str(row.download_id)
    present = _already_present(row, target_dir)
    if present is not None:
        return present
    if dry_run:
        return [FileResult(row, Outcome.DRY_RUN)]

    try:
        with fetcher.open(f"{PGP_BASE}/user_file/download/{row.download_id}", follow=False) as resp:
            location = resp.headers.get("location")
            if resp.is_redirect and location:
                pass  # handled below, after this response is closed
            elif resp.status_code == 200 and not resp.headers.get("content-type", "").startswith("text/html"):
                name = _attachment_name(resp) or f"{row.download_id}.bin"
                dest = target_dir / _safe(name)
                dest.parent.mkdir(parents=True, exist_ok=True)
                return [_finish(row, name, dest, _write_body(resp, dest), row.listed_size)]
            else:
                return [FileResult(row, Outcome.UNEXPECTED_RESPONSE,
                                   reason=f"HTTP {resp.status_code} {resp.headers.get('content-type', '')}".strip())]
    except httpx.HTTPError as err:
        return [FileResult(row, Outcome.HTTP_ERROR, reason=f"{type(err).__name__}: {err}")]

    # The redirect lands either on a collection's directory page (`…/_/`) or straight on a file in it.
    try:
        with fetcher.open(location, follow=True) as page:
            if page.status_code in (401, 403):
                return [FileResult(row, Outcome.NOT_PUBLIC, reason=f"redirect target HTTP {page.status_code}")]
            if page.status_code != 200:
                return [FileResult(row, Outcome.HTTP_ERROR, reason=f"redirect target HTTP {page.status_code}")]
            is_listing = page.url.path.endswith("/") and page.headers.get("content-type", "").startswith("text/html")
            if not is_listing:
                name = _attachment_name(page) or unquote(page.url.path.rsplit("/", 1)[-1]) or f"{row.download_id}.bin"
                dest = target_dir / _safe(name)
                dest.parent.mkdir(parents=True, exist_ok=True)
                return [_finish(row, name, dest, _write_body(page, dest), row.listed_size)]
            page.read()
            listing_html, listing_url = page.text, str(page.url)
    except httpx.HTTPError as err:
        return [FileResult(row, Outcome.HTTP_ERROR, reason=f"{type(err).__name__}: {err}")]
    files = parse_collection(listing_html, listing_url)
    if not files:
        return [FileResult(row, Outcome.EMPTY_COLLECTION, reason="collection lists no files")]

    results: list[FileResult] = []
    for f in files:
        dest = target_dir / _safe(f.name)
        dest.parent.mkdir(parents=True, exist_ok=True)
        try:
            with fetcher.open(f.url, follow=True) as resp:
                resp.raise_for_status()
                written = _write_body(resp, dest)
        except httpx.HTTPError as err:
            results.append(FileResult(row, Outcome.HTTP_ERROR, f.name, reason=f"{type(err).__name__}: {err}"))
            continue
        results.append(_finish(row, f.name, dest, written, f.size))
    return results


def _manifest_key(rec: dict) -> tuple[int, str]:
    return (rec["row"]["download_id"], rec.get("file_name") or "")


def write_manifest(path: Path, results: list[FileResult]) -> None:
    """Merge this run's results into the manifest: newest record per (download id, file) wins, sorted."""
    merged: dict[tuple[int, str], dict] = {}
    if path.exists():
        for line in path.read_text().splitlines():
            if line.strip():
                rec = json.loads(line)
                merged[_manifest_key(rec)] = rec
    for r in results:
        if r.outcome is Outcome.DRY_RUN:
            continue
        rec = asdict(r)
        rec["outcome"] = r.outcome.value
        merged[_manifest_key(rec)] = rec
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text("".join(json.dumps(merged[k], sort_keys=True) + "\n" for k in sorted(merged)))
    tmp.replace(path)


app = typer.Typer(add_completion=False, help=__doc__.split("\n\n")[0])


@app.command()
def main(
    data_type: str = typer.Option("23andMe", help="PGP data type, as the listing's filter spells it."),
    out: Path | None = typer.Option(None, help="Download root. Default: $JUST_VCF_SAMPLES_DIR/pgp-hms, else data/input/pgp-hms."),
    limit: int | None = typer.Option(None, help="Process at most this many listing rows (lowest download ids first)."),
    download_id: list[int] = typer.Option([], "--id", help="Only these download ids. Repeatable."),
    max_bytes: int = typer.Option(DEFAULT_MAX_BYTES, help="Skip rows whose listed size exceeds this."),
    delay: float = typer.Option(DEFAULT_DELAY_S, help="Seconds between requests to one host (robots.txt asks for 5)."),
    dry_run: bool = typer.Option(False, "--dry-run", help="Parse the listing and report; download nothing."),
) -> None:
    load_dotenv()
    if out is None:
        samples = os.environ.get("JUST_VCF_SAMPLES_DIR", "")
        out = Path(samples) / "pgp-hms" if samples else FALLBACK_OUT
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    fetcher = Fetcher(delay)
    listing_url = f"{PGP_BASE}{LISTING_PATH}?data_type={quote(data_type)}"
    resp = fetcher.get(listing_url)
    resp.raise_for_status()
    rows, dropped = parse_listing(resp.text)
    for reason, n in sorted(dropped.items()):
        log.warning("listing: dropped %d row(s): %s", n, reason)
    off_type = Counter(r.data_type for r in rows if r.data_type != data_type)
    for other, n in sorted(off_type.items()):
        log.warning("listing: %d row(s) carry data type %r, not %r; kept", n, other, data_type)
    if download_id:
        wanted = set(download_id)
        missing = sorted(wanted - {r.download_id for r in rows})
        if missing:
            raise typer.BadParameter(f"not in the {data_type} listing: {missing}", param_hint="--id")
        rows = [r for r in rows if r.download_id in wanted]
    selected = rows[:limit] if limit is not None else rows
    log.info("listing: %d row(s), %d selected, %.2f GB listed in selection",
             len(rows), len(selected), sum(r.listed_size for r in selected) / 1e9)

    results: list[FileResult] = []
    manifest = out / "manifest.jsonl"
    for i, row in enumerate(selected, 1):
        row_results = fetch_row(fetcher, row, out, max_bytes=max_bytes, dry_run=dry_run)
        results.extend(row_results)
        for r in row_results:
            if not r and r.outcome not in (Outcome.SKIPPED_TOO_LARGE, Outcome.DRY_RUN):
                log.warning("row %d (%s): %s — %s", row.download_id, row.participant, r.outcome.value, r.reason)
        if not dry_run and i % 10 == 0:
            write_manifest(manifest, results)  # checkpoint, so an interrupted run loses at most ten rows
    if not dry_run:
        write_manifest(manifest, results)

    by_outcome = Counter(r.outcome.value for r in results)
    typer.echo(" · ".join(f"{k}: {v}" for k, v in sorted(by_outcome.items())))
    if not dry_run:
        typer.echo(f"manifest: {manifest}")
    if any(r.outcome in (Outcome.HTTP_ERROR, Outcome.SIZE_MISMATCH, Outcome.UNEXPECTED_RESPONSE) for r in results):
        raise typer.Exit(code=1)


if __name__ == "__main__":
    app()
