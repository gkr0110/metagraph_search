#!/usr/bin/env python3
"""
Galaxy tool wrapper for MetaGraph Online API.
https://metagraph.ethz.ch

Submits FASTA sequences to the MetaGraph async REST API and retrieves
results as a tab-separated file.

API flow
--------
  1. POST /search         → search_id
  2. GET  /search/{id}/status  → poll until "completed"
  3. GET  /search/{id}/results → download and format hits
"""

import argparse
import csv
import json
import sys
import time

try:
    import requests
except ImportError:
    sys.exit("ERROR: 'requests' library is required. Install with: pip install requests")

API_BASE = "https://metagraph.ethz.ch:8081"

FIELDNAMES = [
    "query", "database", "label",
    "normalized_score", "kmer_count",
    "coordinates", "organism", "taxid",
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def eprint(*args, **kwargs):
    """Print to stderr."""
    print(*args, file=sys.stderr, **kwargs)


def parse_fasta(filepath):
    """Parse FASTA file; return list of (name, sequence) tuples."""
    records = []
    name, parts = None, []
    with open(filepath) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(parts)))
                # Use first word of the FASTA header as the sequence name
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line.upper())
    if name is not None:
        records.append((name, "".join(parts)))
    return records


# ---------------------------------------------------------------------------
# API calls
# ---------------------------------------------------------------------------

def submit_search(session, sequence, databases, params, sequence_type,
                  search_name=None):
    """POST a new search request; return search_id string."""
    queries = [{"db": db, "params": dict(params)} for db in databases]
    payload = {
        "sequence": sequence,
        "queries": queries,
        "sequence_type": sequence_type,
    }
    if search_name:
        payload["search_name"] = search_name[:255]

    url = f"{API_BASE}/search"
    resp = session.post(url, json=payload, timeout=60)
    resp.raise_for_status()
    data = resp.json()
    return data["search_id"]


def poll_until_done(session, search_id, poll_interval, max_wait):
    """Poll the status endpoint; return True when completed, False on error/timeout.

    The MetaGraph async API uses the status value "done" (not "completed") to
    indicate a successfully finished search.  Both spellings are accepted here
    for forward-compatibility.
    """
    url = f"{API_BASE}/search/{search_id}/status"
    deadline = time.time() + max_wait
    start = time.time()
    while time.time() < deadline:
        resp = session.get(url, timeout=60)
        resp.raise_for_status()
        data = resp.json()
        status = data.get("status", "unknown").lower()

        # "done" is the terminal success state returned by the current API.
        # "completed" is kept for forward-compatibility with any future change.
        if status in ("done", "completed"):
            return True
        if status in ("failed", "error", "cancelled"):
            eprint(f"  [WARN] Search {search_id} ended with status "
                   f"'{status}': {data.get('message', '')}")
            return False

        elapsed = int(time.time() - start)
        progress = data.get("progress", "")
        eprint(f"  [{elapsed:5d}s] {search_id}: {status}"
               + (f" ({progress})" if progress else ""))
        time.sleep(poll_interval)

    eprint(f"  [TIMEOUT] Gave up waiting for {search_id} after {max_wait}s")
    return False


def fetch_results(session, search_id):
    """GET results JSON for a completed search."""
    url = f"{API_BASE}/search/{search_id}/results"
    # Use a generous timeout — large result sets (e.g. UniParc) can be slow to
    # serialize on the server side.
    resp = session.get(url, timeout=300)
    resp.raise_for_status()
    return resp.json()


# ---------------------------------------------------------------------------
# Result flattening
# ---------------------------------------------------------------------------

def extract_rows(query_name, results_json, debug=False):
    """
    Flatten the MetaGraph async API results JSON into a list of row dicts.

    The API (v0.9+) returns a top-level ``results`` key containing a flat list
    of hit objects.  Each hit carries nested ``core_metrics``, ``payload``,
    ``enrichment``, and ``table_data`` sub-objects.

    Example structure (abbreviated):
    {
      "results": [
        {
          "sequence_id": "UPI0000129DB5",
          "score": 380.0,
          "core_metrics": {
              "database": "uniparc",
              "normalized_score": 1.0,
              "kmer_count": 380,
              "query_sequence_id": "s1"
          },
          "kmer_coords": "0-0-379",
          "enrichment": {
              "organism": "Mycoplasma capricolum",
              "taxid": 2095
          }
        }, ...
      ],
      "metadata": {"resultCount": 5},
      "complete": true
    }
    """
    if debug:
        eprint("DEBUG raw JSON (first 4000 chars):",
               json.dumps(results_json, indent=2)[:4000])

    rows = []

    # Unwrap top-level container
    if isinstance(results_json, list):
        hit_list = results_json
    elif isinstance(results_json, dict):
        hit_list = results_json.get("results") or []
    else:
        eprint("  [WARN] Unexpected results type:", type(results_json))
        return rows

    if debug:
        eprint(f"DEBUG: {len(hit_list)} raw hit(s) in results list")

    for hit in hit_list:
        if not isinstance(hit, dict):
            continue

        core       = hit.get("core_metrics") or {}
        payload    = hit.get("payload")      or {}
        enrichment = hit.get("enrichment")   or {}
        table_data = hit.get("table_data")   or {}

        db = (core.get("database")
              or payload.get("db")
              or table_data.get("database")
              or "unknown")

        label = (hit.get("sequence_id")
                 or core.get("sequence_id")
                 or table_data.get("sequence_id", ""))

        normalized_score = (core.get("normalized_score")
                            or table_data.get("normalized_score", ""))

        kmer_count = (core.get("kmer_count")
                      or table_data.get("kmer_count", ""))

        coords = hit.get("kmer_coords", "")

        organism = (enrichment.get("organism")
                    or table_data.get("organism", ""))

        taxid = enrichment.get("taxid", "")

        rows.append({
            "query":            query_name,
            "database":         db,
            "label":            label,
            "normalized_score": normalized_score,
            "kmer_count":       kmer_count,
            "coordinates":      coords,
            "organism":         organism,
            "taxid":            taxid,
        })

    return rows


# ---------------------------------------------------------------------------
# Mock response (used by --mock for planemo testing without network)
# ---------------------------------------------------------------------------

MOCK_RESULTS = {
    "results": [
        {
            "sequence_id": "NC_000913.3",
            "score": 310.0,
            "core_metrics": {
                "database": "refseq33m",
                "normalized_score": 0.98,
                "kmer_count": 310,
                "query_sequence_id": "s1",
            },
            "kmer_coords": "0-0-309",
            "enrichment": {
                "organism": "Escherichia coli str. K-12 substr. MG1655",
                "taxid": 511145,
            },
            "table_data": {},
            "payload": {"db": "refseq33m"},
        },
        {
            "sequence_id": "NC_002695.2",
            "score": 265.0,
            "core_metrics": {
                "database": "refseq33m",
                "normalized_score": 0.85,
                "kmer_count": 265,
                "query_sequence_id": "s1",
            },
            "kmer_coords": "0-0-264",
            "enrichment": {
                "organism": "Escherichia coli O157:H7 str. Sakai",
                "taxid": 386585,
            },
            "table_data": {},
            "payload": {"db": "refseq33m"},
        },
    ],
    "metadata": {"resultCount": 2},
    "complete": True,
    "search_complete": True,
}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Submit FASTA sequences to MetaGraph Online and write hits as TSV"
    )
    parser.add_argument("--input", required=True,
                        help="Input FASTA file")
    parser.add_argument("--databases", required=True,
                        help="Comma-separated list of MetaGraph database IDs")
    parser.add_argument("--mode", choices=["exact", "align"], default="exact",
                        help="Search mode: exact k-mer match or alignment")
    parser.add_argument("--sequence-type",
                        choices=["nucleotide", "amino_acid"],
                        default="nucleotide",
                        help="Molecule type of the input sequences")
    parser.add_argument("--discovery-threshold", type=float, default=0.5,
                        help="Minimum k-mer match fraction (0.0–1.0)")
    parser.add_argument("--top-labels", type=int, default=100,
                        help="Maximum matched accessions per database")
    parser.add_argument("--min-exact-match", type=float, default=0.0,
                        help="Alignment: min fraction covered by seeds (0.0–1.0)")
    parser.add_argument("--max-alternative-alignments", type=int, default=1,
                        help="Alignment: number of alternative alignments per sequence")
    parser.add_argument("--output", required=True,
                        help="Output TSV file")
    parser.add_argument("--poll-interval", type=int, default=10,
                        help="Seconds between status polls")
    parser.add_argument("--max-wait", type=int, default=1800,
                        help="Maximum seconds to wait per search")
    parser.add_argument("--debug", action="store_true",
                        help="Print raw API JSON to stderr for debugging")
    parser.add_argument("--mock", action="store_true",
                        help="Use a built-in mock response (for testing without network)")
    args = parser.parse_args()

    # Parse databases
    databases = [d.strip() for d in args.databases.split(",") if d.strip()]
    if not databases:
        sys.exit("ERROR: No databases specified (--databases is empty)")

    # Build per-query params dict
    params = {
        "discovery_threshold": args.discovery_threshold,
        "top_labels": args.top_labels,
    }
    if args.mode == "align":
        params["align"] = True
        params["min_exact_match"] = args.min_exact_match
        params["max_alternative_alignments"] = args.max_alternative_alignments

    # Parse FASTA
    eprint("Parsing FASTA input...")
    sequences = parse_fasta(args.input)
    if not sequences:
        sys.exit("ERROR: No sequences found in input FASTA file")

    eprint(f"Found {len(sequences)} sequence(s) to search against "
           f"{len(databases)} database(s): {', '.join(databases)}")

    # HTTP session
    session = requests.Session()
    session.headers.update({
        "Content-Type": "application/json",
        "Accept": "application/json",
    })

    all_rows = []

    for idx, (name, seq) in enumerate(sequences, 1):
        eprint(f"\n[{idx}/{len(sequences)}] Sequence: {name} ({len(seq)} bp)")

        # ----- Mock mode (no network) -----
        if args.mock:
            eprint("  [MOCK] Returning built-in mock results")
            rows = extract_rows(name, MOCK_RESULTS, debug=args.debug)
            all_rows.extend(rows)
            eprint(f"  [MOCK] {len(rows)} hit(s)")
            continue

        # ----- Live API -----
        try:
            eprint("  Submitting search to MetaGraph API...")
            search_id = submit_search(
                session, seq, databases, params,
                args.sequence_type,
                search_name=f"Galaxy:{name[:200]}",
            )
            eprint(f"  search_id = {search_id}")

            eprint("  Polling for completion (this may take several minutes)...")
            ok = poll_until_done(
                session, search_id,
                args.poll_interval, args.max_wait
            )
            if not ok:
                eprint(f"  Skipping {name}: search did not complete")
                continue

            eprint("  Fetching results...")
            results_json = fetch_results(session, search_id)
            rows = extract_rows(name, results_json, debug=args.debug)
            all_rows.extend(rows)
            eprint(f"  {len(rows)} hit(s)")

        except requests.exceptions.HTTPError as exc:
            eprint(f"  [HTTP ERROR] {exc}")
            if args.debug and exc.response is not None:
                eprint(f"  Response body: {exc.response.text[:1000]}")
        except requests.exceptions.ConnectionError as exc:
            eprint(f"  [CONNECTION ERROR] {exc}")
        except requests.exceptions.RequestException as exc:
            eprint(f"  [REQUEST ERROR] {exc}")

    # Write output TSV
    eprint(f"\nWriting {len(all_rows)} total hit(s) to {args.output}")
    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=FIELDNAMES, delimiter="\t",
            extrasaction="ignore", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(all_rows)

    eprint("Done.")


if __name__ == "__main__":
    main()
