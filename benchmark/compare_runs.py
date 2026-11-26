#!/usr/bin/env python3
"""
Stage-by-stage comparison of two annotation runs.

Each run is either an ORIGINAL (Snakemake, mgg_annotation@prodigal-gv) output
directory or a NEW (Nextflow) output directory produced with
PUBLISH_INTERMEDIATES: true.

    python benchmark/compare_runs.py \
        --a ORIG_OUTPUT_DIR --a-layout orig --a-label original \
        --b NEW_OUTPUT_DIR  --b-layout new  --b-label nextflow \
        --out comparison.md

Requires pandas and biopython (both are in the pipeline's Docker image):

    docker run --rm -u $(id -u):$(id -g) -v "$PWD":"$PWD" -w "$PWD" \
        phage_annotation:latest python benchmark/compare_runs.py ...

Stages:
  1. Input & ORFs   - tables must be identical (after row sorting)
  2. Clustering     - partition agreement (PC ids are not stable; matched via proteins)
  3. Hits           - per-PC top-hit agreement and score correlation, per database
  4. Annotation     - per-protein function agreement, with a likely cause for each mismatch
  5. GenBank        - CDS coordinates and translations
"""

import argparse
import math
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd
from Bio import SeqIO

DBS = ["PHROGS", "PFAM", "ECOD"]

# Annotation columns compared per protein: logical name -> column per layout
FUNCTION_COLS = {
    "PHROGS1": {"orig": "PHROGS1_function", "new": "PHROGS1_function"},
    "PHROGS2": {"orig": "PHROGS2_function", "new": "PHROGS2_function"},
    "PFAM":    {"orig": "PFAM_function",    "new": "PFAM1_function"},
    "ECOD":    {"orig": "ECOD_function",    "new": "ECOD1_function"},
}

ORF_TABLES = [
    # (label, relative path in the original layout, sort keys)
    ("prodigal.csv", "2_ORF_PREDICTION/1_GENE_CALLING/1_PRODIGALGV/prodigal.csv", None),
    ("unfiltered-glimmer.csv", "2_ORF_PREDICTION/1_GENE_CALLING/2_GLIMMER/unfiltered-glimmer.csv", None),
    ("glimmer.csv", "2_ORF_PREDICTION/1_GENE_CALLING/2_GLIMMER/glimmer.csv", None),
    ("removed.csv", "2_ORF_PREDICTION/1_GENE_CALLING/2_GLIMMER/removed.csv", None),
    ("orfs.csv", "2_ORF_PREDICTION/2_PROCESSING/orfs.csv", None),
    ("ambigous.csv", "2_ORF_PREDICTION/2_PROCESSING/ambigous.csv", None),
    ("confident_orfs.csv", "2_ORF_PREDICTION/confident_orfs.csv", ["proteinID"]),
]


# ---------------------------------------------------------------------------
# Run layout
# ---------------------------------------------------------------------------
class Run:
    def __init__(self, root, layout, label):
        self.root = Path(root)
        self.layout = layout
        self.label = label
        if not self.root.is_dir():
            sys.exit(f"Not a directory: {self.root}")

    def orf_path(self, orig_rel):
        """ORF-stage and protein files: NEW publishes them under intermediates/ with the original layout."""
        if self.layout == "orig":
            return self.root / orig_rel
        return self.root / "intermediates" / orig_rel

    @property
    def metadata(self):
        return self.orf_path("2_ORF_PREDICTION/metadata.tsv")

    @property
    def proteins(self):
        return self.orf_path("3_ANNOTATION/proteins.fasta")

    @property
    def raw_pcs(self):
        if self.layout == "orig":
            hits = sorted(self.root.glob("3_ANNOTATION/1_CLUSTERING_*/raw_PCs.tsv"))
            return hits[0] if hits else self.root / "3_ANNOTATION/1_CLUSTERING/raw_PCs.tsv"
        return self.root / "intermediates/3_ANNOTATION/1_CLUSTERING/raw_PCs.tsv"

    @property
    def pcs2proteins(self):
        if self.layout == "orig":
            return self.root / "3_ANNOTATION/PCs2proteins.tsv"
        return self.root / "PCs2proteins.tsv"

    @property
    def search(self):
        return self.root / "search.tsv"

    @property
    def filtered_hits(self):
        """Hits that passed the hit filter (NEW only; ORIGINAL had no such table)."""
        return self.root / "report.tsv" if self.layout == "new" else None

    @property
    def annotation(self):
        if self.layout == "orig":
            return self.root / "annotation.tsv"
        full = self.root / "intermediates/annotation_full.tsv"
        return full if full.exists() else self.root / "annotation.tsv"

    @property
    def genbank_dir(self):
        return self.root / ("4_GENBANK" if self.layout == "orig" else "genbanks")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
class Report:
    def __init__(self):
        self.lines = []
        self.verdicts = []

    def h(self, text, level=2):
        self.lines += ["", "#" * level + " " + text, ""]

    def p(self, text=""):
        self.lines.append(text)

    def table(self, header, rows):
        self.lines.append("| " + " | ".join(header) + " |")
        self.lines.append("|" + "---|" * len(header))
        for r in rows:
            self.lines.append("| " + " | ".join(str(x) for x in r) + " |")
        self.lines.append("")

    def verdict(self, stage, status, detail):
        self.verdicts.append((stage, status, detail))

    def render(self, title):
        out = [f"# {title}", "", "## Summary", ""]
        out.append("| Stage | Result | Detail |")
        out.append("|---|---|---|")
        for stage, status, detail in self.verdicts:
            out.append(f"| {stage} | {status} | {detail} |")
        return "\n".join(out + self.lines) + "\n"


def missing(paths):
    return [str(p) for p in paths if p is None or not Path(p).exists()]


def read_str_table(path):
    sep = "\t" if str(path).endswith(".tsv") else ","
    return pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False)


def pct(n, d):
    return f"{100 * n / d:.1f}%" if d else "n/a"


def pearson(xs, ys):
    n = len(xs)
    if n < 3:
        return float("nan")
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    if sxx == 0 or syy == 0:
        return float("nan")
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / math.sqrt(sxx * syy)


def comb2(n):
    return n * (n - 1) // 2


def adjusted_rand_index(labels_a, labels_b):
    """ARI for two labelings of the same items (dicts item -> cluster label)."""
    items = sorted(set(labels_a) & set(labels_b))
    n = len(items)
    if n < 2:
        return float("nan")
    contingency = Counter((labels_a[i], labels_b[i]) for i in items)
    a_sizes = Counter(labels_a[i] for i in items)
    b_sizes = Counter(labels_b[i] for i in items)
    sum_comb = sum(comb2(v) for v in contingency.values())
    sum_a = sum(comb2(v) for v in a_sizes.values())
    sum_b = sum(comb2(v) for v in b_sizes.values())
    expected = sum_a * sum_b / comb2(n)
    max_index = (sum_a + sum_b) / 2
    if max_index == expected:
        return 1.0
    return (sum_comb - expected) / (max_index - expected)


# ---------------------------------------------------------------------------
# Stage 1: input & ORFs
# ---------------------------------------------------------------------------
def compare_tables(path_a, path_b, sort_keys=None):
    """Return (identical, detail) for two CSV/TSV tables compared as strings, order-insensitive."""
    a, b = read_str_table(path_a), read_str_table(path_b)
    if list(a.columns) != list(b.columns):
        return False, f"columns differ: {list(a.columns)} vs {list(b.columns)}"
    keys = sort_keys or list(a.columns)
    a = a.sort_values(keys + [c for c in a.columns if c not in keys], kind="mergesort").reset_index(drop=True)
    b = b.sort_values(keys + [c for c in b.columns if c not in keys], kind="mergesort").reset_index(drop=True)
    if len(a) != len(b):
        ra = set(map(tuple, a.values.tolist()))
        rb = set(map(tuple, b.values.tolist()))
        return False, f"{len(a)} vs {len(b)} rows; {len(ra - rb)} only in A, {len(rb - ra)} only in B"
    diff = (a != b).any(axis=1)
    if diff.any():
        cols = [c for c in a.columns if (a[c] != b[c]).any()]
        return False, f"{int(diff.sum())} of {len(a)} rows differ (columns: {', '.join(cols)})"
    return True, f"{len(a)} rows identical"


def stage_orfs(A, B, rep):
    rep.h("1. Input & ORF prediction")
    rows, ok_all = [], True

    # metadata: contig set and lengths must match; other columns are informational
    if not missing([A.metadata, B.metadata]):
        ma, mb = read_str_table(A.metadata), read_str_table(B.metadata)
        la = dict(zip(ma["contigID"], ma["contig_len [bp]"]))
        lb = dict(zip(mb["contigID"], mb["contig_len [bp]"]))
        same = la == lb
        ok_all &= same
        detail = f"{len(la)} contigs, lengths identical" if same else (
            f"only in A: {sorted(set(la) - set(lb))}; only in B: {sorted(set(lb) - set(la))}; "
            f"length mismatches: {[k for k in set(la) & set(lb) if la[k] != lb[k]]}")
        rows.append(("metadata.tsv (contigs, lengths)", "identical" if same else "DIFFERENT", detail))
        other = []
        for col in ("file_name", "contig_description"):
            if col in ma.columns and col in mb.columns:
                da = dict(zip(ma["contigID"], ma[col]))
                db = dict(zip(mb["contigID"], mb[col]))
                n = sum(1 for k in set(da) & set(db) if da[k] != db[k])
                if n:
                    other.append(f"{col}: {n} differ")
        if other:
            rows.append(("metadata.tsv (other columns)", "informational", "; ".join(other)))
    else:
        rows.append(("metadata.tsv", "missing", ", ".join(missing([A.metadata, B.metadata]))))

    for label, rel, keys in ORF_TABLES:
        pa, pb = A.orf_path(rel), B.orf_path(rel)
        miss = missing([pa, pb])
        if miss:
            rows.append((label, "missing", ", ".join(miss)))
            continue
        same, detail = compare_tables(pa, pb, keys)
        ok_all &= same
        rows.append((label, "identical" if same else "DIFFERENT", detail))

    if not missing([A.proteins, B.proteins]):
        sa = {r.id: str(r.seq) for r in SeqIO.parse(A.proteins, "fasta")}
        sb = {r.id: str(r.seq) for r in SeqIO.parse(B.proteins, "fasta")}
        same = sa == sb
        ok_all &= same
        detail = f"{len(sa)} proteins identical" if same else (
            f"{len(set(sa) - set(sb))} only in A, {len(set(sb) - set(sa))} only in B, "
            f"{sum(1 for k in set(sa) & set(sb) if sa[k] != sb[k])} sequence mismatches")
        rows.append(("proteins.fasta (id -> sequence)", "identical" if same else "DIFFERENT", detail))
    else:
        rows.append(("proteins.fasta", "missing", ", ".join(missing([A.proteins, B.proteins]))))

    rep.table(["File", "Result", "Detail"], rows)
    n_missing = sum(1 for r in rows if r[1] == "missing")
    status = "PASS" if ok_all and not n_missing else ("FAIL" if not ok_all else "INCOMPLETE")
    rep.verdict("1. Input & ORFs", status, f"{sum(1 for r in rows if r[1] == 'identical')} identical, "
                f"{sum(1 for r in rows if r[1] == 'DIFFERENT')} different, {n_missing} missing")


# ---------------------------------------------------------------------------
# Stage 2: clustering
# ---------------------------------------------------------------------------
def load_clusters(path):
    df = pd.read_csv(path, sep="\t", header=None, names=["repr", "proteinID"], dtype=str)
    return dict(zip(df["proteinID"], df["repr"]))


def stage_clustering(A, B, rep, thresholds):
    rep.h("2. Clustering (MMseqs2)")
    miss = missing([A.raw_pcs, B.raw_pcs])
    if miss:
        rep.p("Missing: " + ", ".join(miss))
        rep.verdict("2. Clustering", "INCOMPLETE", "raw_PCs.tsv missing")
        return None, None

    ca, cb = load_clusters(A.raw_pcs), load_clusters(B.raw_pcs)
    common = set(ca) & set(cb)
    members_a, members_b = defaultdict(set), defaultdict(set)
    for p, r in ca.items():
        members_a[r].add(p)
    for p, r in cb.items():
        members_b[r].add(p)

    same_repr = sum(1 for p in common if ca[p] == cb[p])
    same_members = sum(1 for p in common if members_a[ca[p]] == members_b[cb[p]])
    ari = adjusted_rand_index(ca, cb)
    identical_clusters = {frozenset(m) for m in members_a.values()} & {frozenset(m) for m in members_b.values()}

    rep.table(["Metric", "Value"], [
        ("Proteins (A / B / common)", f"{len(ca)} / {len(cb)} / {len(common)}"),
        ("Clusters (A / B)", f"{len(members_a)} / {len(members_b)}"),
        ("Identical clusters (same member set)", f"{len(identical_clusters)} ({pct(len(identical_clusters), max(len(members_a), 1))} of A)"),
        ("Proteins with identical cluster members", f"{same_members} ({pct(same_members, len(common))})"),
        ("Proteins with the same representative", f"{same_repr} ({pct(same_repr, len(common))})"),
        ("Adjusted Rand index", f"{ari:.4f}"),
    ])
    rep.p("Expected: very high but not necessarily identical (MMseqs2 13.45111 -> 16.747c6).")

    frac = same_members / len(common) if common else 0
    status = "PASS" if frac == 1 else ("PASS (explained)" if frac >= thresholds["clustering"] else "REVIEW")
    rep.verdict("2. Clustering", status, f"{pct(same_members, len(common))} of proteins in identical clusters, ARI {ari:.3f}")
    return ca, cb


# ---------------------------------------------------------------------------
# Stage 3: hits
# ---------------------------------------------------------------------------
def load_pc_repr(path):
    """PC id -> representative protein (for matching PCs across runs)."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    return dict(zip(df["PC"], df["repr"]))


def load_hits(run):
    df = pd.read_csv(run.search, sep="\t", low_memory=False)
    df = df[df["db"].isin(DBS)].copy()
    pc2repr = load_pc_repr(run.pcs2proteins)
    df["repr"] = df["query"].map(pc2repr)
    return df


def top_hits(df):
    """(repr, db) -> target of the max-bits hit."""
    idx = df.groupby(["repr", "db"])["bits"].idxmax()
    best = df.loc[idx]
    return {(r, d): t for r, d, t in zip(best["repr"], best["db"], best["target"])}


def stage_hits(A, B, rep, thresholds):
    rep.h("3. HH-suite hits (search.tsv)")
    miss = missing([A.search, B.search, A.pcs2proteins, B.pcs2proteins])
    if miss:
        rep.p("Missing: " + ", ".join(miss))
        rep.verdict("3. Hits", "INCOMPLETE", "search.tsv / PCs2proteins.tsv missing")
        return None

    ha, hb = load_hits(A), load_hits(B)
    common_reprs = set(load_pc_repr(A.pcs2proteins).values()) & set(load_pc_repr(B.pcs2proteins).values())
    ha, hb = ha[ha["repr"].isin(common_reprs)], hb[hb["repr"].isin(common_reprs)]
    ta, tb = top_hits(ha), top_hits(hb)

    rows, worst = [], 1.0
    for db in DBS:
        keys_a = {k for k in ta if k[1] == db}
        keys_b = {k for k in tb if k[1] == db}
        both = keys_a & keys_b
        agree = sum(1 for k in both if ta[k] == tb[k])
        # score agreement on (PC, target) pairs found by both runs, best hit per pair
        sa = ha[ha["db"] == db].groupby(["repr", "target"]).agg(prob=("prob", "max"), bits=("bits", "max"))
        sb = hb[hb["db"] == db].groupby(["repr", "target"]).agg(prob=("prob", "max"), bits=("bits", "max"))
        j = sa.join(sb, lsuffix="_a", rsuffix="_b", how="inner")
        r_bits = pearson(j["bits_a"].tolist(), j["bits_b"].tolist())
        med_dprob = (j["prob_a"] - j["prob_b"]).abs().median() if len(j) else float("nan")
        union_pairs = len(sa.index.union(sb.index))
        rows.append((db, len(keys_a), len(keys_b), f"{agree}/{len(both)} ({pct(agree, len(both))})",
                     f"{len(j)}/{union_pairs} ({pct(len(j), union_pairs)})", f"{r_bits:.4f}", f"{med_dprob:.4f}"))
        if both:
            worst = min(worst, agree / len(both))

    rep.p(f"PCs matched across runs by representative protein: {len(common_reprs)}.")
    rep.table(["DB", "PCs with hits (A)", "PCs with hits (B)", "Same top hit", "Shared (PC, target) pairs",
               "Pearson r (bits)", "Median |Δprob|"], rows)
    rep.p("Small score differences are expected from the MSA format change (true A3M vs. aligned FASTA). "
          "Larger differences point to different database files or search settings.")
    status = "PASS" if worst == 1 else ("PASS (explained)" if worst >= thresholds["hits"] else "REVIEW")
    rep.verdict("3. Hits", status, f"lowest top-hit agreement across DBs: {worst:.1%}")
    return ha, hb, ta, tb


# ---------------------------------------------------------------------------
# Stage 4: annotation
# ---------------------------------------------------------------------------
def norm_function(value):
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "-"
    s = str(value).strip()
    if s in ("", "-", "nan", "0.0"):
        return "-"
    # NEW ECOD strings carry an "X: <x-group> " prefix
    s = re.sub(r"^X: \S+ ", "", s)
    return s


def functions_match(fa, fb):
    """ORIGINAL strings carry a 5-character confidence prefix; compare by suffix."""
    na, nb = norm_function(fa), norm_function(fb)
    if na == "-" or nb == "-":
        return na == nb
    return na == nb or na.endswith(nb) or nb.endswith(na)


def stage_annotation(A, B, rep, thresholds, clusters, hits):
    rep.h("4. Annotation (per protein)")
    miss = missing([A.annotation, B.annotation])
    if miss:
        rep.p("Missing: " + ", ".join(miss))
        rep.verdict("4. Annotation", "INCOMPLETE", "annotation table missing")
        return

    aa = pd.read_csv(A.annotation, sep="\t", dtype=str, keep_default_na=False).set_index("proteinID")
    ab = pd.read_csv(B.annotation, sep="\t", dtype=str, keep_default_na=False).set_index("proteinID")
    common = aa.index.intersection(ab.index)
    ca, cb = clusters if clusters else (None, None)
    ta, tb = (hits[2], hits[3]) if hits else (None, None)
    filtered_b = None
    if B.filtered_hits is not None and B.filtered_hits.exists():
        fb = pd.read_csv(B.filtered_hits, sep="\t", usecols=["query", "target", "db"], dtype=str)
        pc2repr_b = load_pc_repr(B.pcs2proteins)
        filtered_b = set(zip(fb["query"].map(pc2repr_b), fb["db"], fb["target"]))

    rows, examples, worst = [], [], 1.0
    all_causes = Counter()
    for logical, cols in FUNCTION_COLS.items():
        col_a, col_b = cols[A.layout], cols[B.layout]
        if col_a not in aa.columns or col_b not in ab.columns:
            rows.append((logical, "column missing", "", "", ""))
            continue
        db = re.sub(r"\d+$", "", logical)
        agree, causes = 0, Counter()
        for pid in common:
            fa, fb_ = aa.at[pid, col_a], ab.at[pid, col_b]
            if functions_match(fa, fb_):
                agree += 1
                continue
            # likely cause, most upstream first
            if ca is None:
                cause = "unknown (no clustering data)"
            elif pid not in ca or pid not in cb or ca[pid] != cb[pid]:
                cause = "clustering (different representative)"
            elif ta is not None and ta.get((ca[pid], db)) != tb.get((cb[pid], db)):
                cause = "search (different top hit)"
            elif filtered_b is not None and ta is not None and ta.get((ca[pid], db)) is not None and \
                    (cb[pid], db, ta.get((ca[pid], db))) not in filtered_b:
                cause = "hit filter (prob/qcov/tcov)"
            else:
                cause = "best-hit selection / formatting"
            causes[cause] += 1
            if len(examples) < 25:
                examples.append((pid, logical, norm_function(fa)[:70], norm_function(fb_)[:70], cause))
        all_causes.update(causes)
        n = len(common)
        rows.append((logical, f"{col_a} vs {col_b}", f"{agree}/{n} ({pct(agree, n)})",
                     "; ".join(f"{k}: {v}" for k, v in causes.most_common()) or "-"))
        if n:
            worst = min(worst, agree / n)

    rep.p(f"Proteins in A / B / common: {len(aa)} / {len(ab)} / {len(common)}. "
          "Strings are compared after removing the original's confidence prefix and the new `X: ...` ECOD prefix.")
    rep.table(["Function", "Columns", "Agreement", "Mismatch causes"], rows)
    if examples:
        rep.p("First mismatches:")
        rep.table(["proteinID", "Function", "A", "B", "Likely cause"], examples)
    # mismatches not traceable to an upstream (intentional) change always need a look
    unexplained = sum(v for k, v in all_causes.items() if k.startswith(("best-hit", "unknown")))
    if worst == 1:
        status = "PASS"
    elif worst >= thresholds["annotation"] and not unexplained:
        status = "PASS (explained)"
    else:
        status = "REVIEW"
    detail = f"lowest per-column agreement: {worst:.1%}"
    if unexplained:
        detail += f"; {unexplained} mismatches without an upstream cause"
    rep.verdict("4. Annotation", status, detail)


# ---------------------------------------------------------------------------
# Stage 5: GenBank
# ---------------------------------------------------------------------------
def load_genbank(directory):
    """record id -> sorted unique list of (start, end, strand, translation) for CDS features."""
    out = {}
    for gb in sorted(Path(directory).glob("*.gb")):
        for rec in SeqIO.parse(gb, "genbank"):
            cds = set()
            for f in rec.features:
                if f.type != "CDS":
                    continue
                tr = f.qualifiers.get("translation", [""])[0]
                cds.add((int(f.location.start), int(f.location.end), f.location.strand, tr))
            out[rec.id] = sorted(cds)
    return out


def stage_genbank(A, B, rep):
    rep.h("5. GenBank")
    miss = missing([A.genbank_dir, B.genbank_dir])
    if miss:
        rep.p("Missing: " + ", ".join(miss))
        rep.verdict("5. GenBank", "INCOMPLETE", "GenBank directory missing")
        return
    ga, gb = load_genbank(A.genbank_dir), load_genbank(B.genbank_dir)
    common = set(ga) & set(gb)
    same = [k for k in common if ga[k] == gb[k]]
    rows = [
        ("Records (A / B / common)", f"{len(ga)} / {len(gb)} / {len(common)}"),
        ("Records with identical CDS set (location, strand, translation)", f"{len(same)} ({pct(len(same), len(common))})"),
        ("CDS (A / B)", f"{sum(len(v) for v in ga.values())} / {sum(len(v) for v in gb.values())}"),
    ]
    rep.table(["Metric", "Value"], rows)
    rep.p("CDS are de-duplicated before comparing: the original writes the last CDS of each contig twice.")
    ok = len(same) == len(common) and set(ga) == set(gb)
    rep.verdict("5. GenBank", "PASS" if ok else "FAIL", f"{len(same)}/{len(common)} records identical")


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--a", required=True, help="Output directory of run A")
    ap.add_argument("--a-layout", choices=["orig", "new"], required=True)
    ap.add_argument("--a-label", default="A")
    ap.add_argument("--b", required=True, help="Output directory of run B")
    ap.add_argument("--b-layout", choices=["orig", "new"], required=True)
    ap.add_argument("--b-label", default="B")
    ap.add_argument("--out", default="comparison.md", help="Markdown report path")
    ap.add_argument("--min-clustering", type=float, default=0.95,
                    help="Min. fraction of proteins in identical clusters to pass as 'explained'")
    ap.add_argument("--min-hits", type=float, default=0.95, help="Min. top-hit agreement per DB")
    ap.add_argument("--min-annotation", type=float, default=0.95, help="Min. per-column function agreement")
    args = ap.parse_args()

    A = Run(args.a, args.a_layout, args.a_label)
    B = Run(args.b, args.b_layout, args.b_label)
    thresholds = {"clustering": args.min_clustering, "hits": args.min_hits, "annotation": args.min_annotation}

    rep = Report()
    rep.p("")
    rep.p(f"- **A** ({A.layout}): `{A.label}` — `{A.root}`")
    rep.p(f"- **B** ({B.layout}): `{B.label}` — `{B.root}`")
    rep.p(f"- Thresholds for 'explained' differences: clustering {args.min_clustering:.0%}, "
          f"hits {args.min_hits:.0%}, annotation {args.min_annotation:.0%}")

    stage_orfs(A, B, rep)
    clusters = stage_clustering(A, B, rep, thresholds)
    hits = stage_hits(A, B, rep, thresholds)
    stage_annotation(A, B, rep, thresholds, clusters if clusters[0] is not None else None, hits)
    stage_genbank(A, B, rep)

    text = rep.render(f"Comparison: {A.label} vs {B.label}")
    Path(args.out).write_text(text)
    print(text)
    return 0 if all(v[1].startswith("PASS") for v in rep.verdicts) else 1


if __name__ == "__main__":
    sys.exit(main())
