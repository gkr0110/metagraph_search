# MetaGraph Sequence Search — Galaxy Tool Wrapper

A Galaxy tool that submits FASTA sequences to the
[MetaGraph Online](https://metagraph.ethz.ch) service at ETH Zurich and returns
matched accessions as a tab-separated table.  An internet connection from the
Galaxy server to `metagraph.ethz.ch` is required.

---

## What is MetaGraph?

MetaGraph is a petabase-scale full-text search engine for biological sequences.
It indexes raw sequencing data from NCBI SRA, ENA, and DDBJ using compact
annotated de Bruijn graphs, enabling rapid query of hundreds of millions of
samples.  Supported query types are nucleotide (DNA/RNA) and amino acid
(protein).

---

## Repository layout

```
metagraph_search/
├── metagraph_search.xml          # Galaxy tool definition (Planemo)
├── metagraph_search.py           # Python wrapper script
├── .shed.yml                     # Galaxy ToolShed metadata
├── conda_env.yml                 # Conda environment for local development
└── test-data/
    ├── test_input.fa             # Short E. coli 16S rRNA test sequence
    └── expected_output.tabular   # Expected TSV for planemo test
```

---

## Installation

### On a Galaxy instance (usegalaxy.ch or self-hosted)

Install from the [Galaxy ToolShed](https://toolshed.g2.bx.psu.edu/) by
searching for **metagraph_search** (owner: `govind`) and requesting
installation through the Galaxy admin panel.

### Dependencies

The tool requires Python ≥ 3.9 and the `requests` library (≥ 2.31).  When
installed via ToolShed, Conda resolves these automatically from the
`conda-forge` and `bioconda` channels.

---

## Usage

### Inputs

| Parameter | Description |
|-----------|-------------|
| **Input sequences (FASTA)** | DNA, RNA, or protein sequences. Standard limit: 10 sequences / 55,000 characters. |
| **Molecule type** | *Nucleotide* for DNA/RNA; *Amino acid* for protein (UniParc only). |
| **Database(s)** | One or more databases (see table below). Multi-select is supported. |
| **Search mode** | *Exact match* (fast, k-mer based) or *Alignment* (sensitive, seed–chain–extend). |
| **Discovery threshold** | Minimum fraction of query k-mers that must match a label (0.0–1.0; default 0.5). |
| **Max matched accessions** | Hard cap on returned labels per database per query (default 100). |
| **Alignment options** | Min seed coverage and max alternative alignments (alignment mode only). |

### Available databases

| Database ID | Description |
|-------------|-------------|
| `refseq33m` | RefSeq — NCBI Reference Sequences (32.9 M accessions; DNA · Align · Coordinates) |
| `atb` | AllTheBacteria — bacterial WGS assemblies (2.76 M genomes; DNA · Coordinates) |
| `gnomad` | GnomAD — human genomic variation (29 chromosomal segments; DNA · Align) |
| `uhgg` | UHGG Catalog — Unified Human Gut Genome catalogue v1 (4,644 genomes; DNA · Align) |
| `uhgg_all` | UHGG All — all UHGG non-redundant genomes (287 K genomes; DNA · Align) |
| `sra-metagut` | SRA MetaGut — human gut metagenome (241 K samples; DNA/RNA · Align) |
| `tara-oceans` | Tara Oceans Genomes — ocean metagenomes with coordinates (34 K genomes; DNA · Align · Coordinates) |
| `tara-assemblies` | Tara Oceans Scaffolds — assembled ocean scaffolds (318 M scaffolds; DNA · Align) |
| `metasub41` | MetaSUB k=41 — urban microbiome (4,220 city samples; DNA · Align) |
| `metasub19` | MetaSUB k=19 — urban microbiome (4,220 city samples; DNA · Align) |
| `sra-microbe` | SRA Microbe — microbial WGS (446 K samples; DNA/RNA · Align) |
| `sra-fungi` | SRA Fungi — fungal sequences (122 K samples; DNA/RNA · Align) |
| `sra-plants` | SRA Plants — plant sequences (532 K samples; DNA/RNA · Align) |
| `sra-metazoa` | SRA Metazoa — animal sequences excl. human (805 K samples; DNA/RNA · Align) |
| `sra-human` | SRA Human — human raw reads (436 K samples; DNA/RNA · Align) |
| `sra-mus_muculus` | SRA Mus musculus — mouse raw reads (58 K samples; DNA/RNA · Align) |
| `sra-logan-chunks` | SRA Logan Contigs — pre-assembled SRA contigs (21.4 M samples; DNA/RNA) |
| `uniparc` | UniParc — UniProt Archive protein sequences (541 M proteins; Amino Acids · Coordinates) |

### Output

A tab-separated file with one row per matched accession:

| Column | Description |
|--------|-------------|
| `query` | FASTA sequence name (first word of the `>` header) |
| `database` | MetaGraph database ID searched |
| `label` | Matched accession ID (e.g. `NC_000913.3`, `UPI0000129DB5`) |
| `normalized_score` | Fraction of query k-mers matched (0.0–1.0) |
| `kmer_count` | Raw count of matching k-mers |
| `coordinates` | K-mer coordinate ranges in the matched sequence (if available) |
| `organism` | Organism name from enrichment metadata |
| `taxid` | NCBI Taxonomy ID of the matched organism |

---

## API details

The tool wraps the MetaGraph async REST API at `https://metagraph.ethz.ch:8081`:

1. **POST** `/search` — submit query, receive `search_id`
2. **GET** `/search/{id}/status` — poll every 10 s until status is `"done"`
3. **GET** `/search/{id}/results` — download and flatten hit list to TSV

The tool polls for up to 1 hour per sequence; large multi-database jobs may
take several minutes.  Results on the MetaGraph server expire after 48 hours
but are downloaded immediately upon completion.

Full API documentation: https://metagraph.ethz.ch/static/docs/api.html

---

## Local development with Planemo

### Set up the environment

```bash
conda env create -f conda_env.yml
conda activate metagraph-galaxy
```

### Lint the tool

```bash
cd metagraph_search
planemo lint metagraph_search.xml
```

### Run tests (mock mode — no network required)

```bash
planemo test metagraph_search.xml
```

The test injects `--mock` via a hidden XML parameter so the Python script
returns a built-in response instead of calling the live API.

### Serve locally in Galaxy

```bash
planemo serve metagraph_search.xml
# Galaxy UI available at http://localhost:9090
```

### Publish to the Galaxy ToolShed

Configure `~/.planemo.yml` with your ToolShed API key:

```yaml
sheds:
  toolshed:
    key: YOUR_TOOLSHED_API_KEY
```

Then create or update the repository:

```bash
# First time
planemo shed_create --shed_target toolshed metagraph_search/

# Subsequent updates
planemo shed_update --shed_target toolshed metagraph_search/
```

The `--shed_target toolshed` flag targets `https://toolshed.g2.bx.psu.edu/`
(the main Galaxy ToolShed).  Use `--shed_target test_toolshed` for
`https://testtoolshed.g2.bx.psu.edu/` during testing.

---

## Citation

If you use MetaGraph in your research, please cite:

> Karasikov et al. (2023) *Indexing and analysing nucleotide archives at
> petabase scale.* Genome Biology.
> https://doi.org/10.1186/s13059-023-02958-x

---

## Links

- MetaGraph Online: https://metagraph.ethz.ch
- API documentation: https://metagraph.ethz.ch/static/docs/api.html
- Galaxy ToolShed entry: https://toolshed.g2.bx.psu.edu/view/govind/metagraph_search
- usegalaxy.ch: https://usegalaxy.ch
