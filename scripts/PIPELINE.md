# MGEfinder Pipeline

Five-step pipeline for detecting and characterizing mobile genetic elements (MGEs)
and insertion sequences (IS) from bacterial whole-genome sequencing data.

Each step is **idempotent** — it skips already-completed samples and only processes
what's still needed. All scripts (steps 2–5) accept a `BASE_DIR` argument pointing
to the batch directory (e.g., `/global/scratch/users/kh36969/mgefinder_batch1`).

## Directory Layout

```
~/mgefinder_project/
├── download_manager/          ← Python module (sample selection + download)
├── is_extractor/              ← Python module (IS extraction + annotation)
└── sample_registry.tsv        ← Central registry of all picked samples

mgefinder_batchN/
├── target_*.tsv               ← Sample selection target list (Step 1)
├── samples/
│   └── SAMEA1234567/
│       ├── 00.reads/          ← Step 1 output (FASTQ)
│       ├── 00.reference/      ← Step 1 output (reference genome)
│       ├── 00.genome/         ← Step 2 output (genome copy + BWA index)
│       ├── 00.assembly/       ← Step 2 output (SKESA assembly)
│       ├── 00.bam/            ← Step 2 output (aligned BAM)
│       ├── 01.find/ … 03.results/  ← Step 3 output (MGEfinder)
│       ├── is_extraction/     ← Step 4 output (IS element JSON)
│       └── is_circle_analysis/ ← Step 5 output (circularized IS TSV)
├── scripts/                   ← All pipeline scripts
├── jobs/                      ← SLURM log files
├── is_extraction_results/     ← Step 4 summary TSV
└── reference_genomes/         ← Shared reference genomes
```

---

## Step 1: Select & Download Samples

**What it does**: Selects a diverse subset of bacterial samples from the ENA master list
(~2.7M candidates), automatically excluding any samples already picked in previous
batches, then downloads their paired-end FASTQ reads and reference genomes.

**Scripts** (all run from `~/mgefinder_project/`):

| Subcommand | Purpose |
|-----------|---------|
| `match` | Find ENA samples with both reads + assembly |
| `build` | Merge metadata into a master list TSV |
| `subset` | Pick N samples (diversity-first), auto-exclude previous batches |
| `download` | Download FASTQ + reference for each selected sample |
| `register` | Retroactively register an existing batch into the registry |
| `status` | Show registry summary (how many samples per batch) |

### Sample Registry

A central registry at `~/mgefinder_project/sample_registry.tsv` tracks every sample
ever picked. The `subset` command reads it automatically so you never pick the same
sample twice across batches.

```bash
cd ~/mgefinder_project

# Check what's already registered
python -m download_manager.cli status

# Register an existing batch retroactively (from a target TSV or plain ID list)
python -m download_manager.cli register \
    --input /path/to/target_samples.tsv \
    --batch-name batch1
```

### Selecting Samples for a New Batch

The ENA master list is pre-built at:
`/global/home/groups/pc_rubinlab/databases/kuang/ena/master_list.tsv`

To rebuild it from fresh ENA data (only needed if the ENA database files are updated):
```bash
cd ~/mgefinder_project

python -m download_manager.cli match \
    --reads /global/home/groups/pc_rubinlab/databases/kuang/ena/ena_all_bacteria_reads.tsv \
    --assembly /global/home/groups/pc_rubinlab/databases/kuang/ena/ena_all_bacteria_assembly.tsv \
    --output matched_samples.txt

python -m download_manager.cli build \
    --matched matched_samples.txt \
    --reads /global/home/groups/pc_rubinlab/databases/kuang/ena/ena_all_bacteria_reads.tsv \
    --assembly /global/home/groups/pc_rubinlab/databases/kuang/ena/ena_all_bacteria_assembly.tsv \
    --output master_list.tsv
```

To select samples (auto-excludes previous batches via registry):
```bash
cd ~/mgefinder_project

# Pick 100K new samples for batch3
python -m download_manager.cli subset \
    --input /global/home/groups/pc_rubinlab/databases/kuang/ena/master_list.tsv \
    --output /global/scratch/users/kh36969/mgefinder_batch3/target_100k.tsv \
    --size 100000 \
    --batch-name batch3 \
    --seed 42

# Preview without registering
python -m download_manager.cli subset \
    --input /global/home/groups/pc_rubinlab/databases/kuang/ena/master_list.tsv \
    --output /tmp/preview.tsv \
    --size 100000 \
    --no-register
```

The selection strategy is diversity-first: one representative per species first,
then random fill to reach the target size.

### Downloading

```bash
cd ~/mgefinder_project

python -m download_manager.cli download \
    --input /global/scratch/users/kh36969/mgefinder_batch3/target_100k.tsv \
    --output /global/scratch/users/kh36969/mgefinder_batch3/samples
```

The download command skips samples that already have files on disk.

**What it produces**:
- `samples/{SAMPLE}/00.reads/{SAMPLE}_R1.fastq.gz`
- `samples/{SAMPLE}/00.reads/{SAMPLE}_R2.fastq.gz`
- `samples/{SAMPLE}/00.reference/genome.fna` (or `.fna.gz`)

**How to check status**:
```bash
# Count samples with reads
ls -1d samples/SAM*/00.reads/*_R1.fastq.gz 2>/dev/null | wc -l

# Count samples with reference genomes
ls -1d samples/SAM*/00.reference/genome.fna* 2>/dev/null | wc -l

# Check registry
python -m download_manager.cli status
```

---

## Step 2: Prepare Samples

**What it does**: For each sample, copies the reference genome, builds a BWA index,
assembles reads with SKESA, and aligns reads to the reference with BWA MEM.

**Script**: `scripts/prepare_samples.sh`

**Prerequisites**: Step 1 complete (reads + reference genome present).

**How to run**:
```bash
# Discover incomplete samples and submit SLURM array jobs
bash scripts/prepare_samples.sh /global/scratch/users/kh36969/mgefinder_batchN

# Dry run (see what would be submitted)
bash scripts/prepare_samples.sh /global/scratch/users/kh36969/mgefinder_batchN --dry-run

# Options: --max-jobs N, --force, --batch-size N
```

**What it produces** (per sample):
- `00.genome/genome.fna` + BWA index files
- `00.assembly/{SAMPLE}.fna` (SKESA contigs)
- `00.bam/{SAMPLE}.genome.bam` + `.bam.bai`

**How to check status**:
```bash
# Count fully prepared samples
ls -1 samples/SAM*/00.bam/*.genome.bam 2>/dev/null | wc -l

# Re-run to see what's still needed
bash scripts/prepare_samples.sh /global/scratch/users/kh36969/mgefinder_batchN --dry-run
```

---

## Step 3: Run MGEfinder

**What it does**: Runs the MGEfinder `denovo` workflow on each prepared sample to
identify mobile genetic elements and their insertion sites.

**Script**: `scripts/run_mgefinder.sh`

**Prerequisites**: Step 2 complete (genome, assembly, BAM all present).

**How to run**:
```bash
# Discover incomplete samples and submit SLURM array jobs
bash scripts/run_mgefinder.sh /global/scratch/users/kh36969/mgefinder_batchN

# Dry run
bash scripts/run_mgefinder.sh /global/scratch/users/kh36969/mgefinder_batchN --dry-run

# Single sample (no SLURM)
bash scripts/run_mgefinder.sh /global/scratch/users/kh36969/mgefinder_batchN --sample SAMEA1234567

# Options: --max-jobs N, --force, --batch-size N
```

**What it produces** (per sample):
- `01.find/` through `03.results/` directories
- `03.results/genome/04.makefasta.genome.all_seqs.fna` (IS sequences)
- `03.results/genome/NO_MGE_DETECTED.txt` (if no MGEs found)

**How to check status**:
```bash
# Count completed (with or without MGEs)
ls -1 samples/SAM*/mgefinder.log 2>/dev/null | wc -l

# Count samples with IS sequences found
find samples -name "04.makefasta.genome.all_seqs.fna" -not -empty 2>/dev/null | wc -l

# Re-run to see what's still needed
bash scripts/run_mgefinder.sh /global/scratch/users/kh36969/mgefinder_batchN --dry-run
```

---

## Step 4: Extract IS Elements

**What it does**: Extracts detailed IS element information from MGEfinder output:

1. **IS element sequence** — from `clusterseq.genome.tsv` (`inferred_seq` column)
2. **Flanking regions** (default 80 bp) — extracted from the **assembly contig** where the
   IS element sits, ensuring `upstream_flank → IS_element → downstream_flank` are
   contiguous on the same contig. Strand is auto-detected so flanking is always in the
   IS element's reading frame.
3. **ORF annotation** — Prodigal predicts coding regions (ORFs) and noncoding regions
   within each IS element sequence.
4. **Insertion site** — reference genome coordinates from `genotype.genome.tsv`.

IS elements detected only via `inferred_database` or `inferred_reference` (i.e., no
assembly contig mapping) will have empty flanking regions, since their full sequence
does not exist contiguously in the assembly.

Runs on a single node with 30 parallel workers (not a SLURM array — IS extraction is
lightweight per sample).

**Scripts**:
- `scripts/run_is_extraction.sh` — SLURM wrapper
- `scripts/run_is_extraction.py` — Python extraction logic
- `is_extractor/` module — `extractor.py` (extraction), `orf_finder.py` (ORF annotation)

**Prerequisites**: Step 2 (assembly) + Step 3 (MGEfinder results). Requires `module load bio/prodigal`.

**How to run**:
```bash
# Submit SLURM job (auto-discovers incomplete samples)
bash scripts/run_is_extraction.sh /global/scratch/users/kh36969/mgefinder_batchN

# Or run the Python script directly (e.g., inside an interactive job)
python scripts/run_is_extraction.py --samples-dir /path/to/samples

# With explicit sample list
python scripts/run_is_extraction.py --samples-dir /path/to/samples --sample-list my_list.txt

# Python options: --output-dir DIR, --workers N, --flank-size N
```

**What it produces**:
- Per-sample: `samples/{SAMPLE}/is_extraction/is_elements.json`
- Summary: `is_extraction_results/is_elements_summary.tsv`

**JSON output structure** (per IS element):
```
is_id, sample, seqid, cluster, group, pair_id, confidence
is_element:
  sequence, length, contig, start, end, strand, method
flanking_upstream:
  sequence, length, contig_start, contig_end, contig
flanking_downstream:
  sequence, length, contig_start, contig_end, contig
insertion_site:
  contig, pos_5p, pos_3p
orf_annotation:
  orfs: [{orf_id, start, end, strand, dna_sequence, protein_sequence, ...}]
  noncoding_regions: [{start, end, length, type, sequence}]
```

**How to check status**:
```bash
# Count completed samples
ls -1 samples/SAM*/is_extraction/is_elements.json 2>/dev/null | wc -l

# Check summary TSV row count
wc -l is_extraction_results/is_elements_summary.tsv
```

### Verifying IS Extraction Results

`scripts/verify_is_extraction.py` validates that the JSON output is self-consistent:

1. **Contiguity check** — `upstream_flank + IS_sequence + downstream_flank` matches
   the assembly contig at the recorded coordinates (handles both + and - strand).
2. **ORF tiling check** — ORFs + noncoding regions tile the IS sequence completely
   with no gaps, and each subsequence matches at its recorded coordinates.

```bash
cd ~/mgefinder_project

# Live test: run extractor on a sample and verify (no existing JSON needed)
python scripts/verify_is_extraction.py --live /path/to/samples/SAMEA1234567 -v

# Verify existing JSON for a sample
python scripts/verify_is_extraction.py /path/to/samples/SAMEA1234567

# Batch verify all samples with JSON files
python scripts/verify_is_extraction.py --batch /path/to/samples --max-samples 50
```

---

## Step 5: Detect Circularized IS Elements

**What it does**: For each sample with IS sequences and FASTQ reads, aligns reads to
an IS-only reference, runs Circle-Map to detect circular DNA, and validates whether
IS elements are self-circularized.

**Script**: `scripts/detect_circularized_is.sh`

**Prerequisites**: Step 3 complete (IS FASTA) + FASTQ reads still available.
Requires `mgefinder` conda environment (has bwa, samtools, Circle-Map).

**How to run**:
```bash
# Discover incomplete samples and submit SLURM array jobs
bash scripts/detect_circularized_is.sh /global/scratch/users/kh36969/mgefinder_batchN

# Dry run
bash scripts/detect_circularized_is.sh /global/scratch/users/kh36969/mgefinder_batchN --dry-run

# Force re-run all
bash scripts/detect_circularized_is.sh /global/scratch/users/kh36969/mgefinder_batchN --force

# Options: --max-jobs N (default 50), --batch-size N (default 1000)
```

**What it produces** (per sample):
- `is_circle_analysis/circularized_is.tsv` (validated results)
- `is_circle_analysis/is_circles.bed` (raw Circle-Map output)
- `is_circle_analysis/is_reference.fna` (IS-only reference used)
- `is_circle_analysis/status.txt` (machine-readable status)

**How to check status**:
```bash
# Count completed samples
ls -1 samples/SAM*/is_circle_analysis/circularized_is.tsv 2>/dev/null | wc -l

# Count samples with confirmed circular IS
grep -l "CONFIRMED_CIRCULAR_IS" samples/SAM*/is_circle_analysis/circularized_is.tsv 2>/dev/null | wc -l

# Re-run to see what's still needed
bash scripts/detect_circularized_is.sh /global/scratch/users/kh36969/mgefinder_batchN --dry-run
```

---

## Running on a New Batch

Full example for batch3 with 100K samples:

```bash
cd ~/mgefinder_project
BATCH=/global/scratch/users/kh36969/mgefinder_batch3

# Step 1a: Select samples (auto-excludes batch1 + batch2 via registry)
python -m download_manager.cli subset \
    --input /global/home/groups/pc_rubinlab/databases/kuang/ena/master_list.tsv \
    --output "$BATCH/target_100k.tsv" \
    --size 100000 \
    --batch-name batch3 \
    --seed 42

# Step 1b: Download
python -m download_manager.cli download \
    --input "$BATCH/target_100k.tsv" \
    --output "$BATCH/samples"

# Copy scripts from an existing batch
cp -r /global/scratch/users/kh36969/mgefinder_batch1/scripts/ "$BATCH/scripts/"

# Step 2: Prepare (assembly + alignment)
bash "$BATCH/scripts/prepare_samples.sh" "$BATCH"

# Step 3: MGEfinder
bash "$BATCH/scripts/run_mgefinder.sh" "$BATCH"

# Steps 4 & 5 can run in parallel:
bash "$BATCH/scripts/run_is_extraction.sh" "$BATCH"            # Step 4
bash "$BATCH/scripts/detect_circularized_is.sh" "$BATCH"       # Step 5
```

Every step is idempotent — if jobs fail, rerun the same command to process
only the incomplete samples.

## SLURM Settings

All scripts use:
- **Partition**: lr6
- **Account**: pc_rubinlab
- **QoS**: lr_normal

Per-step resources:
| Step | CPUs | Memory | Time | Pattern |
|------|------|--------|------|---------|
| 2. prepare_samples | 32 | full node | 24h | Array job (1 sample/task) |
| 3. run_mgefinder | 8 | 32 GB | 24h | Array job (1 sample/task) |
| 4. run_is_extraction | 32 (exclusive) | full node | 6h | Single node, 30 Python workers |
| 5. detect_circularized_is | 8 | default | 4h | Array job (1 sample/task) |

---

## Archiving a Batch

**What it does**: Creates a compressed tar.gz of the batch directory, excluding
large regenerable files (FASTQ reads, BAM alignments, BWA indexes). This preserves
all results (MGEfinder output, IS extractions, circle analysis) at a fraction of the
original size.

**Script**: `scripts/archive_batch.sh`

**How to run**:
```bash
# Archive with default output path ({parent_dir}/mgefinder_batchN_results.tar.gz)
bash scripts/archive_batch.sh /global/scratch/users/kh36969/mgefinder_batch1

# Custom output path
bash scripts/archive_batch.sh /global/scratch/users/kh36969/mgefinder_batch1 /path/to/batch1_results.tar.gz
```

**Excluded file types**:
- `*.fastq.gz` — raw reads (~1-3 GB per sample)
- `*.bam`, `*.bam.bai` — alignments
- `*.bwt`, `*.pac`, `*.ann`, `*.amb`, `*.sa` — BWA indexes (regenerable)

**To extract**:
```bash
tar xzf mgefinder_batch1_results.tar.gz -C /path/to/destination
```
