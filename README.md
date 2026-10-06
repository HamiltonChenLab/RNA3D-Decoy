# RNA3D-Decoy

RNA3D-Decoy provides target split lists and utilities for preparing RNA reference structures and evaluating predicted structures (decoys). The scoring scripts compare predictions with experimental reference structures using RMSD, lDDT, and, for RNA–protein complexes, DockQ.

The typical workflow is to obtain reference structures, generate predictions with an external prediction method, organize the structure files, and run the scoring scripts. Prediction generation, pretrained models, and pre-generated decoy coordinates are not included in this checkout.

## Repository contents

| File or directory | Purpose |
| --- | --- |
| `download_rna.py` | Download and filter experimental RNA structures from RCSB PDB. |
| `prediction_rmsd.py` | Compute fitted global RMSD and residue-level deviations. |
| `prediction_lddt.py` | Compute global and residue-level lDDT. |
| `prediction_dockq.py` | Compute interface quality scores for RNA–protein complexes. |
| `splits_result/train.txt` | Training target list. |
| `splits_result/val.txt` | Validation target list. |
| `splits_result/test.txt` | Test target list. |

**Current setup requirements:** the `level1`, `level2`, and `level3` download modes depend on helper modules that are not included in this checkout. The `strict` download mode is implemented in `download_rna.py` itself. The lDDT script also requires a path change before use outside the original cluster. Both are explained below.

## Installation

Clone the repository and run the examples from its root directory:

```bash
git clone https://github.com/HamiltonChenLab/RNA3D-Decoy.git
cd RNA3D-Decoy
```

Use Python 3.10 or newer. RMSD and lDDT require **OpenStructure**, including its `ost` Python module and compound library. DockQ scoring requires **DockQ v2** and **Biopython**. The strict downloader uses only the Python standard library and requires internet access.

For a platform supported by the OpenStructure Conda package, create an environment:

```bash
conda create -n rna3d-decoy -c conda-forge -c bioconda python=3.11 openstructure pip
conda activate rna3d-decoy
python -m pip install "DockQ>=2,<3" biopython
```

See the official [OpenStructure installation instructions](https://openstructure.org/install) for Conda and container options, and the [DockQ installation instructions](https://github.com/wallnerlab/DockQ) for its supported setup. If the Conda package is unavailable for your platform, use one of OpenStructure's documented container environments.

Check that the active environment can load the dependencies:

```bash
python -c "from ost import io; from ost.mol import alg; from ost.mol.alg.lddt import lDDTScorer; from Bio import PDB; print('Imports OK')"
DockQ --help
```

The repository does not currently include a pinned environment or dependency lock file.

## Dataset levels and splits

The downloader describes three target categories:

| Level | Target type |
| --- | --- |
| L1 | Clean single-chain RNA without metal ions. |
| L2 | Clean single-chain RNA with metal ions. |
| L3 | One RNA polymer chain and one protein polymer chain. |

The supplied split files contain these numbers of target entries:

| Split | L1 | L2 | L3 | Total |
| --- | ---: | ---: | ---: | ---: |
| Training | 338 | 93 | 115 | 546 |
| Validation | 42 | 12 | 14 | 68 |
| Test | 42 | 12 | 15 | 69 |
| Total | 422 | 117 | 144 | 683 |

Each line contains a level label and a target identifier, for example:

```text
L1. 1A60
L2. 1KXK
L3. 1A1T_complex
```

These are target lists, not coordinate files. Keep the decoys belonging to each target in its assigned split when constructing a training or evaluation dataset. The download and scoring scripts do **not** automatically read these split files or restrict processing to them; select the desired inputs yourself. A fresh RCSB query is not guaranteed to reproduce the supplied lists.

## Download reference structures

### Strict mode available in this checkout

The downloader is configured through environment variables, not command-line arguments. Always set `RCSB_RNA_OUTDIR` to a writable path: its default points to the original author's filesystem.

For a small download of up to five accepted targets:

```bash
RCSB_RNA_MODE=strict \
RCSB_RNA_OUTDIR="$PWD/data/strict" \
RCSB_RNA_END=2026-06-23 \
RCSB_RNA_LIMIT=5 \
python download_rna.py
```

Set `RCSB_RNA_LIMIT=0` to remove the accepted-target limit. Strict mode searches experimental RNA-only entries, excludes nonpolymer and branched entities, and requires a single canonical A/C/G/U RNA chain with coordinates covering all sequence positions. It also checks agreement between the FASTA and mmCIF sequences.

| Variable | Meaning and default |
| --- | --- |
| `RCSB_RNA_MODE` | `strict`, `level1`, `level2`, or `level3`; default: `level1`. |
| `RCSB_RNA_OUTDIR` | Output directory; explicitly override the machine-specific default. |
| `RCSB_RNA_START` | Earliest initial release date; default: `1990-01-01`. |
| `RCSB_RNA_END` | Latest initial release date; defaults: `2030-01-01` for strict/level1, `2026-06-23` for level2, `2026-06-26` for level3. |
| `RCSB_RNA_MAX_RESOLUTION` | Maximum resolution in Å; default: `4.0` for strict/level1, disabled for level2/level3. Use `none` to disable. |
| `RCSB_RNA_LIMIT` | In strict mode, maximum accepted targets; default: `0` (unlimited). Passed to the helper pipeline in other modes. |
| `RCSB_RNA_TIMEOUT` | Network timeout in seconds; default: `45`. |
| `RCSB_RNA_MAIN_MIN`, `RCSB_RNA_MAIN_MAX` | Strict-mode main length range; default: 30–500 nucleotides, inclusive. |
| `RCSB_RNA_EXTENSION_MIN`, `RCSB_RNA_EXTENSION_MAX` | Strict-mode extension length range; default: 500–1000 nucleotides, inclusive. |

Strict mode creates:

```text
data/strict/
├── all/
│   ├── fasta_rna/       # RNA FASTA sequences
│   ├── cif_rna/         # Filtered RNA/allowed-metal atom records
│   ├── cif_native/      # Full downloaded mmCIF structures
│   └── af3_inputs/      # AlphaFold 3 input JSON files
├── main_30_500/         # Same subdirectories for the main length range
├── extension_500_1000/  # Same subdirectories for the extension range
├── pdb_rna_only.txt
├── pdb_rna_main_30_500.txt
├── pdb_rna_extension_500_1000.txt
├── selected_rna_sequences.tsv
└── skipped_rna_sequences.tsv
```

Files use the PDB ID as their basename, such as `1A60.cif`. A 500-nucleotide target belongs to both default length ranges. The `cif_rna` writer can retain allowed metal atoms, but the strict selection rules exclude nonpolymer entities. The JSON files prepare prediction inputs; this script does not run AlphaFold 3.

To score predictions against the flat strict reference directory, use `--native-dir "$PWD/data/strict/all/cif_native" --native-level 1` with RMSD or lDDT. This assigns a level for matching; it does not establish membership in the supplied L1 split.

### Level 1/2/3 modes require additional modules

The following imports must be available beside `download_rna.py` or on `PYTHONPATH`:

| Mode | Required module, currently absent |
| --- | --- |
| `level1`, `level2` | `build_single_chain_rna_datasets.py` |
| `level3` | `build_single_rna_single_protein_targets.py` |

Without these modules, the default `python download_rna.py` command fails with `ModuleNotFoundError`. Obtain the matching helper implementations before using these modes.

Both `level1` and `level2` call the same helper pipeline; the wrapper advertises outputs for both levels. The mode changes defaults such as the date and resolution filters rather than selecting a separate level-specific pipeline. The wrapper exposes `RCSB_RNA_MAX_LENGTH` (default 1000) for this pipeline.

For level 3, the wrapper exposes `RCSB_RNA_MAX_RNA_LENGTH` (default 1000), `RCSB_RNA_CONTACT_CUTOFF` (default 6.0 Å), `RCSB_RNA_REQUIRE_CONTACT` (default true), and `RCSB_RNA_ALLOW_ORGANIC_LIGAND` (default false). Its complete filtering and output behavior depends on the missing helper module.

## Organize prediction and reference files

A **native** is the reference structure; a **prediction** is the decoy to evaluate. The prepared dataset on **h200** is stored at:

```text
/home/schen3/CASP17/RNA3D_decoys/
```

The dataset is organized **by level first, then by prediction method**. The directory layout and examples below were checked against the server files.

```text
RNA3D_decoys/
├── README.md
├── level1/
│   ├── native/
│   ├── fasta/
│   ├── labels/
│   ├── af3/
│   ├── boltz/
│   ├── chai/
│   ├── esmfold2/
│   └── protenix/
├── level2/             # Same eight subdirectories
└── level3/             # Same eight subdirectories
```

| Level directory | Native structures |
| --- | ---: |
| `level1/native/` | 422 |
| `level2/native/` | 117 |
| `level3/native/` | 144 |

### Reference structures and sequences

Native structures use the four-character PDB ID as their filename at every level. In particular, level-3 native filenames do **not** have the `_complex` suffix used by their predictions.

```text
level1/native/17RA.cif
level2/native/1JTW.cif
level3/native/11AO.cif
```

Level-1 and level-2 FASTA files are stored directly under `fasta/`. Level 3 separates RNA, protein, and combined complex sequences:

```text
level1/fasta/17RA.fasta
level2/fasta/1JTW.fasta
level3/fasta/rna/11AO_rna.fasta
level3/fasta/protein/11AO_protein.fasta
level3/fasta/complex/11AO_complex.fasta
```

### Prediction structures

Each level contains outputs from `af3`, `boltz`, `chai`, `esmfold2`, and `protenix`. These are stored prediction files; the scoring scripts do not run the prediction methods.

**AF3** uses a target directory named `<id>_001_<id>_dataset` for levels 1 and 2, and `<id>_complex` for level 3. Sampled predictions are nested by seed and sample:

```text
level1/af3/17RA_001_17RA_dataset/
├── 17RA_001_17RA_dataset_model.cif
├── 17RA_001_17RA_dataset_ranking_scores.csv
├── seed-1_sample-0/17RA_AF3_seed1_model0.cif
└── ... other seed/sample directories and confidence JSON files

level2/af3/1JTW_001_1JTW_dataset/seed-1_sample-0/1JTW_AF3_seed1_model0.cif
level3/af3/11AO_complex/seed-1_sample-0/11AO_complex_AF3_seed1_model0.cif
```

AF3 has 30 sampled CIFs per target (six seeds and five samples), plus a top-level `*_model.cif`. A recursive scan includes that additional file too; select the `seed-*/` outputs when you want just the 30 sampled decoys.

The other methods use the following patterns. Here, `<label>` is the PDB ID for levels 1/2 and `<id>_complex` for level 3; all paths are relative to `level<N>/`.

| Method | Structure path pattern |
| --- | --- |
| Boltz | `boltz/<label>/<label>_boltz_<recycles>_<sampling_steps>_<seed>_<sample>.cif` |
| Chai | `chai/<label>/<label>_chai_<recycles>_<diffusion_steps>_<seed>_<sample>.cif` |
| ESMFold2 | `esmfold2/<label>_ESM_<model_num>_<loops>_<sampling_steps>.cif` |
| Protenix | `protenix/<label>/seed_<seed>/predictions/<label>_sample_<sample>.cif` |

For example:

```text
level1/boltz/17RA/17RA_boltz_3_150_1_0.cif
level1/chai/17RA/17RA_chai_3_150_1_0.cif
level1/esmfold2/17RA_ESM_2_15_100.cif
level1/protenix/17RA/seed_1/predictions/17RA_sample_0.cif
level3/boltz/11AO_complex/11AO_complex_boltz_3_150_1_0.cif
```

ESMFold2 structures are directly inside the method directory; the other methods use target subdirectories. AF3 and Protenix also include confidence JSON files, which the structure scorers ignore.

### Existing score labels

Each level contains these five CSV files, where `<N>` is the level number:

```text
level<N>/labels/
├── prediction_tmscores_level<N>.csv
├── prediction_lddt_summary_level<N>.csv
├── prediction_lddt_local_level<N>.csv
├── prediction_rmsd_summary_level<N>.csv
└── prediction_rmsd_local_residue_level<N>.csv
```

The summary files contain prediction-level metric rows, with separate atom selections where applicable; the local files contain residue-level rows. The CSVs include `Model`, `level`, `target_id`, `native_path`, and `prediction_path`.

Existing label rows retain paths from the original `/hpc/netapp/casp17/CASP17/...` calculation layout. Resolve those entries against the level-first layout above when accessing the packaged structures. TM-score labels are provided on the server, although this repository does not include their scoring script. No DockQ CSV is currently present in these `labels/` directories.

### Use this layout with the current scoring scripts

The server dataset and the scripts' batch discovery layouts differ. Simply setting `--prediction-root` to `RNA3D_decoys/` does **not** make the current batch scorers discover `level<N>/<model>/`.

For a single structure pair, RMSD and DockQ accept the actual server paths directly. For example, run this on h200 from the code repository:

```bash
data_root=/home/schen3/CASP17/RNA3D_decoys
python prediction_rmsd.py \
  --native "$data_root/level1/native/17RA.cif" \
  --prediction "$data_root/level1/af3/17RA_001_17RA_dataset/seed-1_sample-0/17RA_AF3_seed1_model0.cif" \
  --model-name af3 \
  --output-dir results/RMSD
```

For batch scoring, the current scripts require a separate working layout with the prediction method first, or a change to their discovery code. The examples in the following sections use `data/native_cif/level<N>/` for references and `generated_decoys/<model>/...` for predictions. Populate those working directories with the structures you want to score; they are not the server dataset's existing directory names.

For levels 1 and 3, batch prediction paths are `<prediction-root>/<model>/level1/` and `<prediction-root>/<model>/level3/`. RMSD and lDDT require these special level-2 paths:

| Model | Level-2 path below the batch prediction root |
| --- | --- |
| `af3` | `af3/level2/fold_output_with_ions/` |
| `boltz` | `boltz/level2_with_ions/` |
| `chai` | `chai/level2_with_ions/` |
| `esmfold2` | `esmfold2/level2_with_grouped_ions/` |
| `protenix` | `protenix/level2_with_ions/` |

These are discovery conventions in the code, not subdirectories of `RNA3D_decoys/level2/`. Ordinary `level2/` folders for the four non-AF3 methods are skipped. The server README and TM-score labels describe level 2 as `rna_multimer`, while the downloader wrapper describes it as single-chain RNA with metals; the existing terminology is inconsistent.

Additional matching rules:

- RMSD/lDDT match by `(level, target_id)`. For a flat server reference directory such as `RNA3D_decoys/level2/native/`, pass both `--native-dir /home/schen3/CASP17/RNA3D_decoys/level2/native` and `--native-level 2` in batch mode.
- PDB IDs are read from the start of structure filenames by RMSD/lDDT. The stored filenames already follow this convention, including `11AO_complex_...`, which maps to native `11AO.cif`.
- CASP-style IDs such as `R1107`, `H1107`, or `T1107` can be inferred from the full path.
- Supported suffixes are `.pdb`, `.cif`, `.mmcif`, `.ent`, `.pdb.gz`, `.cif.gz`, and `.mmcif.gz`; searches are recursive, but the scripts' directory walks do not follow nested directory symlinks.
- If multiple references match the same target and level, the prediction is scored against each reference.
- lDDT still requires the path configuration described in its section below; it has no direct-pair mode.

## Compute RMSD

For one reference/prediction pair:

```bash
python prediction_rmsd.py \
  --native data/native_cif/level1/1A60.cif \
  --prediction generated_decoys/af3/level1/1A60_model_0.cif \
  --model-name af3 \
  --output-dir results/RMSD
```

For a batch:

```bash
python prediction_rmsd.py \
  --models af3 boltz \
  --levels 1 2 3 \
  --native-dir "$PWD/data/native_cif" \
  --prediction-root "$PWD/generated_decoys" \
  --output-dir "$PWD/results/RMSD"
```

The output directory is created automatically. The two CSV files are:

- `prediction_rmsd_summary.csv`: global scores, matching method, matched atom/residue counts, input paths, and errors.
- `prediction_rmsd_local_residue.csv`: per-residue deviations and the corresponding prediction/reference residue identifiers.

RMSD is reported in Å; lower values indicate closer agreement. Each atom selection receives a global least-squares fit. The local values describe residue deviations **after that shared fit**, not independently fitted residue scores.

Useful options:

| Option | Behavior |
| --- | --- |
| `--match number` | Match residues by number; the default. |
| `--match index` | Match residues by order. Use only when corresponding residue order is known to agree. |
| `--local-atom-selection all_atom` | Write local all-atom rows; the default. |
| `--local-atom-selection all` | Write local rows for every computed selection. |
| `--local-atom-selection hybrid_C4_CA` | Write local rows from one shared RNA C4′/protein Cα fit. |
| `--summary-csv PATH`, `--local-csv PATH` | Override individual output filenames. |

For levels other than 2, the script attempts `all_atom`, `hybrid_C4_CA`, `C4_prime`, `P_atom`, and `CA_atom`, skipping selections with no atoms. At level 2, it computes only RNA `all_atom`, `C4_prime`, and `P_atom` selections over residue names A/C/G/I/U, excluding metal ions.

In direct mode, level is inferred from the prediction path. A level-2 prediction therefore needs a `level2` path component to activate the level-2 RNA-only selections; `--levels` filters tasks and does not assign a level.

## Compute lDDT

**Configure the root first.** The current `prediction_lddt.py` contains:

```python
ROOT = Path("/hpc/netapp/casp17/CASP17").resolve()
```

For the repository layout above, replace that line in your local copy with:

```python
ROOT = Path(__file__).resolve().parent
```

Alternatively, set it to an absolute directory containing `generated_decoys/`. The script reads predictions from `ROOT/generated_decoys/` and writes results to `ROOT/lDDT/`. It currently has no `--prediction-root`, `--output-dir`, or direct `--native/--prediction` options, and it does not use the `CASP17_ROOT` environment variable.

After configuring that path:

```bash
python prediction_lddt.py \
  --models af3 boltz \
  --levels 1 2 3 \
  --native-dir "$PWD/data/native_cif"
```

For a flat native directory, add `--native-level N`. Outputs are:

- `lDDT/prediction_lddt_summary.csv`: global lDDT, mean defined local lDDT, score/atom counts, input paths, and errors.
- `lDDT/prediction_lddt_local.csv`: residue-level scores for each computed atom selection.

The script attempts `all_atom`, `hybrid_C4_CA`, `C4_prime`, `P_atom`, and `CA_atom` at every level. Unlike the RMSD script, it has no special RNA-only filter for level 2.

lDDT measures preservation of local interatomic distances without a global structural fit; higher values indicate better agreement. The wrapper writes OpenStructure's scores directly, on a 0–1 scale, and averages defined local values separately. It uses automatic chain matching and the scorer's default residue-number mapping, with residue-name checking disabled. Ensure residue numbering corresponds between each prediction and reference; chain matching does not establish sequence alignment. See the [OpenStructure lDDT API](https://openstructure.org/docs/2.12/mol/alg/lddt/) for scoring details.

## Compute DockQ for RNA–protein complexes

For one complex pair, create the CSV's parent directory before running the script:

```bash
mkdir -p results/DockQ
python prediction_dockq.py \
  --native data/native_cif/level3/1A1T_complex.cif \
  --prediction generated_decoys/af3/level3/1A1T_complex_model_0.cif \
  --model-name af3 \
  --output-csv results/DockQ/prediction_dockq.csv \
  --dockq-output-dir results/DockQ/raw
```

For a batch:

```bash
mkdir -p results/DockQ
python prediction_dockq.py \
  --models af3 boltz \
  --level 3 \
  --native-dir "$PWD/data/native_cif/level3" \
  --prediction-root "$PWD/generated_decoys" \
  --workers 4 \
  --dockq-n-cpu 1 \
  --output-csv "$PWD/results/DockQ/prediction_dockq.csv" \
  --dockq-output-dir "$PWD/results/DockQ/raw"
```

DockQ scans `<prediction-root>/<model>/level3/` by default and matches references by target ID. Unlike RMSD/lDDT, its `--models` option accepts arbitrary folder names.

By default, the wrapper reads the first structure model, selects the longest recognized protein chain and the longest recognized RNA chain, renumbers their retained polymer residues, and writes temporary structures with protein chain A and RNA chain B. It then uses fixed mapping `AB:AB`. This is intended for one-RNA/one-protein targets; extra chains are not included in the normalized inputs.

Outputs include `prediction_dockq.csv` and raw JSON files under `raw/<model>/level3/<target>/`. The CSV records DockQ, available global/interface metrics, chain mapping, input paths, and errors. Higher DockQ values indicate better interface agreement. For structures with multiple interfaces passed through raw mode, the detailed interface columns come from the highest-DockQ interface; inspect the JSON for all interfaces.

Useful options:

- `--max-tasks 5`: process a small number of pairs.
- `--force`: recompute existing results, including cached JSON.
- `--dockq-bin /path/to/DockQ`: choose a specific executable.
- `--raw-dockq-inputs`: pass the original structures to DockQ.
- `--raw-dockq-inputs --mapping AB:HL`: map original prediction chains A/B to native chains H/L.
- `--no-align`: use residue numbering instead of sequence alignment. With default normalization, this refers to the newly assigned residue numbers.
- `--shard-index 0 --shard-count 4`: process one of four task shards; run indices 0–3 to cover the full task list.

`--workers` controls concurrent structure pairs; `--dockq-n-cpu` controls CPUs requested by each DockQ call. The wrapper uses Unix file locking and should be run in a Unix-like environment.

## Re-running and troubleshooting

All summary CSVs include an `error` column. Inspect it even if a script finishes normally, since individual scoring failures are recorded while processing continues.

| Symptom | What to check |
| --- | --- |
| `No module named 'ost'` | Activate an OpenStructure environment. RMSD/lDDT import it before parsing arguments, so even `--help` needs this dependency. |
| Missing `build_single_...` module | Obtain the helper modules for level1/2/3 downloading, or explicitly use strict mode. |
| Zero tasks or skipped directories | Check the prediction root, model folders, level filters, and special level-2 paths. |
| `No matching native found` | Check filename target IDs, `levelN` directories, and `--native-level` for flat native directories. |
| DockQ cannot open its output CSV | Create the parent directory of `--output-csv` before running. |
| DockQ executable not found | Activate the installation environment or set `--dockq-bin`. |
| Unexpected residue matches or missing scores | Inspect chain content, residue numbering, atom names, and the summary's matched atom/residue counts. |

RMSD and lDDT merge new CSV rows by their identifying fields and preserve existing successful rows with the same keys. To recompute after changing coordinates or settings, use a fresh output location (or move the earlier outputs aside). Run jobs writing the same RMSD/lDDT CSVs sequentially because their merge operations do not use file locks.

DockQ skips existing successful CSV entries, can reuse cached JSON, and locks CSV merges for concurrent workers or shards. Use `--force` to refresh existing scores.

RMSD and DockQ default to the **parent of the repository directory** as `ROOT`, unless `CASP17_ROOT` is set. Their inherited native/prediction paths refer to the original project layout. The explicit paths in the examples avoid those assumptions.

For all available scoring options, after installing the required dependencies:

```bash
python prediction_rmsd.py --help
python prediction_lddt.py --help
python prediction_dockq.py --help
```
