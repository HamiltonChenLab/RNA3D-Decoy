#!/usr/bin/env python3

from pathlib import Path
import argparse
import concurrent.futures
import copy
import csv
import fcntl
import gzip
import io
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import threading


ROOT = Path(os.environ.get("CASP17_ROOT", Path(__file__).resolve().parents[1])).resolve()

NATIVE_DIR = Path(
    os.environ.get(
        "CASP17_DOCKQ_NATIVE_DIR",
        ROOT / "DProQA_RNA" / "data" / "rna_protein" / "native_cif",
    )
).resolve()
GENERATED_DECOYS = Path(
    os.environ.get("CASP17_GENERATED_DECOYS", ROOT / "generated_decoys")
).resolve()

OUTPUT_DIR = ROOT / "DockQ"
CSV_OUTPUT = OUTPUT_DIR / "prediction_dockq.csv"
DOCKQ_OUTPUT_DIR = OUTPUT_DIR / "prediction_dockq_outputs"

STRUCTURE_SUFFIXES = (".pdb", ".cif", ".mmcif", ".ent", ".pdb.gz", ".cif.gz", ".mmcif.gz")

CASP_TARGET_RE = re.compile(r"([RHT]\d{4}(?:v\d+)?)", re.IGNORECASE)
PDB_TARGET_RE = re.compile(r"(?<![A-Za-z0-9])([0-9][A-Za-z0-9]{3})(?:_complex)?", re.IGNORECASE)

PROTEIN_RESIDUES = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
    "MSE",
    "CME",
}
RNA_RESIDUES = {
    "A",
    "C",
    "G",
    "U",
    "I",
    "ADE",
    "CYT",
    "GUA",
    "URA",
    "URI",
    "PSU",
}
RNA_ATOM_MARKERS = {"P", "OP1", "OP2", "O2'", "O3'", "O4'", "O5'", "C1'", "C2'", "C3'", "C4'", "C5'"}
NORMALIZED_DOCKQ_MAPPING = "AB:AB"

CSV_COLUMNS = [
    "Model",
    "level",
    "target_id",
    "native_path",
    "prediction_path",
    "DockQ",
    "GlobalDockQ",
    "iRMSD",
    "LRMSD",
    "fnat",
    "fnonnat",
    "F1",
    "clashes",
    "nat_correct",
    "nat_total",
    "nonnat_count",
    "model_total",
    "interface",
    "chain1",
    "chain2",
    "class1",
    "class2",
    "best_mapping",
    "chain_map",
    "dockq_json",
    "error",
]

CSV_KEY_COLUMNS = ["Model", "level", "target_id", "native_path", "prediction_path"]


def parse_model_names(raw_model_names):
    model_names = []
    for raw_name in raw_model_names:
        model_names.extend(
            name.strip().lower()
            for name in re.split(r"[,\s]+", raw_name)
            if name.strip()
        )
    return model_names


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute DockQ interface scores for level 3 RNA-protein decoys."
    )
    parser.add_argument(
        "--models",
        nargs="+",
        default=None,
        help=(
            "Prediction model folders to process. Example: --models af3, or "
            "--models boltz chai. Default: CASP17_DOCKQ_MODELS or af3."
        ),
    )
    parser.add_argument(
        "--native-dir",
        type=Path,
        default=NATIVE_DIR,
        help=f"Level 3 native RNA-protein complex directory. Default: {NATIVE_DIR}",
    )
    parser.add_argument(
        "--prediction-root",
        type=Path,
        default=GENERATED_DECOYS,
        help=f"Prediction root containing model folders. Default: {GENERATED_DECOYS}",
    )
    parser.add_argument(
        "--level",
        type=int,
        default=3,
        help="Dataset level to score. Default: 3.",
    )
    parser.add_argument(
        "--dockq-bin",
        default=os.environ.get("DOCKQ_BIN", "DockQ"),
        help="DockQ executable path/name. Default: DOCKQ_BIN or DockQ.",
    )
    parser.add_argument(
        "--dockq-n-cpu",
        type=int,
        default=int(os.environ.get("CASP17_DOCKQ_N_CPU", "1")),
        help="Value passed to DockQ --n_cpu for each structure. Default: 1.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=int(os.environ.get("CASP17_DOCKQ_WORKERS", "1")),
        help="Number of DockQ tasks to run concurrently in this process. Default: 1.",
    )
    parser.add_argument(
        "--shard-index",
        type=int,
        default=(
            int(os.environ["CASP17_DOCKQ_SHARD_INDEX"])
            if "CASP17_DOCKQ_SHARD_INDEX" in os.environ
            else None
        ),
        help="Zero-based shard index for Slurm arrays. Requires --shard-count.",
    )
    parser.add_argument(
        "--shard-count",
        type=int,
        default=(
            int(os.environ["CASP17_DOCKQ_SHARD_COUNT"])
            if "CASP17_DOCKQ_SHARD_COUNT" in os.environ
            else None
        ),
        help="Total number of shards for Slurm arrays. Requires --shard-index.",
    )
    parser.add_argument(
        "--no-align",
        action="store_true",
        help="Pass DockQ --no_align to match residues by numbering instead of sequence alignment.",
    )
    parser.add_argument(
        "--mapping",
        default=None,
        help="Optional fixed DockQ chain mapping, e.g. AB:HL.",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=CSV_OUTPUT,
        help=f"Mergeable output CSV. Default: {CSV_OUTPUT}",
    )
    parser.add_argument(
        "--dockq-output-dir",
        type=Path,
        default=DOCKQ_OUTPUT_DIR,
        help=f"Directory for raw DockQ JSON files. Default: {DOCKQ_OUTPUT_DIR}",
    )
    parser.add_argument(
        "--native",
        type=Path,
        help="Optional single native/reference complex. Requires --prediction.",
    )
    parser.add_argument(
        "--prediction",
        type=Path,
        help="Optional single prediction/model complex. Requires --native.",
    )
    parser.add_argument(
        "--model-name",
        default="direct",
        help="Model label for --native/--prediction direct mode.",
    )
    parser.add_argument(
        "--max-tasks",
        type=int,
        default=None,
        help="Optional number of prediction/native pairs to process for smoke tests.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Recompute rows even when output CSV already has a successful score.",
    )
    parser.add_argument(
        "--raw-dockq-inputs",
        action="store_true",
        help=(
            "Pass original structures directly to DockQ. By default RNA-protein "
            "pairs are normalized to temporary single-model protein=A/RNA=B PDBs "
            "and scored with a fixed AB:AB mapping."
        ),
    )

    args = parser.parse_args()

    if (args.native is None) != (args.prediction is None):
        parser.error("--native and --prediction must be provided together.")

    if (args.shard_index is None) != (args.shard_count is None):
        parser.error("--shard-index and --shard-count must be provided together.")
    if args.shard_count is not None and args.shard_count < 1:
        parser.error("--shard-count must be at least 1.")
    if args.shard_index is not None and not 0 <= args.shard_index < args.shard_count:
        parser.error("--shard-index must satisfy 0 <= shard-index < shard-count.")
    if args.workers < 1:
        parser.error("--workers must be at least 1.")

    env_models = os.environ.get("CASP17_DOCKQ_MODELS")
    raw_model_names = args.models or ([env_models] if env_models else ["af3"])
    args.models = parse_model_names(raw_model_names)

    return args


def is_structure_file(path):
    return path.name.lower().endswith(STRUCTURE_SUFFIXES)


def strip_structure_suffix(name):
    lower = name.lower()
    for suffix in STRUCTURE_SUFFIXES:
        if lower.endswith(suffix):
            return name[: -len(suffix)]
    return Path(name).stem


def infer_target_id(path):
    text = str(path)

    casp_match = CASP_TARGET_RE.search(text)
    if casp_match:
        return casp_match.group(1).upper()

    for part in [path.name, *reversed(path.parts)]:
        pdb_match = PDB_TARGET_RE.search(part)
        if pdb_match:
            return pdb_match.group(1).upper()

    return None


def iter_structure_files(root):
    if not root.exists():
        print(f"Skipping missing directory: {root}", file=sys.stderr)
        return

    for dirpath, dirnames, filenames in os.walk(root):
        dirnames.sort()
        filenames.sort()

        for filename in filenames:
            path = Path(dirpath) / filename
            if is_structure_file(path):
                yield path


def build_native_index(native_dir):
    native_index = {}
    for native_path in iter_structure_files(native_dir):
        target_id = infer_target_id(native_path)
        if target_id is None:
            continue
        native_index.setdefault(target_id, []).append(native_path)
    return native_index


def batch_prediction_sets(prediction_root, model_names):
    return {model_name: prediction_root / model_name for model_name in model_names}


def collect_batch_tasks(args):
    native_index = build_native_index(args.native_dir)
    prediction_sets = batch_prediction_sets(args.prediction_root, args.models)
    tasks = []

    for model_name, prediction_dir in prediction_sets.items():
        level_dir = prediction_dir / f"level{args.level}"
        for prediction_path in iter_structure_files(level_dir):
            target_id = infer_target_id(prediction_path)

            if target_id is None:
                tasks.append((model_name, args.level, None, prediction_path, None))
                continue

            matched_natives = native_index.get(target_id, [])
            if not matched_natives:
                tasks.append((model_name, args.level, target_id, prediction_path, None))
            else:
                for native_path in matched_natives:
                    tasks.append((model_name, args.level, target_id, prediction_path, native_path))

    return tasks


def collect_direct_task(args):
    prediction_path = args.prediction.resolve()
    native_path = args.native.resolve()
    target_id = infer_target_id(prediction_path) or infer_target_id(native_path)
    return [(args.model_name, args.level, target_id, prediction_path, native_path)]


def shard_tasks(tasks, shard_index, shard_count):
    if shard_index is None:
        return tasks
    return [task for index, task in enumerate(tasks) if index % shard_count == shard_index]


def fmt(value):
    if value is None:
        return ""
    if isinstance(value, bool):
        return str(value).lower()
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        return f"{float(value):.8f}"
    return str(value)


def safe_filename_label(value):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("_") or "unknown"


def prediction_output_label(model_name, level, target_id, prediction_path, prediction_root):
    prediction_dir = prediction_root / model_name

    try:
        relative_path = prediction_path.relative_to(prediction_dir)
    except ValueError:
        relative_path = Path(prediction_path.name)

    relative_path = Path(strip_structure_suffix(str(relative_path)))
    parts = list(relative_path.parts)

    if parts and re.fullmatch(rf"level{level}", parts[0], re.IGNORECASE):
        parts = parts[1:]

    if target_id and parts and parts[0].upper() in {str(target_id).upper(), f"{str(target_id).upper()}_COMPLEX"}:
        parts = parts[1:]

    return safe_filename_label("__".join(parts) if parts else prediction_path.stem)


def dockq_json_path(args, model_name, level, target_id, prediction_path):
    target_label = safe_filename_label(target_id or "unknown")
    output_label = prediction_output_label(
        model_name,
        level,
        target_id,
        prediction_path,
        args.prediction_root,
    )
    return (
        args.dockq_output_dir
        / model_name
        / f"level{level}"
        / target_label
        / f"{target_label}_vs_{output_label}_dockq.json"
    )


def is_cif_path(path):
    lower = path.name.lower()
    return lower.endswith((".cif", ".mmcif", ".cif.gz", ".mmcif.gz"))


def open_structure_text(path):
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open()


def ensure_mmcif_atom_site_occupancy(text):
    if "_atom_site.occupancy" in text:
        return text

    lines = text.splitlines()
    output = []
    changed = False
    index = 0

    while index < len(lines):
        line = lines[index]
        if line.strip() != "loop_":
            output.append(line)
            index += 1
            continue

        output.append(line)
        index += 1
        headers = []
        while index < len(lines) and lines[index].lstrip().startswith("_"):
            headers.append(lines[index])
            index += 1

        header_names = [header.strip().split()[0] for header in headers]
        is_atom_site_loop = headers and all(name.startswith("_atom_site.") for name in header_names)
        if not is_atom_site_loop or "_atom_site.occupancy" in header_names:
            output.extend(headers)
            continue

        try:
            occupancy_index = header_names.index("_atom_site.B_iso_or_equiv") + 1
        except ValueError:
            occupancy_index = len(headers)

        for header_index, header in enumerate(headers):
            if header_index == occupancy_index:
                output.append("_atom_site.occupancy")
            output.append(header)
        if occupancy_index == len(headers):
            output.append("_atom_site.occupancy")

        while index < len(lines):
            data_line = lines[index]
            stripped = data_line.strip()
            if stripped == "#" or stripped == "loop_" or stripped.startswith("_") or stripped.startswith("data_"):
                break

            fields = data_line.split()
            if len(fields) == len(headers):
                fields.insert(occupancy_index, "1.0")
                output.append(" ".join(fields))
                changed = True
            else:
                output.append(data_line)
            index += 1

    repaired = "\n".join(output)
    if text.endswith("\n"):
        repaired += "\n"
    return repaired if changed else text


def load_biopython_model(path):
    try:
        from Bio.PDB import MMCIFParser, PDBParser
    except ImportError as exc:
        raise RuntimeError("Biopython is required to normalize DockQ inputs") from exc

    parser = MMCIFParser(QUIET=True) if is_cif_path(path) else PDBParser(QUIET=True)
    if is_cif_path(path):
        with open_structure_text(path) as handle:
            cif_text = ensure_mmcif_atom_site_occupancy(handle.read())
        structure = parser.get_structure("dockq_input", io.StringIO(cif_text))
    else:
        with open_structure_text(path) as handle:
            structure = parser.get_structure("dockq_input", handle)
    return next(structure.get_models())


def is_polymer_residue(residue):
    return residue.id[0] == " "


def residue_kind(residue):
    if not is_polymer_residue(residue):
        return None

    resname = residue.get_resname().strip().upper()
    if resname in PROTEIN_RESIDUES:
        return "protein"
    if resname in RNA_RESIDUES:
        return "rna"

    atom_ids = {atom.id for atom in residue.get_atoms()}
    if atom_ids & RNA_ATOM_MARKERS and ("P" in atom_ids or "C1'" in atom_ids):
        return "rna"

    if {"N", "CA", "C"} <= atom_ids:
        return "protein"

    return None


def classified_polymer_residues(chain, kind):
    return [residue for residue in chain.get_residues() if residue_kind(residue) == kind]


def select_chain(model, kind):
    candidates = []
    for chain in model:
        residues = classified_polymer_residues(chain, kind)
        if residues:
            candidates.append((len(residues), chain.id, chain, residues))

    if not candidates:
        raise RuntimeError(f"Could not find a {kind} polymer chain")

    candidates.sort(key=lambda item: (-item[0], item[1]))
    _, _, chain, residues = candidates[0]
    return chain, residues


def clone_residue(residue, residue_number):
    cloned = copy.deepcopy(residue)
    cloned.detach_parent()
    cloned.id = (" ", residue_number, " ")
    cloned.segid = " "
    return cloned


def write_normalized_rna_protein_pdb(source_path, output_path):
    try:
        from Bio.PDB import PDBIO
        from Bio.PDB.Chain import Chain
        from Bio.PDB.Model import Model
        from Bio.PDB.Structure import Structure
    except ImportError as exc:
        raise RuntimeError("Biopython is required to normalize DockQ inputs") from exc

    source_model = load_biopython_model(source_path)
    protein_chain, protein_residues = select_chain(source_model, "protein")
    rna_chain, rna_residues = select_chain(source_model, "rna")

    if protein_chain.id == rna_chain.id:
        raise RuntimeError(f"Protein and RNA were both selected from chain {protein_chain.id}")

    structure = Structure("normalized")
    model = Model(0)
    structure.add(model)

    for output_chain_id, residues in (("A", protein_residues), ("B", rna_residues)):
        output_chain = Chain(output_chain_id)
        model.add(output_chain)
        for residue_number, residue in enumerate(residues, start=1):
            output_chain.add(clone_residue(residue, residue_number))

    io = PDBIO()
    io.set_structure(structure)
    io.save(str(output_path))

    return {
        "protein_chain": protein_chain.id,
        "rna_chain": rna_chain.id,
        "protein_residues": len(protein_residues),
        "rna_residues": len(rna_residues),
    }


def prepare_dockq_inputs(args, prediction_path, native_path, temp_dir):
    if args.raw_dockq_inputs:
        return prediction_path, native_path, args.mapping

    normalized_prediction = Path(temp_dir) / "prediction_normalized.pdb"
    normalized_native = Path(temp_dir) / "native_normalized.pdb"
    write_normalized_rna_protein_pdb(prediction_path, normalized_prediction)
    write_normalized_rna_protein_pdb(native_path, normalized_native)

    mapping = args.mapping or NORMALIZED_DOCKQ_MAPPING
    return normalized_prediction, normalized_native, mapping


def run_dockq(args, prediction_path, native_path, output_json):
    output_json.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="dockq_") as temp_dir:
        temp_json = Path(temp_dir) / "dockq.json"
        dockq_prediction_path, dockq_native_path, mapping = prepare_dockq_inputs(
            args,
            prediction_path,
            native_path,
            temp_dir,
        )
        command = [
            args.dockq_bin,
            str(dockq_prediction_path),
            str(dockq_native_path),
            "--short",
            "--json",
            str(temp_json),
            "--n_cpu",
            str(args.dockq_n_cpu),
        ]
        if args.no_align:
            command.append("--no_align")
        if mapping:
            command.extend(["--mapping", mapping])

        completed = subprocess.run(
            command,
            check=False,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

        if completed.returncode != 0:
            message = completed.stderr.strip() or completed.stdout.strip() or f"DockQ exited {completed.returncode}"
            raise RuntimeError(message)

        if not temp_json.exists():
            raise RuntimeError("DockQ finished but did not write JSON output")

        with temp_json.open() as handle:
            dockq_data = json.load(handle)
        dockq_data["model"] = str(prediction_path)
        dockq_data["native"] = str(native_path)
        with temp_json.open("w") as handle:
            json.dump(dockq_data, handle)

        temp_output_json = output_json.with_name(
            f".{output_json.name}.{os.getpid()}.{threading.get_ident()}.tmp"
        )
        try:
            shutil.copyfile(temp_json, temp_output_json)
            os.replace(temp_output_json, output_json)
        finally:
            temp_output_json.unlink(missing_ok=True)


def best_interface_metrics(dockq_data):
    best_result = dockq_data.get("best_result") or {}
    if not isinstance(best_result, dict) or not best_result:
        return "", {}

    best_interface = None
    best_metrics = None

    for interface, metrics in best_result.items():
        if not isinstance(metrics, dict):
            continue
        if best_metrics is None or float(metrics.get("DockQ", -1.0)) > float(best_metrics.get("DockQ", -1.0)):
            best_interface = interface
            best_metrics = metrics

    return best_interface or "", best_metrics or {}


def parse_dockq_json(path):
    with path.open() as handle:
        data = json.load(handle)

    interface, metrics = best_interface_metrics(data)
    chain_map = metrics.get("chain_map") or data.get("best_mapping") or {}
    if isinstance(chain_map, dict):
        chain_map_text = json.dumps(chain_map, sort_keys=True)
    else:
        chain_map_text = fmt(chain_map)

    dockq_value = data.get("best_dockq", metrics.get("DockQ"))

    return {
        "DockQ": fmt(dockq_value),
        "GlobalDockQ": fmt(data.get("GlobalDockQ")),
        "iRMSD": fmt(metrics.get("iRMSD")),
        "LRMSD": fmt(metrics.get("LRMSD")),
        "fnat": fmt(metrics.get("fnat")),
        "fnonnat": fmt(metrics.get("fnonnat")),
        "F1": fmt(metrics.get("F1")),
        "clashes": fmt(metrics.get("clashes")),
        "nat_correct": fmt(metrics.get("nat_correct")),
        "nat_total": fmt(metrics.get("nat_total")),
        "nonnat_count": fmt(metrics.get("nonnat_count")),
        "model_total": fmt(metrics.get("model_total")),
        "interface": interface,
        "chain1": fmt(metrics.get("chain1")),
        "chain2": fmt(metrics.get("chain2")),
        "class1": fmt(metrics.get("class1")),
        "class2": fmt(metrics.get("class2")),
        "best_mapping": fmt(data.get("best_mapping_str")),
        "chain_map": chain_map_text,
    }


def empty_metric_fields():
    return {
        "DockQ": "",
        "GlobalDockQ": "",
        "iRMSD": "",
        "LRMSD": "",
        "fnat": "",
        "fnonnat": "",
        "F1": "",
        "clashes": "",
        "nat_correct": "",
        "nat_total": "",
        "nonnat_count": "",
        "model_total": "",
        "interface": "",
        "chain1": "",
        "chain2": "",
        "class1": "",
        "class2": "",
        "best_mapping": "",
        "chain_map": "",
    }


def task_key(model_name, level, target_id, prediction_path, native_path):
    return (
        str(model_name),
        str(level or ""),
        str(target_id or ""),
        str(native_path.resolve()) if native_path else "",
        str(prediction_path.resolve()),
    )


def row_key(row):
    return tuple(row.get(field, "") for field in CSV_KEY_COLUMNS)


def existing_success_keys(output_csv):
    success_keys = set()
    if not output_csv.exists() or output_csv.stat().st_size == 0:
        return success_keys

    with output_csv.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if not row.get("error"):
                success_keys.add(row_key(row))

    return success_keys


def existing_dockq_metrics(output_json):
    if not output_json.exists() or output_json.stat().st_size == 0:
        return None

    try:
        metrics = parse_dockq_json(output_json)
    except (OSError, json.JSONDecodeError, ValueError, TypeError):
        return None

    if not metrics.get("DockQ") and not metrics.get("GlobalDockQ"):
        return None

    return metrics


def process_task(args, model_name, level, target_id, prediction_path, native_path):
    base_row = {
        "Model": model_name,
        "level": level,
        "target_id": target_id or "",
        "native_path": str(native_path.resolve()) if native_path else "",
        "prediction_path": str(prediction_path.resolve()),
        "dockq_json": "",
        "error": "",
    }

    if native_path is None:
        return {
            **base_row,
            **empty_metric_fields(),
            "error": "No matching level 3 native RNA-protein complex found",
        }

    output_json = dockq_json_path(args, model_name, level, target_id, prediction_path)
    try:
        metrics = None if args.force else existing_dockq_metrics(output_json)
        if metrics is None:
            run_dockq(args, prediction_path.resolve(), native_path.resolve(), output_json)
            metrics = parse_dockq_json(output_json)
        return {
            **base_row,
            **metrics,
            "dockq_json": str(output_json.resolve()),
            "error": "",
        }
    except Exception as exc:
        return {
            **base_row,
            **empty_metric_fields(),
            "dockq_json": str(output_json.resolve()),
            "error": str(exc),
        }


def temp_csv_path(output_csv, args):
    run_parts = ["run"]
    if os.environ.get("SLURM_JOB_ID"):
        run_parts.append(os.environ["SLURM_JOB_ID"])
    if args.shard_index is not None:
        run_parts.append(f"shard_{args.shard_index}_of_{args.shard_count}")
    run_parts.append(str(os.getpid()))
    run_label = safe_filename_label("_".join(run_parts))
    return output_csv.with_name(f"{output_csv.stem}.{run_label}{output_csv.suffix}")


def merge_temp_csv(run_csv, output_csv, force=False):
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    lock_path = output_csv.with_suffix(output_csv.suffix + ".lock")

    with lock_path.open("w") as lock_handle:
        fcntl.flock(lock_handle, fcntl.LOCK_EX)

        if not output_csv.exists() or output_csv.stat().st_size == 0:
            run_csv.replace(output_csv)
            return

        existing_rows = []
        existing_by_key = {}
        with output_csv.open(newline="") as output_handle:
            for row in csv.DictReader(output_handle):
                normalized = {field: row.get(field, "") for field in CSV_COLUMNS}
                key = row_key(normalized)
                if key not in existing_by_key:
                    existing_by_key[key] = normalized
                    existing_rows.append(normalized)

        with run_csv.open(newline="") as run_handle:
            for row in csv.DictReader(run_handle):
                normalized = {field: row.get(field, "") for field in CSV_COLUMNS}
                key = row_key(normalized)
                existing = existing_by_key.get(key)

                if existing is None:
                    existing_by_key[key] = normalized
                    existing_rows.append(normalized)
                    continue

                existing_has_error = bool(existing.get("error"))
                new_has_score = not bool(normalized.get("error"))
                if force or (existing_has_error and new_has_score):
                    existing.update(normalized)

        temp_output = output_csv.with_name(f"{output_csv.stem}.merge_{os.getpid()}{output_csv.suffix}")
        with temp_output.open("w", newline="") as output_handle:
            writer = csv.DictWriter(output_handle, fieldnames=CSV_COLUMNS)
            writer.writeheader()
            writer.writerows(existing_rows)

        temp_output.replace(output_csv)
        run_csv.unlink(missing_ok=True)


def iter_processable_tasks(args, tasks, success_keys):
    skipped = 0
    selected = []

    for task in tasks:
        model_name, level, target_id, prediction_path, native_path = task
        key = task_key(model_name, level, target_id, prediction_path, native_path)
        if key in success_keys and not args.force:
            skipped += 1
            continue

        if args.max_tasks is not None and len(selected) >= args.max_tasks:
            break

        selected.append(task)

    return selected, skipped


def process_tasks(args, tasks, writer):
    processed = 0

    if args.workers == 1:
        for task in tasks:
            row = process_task(args, *task)
            writer.writerow(row)
            processed += 1

            if processed % 25 == 0:
                print(f"Processed {processed} DockQ tasks", flush=True)

        return processed

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = [executor.submit(process_task, args, *task) for task in tasks]
        for future in concurrent.futures.as_completed(futures):
            row = future.result()
            writer.writerow(row)
            processed += 1

            if processed % 25 == 0:
                print(f"Processed {processed} DockQ tasks", flush=True)

    return processed


def main():
    args = parse_args()
    args.native_dir = args.native_dir.resolve()
    args.prediction_root = args.prediction_root.resolve()
    args.output_csv = args.output_csv.resolve()
    args.dockq_output_dir = args.dockq_output_dir.resolve()

    if args.native and args.prediction:
        tasks = collect_direct_task(args)
    else:
        tasks = collect_batch_tasks(args)

    discovered_tasks = len(tasks)
    tasks = shard_tasks(tasks, args.shard_index, args.shard_count)
    success_keys = existing_success_keys(args.output_csv)
    tasks_to_process, skipped = iter_processable_tasks(args, tasks, success_keys)
    run_csv = temp_csv_path(args.output_csv, args)

    with run_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        processed = process_tasks(args, tasks_to_process, writer)

    merge_temp_csv(run_csv, args.output_csv, force=args.force)

    print(f"Discovered {discovered_tasks} DockQ tasks")
    if args.shard_index is not None:
        print(f"Shard {args.shard_index + 1}/{args.shard_count}: {len(tasks)} tasks")
    print(f"Workers: {args.workers}")
    print(f"Processed {processed} DockQ tasks")
    print(f"Skipped {skipped} existing successful tasks")
    print(f"Wrote/merged: {args.output_csv}")


if __name__ == "__main__":
    main()
