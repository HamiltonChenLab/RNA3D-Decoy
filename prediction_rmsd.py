#!/usr/bin/env python3

from pathlib import Path
import argparse
import csv
import math
import os
import re
import sys

from ost import io
from ost.mol import alg


ROOT = Path(os.environ.get("CASP17_ROOT", Path(__file__).resolve().parents[1])).resolve()

NATIVE_DIRS = [
    Path(os.environ.get("CASP17_NATIVE_DIR", ROOT / "DProQA_RNA" / "data" / "native_cif")).resolve(),
    (ROOT / "DProQA_RNA" / "data" / "rna_protein" / "native_cif").resolve(),
]
GENERATED_DECOYS = Path(
    os.environ.get("CASP17_GENERATED_DECOYS", ROOT / "generated_decoys")
).resolve()

PREDICTION_SETS = {
    "af3": GENERATED_DECOYS / "af3",
    "boltz": GENERATED_DECOYS / "boltz",
    "chai": GENERATED_DECOYS / "chai",
    "esmfold2": GENERATED_DECOYS / "esmfold2",
    "protenix": GENERATED_DECOYS / "protenix",
}

AF3_LEVEL2_FOLDER = "fold_output_with_ions"
BOLTZ_LEVEL2_FOLDER = "level2_with_ions"
CHAI_LEVEL2_FOLDER = "level2_with_ions"
ESMFOLD2_LEVEL2_FOLDER = "level2_with_grouped_ions"
PROTENIX_LEVEL2_FOLDER = "level2_with_ions"

OUTPUT_DIR = ROOT / "RMSD"
SUMMARY_CSV = OUTPUT_DIR / "prediction_rmsd_summary.csv"
LOCAL_CSV = OUTPUT_DIR / "prediction_rmsd_local_residue.csv"

STRUCTURE_SUFFIXES = (".pdb", ".cif", ".mmcif", ".ent", ".pdb.gz", ".cif.gz", ".mmcif.gz")

CASP_TARGET_RE = re.compile(r"([RHT]\d{4}(?:v\d+)?)", re.IGNORECASE)
PDB_TARGET_RE = re.compile(r"^([0-9][A-Za-z0-9]{3})")

HYBRID_C4_CA_SELECTION = "hybrid_C4_CA"
RNA_RESIDUE_NAMES = {"A", "C", "G", "I", "U"}
RNA_RESIDUE_QUERY = "(" + " or ".join(
    f"rname={residue_name}" for residue_name in sorted(RNA_RESIDUE_NAMES)
) + ")"
PROTEIN_CA_QUERY = "aname=CA and peptide=true"
HYBRID_C4_CA_QUERY = f"aname=\"C4'\" or ({PROTEIN_CA_QUERY})"

SELECTIONS = [
    ("all_atom", "ALL", "all", None),
    (HYBRID_C4_CA_SELECTION, "C4'/CA", "all", HYBRID_C4_CA_QUERY),
    ("C4_prime", "C4'", "C4'", None),
    ("P_atom", "P", "P", None),
    ("CA_atom", "CA", "all", PROTEIN_CA_QUERY),
]

LEVEL2_RNA_SELECTIONS = [
    ("all_atom", "ALL", "all", RNA_RESIDUE_QUERY),
    ("C4_prime", "C4'", "C4'", RNA_RESIDUE_QUERY),
    ("P_atom", "P", "P", RNA_RESIDUE_QUERY),
]
LOCAL_ATOM_SELECTION_CHOICES = [
    "all",
    *(selection[0] for selection in SELECTIONS),
]

MATCHERS = {
    "number": alg.MatchResidueByNum,
    "index": alg.MatchResidueByIdx,
}

SUMMARY_FIELDS = [
    "Model",
    "level",
    "target_id",
    "native_path",
    "prediction_path",
    "atom_selection",
    "match_method",
    "global_RMSD",
    "RMSD_from_local_distances",
    "n_atom_pairs",
    "n_residue_pairs",
    "ncycles",
    "error",
]

LOCAL_FIELDS = [
    "Model",
    "level",
    "target_id",
    "native_path",
    "prediction_path",
    "atom_selection",
    "atom_name",
    "local_index",
    "chain_id",
    "residue_number",
    "residue_name",
    "native_chain_id",
    "native_residue_number",
    "native_residue_name",
    "local_RMSD",
    "n_atom_pairs",
]

SUMMARY_KEY_FIELDS = [
    "Model",
    "level",
    "target_id",
    "native_path",
    "prediction_path",
    "atom_selection",
    "match_method",
]

LOCAL_KEY_FIELDS = [
    "Model",
    "level",
    "target_id",
    "native_path",
    "prediction_path",
    "atom_selection",
    "atom_name",
    "local_index",
    "chain_id",
    "residue_number",
    "residue_name",
    "native_chain_id",
    "native_residue_number",
    "native_residue_name",
]


def parse_model_names(raw_model_names):
    model_names = []

    for raw_name in raw_model_names:
        model_names.extend(
            name.strip().lower()
            for name in re.split(r"[,\s]+", raw_name)
            if name.strip()
        )

    return model_names


def parse_level_names(raw_level_names):
    levels = []

    for raw_name in raw_level_names:
        levels.extend(
            int(name)
            for name in re.split(r"[,\s]+", raw_name)
            if name.strip()
        )

    return levels


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Compute fitted global RMSD and residue-level local RMSD for all "
            "atoms, RNA C4' atoms, and RNA P atoms."
        )
    )

    parser.add_argument(
        "--native",
        type=Path,
        help="Optional single native/reference structure. Requires --prediction.",
    )
    parser.add_argument(
        "--prediction",
        type=Path,
        help="Optional single prediction/model structure. Requires --native.",
    )
    parser.add_argument(
        "--model-name",
        default="direct",
        help="Model label for --native/--prediction direct mode.",
    )
    parser.add_argument(
        "--native-dir",
        type=Path,
        default=None,
        help=f"Batch native directory. Default: {NATIVE_DIRS[0]} and {NATIVE_DIRS[1]}",
    )
    parser.add_argument(
        "--native-level",
        type=int,
        default=None,
        help=(
            "Assign this level to natives when --native-dir does not contain "
            "levelN folders. Useful for flat level-specific native directories."
        ),
    )
    parser.add_argument(
        "--prediction-root",
        type=Path,
        default=GENERATED_DECOYS,
        help=(
            "Batch prediction root containing model folders such as af3/, "
            "boltz/, chai/, esmfold2/, and protenix/. "
            f"Default: {GENERATED_DECOYS}"
        ),
    )
    parser.add_argument(
        "--models",
        nargs="+",
        default=None,
        help=(
            "Batch model folders to process. Example: --models af3, or "
            "--models esmfold2 protenix. Default: af3."
        ),
    )
    parser.add_argument(
        "--levels",
        nargs="+",
        default=None,
        help=(
            "Prediction levels to process. Example: --levels 3 or "
            "--levels 1 2. Default: all levels found."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=OUTPUT_DIR,
        help=f"Output directory. Default: {OUTPUT_DIR}",
    )
    parser.add_argument(
        "--summary-csv",
        type=Path,
        default=None,
        help="Optional explicit summary CSV path.",
    )
    parser.add_argument(
        "--local-csv",
        type=Path,
        default=None,
        help="Optional explicit residue-level local RMSD CSV path.",
    )
    parser.add_argument(
        "--local-atom-selection",
        choices=LOCAL_ATOM_SELECTION_CHOICES,
        default="all_atom",
        help=(
            "Atom selection to write to the residue-level local CSV. Use 'all' "
            "to write every selection, or 'hybrid_C4_CA' to write rows from "
            "one shared RNA C4' plus protein CA fit for mixed RNA-protein complexes. "
            "Default: all_atom."
        ),
    )
    parser.add_argument(
        "--match",
        choices=sorted(MATCHERS),
        default="number",
        help="Residue matching method before atom pairing. Default: number.",
    )

    args = parser.parse_args()

    if (args.native is None) != (args.prediction is None):
        parser.error("--native and --prediction must be provided together.")

    env_models = os.environ.get("CASP17_RMSD_MODELS")
    raw_model_names = args.models or ([env_models] if env_models else ["af3"])
    args.models = parse_model_names(raw_model_names)

    env_levels = os.environ.get("CASP17_RMSD_LEVELS")
    raw_level_names = args.levels or ([env_levels] if env_levels else [])
    args.levels = parse_level_names(raw_level_names) if raw_level_names else None

    env_native_dir = os.environ.get("CASP17_RMSD_NATIVE_DIR")
    args.native_dir = Path(args.native_dir or env_native_dir).resolve() if (args.native_dir or env_native_dir) else None

    env_native_level = os.environ.get("CASP17_RMSD_NATIVE_LEVEL")
    if args.native_level is None and env_native_level:
        args.native_level = int(env_native_level)

    unknown_models = sorted(set(args.models) - set(PREDICTION_SETS))
    if unknown_models:
        parser.error(
            f"Unknown model(s): {', '.join(unknown_models)}. "
            f"Known models: {', '.join(PREDICTION_SETS)}"
        )

    return args


def is_structure_file(path):
    return path.name.lower().endswith(STRUCTURE_SUFFIXES)


def infer_target_id(path):
    text = str(path)

    casp_match = CASP_TARGET_RE.search(text)
    if casp_match:
        return casp_match.group(1).upper()

    pdb_match = PDB_TARGET_RE.search(path.name)
    if pdb_match:
        return pdb_match.group(1).upper()

    return None


def infer_level(path, root_dir):
    try:
        rel_parts = path.relative_to(root_dir).parts
    except ValueError:
        rel_parts = path.parts

    for part in rel_parts:
        match = re.fullmatch(r"level(\d+)", part, re.IGNORECASE)
        if match:
            return int(match.group(1))

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


def build_native_index(native_dirs, native_level=None):
    native_index = {}

    if isinstance(native_dirs, Path):
        native_dirs = [native_dirs]

    for native_dir in native_dirs:
        for native_path in iter_structure_files(native_dir):
            target_id = infer_target_id(native_path)
            level = infer_level(native_path, native_dir)

            if level is None and native_level is not None:
                level = native_level

            if target_id is not None and level is not None:
                native_index.setdefault((level, target_id), []).append(native_path)

    return native_index


def normalize_blank_chain_names(entity):
    used = {
        ch.GetName().strip()
        for ch in entity.chains
        if ch.GetName().strip() not in {"", "'", '"'}
    }

    candidates = list("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789")
    editor = entity.EditXCS()

    for chain in list(entity.chains):
        name = chain.GetName().strip()
        if name not in {"", "'", '"'}:
            continue

        new_name = next(c for c in candidates if c not in used)
        used.add(new_name)
        editor.RenameChain(chain, new_name)

    return entity


def load_structure(path):
    entity = io.LoadEntity(str(path))
    normalize_blank_chain_names(entity)
    return entity


def fmt(value):
    if value is None:
        return ""
    return f"{float(value):.8f}"


def residue_number(residue):
    number = residue.GetNumber()
    try:
        return str(number.GetNum())
    except Exception:
        return str(number)


def atom_xyz(atom):
    pos = atom.GetPos()
    return float(pos[0]), float(pos[1]), float(pos[2])


def atom_distance(prediction_atom, native_atom):
    pred_x, pred_y, pred_z = atom_xyz(prediction_atom)
    native_x, native_y, native_z = atom_xyz(native_atom)

    dx = pred_x - native_x
    dy = pred_y - native_y
    dz = pred_z - native_z
    squared_distance = dx * dx + dy * dy + dz * dz

    return math.sqrt(squared_distance), squared_distance


def empty_view(entity):
    try:
        return entity.CreateEmptyView()
    except AttributeError:
        return entity.handle.CreateEmptyView()


def ordered_chain_view(entity, chains):
    view = empty_view(entity)

    for chain in chains:
        view.AddChain(chain)
        for residue in chain.residues:
            view.AddResidue(residue)
            for atom in residue.atoms:
                view.AddAtom(atom)

    return view


def chain_pair_score(prediction_chain, native_chain, matcher):
    try:
        prediction_view, native_view = matcher(
            ordered_chain_view(prediction_chain.GetEntity(), [prediction_chain]),
            ordered_chain_view(native_chain.GetEntity(), [native_chain]),
            "all",
        )
    except Exception:
        return 0

    return prediction_view.atom_count


def best_chain_pairs(prediction, native, match_method):
    matcher = MATCHERS[match_method]
    candidates = []

    for prediction_index, prediction_chain in enumerate(prediction.chains):
        for native_index, native_chain in enumerate(native.chains):
            score = chain_pair_score(prediction_chain, native_chain, matcher)
            if score <= 0:
                continue

            residue_count_delta = abs(
                prediction_chain.residue_count - native_chain.residue_count
            )
            candidates.append((
                -score,
                residue_count_delta,
                native_index,
                prediction_index,
                prediction_chain,
                native_chain,
            ))

    candidates.sort()
    used_prediction_chains = set()
    used_native_chains = set()
    pairs = []

    for (
        _score,
        _residue_count_delta,
        native_index,
        prediction_index,
        prediction_chain,
        native_chain,
    ) in candidates:
        if prediction_index in used_prediction_chains or native_index in used_native_chains:
            continue

        used_prediction_chains.add(prediction_index)
        used_native_chains.add(native_index)
        pairs.append((native_index, prediction_chain, native_chain))

    pairs.sort(key=lambda item: item[0])

    return [(prediction_chain, native_chain) for _, prediction_chain, native_chain in pairs]


def chain_ordered_views(prediction, native, match_method):
    pairs = best_chain_pairs(prediction, native, match_method)

    if not pairs:
        return prediction, native

    prediction_chains = [prediction_chain for prediction_chain, _native_chain in pairs]
    native_chains = [native_chain for _prediction_chain, native_chain in pairs]

    return (
        ordered_chain_view(prediction, prediction_chains),
        ordered_chain_view(native, native_chains),
    )


def atom_row_metadata(prefix, atom):
    residue = atom.GetResidue()
    chain = residue.GetChain()
    x, y, z = atom_xyz(atom)

    return {
        f"{prefix}_chain_id": chain.GetName(),
        f"{prefix}_residue_number": residue_number(residue),
        f"{prefix}_residue_name": residue.GetName(),
        f"{prefix}_atom_name": atom.GetName(),
        f"{prefix}_x": fmt(x),
        f"{prefix}_y": fmt(y),
        f"{prefix}_z": fmt(z),
    }


def residue_pair_key(prediction_atom, native_atom):
    prediction_residue = prediction_atom.GetResidue()
    prediction_chain = prediction_residue.GetChain()
    native_residue = native_atom.GetResidue()
    native_chain = native_residue.GetChain()

    return (
        prediction_chain.GetName(),
        residue_number(prediction_residue),
        prediction_residue.GetName(),
        native_chain.GetName(),
        residue_number(native_residue),
        native_residue.GetName(),
    )


def start_residue_group(prediction_atom, native_atom):
    (
        prediction_chain_id,
        prediction_residue_number,
        prediction_residue_name,
        native_chain_id,
        native_residue_number,
        native_residue_name,
    ) = residue_pair_key(prediction_atom, native_atom)

    return {
        "chain_id": prediction_chain_id,
        "residue_number": prediction_residue_number,
        "residue_name": prediction_residue_name,
        "native_chain_id": native_chain_id,
        "native_residue_number": native_residue_number,
        "native_residue_name": native_residue_name,
        "sum_squared_distance": 0.0,
        "n_atom_pairs": 0,
    }


def add_atom_to_residue_group(group, squared_distance):
    group["sum_squared_distance"] += squared_distance
    group["n_atom_pairs"] += 1


def residue_group_to_row(group, atom_selection, atom_name, local_index):
    n_atom_pairs = group["n_atom_pairs"]
    local_rmsd = math.sqrt(group["sum_squared_distance"] / n_atom_pairs)

    return {
        "atom_selection": atom_selection,
        "atom_name": atom_name,
        "local_index": local_index,
        "chain_id": group["chain_id"],
        "residue_number": group["residue_number"],
        "residue_name": group["residue_name"],
        "native_chain_id": group["native_chain_id"],
        "native_residue_number": group["native_residue_number"],
        "native_residue_name": group["native_residue_name"],
        "local_RMSD": fmt(local_rmsd),
        "n_atom_pairs": n_atom_pairs,
    }


def selections_for_level(level):
    if level == 2:
        return LEVEL2_RNA_SELECTIONS

    return SELECTIONS


def score_selection(
    native,
    prediction,
    atom_selection,
    atom_name,
    atoms,
    match_method,
    selection_query=None,
):
    matcher = MATCHERS[match_method]

    if selection_query is not None:
        native = native.Select(selection_query)
        prediction = prediction.Select(selection_query)

    if native.atom_count == 0 or prediction.atom_count == 0:
        return None, None

    prediction, native = chain_ordered_views(prediction, native, match_method)
    prediction_view, native_view = matcher(prediction, native, atoms)

    if prediction_view.atom_count == 0:
        return None, None

    if prediction_view.atom_count != native_view.atom_count:
        raise RuntimeError(
            f"Matched atom count differs for {atom_selection}: "
            f"prediction={prediction_view.atom_count}, native={native_view.atom_count}"
        )

    superposition = alg.SuperposeSVD(prediction_view, native_view, True)

    local_rows = []
    residue_groups = {}
    sum_squared_distances = 0.0
    n_atom_pairs = 0

    for prediction_atom, native_atom in zip(
        prediction_view.GetAtomList(),
        native_view.GetAtomList(),
    ):
        _, squared_distance = atom_distance(prediction_atom, native_atom)
        sum_squared_distances += squared_distance
        n_atom_pairs += 1

        key = residue_pair_key(prediction_atom, native_atom)
        if key not in residue_groups:
            residue_groups[key] = start_residue_group(prediction_atom, native_atom)
        add_atom_to_residue_group(residue_groups[key], squared_distance)

    for local_index, group in enumerate(residue_groups.values(), start=1):
        local_rows.append(
            residue_group_to_row(group, atom_selection, atom_name, local_index)
        )

    rmsd_from_distances = math.sqrt(sum_squared_distances / n_atom_pairs)

    summary = {
        "atom_selection": atom_selection,
        "match_method": match_method,
        "global_RMSD": fmt(superposition.rmsd),
        "RMSD_from_local_distances": fmt(rmsd_from_distances),
        "n_atom_pairs": n_atom_pairs,
        "n_residue_pairs": prediction_view.residue_count,
        "ncycles": superposition.ncycles,
    }

    return summary, local_rows


def batch_prediction_sets(prediction_root):
    return {
        "af3": prediction_root / "af3",
        "boltz": prediction_root / "boltz",
        "chai": prediction_root / "chai",
        "esmfold2": prediction_root / "esmfold2",
        "protenix": prediction_root / "protenix",
    }


def collect_batch_tasks(native_dir, prediction_sets, native_level=None, levels=None):
    native_dirs = [native_dir] if native_dir else NATIVE_DIRS
    native_index = build_native_index(native_dirs, native_level)
    selected_levels = set(levels) if levels is not None else None
    tasks = []

    for model_name, prediction_dir in prediction_sets.items():
        for prediction_path in iter_structure_files(prediction_dir):
            target_id = infer_target_id(prediction_path)
            level = infer_level(prediction_path, prediction_dir)
            rel_parts = prediction_path.relative_to(prediction_dir).parts

            if model_name == "boltz":
                if rel_parts and rel_parts[0] == BOLTZ_LEVEL2_FOLDER:
                    level = 2
                elif level == 2:
                    continue

            if model_name == "chai":
                if rel_parts and rel_parts[0] == CHAI_LEVEL2_FOLDER:
                    level = 2
                elif level == 2:
                    continue

            if model_name == "esmfold2":
                if rel_parts and rel_parts[0] == ESMFOLD2_LEVEL2_FOLDER:
                    level = 2
                elif level == 2:
                    continue

            if model_name == "protenix":
                if rel_parts and rel_parts[0] == PROTENIX_LEVEL2_FOLDER:
                    level = 2
                elif level == 2:
                    continue

            if selected_levels is not None and level not in selected_levels:
                continue

            if model_name == "af3" and level == 2:
                if len(rel_parts) < 2 or rel_parts[1] != AF3_LEVEL2_FOLDER:
                    continue

            if target_id is None or level is None:
                tasks.append((model_name, level, target_id, prediction_path, None))
                continue

            matched_natives = native_index.get((level, target_id), [])

            if not matched_natives:
                tasks.append((model_name, level, target_id, prediction_path, None))
            else:
                for native_path in matched_natives:
                    tasks.append((model_name, level, target_id, prediction_path, native_path))

    return tasks


def collect_direct_task(args):
    prediction_path = args.prediction.resolve()
    native_path = args.native.resolve()
    target_id = infer_target_id(prediction_path) or infer_target_id(native_path)
    level = infer_level(prediction_path, Path(prediction_path.anchor))

    if args.levels is not None and level not in set(args.levels):
        return []

    return [(args.model_name, level, target_id, prediction_path, native_path)]


def write_missing_native(summary_writer, model_name, level, target_id, prediction_path, match_method):
    summary_writer.writerow({
        "Model": model_name,
        "level": level,
        "target_id": target_id,
        "native_path": "",
        "prediction_path": str(prediction_path.resolve()),
        "atom_selection": "",
        "match_method": match_method,
        "global_RMSD": "",
        "RMSD_from_local_distances": "",
        "n_atom_pairs": "",
        "n_residue_pairs": "",
        "ncycles": "",
        "error": "No matching native found",
    })


def should_write_local_rows(local_atom_selection, atom_selection):
    if local_atom_selection == "all":
        return True

    return atom_selection == local_atom_selection


def temp_csv_path(output_csv):
    return output_csv.with_name(f"{output_csv.stem}.run_{os.getpid()}{output_csv.suffix}")


def merge_temp_csv_path(output_csv):
    return output_csv.with_name(f"{output_csv.stem}.merge_{os.getpid()}{output_csv.suffix}")


def csv_key(row, key_fields):
    return tuple(row.get(field, "") for field in key_fields)


def merge_csv(run_csv, output_csv, fieldnames, key_fields, replace_error_rows=False):
    if not output_csv.exists() or output_csv.stat().st_size == 0:
        run_csv.replace(output_csv)
        return

    existing_keys = set()
    existing_rows = []
    with output_csv.open(newline="") as output_handle:
        reader = csv.DictReader(output_handle)
        for row in reader:
            existing_rows.append(row)
            existing_keys.add(csv_key(row, key_fields))

    if replace_error_rows:
        with run_csv.open(newline="") as run_handle:
            run_rows = list(csv.DictReader(run_handle))

        run_rows_by_key = {csv_key(row, key_fields): row for row in run_rows}
        output_rows = []

        for row in existing_rows:
            key = csv_key(row, key_fields)
            new_row = run_rows_by_key.get(key)
            if row.get("error") and new_row is not None and not new_row.get("error"):
                output_rows.append({field: new_row.get(field, "") for field in fieldnames})
            else:
                output_rows.append({field: row.get(field, "") for field in fieldnames})

        for row in run_rows:
            key = csv_key(row, key_fields)
            if key not in existing_keys:
                output_rows.append({field: row.get(field, "") for field in fieldnames})
                existing_keys.add(key)

        temp_output_csv = merge_temp_csv_path(output_csv)
        with temp_output_csv.open("w", newline="") as output_handle:
            writer = csv.DictWriter(output_handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(output_rows)

        temp_output_csv.replace(output_csv)
        run_csv.unlink()
        return

    with run_csv.open(newline="") as run_handle, output_csv.open("a", newline="") as output_handle:
        reader = csv.DictReader(run_handle)
        writer = csv.DictWriter(output_handle, fieldnames=fieldnames)

        for row in reader:
            key = csv_key(row, key_fields)
            if key in existing_keys:
                continue

            writer.writerow({field: row.get(field, "") for field in fieldnames})
            existing_keys.add(key)

    run_csv.unlink()


def main():
    args = parse_args()

    output_dir = args.output_dir.resolve()
    summary_csv = (args.summary_csv or output_dir / SUMMARY_CSV.name).resolve()
    local_csv = (args.local_csv or output_dir / LOCAL_CSV.name).resolve()

    output_dir.mkdir(parents=True, exist_ok=True)
    summary_csv.parent.mkdir(parents=True, exist_ok=True)
    local_csv.parent.mkdir(parents=True, exist_ok=True)

    summary_run_csv = temp_csv_path(summary_csv)
    local_run_csv = temp_csv_path(local_csv)

    if args.native is not None:
        tasks = collect_direct_task(args)
        print(f"Processing direct model: {args.model_name}")
    else:
        prediction_sets = batch_prediction_sets(args.prediction_root.resolve())
        selected_prediction_sets = {
            model_name: prediction_sets[model_name]
            for model_name in args.models
        }
        tasks = collect_batch_tasks(
            args.native_dir,
            selected_prediction_sets,
            native_level=args.native_level,
            levels=args.levels,
        )
        print(f"Processing models: {', '.join(args.models)}")
        if args.levels is not None:
            print(f"Processing levels: {', '.join(str(level) for level in args.levels)}")
        print(f"Native directories: {args.native_dir or NATIVE_DIRS}")
        if args.native_level is not None:
            print(f"Flat native directory level: {args.native_level}")

    print(f"Found {len(tasks)} prediction tasks.")

    native_cache = {}

    with summary_run_csv.open("w", newline="") as summary_handle, \
         local_run_csv.open("w", newline="") as local_handle:

        summary_writer = csv.DictWriter(summary_handle, fieldnames=SUMMARY_FIELDS)
        local_writer = csv.DictWriter(local_handle, fieldnames=LOCAL_FIELDS)

        summary_writer.writeheader()
        local_writer.writeheader()

        for model_name, level, target_id, prediction_path, native_path in tasks:
            abs_prediction_path = prediction_path.resolve()

            if native_path is None:
                write_missing_native(
                    summary_writer,
                    model_name,
                    level,
                    target_id,
                    abs_prediction_path,
                    args.match,
                )
                continue

            abs_native_path = native_path.resolve()

            try:
                if str(abs_native_path) not in native_cache:
                    native_cache[str(abs_native_path)] = load_structure(abs_native_path)

                native = native_cache[str(abs_native_path)]

                for atom_selection, atom_name, atoms, selection_query in selections_for_level(level):
                    try:
                        prediction = load_structure(abs_prediction_path)
                        summary, local_rows = score_selection(
                            native,
                            prediction,
                            atom_selection,
                            atom_name,
                            atoms,
                            args.match,
                            selection_query,
                        )

                        if summary is None:
                            continue

                        summary_writer.writerow({
                            "Model": model_name,
                            "level": level,
                            "target_id": target_id,
                            "native_path": str(abs_native_path),
                            "prediction_path": str(abs_prediction_path),
                            **summary,
                            "error": "",
                        })

                        if should_write_local_rows(args.local_atom_selection, atom_selection):
                            for row in local_rows:
                                local_writer.writerow({
                                    "Model": model_name,
                                    "level": level,
                                    "target_id": target_id,
                                    "native_path": str(abs_native_path),
                                    "prediction_path": str(abs_prediction_path),
                                    **row,
                                })

                    except Exception as exc:
                        summary_writer.writerow({
                            "Model": model_name,
                            "level": level,
                            "target_id": target_id,
                            "native_path": str(abs_native_path),
                            "prediction_path": str(abs_prediction_path),
                            "atom_selection": atom_selection,
                            "match_method": args.match,
                            "global_RMSD": "",
                            "RMSD_from_local_distances": "",
                            "n_atom_pairs": "",
                            "n_residue_pairs": "",
                            "ncycles": "",
                            "error": str(exc),
                        })

            except Exception as exc:
                summary_writer.writerow({
                    "Model": model_name,
                    "level": level,
                    "target_id": target_id,
                    "native_path": str(abs_native_path),
                    "prediction_path": str(abs_prediction_path),
                    "atom_selection": "",
                    "match_method": args.match,
                    "global_RMSD": "",
                    "RMSD_from_local_distances": "",
                    "n_atom_pairs": "",
                    "n_residue_pairs": "",
                    "ncycles": "",
                    "error": str(exc),
                })

    merge_csv(
        summary_run_csv,
        summary_csv,
        SUMMARY_FIELDS,
        SUMMARY_KEY_FIELDS,
        replace_error_rows=True,
    )
    merge_csv(local_run_csv, local_csv, LOCAL_FIELDS, LOCAL_KEY_FIELDS)

    print(f"Merged summary CSV into: {summary_csv}")
    print(f"Merged local CSV into: {local_csv}")


if __name__ == "__main__":
    main()
