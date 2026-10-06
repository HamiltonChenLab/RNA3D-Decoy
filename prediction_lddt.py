#!/usr/bin/env python3

from pathlib import Path
import argparse
import csv
import os
import re
import sys

from ost import io
from ost.mol.alg.lddt import lDDTScorer


ROOT = Path("/hpc/netapp/casp17/CASP17").resolve()

NATIVE_DIRS = [
    ROOT / "DProQA_RNA" / "data" / "native_cif",
    ROOT / "DProQA_RNA" / "data" / "rna_protein" / "native_cif",
]

GENERATED_DECOYS = ROOT / "generated_decoys"

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

OUTPUT_DIR = ROOT / "lDDT"
SUMMARY_CSV = OUTPUT_DIR / "prediction_lddt_summary.csv"
LOCAL_CSV = OUTPUT_DIR / "prediction_lddt_local.csv"

STRUCTURE_SUFFIXES = (".pdb", ".cif", ".mmcif", ".ent", ".pdb.gz", ".cif.gz", ".mmcif.gz")

CASP_TARGET_RE = re.compile(r"([RHT]\d{4}(?:v\d+)?)", re.IGNORECASE)
PDB_TARGET_RE = re.compile(r"^([0-9][A-Za-z0-9]{3})")

RNA_RESIDUE_NAMES = {"A", "C", "G", "I", "U"}
PROTEIN_CA_QUERY = "aname=CA and peptide=true"

SELECTIONS = [
    ("all_atom", "ALL", None),
    ("hybrid_C4_CA", "C4'/CA", f"aname=\"C4'\" or ({PROTEIN_CA_QUERY})"),
    ("C4_prime", "C4'", "aname=\"C4'\""),
    ("P_atom", "P", "aname=P"),
    ("CA_atom", "CA", PROTEIN_CA_QUERY),
]

SUMMARY_FIELDS = [
    "Model",
    "level",
    "target_id",
    "native_path",
    "prediction_path",
    "atom_selection",
    "global_lDDT",
    "average_local_lDDT",
    "n_local_scores",
    "n_native_atoms",
    "n_prediction_atoms",
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
    "local_lDDT",
]

SUMMARY_KEY_FIELDS = [
    "Model",
    "level",
    "target_id",
    "native_path",
    "prediction_path",
    "atom_selection",
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
        description="Compute lDDT scores for selected prediction model folders."
    )
    parser.add_argument(
        "--models",
        nargs="+",
        default=None,
        help=(
            "Model folders to process. Example: --models af3, or "
            "--models esmfold2 protenix. Default: af3."
        ),
    )
    parser.add_argument(
        "--levels",
        nargs="+",
        default=None,
        help="Prediction levels to process. Example: --levels 3 or --levels 1 2. Default: all.",
    )
    parser.add_argument(
        "--native-dir",
        type=Path,
        default=None,
        help=f"Native/reference structure directory. Default: {NATIVE_DIRS[0]} and {NATIVE_DIRS[1]}",
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

    args = parser.parse_args()
    env_models = os.environ.get("CASP17_LDDT_MODELS")
    raw_model_names = args.models or ([env_models] if env_models else ["af3"])
    args.models = parse_model_names(raw_model_names)

    env_levels = os.environ.get("CASP17_LDDT_LEVELS")
    raw_level_names = args.levels or ([env_levels] if env_levels else [])
    args.levels = parse_level_names(raw_level_names) if raw_level_names else None

    env_native_dir = os.environ.get("CASP17_LDDT_NATIVE_DIR")
    args.native_dir = Path(args.native_dir or env_native_dir).resolve() if (args.native_dir or env_native_dir) else None

    env_native_level = os.environ.get("CASP17_LDDT_NATIVE_LEVEL")
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


def infer_level(path, prediction_dir):
    rel_parts = path.relative_to(prediction_dir).parts

    for part in rel_parts:
        match = re.fullmatch(r"level(\d+)", part, re.IGNORECASE)
        if match:
            return int(match.group(1))

    return None


def iter_structure_files(root):
    if not root.exists():
        print(f"Skipping missing directory: {root}", file=sys.stderr)
        return

    for dirpath, _, filenames in os.walk(root):
        for filename in sorted(filenames):
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


def mean_defined(values):
    values = [float(v) for v in values if v is not None]
    if not values:
        return None
    return sum(values) / len(values)


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


def active_chains(entity):
    return [chain for chain in entity.chains if chain.GetAtomCount() > 0]


def chain_names(entity):
    return [chain.GetName() for chain in active_chains(entity)]


def chain_residue_names(chain):
    return {
        residue.GetName().strip().upper()
        for residue in chain.residues
    }


def is_rna_chain(chain):
    return bool(chain_residue_names(chain) & RNA_RESIDUE_NAMES)


def chain_kind(chain):
    return "rna" if is_rna_chain(chain) else "other"


def chain_residue_count(chain):
    return sum(1 for _ in chain.residues)


def chain_profile(chain):
    return {
        "name": chain.GetName(),
        "kind": chain_kind(chain),
        "atom_count": chain.GetAtomCount(),
        "residue_count": chain_residue_count(chain),
    }


def chain_match_score(prediction_profile, native_profile):
    if prediction_profile["kind"] != native_profile["kind"]:
        return None

    return (
        abs(prediction_profile["residue_count"] - native_profile["residue_count"]),
        abs(prediction_profile["atom_count"] - native_profile["atom_count"]),
        0 if prediction_profile["name"] == native_profile["name"] else 1,
        native_profile["name"],
    )


def is_ambiguous_match(candidates):
    if len(candidates) < 2:
        return False

    best_score, _ = candidates[0]
    next_score, _ = candidates[1]
    return best_score[:3] == next_score[:3]


def infer_profile_chain_mapping(native_view, prediction_view):
    native_profiles = [chain_profile(chain) for chain in active_chains(native_view)]
    prediction_profiles = [
        chain_profile(chain)
        for chain in active_chains(prediction_view)
    ]
    mapping = {}
    used_native_chains = set()

    def compatible_native_count(prediction_profile):
        return sum(
            chain_match_score(prediction_profile, native_profile) is not None
            for native_profile in native_profiles
        )

    for prediction_profile in sorted(
        prediction_profiles,
        key=lambda profile: (
            compatible_native_count(profile),
            -profile["residue_count"],
            profile["name"],
        ),
    ):
        candidates = []

        for native_profile in native_profiles:
            if native_profile["name"] in used_native_chains:
                continue

            score = chain_match_score(prediction_profile, native_profile)
            if score is None:
                continue

            candidates.append((score, native_profile))

        if not candidates:
            continue

        candidates.sort(key=lambda candidate: candidate[0])
        if is_ambiguous_match(candidates):
            continue

        _, native_profile = candidates[0]
        mapping[prediction_profile["name"]] = native_profile["name"]
        used_native_chains.add(native_profile["name"])

    return mapping


def rna_chain_names(entity):
    return [
        chain.GetName()
        for chain in active_chains(entity)
        if is_rna_chain(chain)
    ]


def infer_chain_mapping(native_view, prediction_view):
    native_chains = chain_names(native_view)
    prediction_chains = chain_names(prediction_view)

    if len(native_chains) <= 1 and len(prediction_chains) <= 1:
        return None

    chain_mapping = infer_profile_chain_mapping(native_view, prediction_view)

    if chain_mapping:
        return chain_mapping

    native_chain_set = set(native_chains)
    chain_mapping = {
        chain_name: chain_name
        for chain_name in prediction_chains
        if chain_name in native_chain_set
    }

    if not chain_mapping:
        native_rna_chains = rna_chain_names(native_view)
        prediction_rna_chains = rna_chain_names(prediction_view)

        if len(native_rna_chains) == 1 and len(prediction_rna_chains) == 1:
            return {prediction_rna_chains[0]: native_rna_chains[0]}

    if not chain_mapping:
        raise RuntimeError(
            "Could not infer chain mapping. "
            f"Native chains: {','.join(native_chains) or 'none'}; "
            f"prediction chains: {','.join(prediction_chains) or 'none'}"
        )

    return chain_mapping


def score_selection(native, prediction, atom_selection, atom_name, atom_query):
    native_view = native if atom_query is None else native.Select(atom_query)
    prediction_view = prediction if atom_query is None else prediction.Select(atom_query)

    if native_view.atom_count == 0 or prediction_view.atom_count == 0:
        return None, None

    chain_mapping = infer_chain_mapping(native_view, prediction_view)
    global_lddt, local_lddt = lDDTScorer(native_view).lDDT(
        prediction_view,
        chain_mapping=chain_mapping,
        check_resnames=False,
    )

    residues = list(prediction_view.residues)

    summary = {
        "atom_selection": atom_selection,
        "global_lDDT": fmt(global_lddt),
        "average_local_lDDT": fmt(mean_defined(local_lddt)),
        "n_local_scores": sum(v is not None for v in local_lddt),
        "n_native_atoms": native_view.atom_count,
        "n_prediction_atoms": prediction_view.atom_count,
    }

    local_rows = []
    for index, score in enumerate(local_lddt, start=1):
        residue = residues[index - 1] if index - 1 < len(residues) else None

        local_rows.append({
            "atom_selection": atom_selection,
            "atom_name": atom_name,
            "local_index": index,
            "chain_id": residue.GetChain().GetName() if residue else "",
            "residue_number": residue_number(residue) if residue else "",
            "residue_name": residue.GetName() if residue else "",
            "local_lDDT": fmt(score),
        })

    return summary, local_rows


def collect_tasks(prediction_sets, native_dir=None, native_level=None, levels=None):
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
                tasks.append((model_name, level, None, prediction_path, None))
                continue

            matched_natives = native_index.get((level, target_id), [])

            if not matched_natives:
                tasks.append((model_name, level, target_id, prediction_path, None))
            else:
                for native_path in matched_natives:
                    tasks.append((model_name, level, target_id, prediction_path, native_path))

    return tasks


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
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    summary_run_csv = temp_csv_path(SUMMARY_CSV)
    local_run_csv = temp_csv_path(LOCAL_CSV)

    prediction_sets = {model_name: PREDICTION_SETS[model_name] for model_name in args.models}
    tasks = collect_tasks(
        prediction_sets,
        native_dir=args.native_dir,
        native_level=args.native_level,
        levels=args.levels,
    )
    print(f"Processing models: {', '.join(args.models)}")
    if args.levels is not None:
        print(f"Processing levels: {', '.join(str(l) for l in args.levels)}")
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
                summary_writer.writerow({
                    "Model": model_name,
                    "level": level,
                    "target_id": target_id,
                    "native_path": "",
                    "prediction_path": str(abs_prediction_path),
                    "atom_selection": "",
                    "global_lDDT": "",
                    "average_local_lDDT": "",
                    "n_local_scores": "",
                    "n_native_atoms": "",
                    "n_prediction_atoms": "",
                    "error": "No matching native found",
                })
                continue

            abs_native_path = native_path.resolve()

            try:
                if str(abs_native_path) not in native_cache:
                    native_cache[str(abs_native_path)] = load_structure(abs_native_path)

                native = native_cache[str(abs_native_path)]
                prediction = load_structure(abs_prediction_path)

                for atom_selection, atom_name, atom_query in SELECTIONS:
                    try:
                        summary, local_rows = score_selection(
                            native,
                            prediction,
                            atom_selection,
                            atom_name,
                            atom_query,
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
                            "global_lDDT": "",
                            "average_local_lDDT": "",
                            "n_local_scores": "",
                            "n_native_atoms": "",
                            "n_prediction_atoms": "",
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
                    "global_lDDT": "",
                    "average_local_lDDT": "",
                    "n_local_scores": "",
                    "n_native_atoms": "",
                    "n_prediction_atoms": "",
                    "error": str(exc),
                })

    merge_csv(
        summary_run_csv,
        SUMMARY_CSV,
        SUMMARY_FIELDS,
        SUMMARY_KEY_FIELDS,
        replace_error_rows=True,
    )
    merge_csv(local_run_csv, LOCAL_CSV, LOCAL_FIELDS, LOCAL_KEY_FIELDS)

    print(f"Saved summary CSV to: {SUMMARY_CSV}")
    print(f"Saved local CSV to: {LOCAL_CSV}")


if __name__ == "__main__":
    main()
