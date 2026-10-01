#!/usr/bin/env python3
from __future__ import annotations

import json
import logging
import os
import re
import shlex
import urllib.error
import urllib.request
from argparse import Namespace
from pathlib import Path


SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
FASTA_URL = "https://www.rcsb.org/fasta/entry/{pdb_id}"
CIF_URL = "https://files.rcsb.org/download/{pdb_id}.cif"
USER_AGENT = "rna-project-filtered-downloader/1.0"

DATASET_MODE = os.environ.get("RCSB_RNA_MODE", "level1").lower()
if DATASET_MODE not in {"strict", "level1", "level2", "level3"}:
    raise ValueError("RCSB_RNA_MODE must be 'strict', 'level1', 'level2', or 'level3'")

START = os.environ.get("RCSB_RNA_START", "1990-01-01")
DEFAULT_END = {
    "level2": "2026-06-23",
    "level3": "2026-06-26",
}.get(DATASET_MODE, "2030-01-01")
END = os.environ.get(
    "RCSB_RNA_END",
    DEFAULT_END,
)

MAIN_MIN = int(os.environ.get(
    "RCSB_RNA_MAIN_MIN",
    "1" if DATASET_MODE == "level3" else "30",
))
MAIN_MAX = int(os.environ.get(
    "RCSB_RNA_MAIN_MAX",
    "1000" if DATASET_MODE == "level3" else "500",
))
EXT_MIN = int(os.environ.get(
    "RCSB_RNA_EXTENSION_MIN",
    "0" if DATASET_MODE == "level3" else "500",
))
EXT_MAX = int(os.environ.get(
    "RCSB_RNA_EXTENSION_MAX",
    "0" if DATASET_MODE == "level3" else "1000",
))


def optional_float_env(name: str, default: float | None) -> float | None:
    value = os.environ.get(name)
    if value is None:
        return default
    value = value.strip().lower()
    if value in {"", "none", "all", "na", "n/a", "off"}:
        return None
    return float(value)


def bool_env(name: str, default: bool) -> bool:
    value = os.environ.get(name)
    if value is None:
        return default
    value = value.strip().lower()
    if value in {"1", "true", "yes", "y", "on"}:
        return True
    if value in {"0", "false", "no", "n", "off"}:
        return False
    raise ValueError(f"{name} must be a boolean value")


MAX_RESOLUTION = optional_float_env(
    "RCSB_RNA_MAX_RESOLUTION",
    None if DATASET_MODE in {"level2", "level3"} else 4.0,
)
REQUIRE_SINGLE_RNA_CHAIN = bool_env(
    "RCSB_RNA_REQUIRE_SINGLE_RNA_CHAIN",
    False,
)

LIMIT = int(os.environ.get("RCSB_RNA_LIMIT", "0"))
TIMEOUT = int(os.environ.get("RCSB_RNA_TIMEOUT", "45"))

DEFAULT_OUTDIR = {
    "strict": "/home/schen3/CASP17/DProQA_RNA/data/strict",
    "level1": "/home/schen3/CASP17/DProQA_RNA/data",
    "level2": "/home/schen3/CASP17/DProQA_RNA/data",
    "level3": "/home/schen3/CASP17/DProQA_RNA/data/rna_protein",
}[DATASET_MODE]

OUTDIR = Path(os.environ.get(
    "RCSB_RNA_OUTDIR",
    DEFAULT_OUTDIR,
))

CANONICAL_RNA = set("ACGU")
WATER_COMP_IDS = {"HOH", "WAT", "DOD"}
ION_COMP_IDS = {
    "AG",
    "AL",
    "AR",
    "AU",
    "BA",
    "BR",
    "CA",
    "CD",
    "CL",
    "CO",
    "CR",
    "CS",
    "CU",
    "EU",
    "F",
    "FE",
    "GA",
    "GD",
    "HG",
    "IOD",
    "IR",
    "K",
    "LA",
    "LI",
    "MG",
    "MN",
    "NA",
    "NI",
    "PB",
    "PT",
    "RB",
    "RU",
    "SE",
    "SM",
    "SR",
    "TB",
    "TL",
    "V",
    "Y",
    "YB",
    "ZN",
}
METAL_COMP_IDS = {
    "AC",
    "AG",
    "AL",
    "AM",
    "AU",
    "BA",
    "BE",
    "BI",
    "BK",
    "CA",
    "CD",
    "CE",
    "CF",
    "CM",
    "CO",
    "CR",
    "CS",
    "CU",
    "DY",
    "ER",
    "ES",
    "EU",
    "FE",
    "FM",
    "GA",
    "GD",
    "HF",
    "HG",
    "HO",
    "IN",
    "IR",
    "K",
    "LA",
    "LI",
    "LR",
    "LU",
    "MD",
    "MG",
    "MN",
    "MO",
    "NA",
    "NB",
    "ND",
    "NI",
    "NO",
    "NP",
    "OS",
    "PA",
    "PB",
    "PD",
    "PM",
    "PR",
    "PT",
    "PU",
    "RB",
    "RE",
    "RH",
    "RU",
    "SC",
    "SM",
    "SN",
    "SR",
    "TA",
    "TB",
    "TC",
    "TH",
    "TI",
    "TL",
    "TM",
    "U",
    "V",
    "W",
    "Y",
    "YB",
    "ZN",
    "ZR",
}


def fetch_text(url: str) -> str | None:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(req, timeout=TIMEOUT) as r:
            return r.read().decode("utf-8")
    except (urllib.error.URLError, TimeoutError):
        return None


def post_json(url: str, payload: dict) -> dict:
    req = urllib.request.Request(
        url,
        data=json.dumps(payload).encode("utf-8"),
        headers={
            "Accept": "application/json",
            "Content-Type": "application/json",
            "User-Agent": USER_AGENT,
        },
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=TIMEOUT) as r:
        if r.status == 204:
            return {}
        return json.loads(r.read().decode("utf-8"))


def terminal(attribute: str, operator: str, value) -> dict:
    return {
        "type": "terminal",
        "service": "text",
        "parameters": {
            "attribute": attribute,
            "operator": operator,
            "value": value,
        },
    }


def search_payload() -> dict:
    nodes = [
        terminal(
            "rcsb_accession_info.initial_release_date",
            "range",
            {
                "from": START,
                "to": END,
                "include_lower": True,
                "include_upper": True,
            },
        ),
        terminal(
            "rcsb_entry_info.structure_determination_methodology",
            "exact_match",
            "experimental",
        ),
    ]

    nodes.extend([
        terminal("rcsb_entry_info.polymer_composition", "exact_match", "RNA"),
        terminal("rcsb_entry_info.na_polymer_entity_types", "exact_match", "RNA (only)"),
        terminal("rcsb_entry_info.polymer_entity_count_RNA", "equals", 1),
        terminal("rcsb_entry_info.polymer_entity_count_DNA", "equals", 0),
        terminal("rcsb_entry_info.polymer_entity_count_protein", "equals", 0),
        terminal("rcsb_entry_info.polymer_entity_count_nucleic_acid_hybrid", "equals", 0),
    ])

    if MAX_RESOLUTION is not None:
        nodes.append(
            terminal(
                "rcsb_entry_info.resolution_combined",
                "less_or_equal",
                MAX_RESOLUTION,
            )
        )

    if DATASET_MODE == "strict":
        nodes.extend([
            terminal("rcsb_entry_info.nonpolymer_entity_count", "equals", 0),
            terminal("rcsb_entry_info.branched_entity_count", "equals", 0),
        ])

    return {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": nodes,
        },
        "return_type": "entry",
        "request_options": {
            "return_all_hits": True,
            "results_verbosity": "compact",
        },
    }


def search_ids() -> list[str]:
    response = post_json(SEARCH_URL, search_payload())
    ids = []
    for hit in response.get("result_set") or []:
        ids.append(hit.upper() if isinstance(hit, str) else hit["identifier"].upper())
    return sorted(set(ids))


def fasta_records(text: str) -> list[tuple[str, str]]:
    records = []
    header = None
    seq = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header:
                records.append((header, "".join(seq).upper()))
            header = line
            seq = []
        else:
            seq.append(re.sub(r"\s+", "", line))
    if header:
        records.append((header, "".join(seq).upper()))
    return records


def fasta_chain_count(header: str) -> int | None:
    parts = header.split("|")
    if len(parts) < 2:
        return None
    chain_field = re.sub(r"\[auth\s+[^\]]*\]", "", parts[1]).strip()
    m = re.match(r"Chains?\s+(.+)", chain_field)
    if not m:
        return None
    chains = [c.strip() for c in m.group(1).split(",") if c.strip()]
    return len(chains) if chains else None


def mmcif_tokens(lines: list[str]) -> list[str]:
    tokens = []
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith(";"):
            block = [line[1:]]
            i += 1
            while i < len(lines) and not lines[i].startswith(";"):
                block.append(lines[i])
                i += 1
            tokens.append("\n".join(block).strip())
            if i < len(lines):
                i += 1
        else:
            tokens.extend(shlex.split(line, comments=True, posix=True))
            i += 1
    return tokens


def mmcif_loop_rows(cif: str, category: str) -> list[dict[str, str]]:
    lines = cif.splitlines()
    prefix = f"_{category}."
    rows = []
    i = 0

    while i < len(lines):
        if lines[i].strip() != "loop_":
            i += 1
            continue

        i += 1
        fields = []
        while i < len(lines) and lines[i].strip().startswith("_"):
            fields.append(lines[i].strip().split()[0])
            i += 1

        data_lines = []
        while i < len(lines):
            s = lines[i].strip()
            if not s:
                i += 1
                continue
            if s == "#" or s == "loop_" or s.startswith(("data_", "save_", "_")):
                break
            data_lines.append(lines[i])
            i += 1

        if not any(f.startswith(prefix) for f in fields):
            continue

        tokens = mmcif_tokens(data_lines)
        width = len(fields)
        for start in range(0, len(tokens), width):
            chunk = tokens[start:start + width]
            if len(chunk) == width:
                rows.append(dict(zip(fields, chunk)))

    return rows


def mmcif_scalar_row(cif: str, category: str) -> dict[str, str] | None:
    lines = cif.splitlines()
    prefix = f"_{category}."
    row = {}
    i = 0

    while i < len(lines):
        s = lines[i].strip()
        if not s.startswith(prefix):
            i += 1
            continue

        parts = s.split(None, 1)
        key = parts[0]

        if len(parts) == 2:
            vals = shlex.split(parts[1], comments=True, posix=True)
            if vals:
                row[key] = vals[0]
            i += 1
            continue

        if i + 1 < len(lines) and lines[i + 1].startswith(";"):
            i += 2
            block = []
            while i < len(lines) and not lines[i].startswith(";"):
                block.append(lines[i])
                i += 1
            row[key] = "\n".join(block).strip()
            if i < len(lines):
                i += 1
            continue

        i += 1

    return row or None


def mmcif_rows(cif: str, category: str) -> list[dict[str, str]]:
    rows = mmcif_loop_rows(cif, category)
    if rows:
        return rows
    row = mmcif_scalar_row(cif, category)
    return [row] if row else []


def clean(value: str | None) -> str | None:
    if value is None:
        return None
    value = value.strip()
    return None if value in {".", "?"} else value


def clean_int(value: str | None) -> int | None:
    value = clean(value)
    if not value:
        return None
    try:
        return int(value)
    except ValueError:
        return None


def low(value: str | None) -> str | None:
    value = clean(value)
    return value.lower() if value else None


def first_clean_value(rows: list[dict[str, str]], key: str) -> str | None:
    for row in rows:
        value = clean(row.get(key))
        if value:
            return value
    return None


def cif_metadata(cif: str) -> dict[str, str]:
    exptl_rows = mmcif_rows(cif, "exptl")
    refine_rows = mmcif_rows(cif, "refine")
    em_3d_rows = mmcif_rows(cif, "em_3d_reconstruction")
    accession_rows = mmcif_rows(cif, "pdbx_database_status")

    method = first_clean_value(exptl_rows, "_exptl.method") or "NA"
    resolution = (
        first_clean_value(refine_rows, "_refine.ls_d_res_high")
        or first_clean_value(em_3d_rows, "_em_3d_reconstruction.resolution")
        or "NA"
    )
    release_date = first_clean_value(
        accession_rows,
        "_pdbx_database_status.recvd_initial_deposition_date",
    ) or "NA"

    return {
        "method": method,
        "resolution": resolution,
        "release_date": release_date,
    }


def has_ion_name(name: str | None) -> bool:
    name = (name or "").lower()
    return " ion" in name or name.endswith("ion")


def nonpolymer_metadata(cif: str) -> dict[str, str]:
    rows = mmcif_rows(cif, "pdbx_entity_nonpoly")
    branch_rows = mmcif_rows(cif, "pdbx_entity_branch")
    chem_rows = mmcif_rows(cif, "chem_comp")
    chem_by_id = {
        (clean(row.get("_chem_comp.id")) or "").upper(): row
        for row in chem_rows
        if clean(row.get("_chem_comp.id"))
    }

    ions = set()
    ligands = set()

    for row in rows:
        comp_id = (clean(row.get("_pdbx_entity_nonpoly.comp_id")) or "").upper()
        if not comp_id:
            continue

        chem = chem_by_id.get(comp_id, {})
        name = (
            clean(row.get("_pdbx_entity_nonpoly.name"))
            or clean(chem.get("_chem_comp.name"))
            or ""
        )

        if comp_id in WATER_COMP_IDS:
            continue
        if comp_id in ION_COMP_IDS or comp_id in METAL_COMP_IDS or has_ion_name(name):
            ions.add(comp_id)
        else:
            ligands.add(comp_id)

    return {
        "has_ions": "yes" if ions else "no",
        "has_ligands": "yes" if ligands else "no",
        "has_nonpolymer_entities": "yes" if rows else "no",
        "has_branched_entities": "yes" if branch_rows else "no",
    }


def is_allowed_metal_atom(tokens: list[str], field_index: dict[str, int]) -> bool:
    comp_ids = []
    for field in ("_atom_site.label_comp_id", "_atom_site.auth_comp_id"):
        idx = field_index.get(field)
        if idx is not None and idx < len(tokens):
            comp_id = clean(tokens[idx])
            if comp_id:
                comp_ids.append(comp_id.upper())

    if comp_ids:
        return any(comp_id in METAL_COMP_IDS for comp_id in comp_ids)

    type_idx = field_index.get("_atom_site.type_symbol")
    if type_idx is None or type_idx >= len(tokens):
        return False

    type_symbol = clean(tokens[type_idx])
    if not type_symbol:
        return False

    return type_symbol.upper() in METAL_COMP_IDS


def strip_non_rna_atoms(cif: str, entity_id: str, chain_id: str) -> tuple[str, int, int, int]:
    lines = cif.splitlines()
    out = []
    kept_rna_atoms = 0
    kept_metal_atoms = 0
    stripped_atoms = 0
    i = 0

    while i < len(lines):
        if lines[i].strip() != "loop_":
            out.append(lines[i])
            i += 1
            continue

        loop_start = i
        i += 1
        fields = []
        while i < len(lines) and lines[i].strip().startswith("_"):
            fields.append(lines[i].strip().split()[0])
            i += 1

        if not fields or not all(field.startswith("_atom_site.") for field in fields):
            out.extend(lines[loop_start:i])
            continue

        out.extend(lines[loop_start:i])
        field_index = {field: idx for idx, field in enumerate(fields)}
        entity_idx = field_index.get("_atom_site.label_entity_id")
        label_asym_idx = field_index.get("_atom_site.label_asym_id")
        auth_asym_idx = field_index.get("_atom_site.auth_asym_id")
        width = len(fields)

        row_lines = []
        while i < len(lines):
            s = lines[i].strip()
            if not row_lines and (
                s == "#"
                or s == "loop_"
                or s.startswith(("data_", "save_", "_"))
            ):
                break
            if not s and not row_lines:
                out.append(lines[i])
                i += 1
                continue

            row_lines.append(lines[i])
            tokens = mmcif_tokens(row_lines)
            if len(tokens) < width:
                i += 1
                continue

            entity_matches = (
                entity_idx is not None
                and entity_idx < len(tokens)
                and clean(tokens[entity_idx]) == entity_id
            )
            label_asym = (
                clean(tokens[label_asym_idx])
                if label_asym_idx is not None and label_asym_idx < len(tokens)
                else None
            )
            auth_asym = (
                clean(tokens[auth_asym_idx])
                if auth_asym_idx is not None and auth_asym_idx < len(tokens)
                else None
            )
            chain_matches = chain_id in {label_asym, auth_asym}

            if label_asym or auth_asym:
                keep_rna = chain_matches and (entity_idx is None or entity_matches)
            else:
                keep_rna = entity_matches

            keep_metal = not keep_rna and is_allowed_metal_atom(tokens, field_index)
            if keep_rna or keep_metal:
                out.extend(row_lines)
                if keep_rna:
                    kept_rna_atoms += 1
                else:
                    kept_metal_atoms += 1
            else:
                stripped_atoms += 1

            row_lines = []
            i += 1

        if row_lines:
            stripped_atoms += 1

    return (
        "\n".join(out) + ("\n" if cif.endswith("\n") else ""),
        kept_rna_atoms,
        kept_metal_atoms,
        stripped_atoms,
    )


def split_chain_ids(value: str | None) -> list[str]:
    return [
        x.strip()
        for x in (clean(value) or "").split(",")
        if x.strip()
    ]


def rna_chain_instances(
    rna_rows: list[dict[str, str]],
    struct_asym_rows: list[dict[str, str]],
) -> list[tuple[str, str]]:
    asym_by_entity: dict[str, list[str]] = {}
    for row in struct_asym_rows:
        entity_id = clean(row.get("_struct_asym.entity_id"))
        asym_id = clean(row.get("_struct_asym.id"))
        if not entity_id or not asym_id:
            continue
        asym_by_entity.setdefault(entity_id, []).append(asym_id)

    instances = []
    for row in rna_rows:
        entity_id = clean(row.get("_entity_poly.entity_id"))
        if not entity_id:
            continue
        chain_ids = asym_by_entity.get(entity_id) or split_chain_ids(
            row.get("_entity_poly.pdbx_strand_id")
        )
        for chain_id in chain_ids or [entity_id]:
            instances.append((entity_id, chain_id))

    return instances


def chain_coordinate_seq_nums(
    atom_rows: list[dict[str, str]],
    entity_id: str,
    chain_id: str,
) -> set[int]:
    seq_nums = set()

    for row in atom_rows:
        label_entity = clean(row.get("_atom_site.label_entity_id"))
        label_asym = clean(row.get("_atom_site.label_asym_id"))
        auth_asym = clean(row.get("_atom_site.auth_asym_id"))

        if label_asym or auth_asym:
            if chain_id not in {label_asym, auth_asym}:
                continue
            if label_entity and label_entity != entity_id:
                continue
        elif label_entity != entity_id:
            continue

        seq_id = clean_int(row.get("_atom_site.label_seq_id"))
        if seq_id is not None:
            seq_nums.add(seq_id)

    return seq_nums


def chain_auth_ids(
    atom_rows: list[dict[str, str]],
    entity_id: str,
    chain_id: str,
) -> list[str]:
    auth_ids = []

    for row in atom_rows:
        label_entity = clean(row.get("_atom_site.label_entity_id"))
        label_asym = clean(row.get("_atom_site.label_asym_id"))
        auth_asym = clean(row.get("_atom_site.auth_asym_id"))

        if label_entity and label_entity != entity_id:
            continue
        if chain_id not in {label_asym, auth_asym}:
            continue
        if auth_asym and auth_asym not in auth_ids:
            auth_ids.append(auth_asym)

    return auth_ids


def analyze_cif(cif: str) -> tuple[list[dict[str, str]], list[str]]:
    entity_rows = mmcif_rows(cif, "entity")
    if entity_rows:
        non_polymer = [
            clean(r.get("_entity.id")) or "?"
            for r in entity_rows
            if low(r.get("_entity.type")) != "polymer"
        ]
        branched = [
            clean(r.get("_entity.id")) or "?"
            for r in entity_rows
            if low(r.get("_entity.type")) == "branched"
        ]
        if DATASET_MODE == "strict" and non_polymer:
            return [], ["has_nonpolymer_or_branched_entity"]
        if DATASET_MODE == "strict" and branched:
            return [], ["has_branched_entity"]

    if DATASET_MODE == "strict" and mmcif_rows(cif, "pdbx_entity_nonpoly"):
        return [], ["has_nonpolymer_entity"]
    if DATASET_MODE == "strict" and mmcif_rows(cif, "pdbx_entity_branch"):
        return [], ["has_branched_entity"]

    poly_rows = mmcif_rows(cif, "entity_poly")
    rna_rows = [
        r for r in poly_rows
        if low(r.get("_entity_poly.type")) == "polyribonucleotide"
    ]

    if not rna_rows:
        return [], ["missing_rna_polymer_entity"]

    struct_asym_rows = mmcif_rows(cif, "struct_asym")
    rna_instances = rna_chain_instances(rna_rows, struct_asym_rows)
    if REQUIRE_SINGLE_RNA_CHAIN and len(rna_instances) != 1:
        return [], [f"not_single_rna_chain_{len(rna_instances)}"]

    if len(poly_rows) != 1 or len(rna_rows) != 1:
        return [], ["not_exactly_one_rna_polymer_entity"]

    atom_rows = mmcif_rows(cif, "atom_site")
    candidates = []
    reasons = []

    for rna in rna_rows:
        entity_id = clean(rna.get("_entity_poly.entity_id"))
        if not entity_id:
            reasons.append("missing_rna_entity_id")
            continue

        if low(rna.get("_entity_poly.nstd_linkage")) == "yes":
            reasons.append(f"{entity_id}:nonstandard_rna_linkage")
            continue
        if low(rna.get("_entity_poly.nstd_monomer")) == "yes":
            reasons.append(f"{entity_id}:modified_rna_monomer")
            continue

        strand_ids = split_chain_ids(rna.get("_entity_poly.pdbx_strand_id"))
        if len(strand_ids) > 1:
            return [], ["multiple_rna_chains"]

        asym_ids = sorted({
            clean(r.get("_struct_asym.id"))
            for r in struct_asym_rows
            if clean(r.get("_struct_asym.entity_id")) == entity_id
        } - {None})

        if asym_ids and len(asym_ids) != 1:
            return [], ["multiple_rna_chain_instances"]

        seq_rows = [
            r for r in mmcif_rows(cif, "entity_poly_seq")
            if clean(r.get("_entity_poly_seq.entity_id")) == entity_id
        ]
        if not seq_rows:
            reasons.append(f"{entity_id}:missing_entity_poly_seq")
            continue

        numbered_seq_rows = []
        for row in seq_rows:
            seq_num = clean_int(row.get("_entity_poly_seq.num"))
            if seq_num is None:
                reasons.append(f"{entity_id}:bad_entity_poly_seq_number")
                numbered_seq_rows = []
                break
            numbered_seq_rows.append((seq_num, row))
        if not numbered_seq_rows:
            continue

        numbered_seq_rows.sort(key=lambda item: item[0])
        seq_nums = [seq_num for seq_num, _ in numbered_seq_rows]
        monomers = [
            (clean(r.get("_entity_poly_seq.mon_id")) or "").upper()
            for _, r in numbered_seq_rows
        ]

        invalid = sorted(set(monomers) - CANONICAL_RNA)
        if invalid:
            reasons.append(f"{entity_id}:noncanonical_monomer:{','.join(invalid)}")
            continue

        chain_ids = asym_ids or strand_ids or ["A"]
        chain_ids = chain_ids[:1]

        expected_seq_nums = set(seq_nums)
        for chain_id in chain_ids:
            observed_seq_nums = chain_coordinate_seq_nums(atom_rows, entity_id, chain_id)
            if not observed_seq_nums:
                reasons.append(f"{entity_id}:{chain_id}:missing_coordinate_atoms")
                continue
            if not expected_seq_nums <= observed_seq_nums:
                reasons.append(f"{entity_id}:{chain_id}:chain_discontinuity")
                continue

            candidates.append({
                "chain_id": chain_id,
                "auth_chain_id": ",".join(
                    chain_auth_ids(atom_rows, entity_id, chain_id) or [chain_id]
                ),
                "entity_id": entity_id,
                "sequence": "".join(monomers),
            })

    if not candidates:
        return [], reasons or ["no_valid_rna_chain"]

    if len(candidates) != 1:
        return [], ["not_exactly_one_valid_rna_chain"]

    return candidates, []


def wrap_fasta(header: str, seq: str, width: int = 80) -> str:
    lines = [header]
    for i in range(0, len(seq), width):
        lines.append(seq[i:i + width])
    return "\n".join(lines) + "\n"


def buckets(length: int) -> list[str]:
    out = []
    if MAIN_MIN <= length <= MAIN_MAX:
        out.append("main")
    if EXT_MIN <= length <= EXT_MAX:
        out.append("extension")
    return out


def record_id(pdb_id: str, chain_id: str) -> str:
    return pdb_id


def compact_reasons(reasons: list[str], limit: int = 5) -> str:
    if len(reasons) <= limit:
        return ";".join(reasons)
    return ";".join(reasons[:limit]) + f";+{len(reasons) - limit}_more"


def make_dirs() -> dict[str, Path]:
    dirs = {
        "all": OUTDIR / "all",
        "main": OUTDIR / f"main_{MAIN_MIN}_{MAIN_MAX}",
        "extension": OUTDIR / f"extension_{EXT_MIN}_{EXT_MAX}",
    }
    OUTDIR.mkdir(parents=True, exist_ok=True)
    for d in dirs.values():
        (d / "fasta_rna").mkdir(parents=True, exist_ok=True)
        (d / "cif_rna").mkdir(parents=True, exist_ok=True)
        (d / "cif_native").mkdir(parents=True, exist_ok=True)
    return dirs


def write_entry(base: Path, pdb_id: str, fasta: str, rna_cif: str, native_cif: str) -> None:
    (base / "fasta_rna" / f"{pdb_id}.fasta").write_text(fasta)
    (base / "cif_rna" / f"{pdb_id}.cif").write_text(rna_cif)
    (base / "cif_native" / f"{pdb_id}.cif").write_text(native_cif)


def write_af3_input(base: Path, pdb_id: str, chain_id: str, seq: str) -> None:
    af3_dir = base / "af3_inputs"
    af3_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "name": pdb_id,
        "sequences": [
            {
                "rna": {
                    "id": chain_id,
                    "sequence": seq,
                }
            }
        ],
        "modelSeeds": [1],
        "dialect": "alphafold3",
        "version": 1,
    }
    (af3_dir / f"{pdb_id}.json").write_text(json.dumps(payload, indent=2) + "\n")


def optional_path_env(name: str, default: str | None) -> str | None:
    value = os.environ.get(name)
    if value is None:
        return default
    value = value.strip()
    if value.lower() in {"", "none", "off"}:
        return None
    return value


def run_level12_mode() -> int:
    import build_single_chain_rna_datasets as level12

    args = Namespace(
        out_dir=str(OUTDIR),
        start_date=START,
        end_date=END,
        max_length=int(os.environ.get("RCSB_RNA_MAX_LENGTH", "1000")),
        max_resolution=MAX_RESOLUTION,
        limit=LIMIT,
        timeout=TIMEOUT,
        progress_every=int(os.environ.get("RCSB_RNA_PROGRESS_EVERY", "100")),
        refresh_search=bool_env("RCSB_RNA_REFRESH_SEARCH", False),
        refresh_metadata=bool_env("RCSB_RNA_REFRESH_METADATA", False),
        redownload=bool_env("RCSB_RNA_REDOWNLOAD", False),
        sanity_only=bool_env("RCSB_RNA_SANITY_ONLY", False),
        skip_sanity_checks=bool_env("RCSB_RNA_SKIP_SANITY_CHECKS", False),
        run_known_example_checks=bool_env("RCSB_RNA_RUN_KNOWN_EXAMPLE_CHECKS", False),
        verbose=bool_env("RCSB_RNA_VERBOSE", False),
    )

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
    )

    print(f"Dataset mode: {DATASET_MODE}")
    print("Level 1/2 definition: clean single-chain RNA targets")
    print("Level 1 targets have no metal ions; Level 2 targets have metal ions")
    print(f"Date range: {args.start_date} to {args.end_date}")
    print(f"Max RNA length: {args.max_length} nt")
    if args.max_resolution is None:
        print("Max resolution: disabled")
    else:
        print(f"Max resolution: {args.max_resolution} Å")
    print(f"Output: {args.out_dir}")
    print("Level 1 FASTA: fasta/level1")
    print("Level 2 FASTA: fasta/level2")
    print("Level 1 CSV: level1_single_chain_rna.csv")
    print("Level 2 CSV: level2_single_chain_rna_metal.csv")

    if not args.skip_sanity_checks:
        level12.run_sanity_checks()
    if args.run_known_example_checks:
        level12.run_known_example_checks(args)
    if args.sanity_only:
        return 0

    summary = level12.run_pipeline(args)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


def run_level3_mode() -> int:
    import build_single_rna_single_protein_targets as level3

    args = Namespace(
        out_dir=str(OUTDIR),
        start_date=START,
        end_date=END,
        max_rna_length=int(os.environ.get("RCSB_RNA_MAX_RNA_LENGTH", str(MAIN_MAX))),
        max_resolution=MAX_RESOLUTION,
        broad_search=bool_env("RCSB_RNA_LEVEL3_BROAD_SEARCH", False),
        contact_cutoff=float(os.environ.get("RCSB_RNA_CONTACT_CUTOFF", "6.0")),
        require_contact=bool_env("RCSB_RNA_REQUIRE_CONTACT", True),
        allow_organic_ligand=bool_env("RCSB_RNA_ALLOW_ORGANIC_LIGAND", False),
        limit=LIMIT,
        timeout=TIMEOUT,
        progress_every=int(os.environ.get("RCSB_RNA_PROGRESS_EVERY", "100")),
        shared_cif_cache=optional_path_env("RCSB_RNA_SHARED_CIF_CACHE", "data/cache/mmcif"),
        refresh_search=bool_env("RCSB_RNA_REFRESH_SEARCH", False),
        refresh_metadata=bool_env("RCSB_RNA_REFRESH_METADATA", False),
        redownload=bool_env("RCSB_RNA_REDOWNLOAD", False),
        sanity_only=bool_env("RCSB_RNA_SANITY_ONLY", False),
        skip_sanity_checks=bool_env("RCSB_RNA_SKIP_SANITY_CHECKS", False),
        verbose=bool_env("RCSB_RNA_VERBOSE", False),
    )

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
    )

    print("Dataset mode: level3")
    print("Level 3 definition: one RNA polymer chain plus one protein polymer chain")
    print(f"Date range: {args.start_date} to {args.end_date}")
    print(f"Max RNA length: {args.max_rna_length} nt")
    if args.max_resolution is None:
        print("Max resolution: disabled")
    else:
        print(f"Max resolution: {args.max_resolution} Å")
    print(f"Require RNA-protein contact: {'yes' if args.require_contact else 'no'}")
    print(f"Allow organic ligand: {'yes' if args.allow_organic_ligand else 'no'}")
    print(f"Output: {args.out_dir}")
    print("FASTA: rna, protein, and combined complex FASTA files")

    if not args.skip_sanity_checks:
        level3.run_sanity_checks()
    if args.sanity_only:
        return 0

    summary = level3.run_pipeline(args)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


def main() -> int:
    if DATASET_MODE in {"level1", "level2"}:
        return run_level12_mode()
    if DATASET_MODE == "level3":
        return run_level3_mode()

    dirs = make_dirs()

    print(f"Dataset mode: {DATASET_MODE}")
    print(f"Date range: {START} to {END}")
    print(f"Main: {MAIN_MIN}-{MAIN_MAX} nt")
    if EXT_MIN <= EXT_MAX and EXT_MAX > 0:
        print(f"Extension: {EXT_MIN}-{EXT_MAX} nt")
    else:
        print("Extension: disabled")
    if MAX_RESOLUTION is None:
        print("Max resolution: disabled")
    else:
        print(f"Max resolution: {MAX_RESOLUTION} Å")
    print(f"Require single RNA chain: {'yes' if REQUIRE_SINGLE_RNA_CHAIN else 'no'}")
    print(f"Output: {OUTDIR}")
    print("CIF: RNA plus metal atom_site in cif_rna; full native mmCIF in cif_native")

    candidate_ids = search_ids()
    print("Search candidates:", len(candidate_ids))

    selected = []
    selected_by_bucket = {"main": [], "extension": []}
    skipped = []

    for pdb_id in candidate_ids:
        if LIMIT and len(selected) >= LIMIT:
            break

        cif_text = fetch_text(CIF_URL.format(pdb_id=pdb_id))

        if cif_text is None:
            skipped.append((pdb_id, "missing_mmcif"))
            print(f"Skipping {pdb_id}: missing mmCIF")
            continue

        infos, reasons = analyze_cif(cif_text)
        if not infos:
            reason = compact_reasons(reasons)
            skipped.append((pdb_id, reason))
            print(f"Skipping {pdb_id}: {reason}")
            continue

        fasta_text = fetch_text(FASTA_URL.format(pdb_id=pdb_id))
        if fasta_text is None:
            skipped.append((pdb_id, "missing_fasta"))
            print(f"Skipping {pdb_id}: missing FASTA")
            continue

        records = fasta_records(fasta_text)
        if len(records) != 1:
            skipped.append((pdb_id, "fasta_not_single_record"))
            print(f"Skipping {pdb_id}: FASTA has {len(records)} records")
            continue

        header, fasta_seq = records[0]
        chain_count = fasta_chain_count(header)
        if chain_count is not None and chain_count != 1:
            skipped.append((pdb_id, "fasta_not_single_chain"))
            print(f"Skipping {pdb_id}: FASTA has multiple chains")
            continue

        if fasta_seq != infos[0]["sequence"]:
            skipped.append((pdb_id, "fasta_cif_sequence_mismatch"))
            print(f"Skipping {pdb_id}: FASTA/mmCIF sequence mismatch")
            continue

        meta = cif_metadata(cif_text)
        nonpoly_meta = nonpolymer_metadata(cif_text)

        for info in infos:
            if LIMIT and len(selected) >= LIMIT:
                break

            seq = info["sequence"]
            chain_id = info["chain_id"]
            auth_chain_id = info["auth_chain_id"]
            entity_id = info["entity_id"]
            rid = record_id(pdb_id, chain_id)

            entry_buckets = buckets(len(seq))
            if not entry_buckets:
                skipped.append((rid, f"length_out_of_range_{len(seq)}"))
                print(f"Skipping {rid}: length {len(seq)} outside requested ranges")
                continue

            rna_cif, rna_atom_count, metal_atom_count, stripped_atom_count = strip_non_rna_atoms(
                cif_text,
                entity_id,
                chain_id,
            )
            fasta_out = wrap_fasta(
                (
                    f">{rid}|pdb {pdb_id}|label_chain {chain_id}|"
                    f"auth_chain {auth_chain_id}|length {len(seq)}"
                ),
                seq,
            )
            write_entry(dirs["all"], rid, fasta_out, rna_cif, cif_text)
            write_af3_input(dirs["all"], rid, chain_id, seq)

            for bucket in entry_buckets:
                selected_by_bucket[bucket].append(rid)
                write_entry(dirs[bucket], rid, fasta_out, rna_cif, cif_text)
                write_af3_input(dirs[bucket], rid, chain_id, seq)

            selected.append((
                rid,
                pdb_id,
                len(seq),
                chain_id,
                auth_chain_id,
                ",".join(entry_buckets),
                meta["method"],
                meta["resolution"],
                meta["release_date"],
                nonpoly_meta["has_ions"],
                nonpoly_meta["has_ligands"],
                nonpoly_meta["has_nonpolymer_entities"],
                nonpoly_meta["has_branched_entities"],
                rna_atom_count,
                metal_atom_count,
                stripped_atom_count,
            ))
            print(
                f"Selected {rid}: {len(seq)} nt, PDB {pdb_id}, chain {chain_id}, "
                f"auth chain {auth_chain_id}, "
                f"{','.join(entry_buckets)}, {meta['method']}, resolution {meta['resolution']}, "
                f"ions {nonpoly_meta['has_ions']}, ligands {nonpoly_meta['has_ligands']}, "
                f"metal atoms {metal_atom_count}, stripped atoms {stripped_atom_count}"
            )

    (OUTDIR / "pdb_rna_only.txt").write_text("".join(f"{x[0]}\n" for x in selected))
    (OUTDIR / f"pdb_rna_main_{MAIN_MIN}_{MAIN_MAX}.txt").write_text(
        "".join(f"{x}\n" for x in selected_by_bucket["main"])
    )
    (OUTDIR / f"pdb_rna_extension_{EXT_MIN}_{EXT_MAX}.txt").write_text(
        "".join(f"{x}\n" for x in selected_by_bucket["extension"])
    )

    with (OUTDIR / "selected_rna_sequences.tsv").open("w") as f:
        f.write(
            "record_id\tpdb_id\tlength\tchain_id\tauth_chain_id\tbucket\tmethod\tresolution\t"
            "release_date\thas_ions\thas_ligands\thas_nonpolymer_entities\t"
            "has_branched_entities\trna_atom_count\tmetal_atom_count\tstripped_atom_count\n"
        )
        for row in selected:
            f.write("\t".join(map(str, row)) + "\n")

    with (OUTDIR / "skipped_rna_sequences.tsv").open("w") as f:
        f.write("pdb_id\treason\n")
        for row in skipped:
            f.write("\t".join(row) + "\n")

    print("Selected total:", len(selected))
    print("Selected main:", len(selected_by_bucket["main"]))
    print("Selected extension:", len(selected_by_bucket["extension"]))
    print("Skipped:", len(skipped))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
