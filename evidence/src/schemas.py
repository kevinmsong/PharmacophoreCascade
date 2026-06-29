"""
Column schemas and DataFrame validators for all pipeline result tables.

Validates that loaded DataFrames have the expected columns and numeric types.
Raises ValueError with clear messages on schema violations.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Dict, List

import pandas as pd

logger = logging.getLogger(__name__)


@dataclass
class ColSpec:
    """Minimal column specification: name + expected pandas dtype category."""
    name: str
    dtype: str  # "numeric", "string", or "any"


@dataclass
class TableSchema:
    """Schema for a result table: required columns + optional columns."""
    name: str
    required: List[ColSpec]
    optional: List[ColSpec] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Schema definitions — derived from actual column names in result files
# ---------------------------------------------------------------------------

# results/screening_full_1M_topological_hashed.csv
STAGE3_AUDIT = TableSchema(
    name="stage3_audit",
    required=[
        ColSpec("stage3_screen_rank", "numeric"),
        ColSpec("zinc_id", "string"),
        ColSpec("smiles", "string"),
        ColSpec("hotspot_weighted_pct", "numeric"),
        ColSpec("pair_hash_overlap_pct", "numeric"),
        ColSpec("cascade_score_pct", "numeric"),
        ColSpec("shortlist_rank", "numeric"),
        ColSpec("fit_rmsd_angstrom", "numeric"),
        ColSpec("mean_pair_distance_error_angstrom", "numeric"),
        ColSpec("weighted_coverage_pct", "numeric"),   # stage-3 geometry coverage
        ColSpec("mw", "numeric"),
        ColSpec("logp", "numeric"),
    ],
    optional=[
        ColSpec("hbd", "numeric"),
        ColSpec("hba", "numeric"),
        ColSpec("hotspot_exact_matches", "numeric"),
        ColSpec("hotspot_compatible_matches", "numeric"),
        ColSpec("hotspot_group_count", "numeric"),
        ColSpec("pair_hash_recall_pct", "numeric"),
        ColSpec("pair_hash_precision_pct", "numeric"),
        ColSpec("rerank_conformers", "numeric"),
        ColSpec("anchor_weighted_pct", "numeric"),
        ColSpec("geometry_penalty_factor", "numeric"),
        ColSpec("reference_feature_span_angstrom", "numeric"),
        ColSpec("ligand_feature_span_angstrom", "numeric"),
        ColSpec("features_matched", "numeric"),
        ColSpec("lipinski", "any"),
        ColSpec("chemistry_flagged", "any"),
        ColSpec("pains_alert", "any"),
        ColSpec("reactive_flags", "any"),
    ],
)

# final top_1000 output
FINAL_RANKED = TableSchema(
    name="final_ranked",
    required=[
        ColSpec("final_rank", "numeric"),
        ColSpec("zinc_id", "string"),
        ColSpec("smiles", "string"),
        ColSpec("native_weighted_coverage_pct", "numeric"),
        ColSpec("native_matched_reference_features", "numeric"),
        ColSpec("native_fit_rmsd_angstrom", "numeric"),
        ColSpec("native_mean_pair_distance_error_angstrom", "numeric"),
        ColSpec("shortlist_rank", "numeric"),
        ColSpec("cascade_score_pct", "numeric"),
        ColSpec("hotspot_weighted_pct", "numeric"),
        ColSpec("pair_hash_overlap_pct", "numeric"),
        ColSpec("screen_weighted_coverage_pct", "numeric"),
        ColSpec("mw", "numeric"),
        ColSpec("logp", "numeric"),
    ],
    optional=[
        ColSpec("stage3_screen_rank", "numeric"),
        ColSpec("native_candidate_pool_rank", "numeric"),
        ColSpec("native_selection_rank", "numeric"),
        ColSpec("pool_source_stage3_rank", "any"),
        ColSpec("pool_source_hotspot_breadth", "any"),
        ColSpec("pool_source_native_support", "any"),
        ColSpec("native_microstate_id", "string"),
        ColSpec("native_matched_reference_residues", "string"),
        ColSpec("native_conformer_id", "any"),
        ColSpec("native_conformer_energy", "numeric"),
        ColSpec("fit_rmsd_angstrom", "numeric"),
        ColSpec("mean_pair_distance_error_angstrom", "numeric"),
        ColSpec("hotspot_exact_matches", "numeric"),
        ColSpec("hotspot_group_count", "numeric"),
        ColSpec("hbd", "numeric"),
        ColSpec("hba", "numeric"),
        ColSpec("lipinski", "any"),
        ColSpec("chemistry_flagged", "any"),
    ],
)

# derived/.../analysis/microstate_native_mapping.csv
NATIVE_SCORED = TableSchema(
    name="native_scored",
    required=[
        ColSpec("ligand_id", "string"),
        ColSpec("microstate_id", "string"),
        ColSpec("native_weighted_coverage_pct", "numeric"),
        ColSpec("matched_reference_features", "numeric"),
        ColSpec("fit_rmsd_angstrom", "numeric"),
        ColSpec("mean_pair_distance_error_angstrom", "numeric"),
        ColSpec("shortlist_rank", "numeric"),
        ColSpec("cascade_score_pct", "numeric"),
        ColSpec("hotspot_weighted_pct", "numeric"),
        ColSpec("pair_hash_overlap_pct", "numeric"),
        ColSpec("mw", "numeric"),
        ColSpec("logp", "numeric"),
    ],
    optional=[
        ColSpec("smiles", "string"),
        ColSpec("canonical_smiles", "string"),
        ColSpec("murcko_scaffold", "string"),
        ColSpec("source_input_rank", "numeric"),
        ColSpec("matched_reference_feature_ids", "string"),
        ColSpec("matched_reference_residues", "string"),
        ColSpec("matched_reference_families", "string"),
        ColSpec("weighted_coverage_pct", "numeric"),   # stage-3 geometry coverage in this table
        ColSpec("rdkit_mw", "numeric"),
        ColSpec("rdkit_logp", "numeric"),
        ColSpec("aromatic_rings", "numeric"),
        ColSpec("scientific_interpretation", "string"),
    ],
)

# derived/.../analysis/microstate_native_feature_matches.csv
FEATURE_MATCHES = TableSchema(
    name="feature_matches",
    required=[
        ColSpec("ligand_id", "string"),
        ColSpec("microstate_id", "string"),
        ColSpec("reference_feature_id", "string"),
        ColSpec("reference_residue_label", "string"),
        ColSpec("reference_family", "string"),
        ColSpec("reference_weight", "numeric"),
        ColSpec("ligand_feature_id", "string"),
        ColSpec("ligand_family", "string"),
    ],
    optional=[
        ColSpec("source_input_rank", "numeric"),
        ColSpec("ligand_atom_ids", "string"),
    ],
)

# derived/.../analysis/native_reference_features.csv
REFERENCE_FEATURES = TableSchema(
    name="reference_features",
    required=[
        ColSpec("feature_id", "string"),
        ColSpec("residue_label", "string"),
        ColSpec("family", "string"),
        ColSpec("weight", "numeric"),
    ],
    optional=[
        ColSpec("atom_names", "string"),
        ColSpec("x", "numeric"),
        ColSpec("y", "numeric"),
        ColSpec("z", "numeric"),
        ColSpec("contact_distance_angstrom", "numeric"),
    ],
)

# derived/.../analysis/native_rank_comparison.csv  (pre-computed by existing pipeline)
RANK_COMPARISON = TableSchema(
    name="rank_comparison",
    required=[
        ColSpec("ligand_id", "string"),
        ColSpec("weighted_coverage_pct", "numeric"),
        ColSpec("native_weighted_coverage_pct", "numeric"),
        ColSpec("screen_rank", "numeric"),
        ColSpec("native_rank", "numeric"),
        ColSpec("rank_shift", "numeric"),
        ColSpec("abs_rank_shift", "numeric"),
    ],
)

# derived/.../analysis/ligand_hotspot_compatibility.csv
HOTSPOT_COMPAT = TableSchema(
    name="hotspot_compat",
    required=[
        ColSpec("ligand_id", "string"),
        ColSpec("microstate_id", "string"),
        ColSpec("receptor_residue_label", "string"),
        ColSpec("receptor_residue_class", "string"),
        ColSpec("matched_native_feature_count", "numeric"),
        ColSpec("matched_native_feature_weight_sum", "numeric"),
    ],
    optional=[
        ColSpec("source_input_rank", "numeric"),
        ColSpec("closest_native_distance_angstrom", "numeric"),
        ColSpec("matched_native_feature_ids", "string"),
        ColSpec("matched_native_residues", "string"),
    ],
)

ALL_SCHEMAS: Dict[str, TableSchema] = {
    s.name: s
    for s in [
        STAGE3_AUDIT,
        FINAL_RANKED,
        NATIVE_SCORED,
        FEATURE_MATCHES,
        REFERENCE_FEATURES,
        RANK_COMPARISON,
        HOTSPOT_COMPAT,
    ]
}


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def validate(df: pd.DataFrame, schema: TableSchema) -> pd.DataFrame:
    """Validate a DataFrame against a TableSchema.

    Raises ValueError listing all missing required columns and type mismatches.
    Returns the validated DataFrame (unchanged).
    """
    errors: List[str] = []

    for spec in schema.required:
        if spec.name not in df.columns:
            errors.append(f"  missing required column: '{spec.name}'")
        elif spec.dtype == "numeric":
            if not pd.api.types.is_numeric_dtype(df[spec.name]):
                errors.append(
                    f"  column '{spec.name}' expected numeric, got {df[spec.name].dtype}"
                )

    if errors:
        raise ValueError(
            f"Schema validation failed for table '{schema.name}':\n"
            + "\n".join(errors)
            + f"\n  Available columns: {list(df.columns)}"
        )

    n_missing_opt = sum(
        1 for spec in schema.optional if spec.name not in df.columns
    )
    if n_missing_opt:
        logger.debug(
            "Table '%s': %d optional columns absent (not an error)",
            schema.name,
            n_missing_opt,
        )

    return df
