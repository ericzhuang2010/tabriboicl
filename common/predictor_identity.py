"""Canonical public identity for the evaluated TabRiboICL predictor."""

from __future__ import annotations

import fnmatch
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any


TABRIBOICL_IDENTITY_RELATIVE_PATH = Path("configs/predictors/tabriboicl.json")
PREDICTOR_REGISTRY_RELATIVE_PATH = Path("configs/predictors.json")
PREDICTOR_IDENTITY_MODULE_RELATIVE_PATH = Path(
    "common/predictor_identity.py"
)


@dataclass(frozen=True)
class PredictorIdentity:
    """Validated reader-facing identity and its machine-readable source."""

    public_name: str
    definition: str
    legacy_family: str
    legacy_result_id_pattern: str
    base_learner_name: str
    base_learner_class: str
    spec_path: Path
    spec: dict[str, Any]

    def matches_internal_id(self, predictor_id: str) -> bool:
        return predictor_id == self.legacy_family or fnmatch.fnmatchcase(
            predictor_id,
            self.legacy_result_id_pattern,
        )

    def display_name_for(self, predictor_id: str) -> str:
        if not self.matches_internal_id(predictor_id):
            raise ValueError(
                f"{predictor_id!r} does not match the canonical "
                f"{self.public_name} identity"
            )
        return self.public_name


def default_repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _required_text(mapping: dict[str, Any], key: str, context: str) -> str:
    value = mapping.get(key)
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{context}.{key} must be a non-empty string")
    return value


def load_predictor_registry(repo_root: Path | None = None) -> dict[str, dict[str, Any]]:
    """Read the shared registry for the eleven benchmark predictors."""
    root = (repo_root or default_repo_root()).expanduser().resolve()
    path = root / PREDICTOR_REGISTRY_RELATIVE_PATH
    spec = json.loads(path.read_text(encoding="utf-8"))
    if spec.get("schema_version") != 1:
        raise ValueError("unsupported predictor registry schema version")
    predictors = spec.get("predictors")
    if not isinstance(predictors, dict) or not predictors:
        raise ValueError("predictor registry must contain a nonempty predictors mapping")
    for name, record in predictors.items():
        if not isinstance(record, dict) or record.get("name") != name:
            raise ValueError(f"invalid predictor registry entry: {name}")
        for key in ("display_name", "entrypoint"):
            _required_text(record, key, name)
        for key in ("tags", "input_views"):
            values = record.get(key)
            if not isinstance(values, list) or not all(isinstance(v, str) for v in values):
                raise ValueError(f"{name}.{key} must be a list of strings")
        if not isinstance(record.get("outputs"), dict):
            raise ValueError(f"{name}.outputs must be a mapping")
        _required_text(record["outputs"], "prediction_schema", f"{name}.outputs")
        if "python_module" in record:
            _required_text(record, "python_module", name)
            args = record.get("default_args")
            if not isinstance(args, list) or not all(isinstance(v, str) for v in args):
                raise ValueError(f"{name}.default_args must be a list of strings")
    return predictors


def load_tabriboicl_identity(repo_root: Path | None = None) -> PredictorIdentity:
    root = (repo_root or default_repo_root()).expanduser().resolve()
    spec_path = root / TABRIBOICL_IDENTITY_RELATIVE_PATH
    if not spec_path.is_file():
        raise FileNotFoundError(
            f"canonical TabRiboICL identity specification is missing: {spec_path}"
        )
    spec = json.loads(spec_path.read_text(encoding="utf-8"))
    if spec.get("schema_version") != 1:
        raise ValueError("unsupported TabRiboICL identity schema version")

    identity = spec.get("identity")
    if not isinstance(identity, dict):
        raise ValueError("identity must be an object")
    legacy = identity.get("legacy_internal_identifiers")
    learner = identity.get("base_learner")
    if not isinstance(legacy, dict) or not isinstance(learner, dict):
        raise ValueError(
            "identity must define legacy_internal_identifiers and base_learner"
        )
    if learner.get("is_synonymous_with_public_predictor") is not False:
        raise ValueError("TabICL must not be declared synonymous with TabRiboICL")

    return PredictorIdentity(
        public_name=_required_text(identity, "public_predictor_name", "identity"),
        definition=_required_text(identity, "definition", "identity"),
        legacy_family=_required_text(
            legacy,
            "family",
            "identity.legacy_internal_identifiers",
        ),
        legacy_result_id_pattern=_required_text(
            legacy,
            "result_id_pattern",
            "identity.legacy_internal_identifiers",
        ),
        base_learner_name=_required_text(
            learner,
            "name",
            "identity.base_learner",
        ),
        base_learner_class=_required_text(
            learner,
            "class",
            "identity.base_learner",
        ),
        spec_path=spec_path,
        spec=spec,
    )


def identity_provenance_record(
    identity: PredictorIdentity,
    repo_root: Path,
) -> dict[str, Any]:
    root = repo_root.expanduser().resolve()
    spec_path = identity.spec_path.resolve()
    return {
        "public_predictor_name": identity.public_name,
        "definition": identity.definition,
        "legacy_internal_identifiers": identity.spec["identity"][
            "legacy_internal_identifiers"
        ],
        "base_learner": identity.spec["identity"]["base_learner"],
        "identity_spec": str(spec_path.relative_to(root)),
        "identity_spec_sha256": sha256_file(spec_path),
        "runtime_provenance": identity.spec["runtime_provenance"],
    }


TABRIBOICL_IDENTITY = load_tabriboicl_identity()
TABRIBOICL_PUBLIC_NAME = TABRIBOICL_IDENTITY.public_name
