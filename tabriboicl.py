"""Shared command interface for candidate preparation, prediction, and evaluation."""
from __future__ import annotations

import argparse
from importlib import import_module
import sys

from common.predictor_identity import load_predictor_registry


COMMANDS = {
    "extract": ("preprocessing.features", "Extract candidate features from aligned reads"),
    "split": ("preprocessing.splits", "Create training and evaluation splits"),
    "baselines": ("tools.offset_baselines.predict", "Run one or both fitted offset baselines"),
    "evaluate": ("common.evaluation", "Evaluate candidate predictions"),
}


def predict(argv: list[str]) -> int:
    registry = load_predictor_registry()
    models = [name for name, record in registry.items() if "python_module" in record]
    parser = argparse.ArgumentParser(
        prog="tabriboicl predict", add_help=False, allow_abbrev=False,
        description="Choose a candidate-table predictor; remaining arguments go to that predictor.",
    )
    parser.add_argument("--model", choices=["tabriboicl", *models, "logreg"], default="tabriboicl")
    args, remaining = parser.parse_known_args(argv)
    name = "tabicl" if args.model == "tabriboicl" else args.model
    if name == "logreg":
        module = "tools.shared.gbdt_candidate_model"
        defaults = ["--model", "logreg"]
    else:
        module = registry[name]["python_module"]
        defaults = registry[name]["default_args"]
    if "--help" in remaining or "-h" in remaining:
        parser.print_help()
        print(f"\nArguments for {args.model}:\n", flush=True)
    return import_module(module).main([*remaining, *defaults])


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    subparsers = parser.add_subparsers(dest="command", required=True)
    for command, (_, help_text) in COMMANDS.items():
        subparsers.add_parser(command, add_help=False, help=help_text)
    subparsers.add_parser("predict", add_help=False, help="Run TabRiboICL or a learned/fitted comparator")
    # Each implementation owns its options and help; import only the requested route.
    args = parser.parse_args(argv[:1])
    if args.command == "predict":
        return predict(argv[1:])
    return import_module(COMMANDS[args.command][0]).main(argv[1:])


if __name__ == "__main__":
    raise SystemExit(main())
