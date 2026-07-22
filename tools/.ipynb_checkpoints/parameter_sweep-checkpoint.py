#!/usr/bin/env python3

"""
Create BESFEM config files for a parameter sweep.

Example:

python parameter_sweep.py \
    --base-config run_config.txt \
    --sweep sweep.json \
    --benchmark-root benchmark_runs
"""

import argparse
import csv
import json
from datetime import datetime
from pathlib import Path


def get_arguments():
    parser = argparse.ArgumentParser(
        description="Generate BESFEM config files for a parameter sweep."
    )

    parser.add_argument(
        "--base-config",
        required=True,
        help="Path to the original BESFEM config file.",
    )

    parser.add_argument(
        "--sweep",
        required=True,
        help="Path to the sweep JSON file.",
    )

    parser.add_argument(
        "--benchmark-root",
        default="benchmark_runs",
        help="Folder where benchmark runs will be created.",
    )

    return parser.parse_args()


def read_sweep_file(sweep_path):
    with open(sweep_path, "r") as file:
        sweep = json.load(file)

    if "parameter" not in sweep:
        raise ValueError("sweep.json must contain 'parameter'.")

    if "values" not in sweep:
        raise ValueError("sweep.json must contain 'values'.")

    if not sweep["values"]:
        raise ValueError("The values list cannot be empty.")

    return sweep


def parameter_exists(config_lines, parameter_name):
    for line in config_lines:
        stripped_line = line.strip()

        if not stripped_line:
            continue

        if stripped_line.startswith("#"):
            continue

        if "=" not in line:
            continue

        key, value = line.split("=", 1)

        if key.strip() == parameter_name:
            return True

    return False


def create_new_config(config_lines, parameter_name, new_value):
    new_lines = []

    for line in config_lines:
        stripped_line = line.strip()

        # Keep blank lines and comments unchanged.
        if not stripped_line or stripped_line.startswith("#"):
            new_lines.append(line)
            continue

        # Keep lines without "=" unchanged.
        if "=" not in line:
            new_lines.append(line)
            continue

        key, old_value = line.split("=", 1)

        if key.strip() == parameter_name:
            new_lines.append(f"{key.strip()} = {new_value}\n")
        else:
            new_lines.append(line)

    return new_lines


def main():
    args = get_arguments()

    base_config_path = Path(args.base_config)
    sweep_path = Path(args.sweep)
    benchmark_root = Path(args.benchmark_root)

    if not base_config_path.exists():
        raise FileNotFoundError(
            f"Base config was not found: {base_config_path}"
        )

    if not sweep_path.exists():
        raise FileNotFoundError(
            f"Sweep file was not found: {sweep_path}"
        )

    # Read the original BESFEM config.
    with open(base_config_path, "r") as file:
        config_lines = file.readlines()

    # Read the parameter and values to sweep.
    sweep = read_sweep_file(sweep_path)

    parameter_name = sweep["parameter"]
    parameter_values = sweep["values"]
    sweep_name = sweep.get("sweep_name", "besfem_sweep")

    if not parameter_exists(config_lines, parameter_name):
        raise ValueError(
            f"Parameter '{parameter_name}' was not found "
            f"in {base_config_path}."
        )

    # Create a unique benchmark folder.
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    benchmark_directory = (
        benchmark_root / f"{timestamp}_{sweep_name}"
    )

    input_directory = benchmark_directory / "input"
    results_directory = benchmark_directory / "results"
    logs_directory = benchmark_directory / "logs"

    input_directory.mkdir(parents=True)
    results_directory.mkdir()
    logs_directory.mkdir()

    manifest_rows = []

    # Create one config for each parameter value.
    for run_number, parameter_value in enumerate(
        parameter_values,
        start=1,
    ):
        run_id = f"run_{run_number:04d}"

        config_filename = f"{run_id}.cfg"
        config_path = input_directory / config_filename

        new_config_lines = create_new_config(
            config_lines,
            parameter_name,
            parameter_value,
        )

        with open(config_path, "w") as file:
            file.writelines(new_config_lines)

        result_directory = results_directory / run_id
        result_directory.mkdir()

        manifest_rows.append(
            {
                "run_id": run_id,
                "parameter": parameter_name,
                "value": parameter_value,
                "config_file": f"input/{config_filename}",
                "result_directory": f"results/{run_id}",
                "log_file": f"logs/{run_id}.log",
                "status": "generated",
            }
        )

        print(
            f"Created {config_filename}: "
            f"{parameter_name} = {parameter_value}"
        )

    # Create manifest.csv.
    manifest_path = benchmark_directory / "manifest.csv"

    field_names = [
        "run_id",
        "parameter",
        "value",
        "config_file",
        "result_directory",
        "log_file",
        "status",
    ]

    with open(manifest_path, "w", newline="") as file:
        writer = csv.DictWriter(
            file,
            fieldnames=field_names,
        )

        writer.writeheader()
        writer.writerows(manifest_rows)

    print()
    print(f"Created {len(parameter_values)} configurations.")
    print(f"Benchmark folder: {benchmark_directory}")
    print(f"Manifest: {manifest_path}")


if __name__ == "__main__":
    main()
    