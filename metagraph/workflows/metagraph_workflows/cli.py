import argparse
import difflib
import importlib
import logging
import os
import re
import shlex
import shutil
import sys
import subprocess
import time
from pathlib import Path
from typing import Iterable, Optional, Dict, Any

import snakemake
import snakemake.io
import snakemake.utils
import yaml

from .workflow_configs import SEQS_FILE_LIST_PATH, SEQS_DIR_PATH, \
    AnnotationLabelsSource, AnnotationFormats

WORKFLOW_ROOT = Path(__file__).parent / 'snakemake'

LOGGING_FORMAT='%(asctime)s - %(levelname)s: %(message)s'

logging.basicConfig(format=LOGGING_FORMAT, level=logging.WARNING)


default_path = Path(WORKFLOW_ROOT / 'default.yml')
COUNT_COMPATIBLE_FORMATS = {
    AnnotationFormats.INT_BRWT,
    AnnotationFormats.ROW_DIFF_INT_BRWT,
    AnnotationFormats.ROW_DIFF_INT_DISK,
}
COORD_COMPATIBLE_FORMATS = {
    AnnotationFormats.BRWT_COORD,
    AnnotationFormats.ROW_DIFF_COORD,
    AnnotationFormats.ROW_DIFF_BRWT_COORD,
    AnnotationFormats.ROW_DIFF_DISK_COORD,
}


def _help_formatter(prog: str):
    return argparse.RawTextHelpFormatter(prog, width=100, max_help_position=34)


def _help_color(text: str, color_code: str) -> str:
    if not sys.stdout.isatty():
        return text
    return f"\033[{color_code}m{text}\033[0m"


def _default_threads_auto() -> int:
    try:
        out = subprocess.check_output(["nproc"], text=True).strip()
        n = int(out)
        if n > 0:
            return n
    except Exception:
        pass
    try:
        out = subprocess.check_output(["sysctl", "-n", "hw.logicalcpu"], text=True).strip()
        n = int(out)
        if n > 0:
            return n
    except Exception:
        pass
    return snakemake.utils.available_cpu_count()


def _format_bytes(size_bytes: int) -> str:
    units = ["B", "KB", "MB", "GB", "TB"]
    val = float(size_bytes)
    for unit in units:
        if val < 1024 or unit == units[-1]:
            return f"{val:.1f}{unit}"
        val /= 1024.0
    return f"{size_bytes}B"


def _format_seconds_human(value: str) -> str:
    try:
        seconds = float(value)
    except (TypeError, ValueError):
        return "-"

    if seconds < 1.0:
        return f"{int(round(seconds * 1000))}ms"
    if seconds < 60.0:
        return f"{seconds:.1f}s"
    mins, secs = divmod(seconds, 60.0)
    if mins < 60:
        return f"{int(mins)}m {int(round(secs))}s"
    hours, mins = divmod(mins, 60.0)
    return f"{int(hours)}h {int(mins)}m {int(round(secs))}s"


def _parse_timing_line(line: str) -> Optional[Dict[str, str]]:
    if "[timing]" not in line:
        return None
    parts = {}
    for token in line.strip().split():
        if "=" in token:
            k, v = token.split("=", 1)
            parts[k] = v
    return parts if parts else None


def _collect_stage_summaries(output_dir: Path) -> Iterable[Dict[str, str]]:
    logs_root = output_dir / "logs"
    if not logs_root.exists():
        return []

    rows = []
    for log_path in logs_root.rglob("*.log"):
        try:
            text = log_path.read_text(errors='replace')
        except OSError:
            continue
        timing = None
        for line in reversed(text.splitlines()):
            parsed = _parse_timing_line(line)
            if parsed:
                timing = parsed
                break
        if not timing:
            continue

        stage_name = str(log_path.relative_to(logs_root)).replace(".log", "")
        rss_kb = timing.get("max_rss_kb", "-")
        try:
            rss_gb = f"{(float(rss_kb) / (1024.0 * 1024.0)):.3f}"
        except (TypeError, ValueError):
            rss_gb = "-"
        rows.append({
            "stage": stage_name,
            "wall": _format_seconds_human(timing.get("wall_sec", "-")),
            "user": _format_seconds_human(timing.get("user_sec", "-")),
            "sys": _format_seconds_human(timing.get("sys_sec", "-")),
            "rss": rss_gb,
            "_mtime": str(log_path.stat().st_mtime),
        })
    rows.sort(key=lambda r: float(r.get("_mtime", "0")))
    for r in rows:
        r.pop("_mtime", None)
    return rows


def _directory_size(path: Path) -> int:
    total = 0
    for p in path.rglob("*"):
        if p.is_file():
            try:
                total += p.stat().st_size
            except OSError:
                pass
    return total


def _print_run_summary(output_dir: Path, dryrun: bool) -> None:
    print("\n=== Workflow summary ===")
    if dryrun:
        print("Mode: dry-run (no commands executed, no runtime metrics).")
        return

    rows = list(_collect_stage_summaries(output_dir))
    if rows:
        print("1) Executed steps (timing + memory):")
        idx_w = max(2, len(str(len(rows))))
        stage_w = max(12, min(56, max(len(r["stage"]) for r in rows)))
        header = (
            f"   {'#':>{idx_w}}  "
            f"{'stage':<{stage_w}}  "
            f"{'wall':>8}  {'user':>8}  {'sys':>8}  {'rss_gb':>8}"
        )
        print(header)
        print("   " + "-" * (len(header) - 3))
        for i, r in enumerate(rows, start=1):
            print(
                f"   {i:>{idx_w}}  "
                f"{r['stage']:<{stage_w}}  "
                f"{r['wall']:>8}  {r['user']:>8}  {r['sys']:>8}  {r['rss']:>8}"
            )
    else:
        print("1) Executed steps: unavailable (no timing entries found).")

    def _is_final_artifact(path: Path) -> bool:
        name = path.name
        return name == "graph.dbg" or (name.startswith("graph") and name.endswith(".annodbg"))

    total_size = _directory_size(output_dir)
    final_artifacts = []
    intermediate_artifacts = []
    for p in output_dir.rglob("*"):
        if not p.is_file():
            continue
        if "/logs/" in str(p):
            continue
        try:
            size = p.stat().st_size
        except OSError:
            continue
        rel = p.relative_to(output_dir)
        if _is_final_artifact(p):
            final_artifacts.append((rel, size))
        else:
            intermediate_artifacts.append((rel, size))

    final_total = sum(sz for _, sz in final_artifacts)
    intermediate_total = sum(sz for _, sz in intermediate_artifacts)
    final_artifacts.sort(key=lambda x: x[1], reverse=True)
    intermediate_artifacts.sort(key=lambda x: x[1], reverse=True)

    print(f"2) Output directory: {output_dir}")
    print(f"3) Total output size: {_format_bytes(total_size)}")
    print(f"4) Final artifacts total: {_format_bytes(final_total)}")
    for rel, sz in final_artifacts:
        print(f"   - {_format_bytes(sz):>8}  {rel}")
    print(f"5) Intermediate artifacts total: {_format_bytes(intermediate_total)}")

    if intermediate_artifacts:
        print("6) Largest intermediate artifacts:")
        top = intermediate_artifacts[:10]
        idx_w = max(2, len(str(len(top))))
        size_vals = [_format_bytes(sz) for _, sz in top]
        size_w = max(8, max(len(s) for s in size_vals))
        print(f"   {'#':>{idx_w}}  {'size':>{size_w}}  path")
        print("   " + "-" * (idx_w + size_w + 8 + 20))
        for i, ((rel, _), size_s) in enumerate(zip(top, size_vals), start=1):
            rel_str = str(rel)
            if sys.stdout.isatty():
                rel_str = f"\033[90m{rel_str}\033[0m"
            print(f"   {i:>{idx_w}}  {size_s:>{size_w}}  {rel_str}")


def _parse_annotation_format_value(value: str) -> AnnotationFormats:
    try:
        return AnnotationFormats(value)
    except ValueError:
        valid_values = [v.value for v in AnnotationFormats]
        suggestion = difflib.get_close_matches(value, valid_values, n=1)
        suggestion_msg = f" Did you mean '{suggestion[0]}'?" if suggestion else ""
        raise ValueError(
            f"Unsupported annotation format '{value}'. "
            f"Valid values: {', '.join(valid_values)}.{suggestion_msg}"
        )


def _format_error(text: str) -> str:
    if sys.stderr.isatty():
        red = "\033[31m"
        pink = "\033[95m"
        reset = "\033[0m"
        lines = text.splitlines() or [text]

        rendered = [f"{red}Error: {lines[0]}{reset}"]
        in_metagraph = False
        for line in lines[1:]:
            stripped = line.strip()
            if stripped == "metagraph output:":
                in_metagraph = True
                rendered.append(f"{red}{line}{reset}")
                continue
            if in_metagraph and stripped.startswith("[timing]"):
                in_metagraph = False
                rendered.append(f"{red}{line}{reset}")
                continue
            if in_metagraph:
                rendered.append(f"{pink}{line}{reset}")
            else:
                rendered.append(f"{red}{line}{reset}")

        return "\n".join(rendered)
    return f"Error: {text}"


def _validate_metagraph_cmd(cmd: str) -> None:
    cmd_parts = shlex.split(cmd)
    if not cmd_parts:
        raise ValueError("--metagraph-cmd is empty. Provide a valid executable path or command name.")

    executable = cmd_parts[0]
    if "/" in executable:
        exe_path = Path(executable).expanduser()
        if not exe_path.exists():
            raise ValueError(
                f"MetaGraph executable not found at '{executable}'. "
                "Provide a valid path via --metagraph-cmd or add 'metagraph' to PATH."
            )
        if not exe_path.is_file() or not os.access(exe_path, os.X_OK):
            raise ValueError(
                f"MetaGraph executable '{executable}' is not executable. "
                "Fix permissions or use --metagraph-cmd with a valid executable."
            )
    else:
        if shutil.which(executable) is None:
            raise ValueError(
                f"MetaGraph executable '{executable}' was not found in PATH. "
                "Use --metagraph-cmd with an absolute path or add it to PATH."
            )


def _extract_first_relevant_error_line(log_text: str) -> Optional[str]:
    lines = log_text.splitlines()
    for line in lines:
        s = line.strip()
        if not s:
            continue
        lowered = s.lower()
        if lowered.startswith("usage:") or "unrecognized option" in lowered or "unknown option" in lowered \
                or lowered.startswith("error:") or "[error]" in lowered:
            return s
    for line in reversed(lines):
        s = line.strip()
        if not s:
            continue
        lowered = s.lower()
        if "failed_to_exec=" in s or "command not found" in lowered or "no such file or directory" in lowered:
            return s
        if "[error]" in lowered or lowered.startswith("error:") or "invalid argument" in lowered:
            return s
    # Fallback: show last meaningful non-timing line to avoid empty/opaque summaries.
    for line in reversed(lines):
        s = line.strip()
        if not s:
            continue
        lowered = s.lower()
        if s.startswith("[timing]") or "[trace]" in lowered or "[debug]" in lowered or "[info]" in lowered:
            continue
        if "\t" in s and "%" in s:
            continue
        return s
    return None


def _latest_snakemake_log_for_run(run_started_at: float) -> Optional[Path]:
    snakemake_log_dir = Path(".snakemake/log")
    if not snakemake_log_dir.exists():
        return None
    candidates = [p for p in snakemake_log_dir.glob("*.snakemake.log") if p.is_file()]
    current_run = [p for p in candidates if p.stat().st_mtime >= run_started_at - 1.0]
    pool = current_run if current_run else candidates
    if not pool:
        return None
    return max(pool, key=lambda p: p.stat().st_mtime)


def _extract_failing_shell_command(cli_output: str) -> Optional[str]:
    lines = cli_output.splitlines()
    shell_idx = None
    for i, line in enumerate(lines):
        if line.strip() == "shell:":
            shell_idx = i
    if shell_idx is None:
        return None

    cmd_lines = []
    for line in lines[shell_idx + 1:]:
        s = line.rstrip()
        if "(command exited with non-zero exit code)" in s:
            break
        if not s.strip():
            continue
        cmd_lines.append(s.strip())
    if not cmd_lines:
        return None
    cmd = " ".join(cmd_lines)
    return cmd if len(cmd) <= 500 else (cmd[:497] + "...")


def _compact_command(cmd: str, max_len: int = 220) -> str:
    # Show the most relevant subcommand when piped through wrappers.
    metagraph_idx = cmd.find("metagraph ")
    if metagraph_idx != -1:
        cmd = cmd[metagraph_idx:]
    cmd = " ".join(cmd.split())
    return cmd if len(cmd) <= max_len else (cmd[:max_len - 3] + "...")


def _extract_stage_log_block(log_text: str, head_lines: int = 20, tail_lines: int = 10) -> Iterable[str]:
    lines = log_text.splitlines()
    max_lines = head_lines + tail_lines
    if len(lines) <= max_lines:
        return lines
    return (
        lines[:head_lines]
        + ["... (truncated; see failing stage log for full output)"]
        + lines[-tail_lines:]
    )


def _extract_timing_exit_code(log_text: str) -> Optional[int]:
    for line in reversed(log_text.splitlines()):
        parsed = _parse_timing_line(line)
        if not parsed:
            continue
        raw = parsed.get("exit_code")
        if raw is None:
            return None
        try:
            return int(raw)
        except ValueError:
            return None
    return None


def _extract_snakemake_failure_line(cli_output: str) -> Optional[str]:
    for line in cli_output.splitlines():
        s = line.strip()
        if not s:
            continue
        if "MissingOutputException" in s or "RuleException" in s or s.endswith("WorkflowError:"):
            return s
    return None


def _extract_failing_stage_log_from_snakemake_output(cli_output: str, output_dir: Path) -> Optional[Path]:
    lines = cli_output.splitlines()
    for line in reversed(lines):
        if "log:" not in line:
            continue
        after = line.split("log:", 1)[1].strip()
        if not after:
            continue
        candidate = after.split(" (", 1)[0].strip().rstrip(",")
        if not candidate:
            continue
        p = Path(candidate)
        if not p.is_absolute():
            p = Path(output_dir).parent / p
        return p
    return None


def _summarize_snakemake_failure_for_run(output_dir: Path, run_started_at: float) -> str:
    details = []
    stage_log = None
    snakemake_log = _latest_snakemake_log_for_run(run_started_at)
    snakemake_output = None
    if snakemake_log and snakemake_log.exists():
        try:
            snakemake_output = snakemake_log.read_text(errors='replace')
        except OSError:
            snakemake_output = None

    if snakemake_output:
        stage_log = _extract_failing_stage_log_from_snakemake_output(snakemake_output, Path(output_dir))

    if stage_log and stage_log.exists():
        try:
            stage_text = stage_log.read_text(errors='replace')
            stage_block = list(_extract_stage_log_block(stage_text))
            stage_exit_code = _extract_timing_exit_code(stage_text)
            root_line = _extract_first_relevant_error_line(stage_text)
            if stage_block and stage_exit_code != 0:
                details.append("metagraph output:")
                details.extend(stage_block)
            if root_line:
                if stage_exit_code == 0:
                    details.append(
                        "Root cause: metagraph command exited successfully, "
                        "but Snakemake reported missing/incorrect output files."
                    )
                else:
                    details.append(f"Root cause: {root_line}")
            details.append(f"Failing stage log: {stage_log}")
        except OSError:
            pass
    if snakemake_output:
        failure_line = _extract_snakemake_failure_line(snakemake_output)
        if failure_line:
            details.append(f"Snakemake failure: {failure_line}")
        cmd = _extract_failing_shell_command(snakemake_output)
        if cmd:
            details.append(f"Failing command: {_compact_command(cmd)}")

    details.append(f"Stage logs: {Path(output_dir) / 'logs'}")
    details.append("Snakemake logs: .snakemake/log/")
    return "Workflow execution failed.\n" + "\n".join(details)


def run_build_workflow(
        output_dir: Path,
        seqs_file_list_path: Optional[Path] = None,
        seqs_dir_path: Optional[Path] = None,
        k: Optional[int] = None,
        base_name: Optional[str] = None,
        build_primary_graph: bool = False,
        annotation_formats: Iterable[AnnotationFormats] = (),
        annotation_labels_source: Optional[AnnotationLabelsSource] = None,
        with_counts: bool = False,
        with_coordinates: bool = False,
        count_width: Optional[int] = None,
        metagraph_cmd: Optional[str] = None,
        threads: Optional[int] = None,
        force: bool = False,
        verbose: bool = False,
        dryrun: bool = False,
        additional_snakemake_args: Optional[Dict[str, Any]] = None
) -> None:
    # TODO: support str argumt?

    snakefile_path = Path(WORKFLOW_ROOT / 'Snakefile')

    with open(default_path, 'r') as f:
        config = yaml.safe_load(f)

    if not seqs_file_list_path and not seqs_dir_path:
        raise ValueError("seqs_file_list_path and seqs_dir_path cannot both be None")

    if seqs_file_list_path:
        config[SEQS_FILE_LIST_PATH] = str(seqs_file_list_path)
    if seqs_dir_path:
        config[SEQS_DIR_PATH] = str(seqs_dir_path)

    config['output_directory'] = str(output_dir)

    config['k'] = k if k else config['k']

    if annotation_labels_source:
        config['annotation_labels_source'] = annotation_labels_source.value

    config['base_name'] = base_name if base_name else config['base_name']
    config['build_primary_graph'] = build_primary_graph

    selected_formats = set(annotation_formats)
    has_count_formats = any(af in COUNT_COMPATIBLE_FORMATS for af in selected_formats)
    has_coord_formats = any(af in COORD_COMPATIBLE_FORMATS for af in selected_formats)
    effective_with_counts = with_counts or has_count_formats
    effective_with_coordinates = with_coordinates or has_coord_formats

    if effective_with_counts and effective_with_coordinates:
        raise ValueError(
            "Count-aware and coordinate-aware modes are mutually exclusive in this workflow. "
            "Choose either count formats (--with-counts / int_* / row_diff_int_*) or "
            "coordinate formats (--with-coords / *_coord)."
        )

    if effective_with_counts and not annotation_formats:
        config['annotation_formats'] = [AnnotationFormats.ROW_DIFF_INT_BRWT.value]
    elif effective_with_coordinates and not annotation_formats:
        config['annotation_formats'] = [AnnotationFormats.ROW_DIFF_BRWT_COORD.value]
    else:
        config['annotation_formats'] = [af.value for af in
                                        annotation_formats] if annotation_formats else config['annotation_formats']

    if effective_with_counts:
        invalid_formats = [af.value for af in annotation_formats if af not in COUNT_COMPATIBLE_FORMATS]
        if invalid_formats:
            raise ValueError(
                "Count-aware mode is enabled (--with-counts or count-capable --annotation-format), "
                "--annotation-format must be one of: "
                + ", ".join(sorted([f.value for f in COUNT_COMPATIBLE_FORMATS]))
                + f". Got: {', '.join(invalid_formats)}"
            )
    if effective_with_coordinates:
        invalid_formats = [af.value for af in annotation_formats if af not in COORD_COMPATIBLE_FORMATS]
        if invalid_formats:
            raise ValueError(
                "Coordinate-aware mode is enabled (--with-coords or *_coord --annotation-format), "
                "--annotation-format must be one of: "
                + ", ".join(sorted([f.value for f in COORD_COMPATIBLE_FORMATS]))
                + f". Got: {', '.join(invalid_formats)}"
            )
    config['with_counts'] = effective_with_counts
    config['with_coordinates'] = effective_with_coordinates
    if count_width is not None and not (2 <= count_width <= 32):
        raise ValueError(f"--count-width must be in range [2, 32], got {count_width}")
    if count_width is not None and not effective_with_counts:
        raise ValueError("--count-width can only be used with count-aware mode (--with-counts or count formats).")
    if count_width is not None:
        config['count_width'] = count_width

    config['metagraph_cmd'] = metagraph_cmd if metagraph_cmd else config['metagraph_cmd']
    if not dryrun:
        _validate_metagraph_cmd(config['metagraph_cmd'])
    config['max_threads'] = threads if threads else _default_threads_auto()

    if verbose:
        importlib.reload(logging)
        logging.basicConfig(format=LOGGING_FORMAT, level=logging.INFO)
        logging.info("Dumping config:")
        for k, v in sorted(config.items(), key=lambda t: t[0]):
            logging.info(f"\t{k}: {v}")

    additional_args = additional_snakemake_args if additional_snakemake_args else {}

    # Build snakemake command
    cmd = ['python', '-m', 'snakemake', '--snakefile', str(snakefile_path)]

    # Add config file
    output_dir_path = Path(output_dir)
    config_file = output_dir_path / 'config.yaml'
    config_file.parent.mkdir(parents=True, exist_ok=True)
    with open(config_file, 'w') as f:
        yaml.dump(config, f)
    cmd.extend(['--configfile', str(config_file)])

    # Add other arguments
    if force:
        cmd.append('--forceall')
    if dryrun:
        cmd.append('--dry-run')
    if threads:
        cmd.extend(['--cores', str(threads)])
    else:
        # Add default cores if not specified
        cmd.extend(['--cores', str(_default_threads_auto())])

    # Add additional arguments
    for key, value in additional_args.items():
        if isinstance(value, bool):
            if value:
                cmd.append(f'--{key}')
        else:
            cmd.extend([f'--{key}', str(value)])

    # Run snakemake
    # Keep stdout/stderr attached so Snakemake preserves colors and rich formatting.
    run_started_at = time.time()
    result = subprocess.run([' '.join(cmd)], shell=True)

    if result.returncode != 0:
        raise RuntimeError(_summarize_snakemake_failure_for_run(Path(output_dir), run_started_at))
    _print_run_summary(Path(output_dir), dryrun)


def setup_build_parser(parser):
    label_sources = [v.value for v in AnnotationLabelsSource]
    count_formats = sorted([f.value for f in COUNT_COMPATIBLE_FORMATS])
    count_formats_display = [_help_color(fmt, "33") for fmt in count_formats]
    coord_formats = sorted([f.value for f in COORD_COMPATIBLE_FORMATS])
    coord_formats_display = [_help_color(fmt, "35") for fmt in coord_formats]
    classic_formats_display = [
        _help_color(fmt, "36")
        for fmt in [
            "row", "bin_rel_wt", "flat", "rbfish", "brwt", "relax.brwt",
            "rb_brwt", "row_diff_brwt", "relax.row_diff_brwt",
        ]
    ]
    with_counts_label = _help_color("--with-counts", "33")
    with_coords_label = _help_color("--with-coords", "35")
    default_base_fmt = _help_color("relax.row_diff_brwt", "36")
    default_count_fmt = _help_color("row_diff_int_brwt", "33")
    default_coord_fmt = _help_color("row_diff_brwt_coord", "35")
    default_count_width = _help_color("8", "33")
    zero_count_width = _help_color("0", "36")

    parser.description = (
        "Build a MetaGraph graph + annotation workflow from a sequence list or directory."
    )
    parser.epilog = (
        "Examples:\n"
        "  metagraph-workflows build --seqs-file-list-path files.txt -k 31 -o out/\n"
        "  metagraph-workflows build --seqs-file-list-path files.txt --with-counts -o out/\n"
        "  metagraph-workflows build --seqs-file-list-path files.txt --with-coords -o out/"
    )

    input_seq_group = parser.add_argument_group('input/output')

    input_seq_group_xor = input_seq_group.add_mutually_exclusive_group(required=True)
    input_seq_group_xor.add_argument('--seqs-file-list-path',
                                     metavar='PATH',
                                     help='Path to a text file with sample paths (one per line) []')
    input_seq_group_xor.add_argument('--seqs-dir-path',
                                     metavar='DIR',
                                     help="Directory containing samples []")
    input_seq_group.add_argument('-o', '--output_dir', type=Path, required=True,
                                 help='Output directory [required]')

    graph = parser.add_argument_group('graph')
    graph.add_argument('-k', type=int, default=31, metavar='K',
                       help='k-mer length [31]')
    graph.add_argument('--base-name', default='graph', metavar='NAME',
                       help='Base output name [graph]')
    graph.add_argument('--primary', dest='build_primary_graph', default=False,
                       action='store_true',
                       help='Build canonical graph first, then derive/build primary graph [False]')

    annotation = parser.add_argument_group('annotation')
    all_formats_help = "\n".join([
        f"    {classic_formats_display[0]}, {classic_formats_display[1]}, {classic_formats_display[2]}, {classic_formats_display[3]}, {classic_formats_display[4]}, {classic_formats_display[5]}, {classic_formats_display[6]},",
        f"             {classic_formats_display[7]}, {classic_formats_display[8]}",
        f"    {count_formats_display[0]}, {count_formats_display[1]}, {count_formats_display[2]}",
        f"    {coord_formats_display[0]}, {coord_formats_display[2]}, {coord_formats_display[1]}, {coord_formats_display[3]}",
    ])
    annotation.add_argument('--anno-source',
                            dest='annotation_labels_source',
                            type=AnnotationLabelsSource,
                            default=AnnotationLabelsSource.SEQUENCE_HEADERS,
                            metavar='SOURCE',
                            help=f"Column label source: {', '.join(label_sources)} [sequence_headers]\n"
                                 "  ")
    annotation.add_argument('--annotation-format', action='append',
                            default=[],
                            metavar='FORMAT',
                            help=f"Annotation format (can be used multiple times).\n"
                                 f"{all_formats_help}\n"
                                 f"  [{default_base_fmt}/{default_count_fmt}/{default_coord_fmt}]\n"
                                 "  ")
    annotation.add_argument('--with-counts', default=False, action='store_true',
                            help=f"Index with k-mer counts [False]\n"
                                 f"  Supported for {', '.join(count_formats_display)}.")
    annotation.add_argument('--count-width', type=int, default=None,
                            metavar='BITS',
                            help=f"Bit width for count values (passed to annotate/transform_anno) [{zero_count_width}/{default_count_width}]\n"
                                 "  ")
    annotation.add_argument('--with-coords', dest='with_coordinates', default=False, action='store_true',
                            help=f"Index with k-mer positions [False]\n"
                                 f"  Supported for {coord_formats_display[0]}, {coord_formats_display[2]}, {coord_formats_display[1]},\n"
                                 f"                     {coord_formats_display[3]}.")

    workflow = parser.add_argument_group('other')
    workflow.add_argument('--threads', type=int, default=None, metavar='N',
                          help='Max cores for Snakemake execution [num_cores]')
    workflow.add_argument('--force', default=False, action='store_true',
                          help='Force re-run all rules [False]')
    workflow.add_argument('--verbose', default=False, action='store_true',
                          help='Print verbose config/runtime logs [False]')
    workflow.add_argument('--dryrun', default=False, action='store_true',
                          help='Render DAG and config only; do not execute rules [False]')
    workflow.add_argument('--metagraph-cmd', type=str, default=None, metavar='CMD',
                          help='Path/command for metagraph executable [metagraph from PATH]')
    workflow.add_argument('--extra-args', dest='additional_snakemake_args', metavar='ARGS', type=str, default='',
                          help='Extra arguments to pass to snakemake [none]\n'
                               '  Example: --extra-args="arg1=val1 arg2=val2"')
    options = parser.add_argument_group('options')
    options.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                         help='Show this help message and exit')

    parser.set_defaults(func=init_build)


def _convert_type(v: str) -> Any:
    if v.lower() == 'true' or v == '1':
        return True
    elif v.lower() == 'false' or v == '0':
        return False

    try:
        return float(v)
    except:
        pass

    return v


def _parse_additional_snakemake_args(arg: str) -> Dict[str, Any]:
    ret = {}
    for a in shlex.split(arg):
        if '=' not in a:
            raise ValueError("ex")

        k, v = a.split('=')
        ret[k] = _convert_type(v)

    return ret


def init_build(args):
    run_build_workflow(
        args.output_dir,
        seqs_file_list_path=args.seqs_file_list_path,
        seqs_dir_path=args.seqs_dir_path,
        k=args.k,
        base_name=args.base_name,
        build_primary_graph=args.build_primary_graph,
        annotation_formats=[_parse_annotation_format_value(af) for af in args.annotation_format],
        annotation_labels_source=args.annotation_labels_source,
        with_counts=args.with_counts,
        with_coordinates=args.with_coordinates,
        count_width=args.count_width,
        metagraph_cmd=args.metagraph_cmd,
        threads=args.threads,
        force=args.force,
        verbose=args.verbose,
        dryrun=args.dryrun,
        additional_snakemake_args=_parse_additional_snakemake_args(
            getattr(args, "additional_snakemake_args", "")
        )
    )


def main(args=tuple(sys.argv[1:])):
    parser = argparse.ArgumentParser(
        description='MetaGraph workflow utilities',
        formatter_class=_help_formatter,
    )

    subparsers = parser.add_subparsers(
        help="Available subcommands",
        required=True,
        dest="command",
        parser_class=argparse.ArgumentParser,
    )

    build_parser = subparsers.add_parser(
        "build",
        help="Build graph + annotation workflow",
        formatter_class=_help_formatter,
        add_help=False,
    )
    setup_build_parser(build_parser)

    parsed_arguments = parser.parse_args(args)

    try:
        if parsed_arguments.func:
            parsed_arguments.func(parsed_arguments)
        else:
            sys.exit("Unknown function call")
    except (ValueError, RuntimeError) as e:
        print(_format_error(str(e)), file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
