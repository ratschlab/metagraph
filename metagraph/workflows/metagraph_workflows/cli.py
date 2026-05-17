import argparse
import difflib
import importlib
import logging
import os
import shlex
import shutil
import sys
import subprocess
import time
from pathlib import Path
from typing import Iterable, Optional, Dict, Any

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


# Colorize a few flag names in --help so users can see at a glance which
# flags belong to count-aware (yellow, 33) vs coord-aware (purple, 35)
# modes. Matches the colors already used for the format-name lists in
# `_add_annotation_args`.
_FLAG_COLORS = {
    '--with-counts': '33',
    '--count-width': '33',
    '--with-coords': '35',
}


class _ColorHelpFormatter(argparse.RawTextHelpFormatter):
    """RawTextHelpFormatter that colorizes specific option strings.

    We wrap the flag in ANSI codes after argparse has already computed
    its column widths, so the alignment of the help text is unaffected
    (ANSI escape codes render at zero width in the terminal).
    """

    def _format_action(self, action):
        text = super()._format_action(action)
        if not sys.stdout.isatty():
            return text
        for flag, code in _FLAG_COLORS.items():
            if flag in action.option_strings:
                return text.replace(flag, f"\033[{code}m{flag}\033[0m", 1)
        return text


def _help_formatter(prog: str):
    return _ColorHelpFormatter(prog, width=120, max_help_position=34)


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
        # Final .dbg / .annodbg / .seqs files live at the top of
        # output_dir; per-column intermediates live in subdirectories
        # (columns.<mode>/, rd_cols.<mode>/). Build sidecars like
        # graph.dbg.pred / .succ / .anchors have a different suffix
        # so they're naturally excluded.
        if path.parent != output_dir:
            return False
        return path.suffix in ('.dbg', '.annodbg', '.seqs')

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


def _apply_annotation_options(config, annotation_formats, annotation_labels_source,
                              with_counts, with_coordinates, count_width):
    """Apply annotation-related fields to config; raise on incompatible options."""
    if annotation_labels_source:
        config['annotation_labels_source'] = annotation_labels_source.value

    selected = set(annotation_formats)
    has_count = any(af in COUNT_COMPATIBLE_FORMATS for af in selected)
    has_coord = any(af in COORD_COMPATIBLE_FORMATS for af in selected)
    effective_counts = with_counts or has_count
    effective_coords = with_coordinates or has_coord

    if effective_counts and effective_coords:
        raise ValueError(
            "Count-aware and coordinate-aware modes are mutually exclusive in this workflow. "
            "Choose either count formats (--with-counts / int_* / row_diff_int_*) or "
            "coordinate formats (--with-coords / *_coord)."
        )

    if effective_counts and not annotation_formats:
        config['annotation_formats'] = [AnnotationFormats.ROW_DIFF_INT_BRWT.value]
    elif effective_coords and not annotation_formats:
        config['annotation_formats'] = [AnnotationFormats.ROW_DIFF_BRWT_COORD.value]
    elif annotation_formats:
        config['annotation_formats'] = [af.value for af in annotation_formats]
    # else: keep whatever default.yml provided.

    if effective_counts:
        invalid = [af.value for af in annotation_formats if af not in COUNT_COMPATIBLE_FORMATS]
        if invalid:
            raise ValueError(
                "Count-aware mode is enabled (--with-counts or count-capable --anno-type), "
                "--anno-type must be one of: "
                + ", ".join(sorted([f.value for f in COUNT_COMPATIBLE_FORMATS]))
                + f". Got: {', '.join(invalid)}"
            )
    if effective_coords:
        invalid = [af.value for af in annotation_formats if af not in COORD_COMPATIBLE_FORMATS]
        if invalid:
            raise ValueError(
                "Coordinate-aware mode is enabled (--with-coords or *_coord --anno-type), "
                "--anno-type must be one of: "
                + ", ".join(sorted([f.value for f in COORD_COMPATIBLE_FORMATS]))
                + f". Got: {', '.join(invalid)}"
            )

    config['with_counts'] = effective_counts
    config['with_coordinates'] = effective_coords
    if count_width is not None:
        if not (2 <= count_width <= 32):
            raise ValueError(f"--count-width must be in range [2, 32], got {count_width}")
        if not effective_counts:
            raise ValueError(
                "--count-width can only be used with count-aware mode (--with-counts or count formats).")
        config['count_width'] = count_width


def _set_samples_config(config, samples, output_dir):
    """Set the seqs config keys from a single samples path.

    Auto-detects: directory -> SEQS_DIR_PATH; regular file ->
    SEQS_FILE_LIST_PATH. Process substitution `<(...)` shows up as a
    FIFO that becomes unreadable once the parent CLI exits, so we
    snapshot its content to ``<output_dir>/samples.txt`` and point the
    config at the snapshot.
    """
    samples_path = Path(samples).expanduser()
    if samples_path.is_dir():
        config[SEQS_DIR_PATH] = str(samples_path)
        return
    if samples_path.is_fifo() or samples_path.is_char_device():
        output_dir_path = Path(output_dir)
        output_dir_path.mkdir(parents=True, exist_ok=True)
        snapshot = output_dir_path / 'samples.txt'
        with open(samples_path, 'r') as src, open(snapshot, 'w') as dst:
            dst.write(src.read())
        config[SEQS_FILE_LIST_PATH] = str(snapshot)
        return
    if samples_path.is_file():
        config[SEQS_FILE_LIST_PATH] = str(samples_path)
        return
    raise ValueError(f"samples path not found: {samples_path}")


def _apply_runtime_options(config, threads, annotate_threads_each, metagraph_cmd,
                           disk_swap_dir, mem_gb, brwt_subsample, dryrun):
    """Apply runtime / resource fields to config; raise on invalid values."""
    config['metagraph_cmd'] = metagraph_cmd or config['metagraph_cmd']
    if not dryrun:
        _validate_metagraph_cmd(config['metagraph_cmd'])
    config['max_threads'] = threads or _default_threads_auto()
    if annotate_threads_each is not None:
        if annotate_threads_each < 1:
            raise ValueError(
                f"--anno-threads-each must be >= 1, got {annotate_threads_each}")
        config['annotate_threads_each'] = annotate_threads_each
    if disk_swap_dir is not None:
        config['tmpdir'] = str(disk_swap_dir)
    if mem_gb is not None:
        if mem_gb <= 0:
            raise ValueError(f"--mem-gb must be > 0, got {mem_gb}")
        config['max_memory_mb'] = int(mem_gb * 1024)
    if brwt_subsample is not None:
        if brwt_subsample < 1000:
            raise ValueError(f"--brwt-subsample must be >= 1000, got {brwt_subsample}")
        config['brwt_linkage_subsample'] = brwt_subsample


def run_workflow(
        output_dir: Path,
        samples: Path,
        *,
        graph: Optional[Path] = None,
        k: Optional[int] = None,
        base_name: Optional[str] = None,
        build_primary_graph: bool = False,
        annotation_formats: Iterable[AnnotationFormats] = (),
        annotation_labels_source: Optional[AnnotationLabelsSource] = None,
        with_counts: bool = False,
        with_coordinates: bool = False,
        count_width: Optional[int] = None,
        annotate_threads_each: Optional[int] = None,
        disk_swap_dir: Optional[Path] = None,
        mem_gb: Optional[float] = None,
        brwt_subsample: Optional[int] = None,
        metagraph_cmd: Optional[str] = None,
        threads: Optional[int] = None,
        force: bool = False,
        verbose: bool = False,
        dryrun: bool = False,
        additional_snakemake_args: Optional[Dict[str, Any]] = None,
) -> None:
    """Run the metagraph-workflows pipeline.

    With ``graph=None`` (default), build a fresh graph + annotation
    from ``samples``. With ``graph`` pointing at an existing .dbg file,
    skip the build pipeline and run annotation-only against that graph
    (the file is symlinked into the output dir as ``<base_name>.dbg``).
    """
    with open(default_path, 'r') as f:
        config = yaml.safe_load(f)

    _set_samples_config(config, samples, output_dir)
    config['output_directory'] = str(output_dir)
    config['base_name'] = base_name or config['base_name']

    if graph is not None:
        graph_path = Path(graph).expanduser().resolve()
        if not graph_path.exists():
            raise ValueError(f"Graph file not found: {graph_path}")
        config['external_graph'] = True
        # Per-sample primarization needs build.smk rules; not available.
        config['primarize_samples_separately'] = False
        config['build_primary_graph'] = False  # irrelevant without build
    else:
        config['k'] = k or config['k']
        config['build_primary_graph'] = build_primary_graph

    _apply_annotation_options(config, annotation_formats, annotation_labels_source,
                              with_counts, with_coordinates, count_width)
    _apply_runtime_options(config, threads, annotate_threads_each, metagraph_cmd,
                           disk_swap_dir, mem_gb, brwt_subsample, dryrun)

    if graph is not None:
        # Symlink the user's graph as <base_name>.dbg so downstream
        # rules see it as a satisfied input and snakemake doesn't try
        # to rebuild it.
        output_dir_path = Path(output_dir)
        output_dir_path.mkdir(parents=True, exist_ok=True)
        target = output_dir_path / f"{config['base_name']}.dbg"
        if target.is_symlink() or target.exists():
            target.unlink()
        target.symlink_to(graph_path)

    _invoke_snakemake(config, output_dir, threads=threads, force=force,
                      dryrun=dryrun, verbose=verbose,
                      additional_snakemake_args=additional_snakemake_args)


def _invoke_snakemake(config, output_dir, threads, force, dryrun, verbose,
                      additional_snakemake_args):
    """Write the merged config to <output_dir>/config.yaml, run snakemake,
    and emit the run summary. Raises RuntimeError on snakemake failure."""
    snakefile_path = Path(WORKFLOW_ROOT / 'Snakefile')
    output_dir_path = Path(output_dir)

    # The Snakefile reads `verbose` from config to decide whether to add
    # `-v` to each metagraph invocation.
    config['verbose'] = verbose

    if verbose:
        importlib.reload(logging)
        logging.basicConfig(format=LOGGING_FORMAT, level=logging.INFO)
        logging.info("Dumping config:")
        for k, v in sorted(config.items(), key=lambda t: t[0]):
            logging.info(f"\t{k}: {v}")

    config_file = output_dir_path / 'config.yaml'
    config_file.parent.mkdir(parents=True, exist_ok=True)
    with open(config_file, 'w') as f:
        yaml.dump(config, f)

    cmd = ['python', '-m', 'snakemake', '--snakefile', str(snakefile_path),
           '--configfile', str(config_file)]
    if force:
        cmd.append('--forceall')
    if dryrun:
        cmd.append('--dry-run')
    cmd.extend(['--cores', str(threads if threads else _default_threads_auto())])

    additional_args = additional_snakemake_args if additional_snakemake_args else {}
    for key, value in additional_args.items():
        if isinstance(value, bool):
            if value:
                cmd.append(f'--{key}')
        else:
            cmd.extend([f'--{key}', str(value)])

    # Keep stdout/stderr attached so Snakemake preserves colors and rich formatting.
    run_started_at = time.time()
    result = subprocess.run([' '.join(cmd)], shell=True)

    if result.returncode != 0:
        raise RuntimeError(_summarize_snakemake_failure_for_run(output_dir_path, run_started_at))
    _print_run_summary(output_dir_path, dryrun)



def _add_seq_input_args(group):
    """Add the positional `samples`, `-o`, and `--base-name` arguments.

    `samples` accepts either a directory (interpreted as a directory of
    sample files) or a regular file (interpreted as a text file listing
    sample paths, one per line). Process substitution `<(...)` is also
    accepted: the CLI snapshots the FIFO contents to a real file inside
    the output dir before invoking snakemake.
    """
    group.add_argument('samples', type=Path, metavar='SAMPLES',
                       help='Either a directory of sample files OR a text file listing sample\n'
                            '  paths (one per line). The type is auto-detected.')
    group.add_argument('-o', dest='output_dir', type=Path, required=True,
                       metavar='DIR',
                       help='Output directory [required]')
    group.add_argument('--base-name', default='graph', metavar='NAME',
                       help='Base output name [graph]')


def _add_annotation_args(annotation):
    """Add the shared annotation argument group (anno-source, anno-type,
    with-counts, count-width, with-coords) with help text that highlights the
    per-mode default formats inline."""
    label_sources = [v.value for v in AnnotationLabelsSource]
    count_formats = sorted([f.value for f in COUNT_COMPATIBLE_FORMATS])
    coord_formats = sorted([f.value for f in COORD_COMPATIBLE_FORMATS])
    count_formats_display = [_help_color(fmt, "33") for fmt in count_formats]
    coord_formats_display = [_help_color(fmt, "35") for fmt in coord_formats]
    plain_fmt_names = [
        "row", "bin_rel_wt", "flat", "rbfish", "brwt", "relax.brwt", "rb_brwt",
    ]
    rd_fmt_names = [
        "row_diff_brwt", "relax.row_diff_brwt",
        "row_diff_flat", "row_diff_sparse", "row_diff_disk",
    ]

    def _with_default(name: str, color: str, default_name: str) -> str:
        colored = _help_color(name, color)
        return f"[{colored}]" if name == default_name else colored

    plain_fmt_help = [_with_default(f, "36", "") for f in plain_fmt_names]
    rd_fmt_help = [_with_default(f, "36", "relax.row_diff_brwt") for f in rd_fmt_names]
    count_fmt_help = [_with_default(f, "33", "row_diff_int_brwt") for f in count_formats]
    coord_fmt_help = [_with_default(f, "35", "row_diff_brwt_coord") for f in coord_formats]
    default_count_width = _help_color("8", "33")
    zero_count_width = _help_color("0", "36")

    all_formats_help = "\n".join([
        f"    {', '.join(plain_fmt_help)}",
        f"    {', '.join(rd_fmt_help)}",
        f"    {', '.join(count_fmt_help)}",
        f"    {', '.join(coord_fmt_help)}",
    ])

    annotation.add_argument('--anno-source',
                            dest='annotation_labels_source',
                            type=AnnotationLabelsSource,
                            default=AnnotationLabelsSource.FILENAME,
                            metavar='SOURCE',
                            help=f"Column label source: {', '.join(label_sources)} [filename]\n"
                                 "  ")
    annotation.add_argument('--anno-type', action='append',
                            dest='annotation_format',
                            default=[],
                            metavar='FORMAT',
                            help=f"Annotation format (can be used multiple times).\n"
                                 f"{all_formats_help}\n"
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
                                 f"  Supported for {coord_formats_display[0]}, {coord_formats_display[2]}, {coord_formats_display[1]}, {coord_formats_display[3]}")


def _add_workflow_args(workflow):
    """Add the shared `other` argument group (threads, disk/mem, force,
    verbose, dryrun, metagraph-cmd, extra-args)."""
    workflow.add_argument('--threads', type=int, default=None, metavar='N',
                          help='Maximum CPU cores to use [num_cores]')
    workflow.add_argument('--disk-swap-dir', dest='disk_swap_dir', type=Path, default=None,
                          metavar='DIR',
                          help='Directory for on-disk buffers; omit to stay in RAM [none]')
    workflow.add_argument('--mem-gb', type=float, default=None,
                          metavar='GB',
                          help='Approximate RAM budget in GB; used to derive --mem-cap-gb for each stage [16]')
    workflow.add_argument('--anno-threads-each', dest='annotate_threads_each',
                          type=int, default=None, metavar='N',
                          help='Threads used to annotate each input file. Parallel columns = ceil(--threads / N);\n'
                               '  raise N to give each column more --mem-cap-gb buffer.\n'
                               '  [8 for binary/counts, 16 for coords]')
    workflow.add_argument('--brwt-subsample', type=int, default=None, metavar='N',
                          help='Number of bits subsampled for distance estimation when clustering BRWT\n'
                               '  columns (passed as --subsample to transform_anno --anno-type *_brwt*). [1000000]')
    workflow.add_argument('--force', default=False, action='store_true',
                          help='Force re-run all rules [False]')
    workflow.add_argument('-v', '--verbose', default=False, action='store_true',
                          help='Print verbose config/runtime logs and pass -v to each\n'
                               '  underlying metagraph invocation [False]')
    workflow.add_argument('--dryrun', default=False, action='store_true',
                          help='Render DAG and config only; do not execute rules [False]')
    workflow.add_argument('--metagraph-cmd', type=str, default=None, metavar='CMD',
                          help='Path/command for metagraph executable [metagraph from PATH]')
    workflow.add_argument('--extra-args', dest='additional_snakemake_args', metavar='ARGS', type=str, default='',
                          help='Extra arguments to pass to snakemake [none]\n'
                               '  Example: --extra-args="arg1=val1 arg2=val2"')


def _add_help_arg(parser):
    options = parser.add_argument_group('options')
    options.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                         help='Show this help message and exit')


def setup_build_parser(parser):
    parser.description = (
        "Build a MetaGraph graph + annotation workflow from a sequence list or directory.\n"
        "\n"
        "Inputs are assumed to be contigs with deduplicated k-mers (one fasta.gz\n"
        "per sample). When building a primary graph (--primary), they must be\n"
        "primary contigs. For --with-coords, inputs must instead be the full,\n"
        "non-deduplicated samples."
    )
    parser.epilog = (
        "Examples:\n"
        "  metagraph-workflows build samples_dir/ -k 31 -o out/\n"
        "  metagraph-workflows build <(ls /data/samples/*.fa) -k 31 -o out/\n"
        "  metagraph-workflows build files.txt --with-counts -o out/\n"
        "  metagraph-workflows build samples_dir/ --with-coords -o out/"
    )

    _add_seq_input_args(parser.add_argument_group('input/output'))

    graph = parser.add_argument_group('graph')
    graph.add_argument('--graph', type=Path, default=None, metavar='PATH',
                       help='Reuse an existing .dbg graph instead of building one from SAMPLES.\n'
                            '  Skips the build pipeline; runs annotation + row-diff transforms only.')
    graph.add_argument('-k', type=int, default=31, metavar='K',
                       help='k-mer length [31]')
    graph.add_argument('--primary', dest='build_primary_graph', default=False,
                       action='store_true',
                       help='Build canonical graph first, then derive/build primary graph [False]')

    _add_annotation_args(parser.add_argument_group('annotation'))
    _add_workflow_args(parser.add_argument_group('other'))
    _add_help_arg(parser)

    parser.set_defaults(func=init_build)


def _convert_type(v: str) -> Any:
    if v.lower() == 'true' or v == '1':
        return True
    if v.lower() == 'false' or v == '0':
        return False
    try:
        return float(v)
    except ValueError:
        return v


def _parse_additional_snakemake_args(arg: str) -> Dict[str, Any]:
    ret = {}
    for a in shlex.split(arg):
        if '=' not in a:
            raise ValueError(
                f"--extra-args expects key=value tokens; got {a!r}. "
                f"Example: --extra-args=\"keep-going=True jobs=4\""
            )
        k, v = a.split('=', 1)
        ret[k] = _convert_type(v)
    return ret


def init_build(args):
    run_workflow(
        output_dir=args.output_dir,
        samples=args.samples,
        graph=args.graph,
        k=args.k,
        base_name=args.base_name,
        build_primary_graph=args.build_primary_graph,
        annotation_formats=[_parse_annotation_format_value(af) for af in args.annotation_format],
        annotation_labels_source=args.annotation_labels_source,
        with_counts=args.with_counts,
        with_coordinates=args.with_coordinates,
        count_width=args.count_width,
        annotate_threads_each=args.annotate_threads_each,
        disk_swap_dir=args.disk_swap_dir,
        mem_gb=args.mem_gb,
        brwt_subsample=args.brwt_subsample,
        metagraph_cmd=args.metagraph_cmd,
        threads=args.threads,
        force=args.force,
        verbose=args.verbose,
        dryrun=args.dryrun,
        additional_snakemake_args=_parse_additional_snakemake_args(args.additional_snakemake_args),
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
