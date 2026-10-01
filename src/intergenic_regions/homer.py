"""Optional, explicit HOMER execution using the extracted FASTA background."""

import shutil
import subprocess
from pathlib import Path
from typing import Any

from intergenic_regions.background import check_sequence_sets
from intergenic_regions.io import output_bundle, read_fasta, write_json


def homer_command(
    *,
    positive: Path,
    negative: Path,
    output: Path,
    executable: str = "findMotifs.pl",
    lengths: tuple[int, ...] = (8, 10, 12),
    threads: int = 1,
    known_motifs: Path | None = None,
    known_only: bool = False,
) -> list[str]:
    """Build HOMER's FASTA-mode argument vector without a shell.

    Args:
        positive: Foreground FASTA.
        negative: Explicit background FASTA.
        output: HOMER results directory.
        executable: HOMER findMotifs.pl command or full path.
        lengths: De novo motif widths.
        threads: Worker count.
        known_motifs: Optional HOMER-formatted known motif database.
        known_only: Skip de novo motif discovery.

    Returns:
        A subprocess argument list retaining literal path characters.

    Raises:
        ValueError: Widths, threads or executable are invalid.
    """
    if (
        threads < 1
        or not lengths
        or any(k < 2 or k > 50 for k in lengths)
        or not executable
    ):
        raise ValueError("Invalid HOMER settings")
    command = [
        executable,
        str(positive.resolve()),
        "fasta",
        str(output.resolve()),
        "-fasta",
        str(negative.resolve()),
        "-len",
        ",".join(str(k) for k in lengths),
        "-p",
        str(threads),
    ]
    if known_motifs is not None:
        command += ["-mknown", str(known_motifs.resolve())]
    if known_only:
        command.append("-nomotif")
    return command


def run_homer(
    *,
    positive: Path,
    negative: Path,
    output: Path,
    executable: str = "findMotifs.pl",
    lengths: tuple[int, ...] = (8, 10, 12),
    threads: int = 1,
    known_motifs: Path | None = None,
    known_only: bool = False,
    timeout: float = 3600,
    dry_run: bool = False,
) -> dict[str, Any]:
    """Run HOMER, or save a clearly labelled command plan without executing.

    Args:
        positive: Foreground FASTA.
        negative: Explicit control FASTA.
        output: New output directory.
        executable: HOMER executable.
        lengths: De novo widths.
        threads: Worker count.
        known_motifs: Optional known-motif file.
        known_only: Disable de novo discovery.
        timeout: Maximum execution time in seconds.
        dry_run: Save the planned arguments only.

    Returns:
        Execution metadata with an explicit completed/planned status.

    Raises:
        ValueError: Inputs or settings are invalid.
        FileNotFoundError: HOMER is unavailable when execution is requested.
        RuntimeError: HOMER fails or times out.
    """
    check_sequence_sets(
        positive=read_fasta(path=positive), negative=read_fasta(path=negative)
    )
    if timeout <= 0:
        raise ValueError("HOMER timeout must be positive")
    if known_motifs is not None and not known_motifs.is_file():
        raise FileNotFoundError(known_motifs)
    resolved = shutil.which(cmd=executable)
    if not dry_run and resolved is None:
        raise FileNotFoundError(
            "HOMER findMotifs.pl is not installed; "
            "install HOMER or use --dry-run"
        )
    with output_bundle(path=output) as stage:
        results = stage / "results"
        command = homer_command(
            positive=positive,
            negative=negative,
            output=results,
            executable=resolved or executable,
            lengths=lengths,
            threads=threads,
            known_motifs=known_motifs,
            known_only=known_only,
        )
        metadata = {
            "status": "planned" if dry_run else "completed",
            "command": command,
            "threads": threads,
            "timeout_seconds": timeout,
        }
        if not dry_run:
            with (stage / "homer.log").open(mode="w", encoding="utf-8") as log:
                try:
                    completed = subprocess.run(
                        args=command,
                        stdin=subprocess.DEVNULL,
                        stdout=log,
                        stderr=subprocess.STDOUT,
                        shell=False,
                        check=False,
                        timeout=timeout,
                    )
                except subprocess.TimeoutExpired as exc:
                    raise RuntimeError(
                        f"HOMER exceeded timeout of {timeout} seconds"
                    ) from exc
                if completed.returncode:
                    tail = (stage / "homer.log").read_text(
                        encoding="utf-8", errors="replace"
                    )[-4000:]
                    raise RuntimeError(
                        f"HOMER exited with {completed.returncode}: {tail}"
                    )
                if not results.is_dir():
                    raise RuntimeError(
                        "HOMER returned success without producing results"
                    )
        # Record durable output paths rather than temporary staging paths.
        metadata["command"] = homer_command(
            positive=positive,
            negative=negative,
            output=output / "results",
            executable=resolved or executable,
            lengths=lengths,
            threads=threads,
            known_motifs=known_motifs,
            known_only=known_only,
        )
        write_json(path=stage / "command.json", data=metadata)
    return metadata
