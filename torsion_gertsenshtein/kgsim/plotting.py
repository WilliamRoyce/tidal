from __future__ import annotations

import pathlib
import shutil

import matplotlib as mpl

mpl.use("Agg")
from matplotlib import animation
from matplotlib.animation import FFMpegWriter, PillowWriter


def choose_writer_and_out(
    snap_count: int, t_end: float, out_path: str | pathlib.Path, fps: int | None = None
) -> tuple[FFMpegWriter | PillowWriter, str, str]:
    """
    Select an appropriate matplotlib animation writer and output filename.

    The caller provides out_path (without extension) and this function will
    append the appropriate extension (".mp4" for ffmpeg, ".gif" for pillow).

    Parameters
    ----------
    snap_count : int
        Number of frames/snapshots available for the animation.
    t_end : float
        Total simulation time (same units as snapshot spacing). Used to compute a
        target frames-per-second (fps) for playback.

    Returns
    -------
    tuple[Any, str, str]
        A tuple containing:
        - writer: an instance of a matplotlib.animation writer (FFMpegWriter or PillowWriter).
        - out: str, path to the output file with extension appended.
        - use_writer: str, identifier of the writer used ("ffmpeg" or "pillow").

    Behavior
    --------
    - Computes fps as the number of snapshots divided by (t_end / 5.0), coerced to an int
      and clamped to at least 1 to avoid zero or fractional fps.
    - Prefers ffmpeg if matplotlib.animation reports "ffmpeg" as available and the
      ffmpeg executable is present on the system; in that case returns an FFMpegWriter
      configured with metadata and a bitrate.
    - Falls back to PillowWriter when ffmpeg is not available.
    - Chooses output filename and a string label for which writer was chosen.
    """
    if fps is None:
        fps = max(1, int(snap_count / max(1.0, t_end / 5.0)))

    p = pathlib.Path(out_path)
    base = str(
        p.with_suffix("")
    )  # strip any existing suffix, we'll append the chosen one
    if animation.writers.is_available("ffmpeg") and shutil.which("ffmpeg") is not None:
        writer = animation.FFMpegWriter(
            fps=fps, metadata={"artist": "kgsim"}, bitrate=2000
        )
        out = base + ".mp4"
        use_writer = "ffmpeg"
    else:
        writer = PillowWriter(fps=fps)
        out = base + ".gif"
        use_writer = "pillow"
    return writer, out, use_writer
