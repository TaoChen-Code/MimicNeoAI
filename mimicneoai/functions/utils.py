import re
from typing import Tuple

_MEM_RE = re.compile(r"^\s*(\d+)\s*([KMGkmg]?)\s*$")


def format_java_heap(mem: str, xms_ratio: float = 0.5) -> Tuple[str, str, str]:
    """
    Normalize a user-provided heap size string and derive JVM -Xmx/-Xms.

    You want:
      - mem  : output (normalized) e.g. "128G"
      - xmx/xms: inputs to your command string, but computed here

    Args:
        mem: Heap size in a compact form like "128G", "128g", "64000M", "800K", or "128".
             If unit omitted, defaults to "G".
        xms_ratio: Xms as a fraction of Xmx. Default 0.5 (i.e., Xms = 50% of Xmx).
                   Must be in (0, 1].

    Returns:
        (mem_norm, xmx, xms) where each is like "128G", "64G".

    Raises:
        ValueError: If mem format is invalid or ratio is out of range.
    """
    if not (0 < xms_ratio <= 1.0):
        raise ValueError(f"xms_ratio must be in (0, 1], got {xms_ratio}")

    s = str(mem)
    m = _MEM_RE.match(s)
    if not m:
        raise ValueError(f"Invalid mem format: {mem!r}. Expected like '128G', '64000M', '128'.")

    mem_num = int(m.group(1))
    mem_unit = (m.group(2) or "G").upper()

    # normalize mem
    mem_norm = f"{mem_num}{mem_unit}"
    xmx = mem_norm

    # compute xms in the same unit (integer)
    xms_num = max(1, int(mem_num * xms_ratio))
    xms = f"{xms_num}{mem_unit}"

    return mem_norm, xmx, xms
# ---- Example usage ----
# mem_norm, xmx, xms = format_java_heap(configure["args"]["mem"], xms_ratio=0.5)
# gatk = f"java -Xms{xms} -Xmx{xmx} -jar {gatk_jar}"