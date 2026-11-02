#!/usr/bin/env python3
"""
M7 ParaView State (.pvsm) Portability Remapper & Verifier (Root CLI Wrapper)
============================================================================
Delegates to reference_dataset/openfoam_case/remap_pvsm_paths.py so the tool
can be invoked from either the M7Multiscale root or reference_dataset/openfoam_case/.

Strictly zero em-dashes and zero en-dashes across all code and output.
"""

import os
import sys

OPENFOAM_CASE_DIR = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "reference_dataset",
    "openfoam_case",
)
if OPENFOAM_CASE_DIR not in sys.path:
    sys.path.insert(0, OPENFOAM_CASE_DIR)

from remap_pvsm_paths import (  # noqa: E402
    REQUIRED_RELATIVE_READER_FILES,
    main,
    remap_pvsm_content,
    remap_pvsm_file,
    sha256_prefixed,
    verify_pvsm_portability,
)

if __name__ == "__main__":
    sys.exit(main())
