#!/usr/bin/env python3
"""
M7 ParaView State (.pvsm) Portability Remapper & Verifier
=========================================================
Strips hardcoded absolute host paths (such as /Users/.../M7Multiscale/...)
from ParaView .pvsm state files and verification reports, replacing them with
portable relative paths (or an optional user-specified target directory).

Strictly zero em-dashes and zero en-dashes across all code and output.
"""

import argparse
import hashlib
import json
import os
import re
import sys

REQUIRED_RELATIVE_READER_FILES = (
    "m7_topological_timeseries.pvd",
    "case.foam",
    "m7_topology_tubes.pvd",
    "m7_topological_timeseries.xdmf",
)


def sha256_prefixed(filepath):
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while True:
            chunk = f.read(65536)
            if not chunk:
                break
            h.update(chunk)
    return "sha256:" + h.hexdigest()


def remap_pvsm_content(xml_text, target_dir=""):
    """
    Replace absolute paths inside <Element ... value="..." /> tags and
    any /Users/.../reference_dataset/openfoam_case/ prefixes with portable
    relative paths (or target_dir/<filename> when target_dir is provided).
    """
    prefix = target_dir.rstrip("/\\") + "/" if target_dir else ""

    # 1. Strip any /Users/.../openfoam_case/ prefix
    remapped = re.sub(
        r"/Users/[^\"'<>\s]+/openfoam_case/([^\"'<>\s]+)",
        lambda m: prefix + m.group(1),
        xml_text,
    )

    # 2. Strip any remaining /Users/.../M7Multiscale/ prefix
    remapped = re.sub(
        r"/Users/[^\"'<>\s]+/M7Multiscale/([^\"'<>\s]+)",
        lambda m: prefix + m.group(1),
        remapped,
    )

    # 3. Strip any generic /Users/.../<reader_file> on known reader filenames
    for reader_fn in REQUIRED_RELATIVE_READER_FILES:
        escaped_fn = re.escape(reader_fn)
        remapped = re.sub(
            r"/Users/[^\"'<>\s]+/" + escaped_fn,
            prefix + reader_fn,
            remapped,
        )

    return remapped


def verify_pvsm_portability(pvsm_path):
    if not os.path.isfile(pvsm_path):
        raise FileNotFoundError(f"ParaView state file not found: {pvsm_path}")
    with open(pvsm_path, "r", encoding="utf-8") as f:
        txt = f.read()

    users_matches = re.findall(r"/Users/[^\"'<>\s]+", txt)
    has_em_or_en_dash = ("\u2014" in txt) or ("\u2013" in txt)
    readers_present = {
        fn: (f'value="{fn}"' in txt or f"/{fn}\"" in txt)
        for fn in REQUIRED_RELATIVE_READER_FILES
    }
    portable_relative = all(f'value="{fn}"' in txt for fn in REQUIRED_RELATIVE_READER_FILES)

    return {
        "pvsm_path": os.path.basename(pvsm_path),
        "zero_users_absolute_paths": len(users_matches) == 0,
        "users_path_occurrences": len(users_matches),
        "portable_relative_readers": portable_relative,
        "readers_present": readers_present,
        "zero_em_en_dashes": not has_em_or_en_dash,
        "pvsm_sha256": sha256_prefixed(pvsm_path),
    }


def remap_pvsm_file(pvsm_path, target_dir="", report_json_path=None):
    with open(pvsm_path, "r", encoding="utf-8") as f:
        orig_txt = f.read()

    new_txt = remap_pvsm_content(orig_txt, target_dir=target_dir)
    if new_txt != orig_txt:
        with open(pvsm_path, "w", encoding="utf-8") as f:
            f.write(new_txt)

    if report_json_path is None:
        cand_report = os.path.join(os.path.dirname(os.path.abspath(pvsm_path)), "m7_paraview_verification_report.json")
        if os.path.isfile(cand_report):
            report_json_path = cand_report

    if report_json_path and os.path.isfile(report_json_path):
        with open(report_json_path, "r", encoding="utf-8") as f:
            rep = json.load(f)
        rep["case_dir"] = target_dir if target_dir else "."
        rep["pvsm_sha256"] = sha256_prefixed(pvsm_path)
        rep["pvsm_portable_relative_paths"] = (target_dir == "")
        if "step_summaries" in rep and isinstance(rep["step_summaries"], list):
            for step in rep["step_summaries"]:
                if "screenshot_path" in step and isinstance(step["screenshot_path"], str):
                    bname = os.path.basename(step["screenshot_path"])
                    step["screenshot_path"] = (
                        os.path.join(target_dir, "renders", bname)
                        if target_dir
                        else f"renders/{bname}"
                    )
        with open(report_json_path, "w", encoding="utf-8") as f:
            json.dump(rep, f, indent=2)
            f.write("\n")

    return verify_pvsm_portability(pvsm_path)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Strip hardcoded /Users/ paths from M7 ParaView .pvsm state files and verify portability."
    )
    default_pvsm = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "m7_topological_glmm_timeseries.pvsm",
    )
    if not os.path.isfile(default_pvsm):
        alt_pvsm = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "reference_dataset",
            "openfoam_case",
            "m7_topological_glmm_timeseries.pvsm",
        )
        if os.path.isfile(alt_pvsm):
            default_pvsm = alt_pvsm

    parser.add_argument("--pvsm", default=default_pvsm, help="Path to .pvsm file")
    parser.add_argument(
        "--target-dir",
        default="",
        help="Optional target directory prefix (leave empty for portable relative paths)",
    )
    parser.add_argument(
        "--check-only",
        action="store_true",
        help="Only verify portability without modifying the file",
    )
    args = parser.parse_args(argv)

    if args.check_only:
        res = verify_pvsm_portability(args.pvsm)
    else:
        res = remap_pvsm_file(args.pvsm, target_dir=args.target_dir)

    print(json.dumps(res, indent=2))
    if not args.target_dir and not res["zero_users_absolute_paths"]:
        sys.exit(1)
    if not res["zero_em_en_dashes"]:
        sys.exit(1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
