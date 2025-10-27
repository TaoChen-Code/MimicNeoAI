# coding=utf-8
"""
Download and extract the MimicNeoAI database (manual operation).
Default target: mimicneoai/database/
Optional: specify a custom path for download and extraction.
If a custom path is used, a symbolic link will be created at mimicneoai/database.
"""

import os
import sys
import tarfile
import argparse
import shutil
from pathlib import Path
import urllib.request
import subprocess
import time

URL = "https://mimicneoai.biostacs.com/database/MimicNeoAI_database_v1.0.tar.gz"
DEFAULT_DIR = Path(__file__).resolve().parent / "database"
ARCHIVE_NAME = "MimicNeoAI_database_v1.0.tar.gz"

def download_file(url: str, dest: Path, retries: int = 3, backoff: float = 1.6):
    """
    Download file to `dest` using wget (preferred), falling back to curl, then urllib.
    - Resumable (-c / -C -) and retriable.
    - Creates parent dirs.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"[INFO] Downloading from {url}")
    # Prefer wget
    if shutil.which("wget"):
        cmd = [
            "wget", "-c",               # resume partial download
            "--tries=3",                # internal retries
            "--timeout=60",
            "--retry-connrefused",
            "--show-progress",
            "-O", str(dest),
            url,
        ]
        for i in range(retries):
            try:
                subprocess.run(cmd, check=True)
                print(f"[OK] Download completed: {dest}")
                return
            except subprocess.CalledProcessError as e:
                if i == retries - 1:
                    raise RuntimeError(f"wget failed: {e}") from e
                sleep_s = backoff ** i
                print(f"[WARN] wget failed (attempt {i+1}/{retries}), retry in {sleep_s:.1f}s...")
                time.sleep(sleep_s)

    # Fallback: curl (also resumable)
    if shutil.which("curl"):
        cmd = [
            "curl", "-fL", "-C", "-",   # -C - resume, -L follow redirects
            "--retry", "3",
            "--retry-delay", "2",
            "--connect-timeout", "30",
            "-A", "Wget/1.21 (compatible; MimicNeoAI/1.0)",  # friendlier UA for strict servers
            "-o", str(dest),
            url,
        ]
        for i in range(retries):
            try:
                subprocess.run(cmd, check=True)
                print(f"[OK] Download completed: {dest}")
                return
            except subprocess.CalledProcessError as e:
                if i == retries - 1:
                    raise RuntimeError(f"curl failed: {e}") from e
                sleep_s = backoff ** i
                print(f"[WARN] curl failed (attempt {i+1}/{retries}), retry in {sleep_s:.1f}s...")
                time.sleep(sleep_s)

    # Last resort: urllib (adds headers to avoid 406)
    req = urllib.request.Request(
        url,
        headers={
            "User-Agent": "MimicNeoAI/1.0 (+wget-fallback)",
            "Accept": "*/*",
        },
    )
    for i in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=60) as resp, open(dest, "wb") as out:
                shutil.copyfileobj(resp, out)
            print(f"[OK] Download completed: {dest}")
            return
        except Exception as e:
            if i == retries - 1:
                raise
            sleep_s = backoff ** i
            print(f"[WARN] urllib failed (attempt {i+1}/{retries}), retry in {sleep_s:.1f}s...")
            time.sleep(sleep_s)



def extract_archive(archive_path: Path, extract_to: Path):
    """Extract .tar.gz archive to a target directory."""
    print(f"[INFO] Extracting {archive_path} → {extract_to}")
    extract_to.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive_path, "r:gz") as tar:
        tar.extractall(path=extract_to)
    print(f"[OK] Extraction completed.")


def create_symlink(target: Path, link_path: Path):
    """Create or update a symbolic link to the extracted database."""
    if link_path.is_symlink():
        print(f"[WARN] Removing existing symlink: {link_path}")
        link_path.unlink()

    elif link_path.exists():
        if any(link_path.iterdir()):
            print(f"[ERROR] '{link_path}' already exists and is not empty. "
                  f"Please remove or rename it before continuing.")
            return
        else:
            print(f"[WARN] Removing empty directory: {link_path}")
            link_path.rmdir()

    # ensure parent dir exists
    link_path.parent.mkdir(parents=True, exist_ok=True)
    link_path.symlink_to(target.resolve(), target_is_directory=True)
    print(f"[OK] Created symbolic link: {link_path} → {target}")



def main():
    parser = argparse.ArgumentParser(
        description="Download and extract MimicNeoAI database manually."
    )
    parser.add_argument(
        "--target-dir",
        type=str,
        default=None,
        help="Custom path to store and extract the database (default: mimicneoai/database)",
    )
    args = parser.parse_args()

    target_dir = Path(args.target_dir).resolve() if args.target_dir else DEFAULT_DIR.resolve()
    archive_path = target_dir / ARCHIVE_NAME

    print(f"[INFO] Target directory: {target_dir}")

    # Step 1: download
    if not archive_path.exists():
        download_file(URL, archive_path)
    else:
        print(f"[SKIP] Archive already exists: {archive_path}")

    # Step 2: extract
    extract_archive(archive_path, target_dir)

    # Step 3: create symlink if custom path is used
    if args.target_dir:
        create_symlink(target_dir, DEFAULT_DIR)

    print("\n[DONE] Database setup complete.")
    print(f"  - Downloaded archive: {archive_path}")
    print(f"  - Extracted to: {target_dir}")
    if args.target_dir:
        print(f"  - Linked: {DEFAULT_DIR} → {target_dir}")


if __name__ == "__main__":
    main()
