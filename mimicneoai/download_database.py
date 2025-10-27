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

URL = "https://mimicneoai.biostacs.com/database/MimicNeoAI_database_v1.0.tar.gz"
DEFAULT_DIR = Path(__file__).resolve().parent / "database"
ARCHIVE_NAME = "MimicNeoAI_database_v1.0.tar.gz"


def download_file(url: str, dest: Path):
    """Download file from URL to destination."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"[INFO] Downloading from {url}")
    with urllib.request.urlopen(url) as response, open(dest, "wb") as out_file:
        shutil.copyfileobj(response, out_file)
    print(f"[OK] Download completed: {dest}")


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
