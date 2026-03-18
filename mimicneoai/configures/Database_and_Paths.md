# Database and Paths

All pipelines use the same reference bundle and `paths.yaml`.

## Automatic Setup

On first run, MimicNeoAI can download and initialize the required database bundle automatically.
`mimicneoai/configures/paths.yaml` is then used to locate tools/containers/references.

## Manual Download (optional)

Use the built-in helper if you want to control download location:

```bash
# Default location (under repository)
python -m mimicneoai.download_database
# or
mimicneoai download_database

# Custom location
python -m mimicneoai.download_database --target-dir /mnt/data/MimicNeoAI_DB
# or
mimicneoai download_database --target-dir /mnt/data/MimicNeoAI_DB
```

If `--target-dir` is used, MimicNeoAI links the extracted database back to the expected project location.

## paths.yaml

Default file:
- `mimicneoai/configures/paths.yaml`

You usually do not need to edit it unless:
- Your tools/containers are installed in custom locations.
- You want to point to custom references.

When needed, pass a custom paths file:

```bash
mimicneoai <pipeline> -c <configure.yaml> -p /path/to/paths.yaml
```
