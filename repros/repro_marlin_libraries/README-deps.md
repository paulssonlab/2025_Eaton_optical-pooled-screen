## Reproducibility: software + OS

This repo uses **micromamba/conda** environments. For each environment YAML in the repository, we record:

- Fully pinned export (includes versions and build strings): `*-conda-env-export.yml`
- Package list (includes versions): `*-conda-list.txt`
- Explicit lockfile (fully pinned, rebuildable): `*-conda-explicit.lock`
- pip freeze (if pip is present): `*-pip-freeze.txt`

### Operating system / platform

- Platform details: `repros/repro_marlin_libraries/platform.txt`

### Environments discovered in subdirectories

| Env YAML path | Conda env name | Exports |
|---|---:|---|
| `./Library_Design/environment.yml` | `marlin_libraries` | ./repros/repro_marlin_libraries/envs/marlin_libraries_environment-conda-env-export.yml<br>./repros/repro_marlin_libraries/envs/marlin_libraries_environment-conda-list.txt<br>./repros/repro_marlin_libraries/envs/marlin_libraries_environment-conda-explicit.lock<br>./repros/repro_marlin_libraries/envs/marlin_libraries_environment-pip-freeze.txt |
