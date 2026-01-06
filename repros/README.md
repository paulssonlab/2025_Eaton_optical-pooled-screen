## Reproducibility: Software + OS

This repo uses **micromamba/conda** environments. For each environment YAML in the repository, we record:

- Fully pinned export (includes versions and build strings): `*-conda-env-export.yml`
- Package list (includes versions): `*-conda-list.txt`
- Explicit lockfile (fully pinned, rebuildable): `*-conda-explicit.lock`
- pip freeze (if pip is present): `*-pip-freeze.txt`

### Operating system / platform

- Platform details: `repros/platform.txt`

### Environments discovered in subdirectories

| Env YAML path | Conda env name | Exports |
|---|---:|---|
| `./Image_Analysis/environment.yml` | `crispri` | ./repros/repro_crispri/envs/crispri_environment-conda-env-export.yml<br>./repros/repro_crispri/envs/crispri_environment-conda-list.txt<br>./repros/repro_crispri/envs/crispri_environment-conda-explicit.lock<br>./repros/repro_crispri/envs/crispri_environment-pip-freeze.txt |
| `./Library_Design/environment.yml` | `marlin_libraries` | ./repros/repro_marlin_libraries/envs/marlin_libraries_environment-conda-env-export.yml<br>./repros/repro_marlin_libraries/envs/marlin_libraries_environment-conda-list.txt<br>./repros/repro_marlin_libraries/envs/marlin_libraries_environment-conda-explicit.lock<br>./repros/repro_marlin_libraries/envs/marlin_libraries_environment-pip-freeze.txt |
| `./Replication_Runout/environment.yml` | `flowcyto` | ./repros/repro_fluocyto/envs/flowcyto_environment-conda-env-export.yml<br>./repros/repro_fluocyto/envs/flowcyto_environment-conda-list.txt<br>./repros/repro_fluocyto/envs/flowcyto_environment-conda-explicit.lock<br>./repros/repro_fluocyto/envs/flowcyto_environment-pip-freeze.txt |
| `./Agar_Pad_Image_Analysis/environment.yml` | `agarpad` | ./repros/repro_agarpad/envs/agarpad_environment-conda-env-export.yml<br>./repros/repro_agarpad/envs/agarpad_environment-conda-list.txt<br>./repros/repro_agarpad/envs/agarpad_environment-conda-explicit.lock<br>./repros/repro_agarpad/envs/agarpad_environment-pip-freeze.txt |
| `./Sequencing/environment.yml` | `nanopore` | ./repros/repro_nanopore/envs/nanopore_environment-conda-env-export.yml<br>./repros/repro_nanopore/envs/nanopore-conda-list.txt<br>./repros/repro_nanopore/envs/nanopore-conda-explicit.lock<br>./repros/repro_nanopore/envs/nanopore-pip-freeze.txt |