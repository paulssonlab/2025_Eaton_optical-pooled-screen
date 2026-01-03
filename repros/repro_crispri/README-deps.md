## Reproducibility: software + OS

This repo uses **micromamba/conda** environments. For each environment YAML in the repository, we record:

- Fully pinned export (includes versions and build strings): `*-conda-env-export.yml`
- Package list (includes versions): `*-conda-list.txt`
- Explicit lockfile (fully pinned, rebuildable): `*-conda-explicit.lock`
- pip freeze (if pip is present): `*-pip-freeze.txt`

### Operating system / platform

- Platform details: `repro/platform.txt`

### Environments discovered in subdirectories

| Env YAML path | Conda env name | Exports |
|---|---:|---|
| `/home/de64/2025_Eaton_optical-pooled-screen/Image_Analysis/environment.yml` | `crispri` | repro/envs/__home__de64__2025_Eaton_optical-pooled-screen__Image_Analysis__environment.yml-conda-env-export.yml<br>repro/envs/__home__de64__2025_Eaton_optical-pooled-screen__Image_Analysis__environment.yml-conda-list.txt<br>repro/envs/__home__de64__2025_Eaton_optical-pooled-screen__Image_Analysis__environment.yml-conda-explicit.lock<br>repro/envs/__home__de64__2025_Eaton_optical-pooled-screen__Image_Analysis__environment.yml-pip-freeze.txt |
