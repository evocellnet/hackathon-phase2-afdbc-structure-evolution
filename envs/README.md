# Environment `conda_envs/main`

## Activate
    source /opt/conda/etc/profile.d/conda.sh && conda activate /scratch/spascare/agentic-ai/workspaces/phase2-afdbc-evolution/conda_envs/main

## Rebuild
    mamba env create -p conda_envs/main -f envs/main.yml      # conda part + pip (foldcomp, torch)
    # torch must come from the CUDA wheel index:
    pip install torch==2.11.0 --index-url https://download.pytorch.org/whl/cu128
Exact pins: `envs/main.explicit.txt` (`mamba create -p conda_envs/main --file envs/main.explicit.txt`, then pip lines above).
Created with: `mamba create -p conda_envs/main -c conda-forge -c bioconda --override-channels python=3.12 <tools...>`; pip: `foldcomp`, `torch` (cu128).
Machine: A100 40GB, driver 580.178.04 (CUDA 13.0 capable); torch wheel cu128 works.

## Versions
| Tool | Version | Tool | Version |
|---|---|---|---|
| foldseek (GPU-capable) | 10.941cd33 | mmseqs2 | 18.8cc5c |
| foldmason | 4.dd3c235 | USalign / TMalign | 20241201 / 20240303 |
| mafft | 7.526 | famsa | 2.4.1 |
| hmmer | 3.4 | iqtree (`iqtree`, `iqtree3`) | 3.1.4 |
| fasttree | 2.2.0 | diamond | 2.2.8 |
| seqkit | 2.14.0 | taxonkit | 0.20.0 |
| csvtk | 0.38.0 | pigz / aria2 / parallel | 2.8 / 1.37.0 / 20260922 |
| python | 3.12 | numpy / scipy / pandas | 2.5.4 / 1.18.1 / 3.1.0 |
| polars / pyarrow / duckdb | 2.1.0 / 25.0.0 / 2.0.0 | statsmodels / scikit-learn | 0.15.0 / 1.9.1 |
| matplotlib / seaborn | 3.11.2 / 0.13.2 | biopython / ete3 / dendropy | 1.88 / 3.1.3 / 5.1.3 |
| networkx / tqdm / joblib | 3.7 / 4.70.1 / 1.6.0 | gemmi / tmtools / foldcomp | 0.7.5 / 0.3.0 / 1.0.0 |
| requests / pyyaml / numba | 2.34.2 / 6.0.3 / 0.68.0 | torch | 2.11.0+cu128 (cuda True) |

## Notes
- IQ-TREE 3 only ships the `iqtree3`/`iqtree` binaries (no `iqtree2`).
- foldseek GPU: `foldseek makepaddedseqdb db dbpad` then `foldseek search q dbpad res tmp --gpu 1` (verified on A100).
- AlphaFold DB files are now `AF-<acc>-F1-model_v6.pdb`; `_v4` returns 404. The EBI server returns 403 to curl's default
  user agent: pass `-A "Mozilla/5.0"`.
- Package cache: `conda_envs/pkgs_cache` (CONDA_PKGS_DIRS), pip cache `conda_envs/pip_cache`.
