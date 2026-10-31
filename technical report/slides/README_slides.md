# Slide deck

- `project_slides.pptx`: 47-slide deck documenting the hackathon, prompt, decision process, method explainers, hypotheses, results, take-home messages and lessons. `project_slides.pdf` is a PDF export.
- Rebuild:
  ```
  npm install pptxgenjs
  NODE_PATH=$PWD/node_modules PPTX_SKILL_DIR=<pptx skill dir> node build_deck.js
  ```
  The theme colours are applied by the pptx skill's `apply_theme.js` if it is available.
- Schematics: `python make_schematics.py` writes them to `assets/`. Project figures are taken from `../figures/`.
- Render for QA: `./render.sh` (headless LibreOffice extracted under `tools/`).
- Placeholders to fill before presenting: presenter, date, hackathon dates and participants, repository link.
