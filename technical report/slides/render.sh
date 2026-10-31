#!/usr/bin/env bash
# Render project_slides.pptx -> project_slides.pdf + render/s-NN.jpg for visual QA.
# Uses LibreOffice 26.8 extracted under tools/libreoffice (no root) + libs from conda_envs/office.
set -e
ROOT=/scratch/spascare/agentic-ai/workspaces/phase2-afdbc-evolution
P=$ROOT/tools/libreoffice/root/opt/libreoffice26.8/program
export LD_LIBRARY_PATH=$ROOT/conda_envs/office/lib
export FONTCONFIG_PATH=$ROOT/conda_envs/office/etc/fonts
SL="$ROOT/technical report/slides"
T=$(mktemp -d /tmp/render_XXXX)
cp "$SL/project_slides.pptx" "$T/"
(cd "$T" && $P/soffice --headless --norestore -env:UserInstallation=file:///tmp/lo_profile_y --convert-to pdf project_slides.pptx >/dev/null 2>&1)
cp "$T/project_slides.pdf" "$SL/"
mkdir -p "$SL/render" && rm -f "$SL/render/"*.jpg
$ROOT/conda_envs/main/bin/pdftoppm -jpeg -r 60 "$SL/project_slides.pdf" "$SL/render/s"
rm -rf "$T"; ls "$SL/render" | wc -l
