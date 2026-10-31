// Build the project slide deck: node build_deck.js  (run inside "technical report/slides")
// Images: ../figures/*.png (project figures) and assets/*.png (schematics, make_schematics.py)
const fs = require("fs");
const path = require("path");
const pptxgen = require("pptxgenjs");

const THEME = {
  name: "Structural clock",
  headFontFace: "Cambria",
  bodyFontFace: "Calibri",
  colors: {
    dk1: "1B262C", lt1: "FFFFFF", dk2: "0B3C49", lt2: "E8F1F2",
    accent1: "1B998B", accent2: "F4A259", accent3: "E76F51", accent4: "8D99AE",
    accent5: "3D6B78", accent6: "C9E4E7", hlink: "1B998B", folHlink: "3D6B78",
  },
};
const HEX = THEME.colors;

const pres = new pptxgen();
pres.layout = "LAYOUT_WIDE"; // 13.33 x 7.5 in
pres.title = "Is protein structure a molecular clock?";
pres.author = "AI-generated (Claude Opus 5.5 PI agent)";
pres.theme = { headFontFace: THEME.headFontFace, bodyFontFace: THEME.bodyFontFace };
const C = pres.SchemeColor;
const W = 13.333, H = 7.5;

// ---------- layouts ----------
pres.defineSlideMaster({
  title: "TITLE_DARK",
  background: { color: HEX.dk2 },
  objects: [
    { placeholder: { options: { name: "title", type: "title", x: 0.8, y: 2.0, w: 11.7, h: 1.9, fontSize: 44, bold: true, color: C.background1, valign: "bottom", align: "left" }, text: "" } },
    { placeholder: { options: { name: "body", type: "body", x: 0.8, y: 4.1, w: 11.7, h: 1.6, fontSize: 20, color: C.background2, valign: "top" }, text: "" } },
  ],
});
pres.defineSlideMaster({
  title: "SECTION",
  background: { color: HEX.dk2 },
  objects: [
    { placeholder: { options: { name: "title", type: "title", x: 2.3, y: 2.7, w: 10.2, h: 1.3, fontSize: 40, bold: true, color: C.background1, valign: "middle", align: "left" }, text: "" } },
    { placeholder: { options: { name: "body", type: "body", x: 2.3, y: 4.0, w: 10.2, h: 1.2, fontSize: 18, color: C.background2, valign: "top" }, text: "" } },
  ],
  slideNumber: { x: 12.4, y: 7.0, w: 0.6, h: 0.3, fontSize: 10, color: C.background2 },
});
pres.defineSlideMaster({
  title: "CONTENT",
  background: { color: HEX.lt1 },
  objects: [
    { placeholder: { options: { name: "title", type: "title", x: 0.5, y: 0.3, w: 12.3, h: 0.85, fontSize: 32, bold: true, color: C.text2, valign: "middle" }, text: "" } },
    { text: { text: "AI hackathon · Phase 2 · structural divergence", options: { x: 0.5, y: 7.05, w: 6, h: 0.3, fontSize: 10, color: C.accent4 } } },
  ],
  slideNumber: { x: 12.4, y: 7.05, w: 0.6, h: 0.3, fontSize: 10, color: C.accent4 },
});

// ---------- helpers ----------
function pngSize(p) {
  const b = fs.readFileSync(p);
  return { w: b.readUInt32BE(16), h: b.readUInt32BE(20) };
}
function img(slide, file, x, y, w, h, alt) {
  const p = path.resolve(file);
  if (!fs.existsSync(p)) {
    slide.addShape(pres.shapes.RECTANGLE, { x, y, w, h, fill: { color: C.background2 }, line: { color: C.accent4, dashType: "dash" } });
    slide.addText(`[image missing: ${path.basename(file)}]`, { x, y, w, h, align: "center", fontSize: 14, color: C.accent4, isTextBox: true });
    return;
  }
  const s = pngSize(p), ar = s.w / s.h;
  let iw = w, ih = w / ar;
  if (ih > h) { ih = h; iw = h * ar; }
  slide.addImage({ path: p, x: x + (w - iw) / 2, y: y + (h - ih) / 2, w: iw, h: ih, altText: alt || path.basename(file) });
}
let sectionName = "";
function content(title, notes) {
  const s = pres.addSlide({ masterName: "CONTENT", sectionTitle: sectionName });
  s.addText(title, { placeholder: "title" });
  if (notes) s.addNotes(notes);
  return s;
}
function section(num, title, sub, notes) {
  sectionName = title;
  pres.addSection({ title });
  const s = pres.addSlide({ masterName: "SECTION", sectionTitle: title });
  s.addShape(pres.shapes.OVAL, { x: 0.9, y: 2.75, w: 1.15, h: 1.15, fill: { color: C.accent2 }, line: { color: C.accent2 } });
  s.addText(String(num), { x: 0.9, y: 2.75, w: 1.15, h: 1.15, align: "center", valign: "middle", fontSize: 36, bold: true, color: C.text2, fontFace: THEME.headFontFace, isTextBox: true });
  s.addText(title, { placeholder: "title" });
  if (sub) s.addText(sub, { placeholder: "body" });
  if (notes) s.addNotes(notes);
  return s;
}
function badge(s, x, y, txt, fill, d = 0.5) {
  s.addShape(pres.shapes.OVAL, { x, y, w: d, h: d, fill: { color: fill || C.accent1 }, line: { color: fill || C.accent1 } });
  s.addText(txt, { x, y, w: d, h: d, align: "center", valign: "middle", fontSize: 16, bold: true, color: C.background1, isTextBox: true, margin: 0 });
}
function card(s, x, y, w, h, head, body, opt = {}) {
  s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y, w, h, rectRadius: 0.12, fill: { color: opt.fill || C.background2 }, line: { color: opt.fill || C.background2 } });
  s.addText(head, { x: x + 0.2, y: y + 0.12, w: w - 0.4, h: 0.45, fontSize: opt.headSize || 18, bold: true, color: opt.headColor || C.text2, isTextBox: true, margin: 0, valign: "top" });
  if (body) s.addText(body, { x: x + 0.2, y: y + 0.62, w: w - 0.4, h: h - 0.75, fontSize: opt.bodySize || 14, color: opt.bodyColor || C.text1, isTextBox: true, margin: 0, valign: "top", paraSpaceAfter: 4 });
}
function bullets(items) {
  return items.map((t, i) => ({ text: t, options: { bullet: true, breakLine: i < items.length - 1 } }));
}
function stat(s, x, y, w, big, label, color) {
  s.addText(big, { x, y, w, h: 0.9, fontSize: 44, bold: true, color: color || C.accent1, fontFace: THEME.headFontFace, isTextBox: true, margin: 0 });
  s.addText(label, { x, y: y + 0.9, w, h: 0.7, fontSize: 14, color: C.text1, isTextBox: true, margin: 0, valign: "top" });
}
function chip(s, x, y, w, txt, fill, color) {
  s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y, w, h: 0.42, rectRadius: 0.2, fill: { color: fill }, line: { color: fill } });
  s.addText(txt, { x, y, w, h: 0.42, align: "center", valign: "middle", fontSize: 13, bold: true, color: color || C.background1, isTextBox: true, margin: 0 });
}
const FIG = (f) => "../figures/" + f;
const AS = (f) => "assets/" + f;

// ================= TITLE =================
sectionName = "Opening";
pres.addSection({ title: "Opening" });
{
  const s = pres.addSlide({ masterName: "TITLE_DARK", sectionTitle: "Opening" });
  s.addText("Is protein structure a molecular clock?", { placeholder: "title" });
  s.addText([
    { text: "Revisiting structural divergence with AlphaFold, orthology and divergence times", options: { breakLine: true } },
    { text: "An AI-led Phase 2 project of the Beltrao-lab agentic AI hackathon · [presenter] · [date]", options: { fontSize: 16, color: C.accent2 } },
  ], { placeholder: "body" });
  s.addNotes("Title slide. Fill in presenter and date. The project was run end-to-end by Claude Opus 5.5 acting as PI, delegating to Opus/Sonnet subagents.");
}
{
  const s = content("Disclosure: every result here was produced by AI", "Human-written disclosure required by the hackathon. Keep on screen briefly.");
  s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: 0.9, y: 1.6, w: 11.5, h: 2.9, rectRadius: 0.15, fill: { color: C.background2 }, line: { color: C.background2 } });
  badge(s, 1.3, 2.0, "!", C.accent3, 0.8);
  s.addText([
    { text: "Exercise in assessing the bioinformatics capabilities of agentic AI (Beltrao-lab hackathon, github.com/evocellnet/ai-hackathon).", options: { breakLine: true } },
    { text: "Executed by Claude Opus 5.5 with Opus/Sonnet 5.5 subagents, within a 48 h compute cap (48 CPUs, 500 GB RAM, 1 A100).", options: { breakLine: true } },
    { text: "Results have not yet been verified by a human. This is not a traditional publication.", options: { bold: true } },
  ], { x: 2.4, y: 1.9, w: 9.6, h: 3.0, fontSize: 18, color: C.text1, isTextBox: true, paraSpaceAfter: 10, valign: "top" });
  s.addText("Repository: [link to project repository]  ·  Report: technical report/report.md", { x: 0.9, y: 4.8, w: 11.5, h: 0.5, fontSize: 14, color: C.accent4, isTextBox: true });
}

// ================= 1 HACKATHON =================
section(1, "The hackathon", "Can an AI agent run a bioinformatics project from prompt to paper?");
{
  const s = content("A week-long test of agentic AI in bioinformatics",
    "From the hackathon README: run by the Beltrao (ETH Zurich), Ryan (UCD), Goncalves (IST) and Dessimoz (UNIL) labs. Claude ran each project with minimal human intervention; humans verify afterwards. Blanks: dates, participants, compute details beyond this project.");
  const labs = [["Beltrao lab", "ETH Zürich"], ["Colm Ryan lab", "UCD"], ["Emanuel Gonçalves lab", "IST"], ["Christophe Dessimoz lab", "UNIL"]];
  labs.forEach((l, i) => card(s, 0.6, 1.45 + i * 1.3, 4.6, 1.1, l[0], l[1], { bodySize: 14 }));
  stat(s, 5.8, 1.5, 3.2, "2 phases", "Phase 1: 'predict X' tools\nPhase 2: open discovery");
  stat(s, 9.3, 1.5, 3.4, "15", "Phase-1 tools listed (pockets, docking, orthology, kinases, MS…)");
  stat(s, 5.8, 3.5, 3.2, "Claude", "Opus 5.5 as PI agent, delegating to Opus/Sonnet subagents", C.accent2);
  stat(s, 9.3, 3.5, 3.4, "[__]", "participants · [dates] · [Euler cluster details]", C.accent4);
  s.addText("Humans verify methodology and benchmarking afterwards; each repository opens with a human-written disclosure.", { x: 5.8, y: 5.6, w: 6.9, h: 0.9, fontSize: 14, italic: true, color: C.text2, isTextBox: true });
}
{
  const s = content("From prediction tools to discovery projects",
    "The attached project-flow diagram. Phase 1 built tools (blue); Phase 2 built models and datasets (green); Phase 3 discovery (orange). This project is a Phase 2 discovery project on structural evolution; it consumed the Phase-1 orthology tool Ancestra and public AFDB-cluster resources.");
  img(s, AS("hackathon_flow.png"), 0.5, 1.3, 8.6, 5.6, "Hackathon project flow diagram");
  card(s, 9.4, 1.4, 3.4, 2.4, "This project", "Phase 2 discovery: does protein 3D structure evolve like a clock?", { fill: C.accent2, headColor: C.text2 });
  card(s, 9.4, 4.0, 3.4, 2.6, "Built on", "Ancestra (Phase-1 orthology tool) + AlphaFold DB clusters (Foldseek) + TimeTree", {});
}

// ================= 2 PROMPT =================
section(2, "The prompt", "One generic research scaffold + one project-specific question");
{
  const s = content("The prompt has two layers",
    "Left: generic instructions shared by all hackathon projects. Right: what is specific to this project. The project ideas in the prompt were explicitly 'human priors, not a plan'.");
  s.addText("Shared by all projects", { x: 0.6, y: 1.35, w: 6, h: 0.5, fontSize: 20, bold: true, color: C.accent1, isTextBox: true, margin: 0 });
  s.addText("Specific to this project", { x: 6.9, y: 1.35, w: 6, h: 0.5, fontSize: 20, bold: true, color: C.accent2, isTextBox: true, margin: 0 });
  const gen = [["PI role", "plan, delegate, review, decide"], ["Design by critique", "literature → design → ≤3 review rounds"], ["Reproducibility", "seeds, DATA.md, git, no hard-coding"], ["Memory", "research tree + lessons file"], ["Verification", "toy data, controls, non-author review"], ["Report", "paper-style, disclaimer, figures"]];
  gen.forEach((g, i) => { const col = i % 2, row = Math.floor(i / 2); card(s, 0.6 + col * 3.05, 1.95 + row * 1.6, 2.9, 1.45, g[0], g[1], { bodySize: 13, headSize: 16 }); });
  const spec = [["Revisit Elofsson 2009", "linear, 3–10× more conserved?"], ["AFDB clusters 2023", "evolution never analysed"], ["Use Ancestra", "orthologs vs paralogs"], ["Neutral clock?", "comment on structure"], ["Ideas (priors)", "LDO, kinases, cytoskeleton"], ["Budget", "48 h, node move at 12 h"]];
  spec.forEach((g, i) => { const col = i % 2, row = Math.floor(i / 2); card(s, 6.9 + col * 3.05, 1.95 + row * 1.6, 2.9, 1.45, g[0], g[1], { bodySize: 13, headSize: 16, fill: "FDEBD9" }); });
}
{
  const s = content("Two goals, five questions",
    "Goal 1 is the confirmatory core (pre-registered endpoints E1–E5). Goal 2 is open-ended; the duplication test E6 doubles as a positive control.");
  card(s, 0.6, 1.4, 6.0, 5.2, "Goal 1 · Is structure a molecular clock?", null, {});
  const g1 = ["Shape: linear or saturating vs sequence?", "Ratio: how much slower than sequence?", "Time: same ratio on young and old branches?", "Lineages: structure-specific rates?"];
  g1.forEach((t, i) => { badge(s, 0.9, 2.15 + i * 1.05, String(i + 1)); s.addText(t, { x: 1.6, y: 2.12 + i * 1.05, w: 4.8, h: 0.6, fontSize: 16, color: C.text1, isTextBox: true, margin: 0, valign: "middle" }); });
  card(s, 6.9, 1.4, 5.9, 5.2, "Goal 2 · AFDB clusters in the light of evolution", null, { fill: "FDEBD9" });
  const g2 = ["Duplication: do paralogs diverge more in structure?", "Do structural clusters match ortholog groups?", "Case studies: cytoskeleton & kinases"];
  g2.forEach((t, i) => { badge(s, 7.2, 2.15 + i * 1.05, String(i + 5), C.accent2); s.addText(t, { x: 7.9, y: 2.12 + i * 1.05, w: 4.7, h: 0.6, fontSize: 16, color: C.text1, isTextBox: true, margin: 0, valign: "middle" }); });
}
{
  const s = content("24 hours of wall clock, two compute nodes",
    "Node 1 (10.9 h): scouting, design, panel, Ancestra, pair metrics, tuning part 1. Node 2 (~12 h): tuning part 2, dry run, single confirmation run, verification, post-hoc review, report. Total compute ~23 h of the 48 h cap.");
  // timeline axis
  s.addShape(pres.shapes.LINE, { x: 0.7, y: 3.6, w: 11.9, h: 0, line: { color: C.accent4, width: 2 } });
  s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: 0.7, y: 3.35, w: 5.4, h: 0.5, rectRadius: 0.1, fill: { color: C.accent1 }, line: { color: C.accent1 } });
  s.addText("Node 1 · 10.9 h", { x: 0.7, y: 3.35, w: 5.4, h: 0.5, align: "center", valign: "middle", fontSize: 14, bold: true, color: C.background1, isTextBox: true, margin: 0 });
  s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: 6.5, y: 3.35, w: 6.1, h: 0.5, rectRadius: 0.1, fill: { color: C.text2 }, line: { color: C.text2 } });
  s.addText("Node 2 · ~12 h", { x: 6.5, y: 3.35, w: 6.1, h: 0.5, align: "center", valign: "middle", fontSize: 14, bold: true, color: C.background1, isTextBox: true, margin: 0 });
  const ms = [
    [0.7, "Scouting", "literature, data, Ancestra, env", true], [1.9, "Design", "3 review rounds", false],
    [3.0, "Ancestra", "180 proteomes, 5 h", true], [4.3, "Pairs", "1.4 M pairs × 3 alignments", false],
    [6.5, "Tuning", "20% of families", true], [7.9, "Freeze → confirm", "once, 80%", false],
    [9.4, "Verify", "re-derive 42/42 · refute", true], [11.0, "Report", "+ review", false]];
  ms.forEach(m => {
    const up = m[3]; const y = up ? 1.5 : 4.2;
    s.addShape(pres.shapes.LINE, { x: m[0] + 0.15, y: up ? 2.85 : 3.85, w: 0, h: 0.5, line: { color: C.accent4, width: 1 } });
    s.addText([{ text: m[1], options: { bold: true, fontSize: 15, color: C.text2, breakLine: true } }, { text: m[2], options: { fontSize: 12, color: C.text1 } }],
      { x: m[0], y, w: 1.9, h: 1.3, isTextBox: true, margin: 0, valign: up ? "bottom" : "top" });
  });
  s.addText("Node move at 20:30 on day 1 (checkpointed); resumed from the research tree on day 2.", { x: 0.7, y: 6.0, w: 11.9, h: 0.5, fontSize: 14, italic: true, color: C.accent4, isTextBox: true });
}
{
  const s = content("An AI research team: one PI agent, many specialists",
    "Roughly 40 subagent tasks across scouting, building, reviewing, verifying, analysing and writing. The PI read written summaries only, keeping its own context lean.");
  img(s, AS("s_agents.png"), 0.5, 1.3, 8.2, 5.6, "Organisation of the agent team");
  stat(s, 9.2, 1.5, 3.6, "≈40", "subagent tasks with written reports");
  stat(s, 9.2, 3.2, 3.6, "13", "independent reviews (design, results, report, retrospective)", C.accent2);
  stat(s, 9.2, 4.9, 3.6, "85+", "git commits: every decision traceable", C.text2);
}

// ================= 3 DECISIONS =================
section(3, "How the method was decided", "Scout → design → independent critique → pilot → freeze");
{
  const s = content("Scouting first: four agents in parallel",
    "Design was written only after four scouting reports. Each changed the design: time had never been the axis; in-cluster pairs are biased; TimeTree is eukaryote-centric; Ancestra orthologs are reliable but duplication timing is not.");
  const sc = [["Literature", "Elofsson 2009 used substitutions, not time → a time test is new", "1"], ["Data", "Cluster pairs are pre-filtered by structure → biased toward plateau", "2"], ["Orthology", "Ancestra: 148/149 Ensembl orthologs; duplication timing only 36% right", "3"], ["Compute", "Foldseek GPU works; dated trees reliable for eukaryotes only", "4"]];
  sc.forEach((c, i) => { card(s, 0.6 + i * 3.1, 1.5, 2.9, 2.6, c[0], c[1], { bodySize: 16 }); badge(s, 0.6 + i * 3.1 + 2.25, 1.6, c[2], C.accent2); });
  s.addText("→ Design choices: eukaryote panel with dated tree · orthologs chosen by sequence only · never trust duplication ages · measure the AlphaFold bias", { x: 0.6, y: 4.5, w: 12.1, h: 1.0, fontSize: 18, bold: true, color: C.text2, isTextBox: true });
}
{
  const s = content("Design by independent critique: 5 → 3 → 1 reviewers",
    "Round 1: external PI, molecular evolution, statistics, algorithms/compute, AI failure modes. Round 2: external PI, statistics, molecular evolution. Round 3: consistency check that produced the executable ANALYSIS_PLAN.");
  const r = [[5, "Round 1"], [3, "Round 2"], [1, "Round 3"]];
  r.forEach((x, i) => {
    for (let k = 0; k < x[0]; k++) s.addShape(pres.shapes.OVAL, { x: 0.8 + i * 1.6 + (k % 3) * 0.42, y: 1.6 + Math.floor(k / 3) * 0.42, w: 0.36, h: 0.36, fill: { color: [C.accent1, C.accent5, C.text2][i] }, line: { color: C.background1 } });
    s.addText(x[1], { x: 0.7 + i * 1.6, y: 2.55, w: 1.5, h: 0.4, fontSize: 14, bold: true, color: C.text2, isTextBox: true, margin: 0 });
  });
  const ch = ["Branch lengths on the dated species tree replace pair medians", "Pre-set decision rules (equivalence margins) replace AIC/BIC", "Positive control (duplication) + power checks", "Tune on 20%, confirm once on 80%; Holm correction", "Noise floors; errors-in-variables correction"];
  ch.forEach((t, i) => { badge(s, 5.8, 1.5 + i * 0.98, "✓", C.accent1, 0.45); s.addText(t, { x: 6.45, y: 1.45 + i * 0.98, w: 6.4, h: 0.6, fontSize: 16, color: C.text1, isTextBox: true, margin: 0, valign: "middle" }); });
  s.addText("Every reviewer had to propose concrete changes within scope, or say none were needed.", { x: 0.6, y: 3.4, w: 4.8, h: 1.5, fontSize: 15, italic: true, color: C.accent4, isTextBox: true });
}
{
  const s = content("Pilots on real data caught what reviewers missed",
    "Running the full pipeline on a 20-species test set, before any endpoint, showed that TM-align (the planned primary alignment) mis-aligns residues of close relatives and gives non-zero distance to identical sequences. Decision D22: switch to Foldseek 3Di+AA alignment; keep TM-align as sensitivity.");
  s.addChart(pres.charts.BAR, [{ name: "median core AA p-distance, human–mouse orthologs", labels: ["TM-align", "Sequence-only", "3Di+AA"], values: [0.091, 0.064, 0.064] }], {
    x: 0.6, y: 1.4, w: 6.6, h: 5.2, barDir: "col", chartColors: [HEX.accent3, HEX.accent1, HEX.accent1],
    showTitle: true, title: "Same orthologs, three alignments", titleFontSize: 16, titleColor: HEX.dk2, titleFontFace: "+mj-lt",
    showValue: true, dataLabelPosition: "outEnd", dataLabelFontSize: 14, dataLabelColor: HEX.dk1, dataLabelFormatCode: "0.000", dataLabelFontFace: "+mn-lt",
    catAxisLabelFontSize: 14, valAxisLabelFontSize: 12, catAxisLabelColor: HEX.dk1, valAxisLabelColor: HEX.accent4, catAxisLabelFontFace: "+mn-lt", valAxisLabelFontFace: "+mn-lt",
    valGridLine: { color: "E5E5E5", size: 0.5 }, catGridLine: { style: "none" }, showLegend: false, valAxisMinVal: 0, valAxisMaxVal: 0.1,
  });
  card(s, 7.6, 1.5, 5.1, 1.6, "8 / 104", "identical-sequence pairs got non-zero AA distance under TM-align", { headSize: 28, fill: "FBE3DC" });
  card(s, 7.6, 3.3, 5.1, 1.6, "Decision D22", "primary alignment → Foldseek 3Di+AA; TM-align and sequence-only kept as brackets", { headSize: 20 });
  card(s, 7.6, 5.1, 5.1, 1.4, "Lesson", "no reviewer predicted it; a small real pilot did", { headSize: 20, fill: "FDEBD9" });
}
{
  const s = content("Pre-registration: tune on 20%, confirm once on 80%",
    "The confirmation script refuses to run unless plan, tuning file and family split match their checksums in the results ledger, the git tree is clean and an explicit flag is set. Every change to the plan is logged with reason and effect (34 entries).");
  img(s, AS("s_prereg.png"), 0.5, 1.3, 8.4, 5.6, "Pre-registration workflow");
  stat(s, 9.3, 1.5, 3.5, "10/10", "automated guards passed before the single confirmation run");
  stat(s, 9.3, 3.2, 3.5, "34", "logged deviations (32 before confirmation)", C.accent2);
  stat(s, 9.3, 4.9, 3.5, "42/42", "numbers re-derived by a non-author agent", C.text2);
}
{
  const s = content("Key decisions and why", "Each row links to DEVIATIONS.md / RESEARCH_TREE.md.");
  const rows = [
    ["Use only proteins whose AlphaFold model sequence = UniProt sequence", "sequence and structure must describe the same molecule"],
    ["Never use Ancestra's duplication ages", "only 36% agreement with Ensembl"],
    ["One Ancestra run on 180 whole proteomes", "per-cluster runs would cost thousands of CPU-hours"],
    ["3Di+AA alignment as primary (D22)", "TM-align failed the identical-sequence check"],
    ["Fix two estimator biases found in simulation (D31, D32)", "both would have produced false positives"],
    ["Keep pre-set margins even when power was low", "never change tests to make results look better"],
    ["Report post-hoc reinterpretation separately, after independent refutation attempt", "a fixed noise floor turned out incomplete"],
  ];
  const tbl = [[{ text: "Decision", options: { bold: true, color: HEX.lt1, fill: { color: HEX.dk2 } } }, { text: "Why", options: { bold: true, color: HEX.lt1, fill: { color: HEX.dk2 } } }]]
    .concat(rows.map((r, i) => [{ text: r[0], options: { fill: { color: i % 2 ? HEX.lt1 : HEX.lt2 } } }, { text: r[1], options: { fill: { color: i % 2 ? HEX.lt1 : HEX.lt2 }, color: HEX.accent5 } }]));
  s.addTable(tbl, { x: 0.6, y: 1.4, w: 12.1, colW: [7.0, 5.1], fontSize: 15, fontFace: "Calibri", border: { type: "solid", color: "FFFFFF", pt: 1 }, rowH: 0.62, valign: "middle" });
}

// ================= 4 METHOD =================
section(4, "The method, step by step", "With plain-language explainers for the tricky parts");
{
  const s = content("Pipeline at a glance", "Species panel → AlphaFold structures → orthology (Ancestra) → family sampling and split → three alignments per pair → amino-acid, 3Di and 3D distances → noise floors → two estimation engines → six pre-registered endpoints.");
  img(s, FIG("fig1_workflow.png"), 0.5, 1.3, 12.3, 5.6, "Workflow");
}
{
  const s = content("Data: 180 eukaryotes, 3.3 M predicted structures", "Panel selected for AlphaFold v4 coverage (≥80% identical sequences), BUSCO completeness, TimeTree dating, and an even spread of split ages; 100 prokaryotes for case studies only.");
  img(s, FIG("fig2_panel.png"), 0.4, 1.3, 8.4, 5.6, "Species panel overview");
  stat(s, 9.2, 1.4, 3.6, "180", "eukaryotic proteomes (66 animals, 40 fungi, 20 plants, 54 protists…)");
  stat(s, 9.2, 3.0, 3.6, "109,598", "ortholog families (HOGs) from Ancestra", C.accent2);
  stat(s, 9.2, 4.6, 3.6, "1.4 M", "protein pairs from 7,000 sampled families", C.text2);
}
{
  const s = content("Explainer: orthologs, paralogs and HOGs", "Orthologs: separated by speciation (same gene in two species). Paralogs: separated by duplication. A HOG gathers all descendants of one ancestral gene. Only orthologs have a split time we can read from the species tree.");
  img(s, AS("s_orthology.png"), 0.5, 1.3, 8.2, 5.6, "Orthologs, paralogs and HOGs");
  card(s, 9.0, 1.5, 3.8, 1.6, "Ortholog", "split by speciation → age known from the species tree", {});
  card(s, 9.0, 3.3, 3.8, 1.6, "Paralog", "split by duplication → age unknown", { fill: "FDEBD9" });
  card(s, 9.0, 5.1, 3.8, 1.5, "Why it matters", "orthologs give a clean clock test", {});
}
{
  const s = content("Explainer: AlphaFold models and 3Di letters", "Foldseek turns each residue's local 3D contact geometry into one of 20 letters (3Di). Proteins become two strings of equal alphabet size, so structure and sequence change can be counted with the same machinery.");
  img(s, AS("s_3di.png"), 0.5, 1.3, 8.2, 5.6, "Amino acids and 3Di letters");
  card(s, 9.0, 1.5, 3.8, 1.7, "20 vs 20", "amino acids vs structural letters: a like-for-like comparison", {});
  card(s, 9.0, 3.4, 3.8, 1.6, "pLDDT ≥ 70", "only confidently predicted residues used", { fill: "FDEBD9" });
  card(s, 9.0, 5.2, 3.8, 1.4, "Illergård 2009", "had to approximate this by grouping amino acids", {});
}
{
  const s = content("Explainer: counting change when changes overlap", "Raw percent difference saturates because sites change more than once. A maximum-likelihood (ML) distance with a substitution model (LG for amino acids, Q.3Di.AF for 3Di) estimates substitutions per site; beyond ~2 it becomes unreliable, so families were restricted to ≤2.");
  img(s, AS("s_saturation.png"), 0.5, 1.3, 8.2, 5.6, "Saturation and ML correction");
  card(s, 9.0, 1.5, 3.8, 2.2, "ML distance", "expected substitutions per site, corrected for multiple hits", {});
  card(s, 9.0, 3.9, 3.8, 2.4, "Validated", "unbiased within 1.6% on simulated data; matches IQ-TREE within 0.004%", { fill: "FDEBD9" });
}
{
  const s = content("Explainer: three ways to compare two structures", "TM-score: global fold match. LDDT: fraction of local inter-residue distances preserved. 3Di distance: substitutions in the structural alphabet. All are measured only over confidently predicted residues.");
  img(s, AS("s_structcompare.png"), 0.5, 1.3, 8.2, 5.6, "TM-score, LDDT, pLDDT");
  stat(s, 9.2, 1.5, 3.6, "0.95", "median TM-score between AlphaFold models of identical sequences: predictions are not identical");
  stat(s, 9.2, 3.6, 3.6, "0.052", "3Di distance at zero sequence change (prediction noise floor)", C.accent3);
}
{
  const s = content("Explainer: two estimation engines", "Pairwise engine: fit structure distance = floor + a·(sequence distance)^c; c = 1 means proportional, c < 1 saturating. Branch engine: place each family's distances on the dated species tree, giving an AA length and a 3Di length for every branch; ratios of sums give k.");
  img(s, AS("s_engines.png"), 0.4, 1.3, 9.0, 5.6, "Pairwise and branch engines");
  card(s, 9.7, 1.5, 3.1, 2.4, "Shape", "curvature c from pairs (E1, E2)", {});
  card(s, 9.7, 4.1, 3.1, 2.5, "Ratio · time · lineage", "branch lengths on the species tree (E3–E5)", { fill: "FDEBD9" });
}
{
  const s = content("Explainer: molecular clock and the neutral null", "Clock: change accumulates in proportion to time. Neutral structural null (N0): structure changes only as a passive by-product of sequence change, one constant ratio per family across lineages and ages. A good fit would show 'proportional', not prove 'neutral'.");
  img(s, AS("s_clock.png"), 0.5, 1.3, 8.2, 5.6, "Molecular clock schematic");
  card(s, 9.0, 1.5, 3.8, 2.2, "Proportional", "constant 3Di/AA ratio across ages and lineages", {});
  card(s, 9.0, 3.9, 3.8, 2.4, "Violation", "ratio drifts with age, or some lineages change structure faster", { fill: "FBE3DC" });
}
{
  const s = content("How this differs from previous approaches", "Commonalities: all relate structural to sequence divergence. Differences: this project selects pairs by orthology only, uses a real time axis, matched 20-state alphabets, and pre-registration at AFDB scale.");
  const H = (t) => ({ text: t, options: { bold: true, color: HEX.lt1, fill: { color: HEX.dk2 } } });
  const rows = [
    ["Illergård et al. 2009", "SCOP domain pairs (PDB)", "substitutions", "core states / RMSD", "~hundreds of families"],
    ["Pascual-García et al. 2019", "superfamilies (PDB)", "substitutions", "structural distances", "clock violations vs function"],
    ["Barrio-Hernandez et al. 2023", "AFDB, 214 M models", "none (cluster LCA)", "Foldseek clustering", "no rate analysis"],
    ["Moi et al. 2025 (FoldTree)", "families (AFDB)", "tree topology", "3Di + AA distances", "phylogenies"],
    ["Warwick Vesztrocy et al. 2025", "paralogs (AFDB)", "branch lengths", "Foldseek LDDT", "least-diverged ortholog"],
    ["This project", "orthologs, 180 eukaryotes", "TimeTree ages + ML", "3Di (Q.3Di.AF), LDDT, TM", "7,000 families, pre-registered"],
  ];
  const tbl = [[H("Study"), H("Pairs"), H("Divergence axis"), H("Structure measure"), H("Scale / focus")]].concat(rows.map((r, i) => r.map(c => ({ text: c, options: { fill: { color: i === rows.length - 1 ? "FDEBD9" : (i % 2 ? HEX.lt1 : HEX.lt2) }, bold: i === rows.length - 1 } }))));
  s.addTable(tbl, { x: 0.6, y: 1.4, w: 12.1, colW: [2.9, 2.5, 2.2, 2.4, 2.1], fontSize: 14, fontFace: "Calibri", border: { type: "solid", color: "FFFFFF", pt: 1 }, rowH: 0.68, valign: "middle" });
}

// ================= 5 HYPOTHESES =================
section(5, "The hypotheses, one by one", "Each with: what it is, prior literature, pros and cons of our test, and the verdict");
function hypo(title, image, what, lit, pros, cons, verdictPre, verdictPost, notes) {
  const s = content(title, notes);
  img(s, image, 0.4, 1.3, 6.3, 4.6, title);
  s.addText(what, { x: 0.5, y: 6.0, w: 6.2, h: 0.9, fontSize: 14, italic: true, color: C.text2, isTextBox: true, margin: 0, valign: "top" });
  s.addText("Prior work", { x: 7.0, y: 1.3, w: 5.8, h: 0.35, fontSize: 15, bold: true, color: C.accent1, isTextBox: true, margin: 0 });
  s.addText(lit, { x: 7.0, y: 1.65, w: 5.8, h: 0.85, fontSize: 13, color: C.text1, isTextBox: true, margin: 0, valign: "top" });
  card(s, 7.0, 2.55, 2.8, 2.05, "+ Pros", null, { headSize: 15 });
  s.addText(bullets(pros), { x: 7.15, y: 2.95, w: 2.55, h: 1.6, fontSize: 14, color: C.text1, isTextBox: true, margin: 0, valign: "top", paraSpaceAfter: 3 });
  card(s, 10.0, 2.55, 2.8, 2.05, "− Cons", null, { headSize: 15, fill: "FBE3DC" });
  s.addText(bullets(cons), { x: 10.15, y: 2.95, w: 2.55, h: 1.6, fontSize: 14, color: C.text1, isTextBox: true, margin: 0, valign: "top", paraSpaceAfter: 3 });
  if (verdictPost) {
    chip(s, 7.0, 4.8, 2.3, "Pre-registered", HEX.dk2);
    s.addText(verdictPre, { x: 9.4, y: 4.75, w: 3.4, h: 0.9, fontSize: 14, color: C.text1, isTextBox: true, margin: 0, valign: "top" });
    chip(s, 7.0, 5.8, 2.3, "Post hoc (verified)", HEX.accent2, HEX.dk2);
    s.addText(verdictPost, { x: 9.4, y: 5.75, w: 3.4, h: 1.1, fontSize: 14, bold: true, color: C.text2, isTextBox: true, margin: 0, valign: "top" });
  } else {
    chip(s, 7.0, 4.8, 2.3, "Exploratory", HEX.accent2, HEX.dk2);
    s.addText(verdictPre.replace(/^Exploratory: /, ""), { x: 9.4, y: 4.75, w: 3.4, h: 1.9, fontSize: 15, bold: true, color: C.text2, isTextBox: true, margin: 0, valign: "top" });
  }
  return s;
}
hypo("H1 · Shape: does structure diverge linearly with sequence?", FIG("figA_shape.png"),
  "3Di (and LDDT) distance vs amino-acid distance across ortholog pairs; curvature c = 1 is linear.",
  "Chothia & Lesk 1986 (exponential vs identity); Koehl & Levitt 2002; Illergård et al. 2009 (linear vs ML distance).",
  ["orthology-selected pairs", "matched 20-state models", "errors-in-variables correction", "3 alignments"],
  ["per-model AlphaFold noise biases toward concave", "mild curvature hard to detect"],
  "Concave: ĉ = 0.72 [0.64–0.80], robust to alignment",
  "Curvature bracketed 0.75–1.0; moves toward 1 in high-confidence models",
  "E1/E2. The pre-registered estimator assumed a fixed noise floor. A per-protein offset (~0.01–0.017) biases it toward concavity (synthetic truth c=1 → 0.75–0.78).");
hypo("H2 · Is structure 3–10× more conserved than sequence?", FIG("figD_ratio_k.png"),
  "k = total amino-acid change / total 3Di change over tree branches.",
  "Illergård et al. 2009: 3–10× after alphabet matching (15–30× with raw states); β faster than α.",
  ["true 20 vs 20 comparison", "branch-level sums on the species tree", "4,891 families"],
  ["specific to 3Di and predicted structures", "terminal branches inflated by noise"],
  "k = 3.3 [3.09–3.57]; α 4.7× vs β 2.6×",
  "Internal branches: k = 4.6; α 6.1× vs β 3.9× (gap persists)",
  "E3. The plain-language reading: for every 4–5 amino-acid substitutions per site, about one residue's local 3D-contact letter changes.");
hypo("H3 · Time: is the structure/sequence ratio constant?", FIG("figE_terminal_internal.png"),
  "Ratio k by branch age (TimeTree), split into terminal (species tips) and internal branches.",
  "Zuckerkandl & Pauling 1965; Kimura 1968 (neutral clock); Pascual-García 2019 vs Moi 2025 (mixed evidence).",
  ["dated species tree", "alternative deep timescales", "terminal vs internal split"],
  ["only 224 young families", "deep dates uncertain", "tip-correction sensitive"],
  "Ratio rises with age (p_Holm 3×10⁻¹²), robust to α and timescale",
  "Age trend = terminal-branch noise artefact; internal ratio flat ≈4.6 from 30 My to >1 Gy",
  "E4. The slope tracks the size of the tip correction almost linearly and vanishes near 0.036; internal branches give −0.047, n.s.");
hypo("H4 · Lineages: is there a neutral structural clock?", AS("s_clock.png"),
  "Do some lineages change structure faster than their sequence change predicts?",
  "Ohta 1973 (nearly neutral); Bromham & Penny 2003 (rate variation); Pascual-García 2019 (structural clock violations).",
  ["sister-lineage contrasts", "spike-in power test", "guards against artefacts"],
  ["power only 24%", "accepting the null practically unreachable"],
  "Artefact-suspect: excess vanishes on internal branches",
  "No robust lineage-specific structural rates; neutrality not testable with this design",
  "E5. Margins were not widened after seeing the low power.");
hypo("H5 · Do AlphaFold models make orthologs look too similar?", AS("s_af_bias.png"),
  "Same ortholog pairs and same residues: AlphaFold vs experimental (PDB) structures.",
  "Jumper et al. 2021; Buel & Walters 2022 and Pak et al. 2023 (AF insensitive to point mutations).",
  ["identical residue sets", "post-training-cutoff structures", "197 post-cutoff pairs"],
  ["few deep pairs", "crystal conditions add own noise"],
  "Exploratory: AF–AF divergence = 0.59× PDB–PDB for 3Di",
  null,
  "X4. Bias is a roughly constant offset, not growing with divergence. A model vs its own crystal structure differs more (0.24) than two ortholog models (0.17).");
hypo("H6 · Do structural clusters hide divergent pairs?", AS("s_selection_bias.png"),
  "Pairs inside AFDB clusters passed a structural filter; ortholog pairs did not.",
  "Barrio-Hernandez et al. 2023 (clustering: 90% coverage, E < 0.01).",
  ["27.7 M cluster pairs vs our orthologs", "same metrics"],
  ["ranges barely overlap (10–35% vs 47–99% identity)", "scoring differs"],
  "Exploratory: in-cluster TM plateaus ≈0.83 = prediction-quality ceiling",
  null,
  "X7 / P-clu. The plateau tracks representative pLDDT and length rather than evolutionary saturation.");
hypo("H7 · Duplication: do paralogs diverge more in structure?", AS("s_ldo.png"),
  "Positive control and least-diverged-ortholog (LDO) test with split alignment columns.",
  "Ortholog conjecture: Nehrt et al. 2011, Altenhoff et al. 2012; LDO: Warwick Vesztrocy et al. 2025.",
  ["design avoids regression to the mean", "3,525 families"],
  ["paralog part not alignment-robust", "effect ~ size of per-model noise"],
  "Detected: paralog excess 0.030; faster copy +0.017",
  "Faster-copy excess halves after pLDDT adjustment (0.011) but stays positive; direction replicates",
  "E6. Direction of the published result replicates; biological size uncertain.");
hypo("H8 · Are AFDB structural clusters ortholog groups?", AS("s_cluster_vs_hog.png"),
  "Map 109,598 ortholog groups onto AFDB structural clusters.",
  "Barrio-Hernandez et al. 2023 (cluster ages by last common ancestor); TED domains (Lau et al. 2024).",
  ["all families", "logistic model of splits"],
  ["Ancestra orthology itself unverified", "domain labels noisy (21–25%)"],
  "Exploratory: 70% of groups split across clusters",
  null,
  "X5. Length ratio dominates (log-odds 4.06 per SD); domain architecture adds nothing once length is in. Cluster ages agree exactly with orthology ages in only 27%.");
hypo("H9 · Case studies: bursts of structural change?", FIG("figI_cytoskeleton_methods.png"),
  "Actin/tubulin folds (with prokaryotic & Asgard homologs) and protein kinases.",
  "Prokaryotic cytoskeleton homologs (MreB, FtsZ), Asgard archaea; kinase evolution (Kannan & Taylor).",
  ["fixed reference trees", "3 alignments", "two estimators"],
  ["per-branch ratios disagree between methods", "deep sequence saturation"],
  "Exploratory: no robust stem bursts",
  null,
  "X6 / X8. Tyrosine-kinase origin shows a weak, consistent excess (88–91st percentile); not a burst.");

// ================= 6 RESULTS =================
section(6, "Results and take-home messages", "What the pre-registered tests said, and what survived scrutiny");
{
  const s = content("Scorecard: pre-registered vs post-hoc reading", "Pre-registered decisions are reported unchanged. Post-hoc readings were checked by an independent agent that tried to refute them.");
  const H = (t) => ({ text: t, options: { bold: true, color: HEX.lt1, fill: { color: HEX.dk2 } } });
  const rows = [
    ["Shape (E1/E2)", "concave, ĉ 0.72", "0.75–1.0: near-proportional"],
    ["Ratio (E3)", "3.3× (α 4.7, β 2.6)", "4.6× on internal branches"],
    ["Time (E4)", "ratio rises with age", "flat 30 My – >1 Gy (tip artefact)"],
    ["Lineages (E5)", "artefact-suspect", "no robust lineage effects"],
    ["Duplication (E6)", "detected", "direction robust; size ~halved"],
    ["AlphaFold bias (X4)", "—", "orthologs look 30–40% more similar"],
    ["Clusters vs orthology (X5)", "—", "70% of groups split, by length"],
  ];
  const tbl = [[H("Question"), H("Pre-registered (confirm set)"), H("Best reading")]].concat(rows.map((r, i) => r.map((c, j) => ({ text: c, options: { fill: { color: j === 2 ? "FDEBD9" : (i % 2 ? HEX.lt1 : HEX.lt2) }, bold: j === 2 } }))));
  s.addTable(tbl, { x: 0.6, y: 1.4, w: 12.1, colW: [3.6, 4.0, 4.5], fontSize: 16, fontFace: "Calibri", border: { type: "solid", color: "FFFFFF", pt: 1 }, rowH: 0.63, valign: "middle" });
}
{
  const s = content("The twist: one hidden noise source explains most signals", "Each AlphaFold model carries a small protein-specific structural offset (~0.01–0.017 3Di substitutions per site, larger at low confidence). Pre-registered floors from identical sequences missed it. It sits on terminal branches: young branches look structurally 'fast', which bends curves and creates an age trend.");
  img(s, FIG("figT_tip_artefact.png"), 0.5, 1.3, 12.3, 4.4, "Tip artefact schematic");
  chip(s, 0.7, 6.0, 3.8, "Bends E1 toward concave", HEX.accent3);
  chip(s, 4.75, 6.0, 3.8, "Creates E4 age trend", HEX.accent3);
  chip(s, 8.8, 6.0, 3.8, "Fakes E5 lineage effects", HEX.accent3);
}
{
  const s = content("What is new, and how it relates to the literature", "Green: confirmed prior claims. Amber: refined. Coral: new observations from this project.");
  card(s, 0.6, 1.4, 3.9, 4.4, "Confirmed", null, { fill: "DDF1EE" });
  s.addText(bullets(["Structure changes more slowly than sequence (Illergård 2009): 3–5×", "β-rich folds change structure faster than α (Illergård 2009)", "Direction of the LDO structural result (Warwick Vesztrocy 2025)"]), { x: 0.8, y: 2.0, w: 3.5, h: 3.7, fontSize: 17, color: C.text1, isTextBox: true, margin: 0, valign: "top", paraSpaceAfter: 10 });
  card(s, 4.7, 1.4, 3.9, 4.4, "Refined", null, { fill: "FDEBD9" });
  s.addText(bullets(["'Linear' holds roughly, now on matched 20-state models and thousands of families", "Structure is proportional to sequence on a real time axis (30 My – >1 Gy, internal branches)", "Structural clock violations not robust once artefacts handled"]), { x: 4.9, y: 2.0, w: 3.5, h: 3.7, fontSize: 17, color: C.text1, isTextBox: true, margin: 0, valign: "top", paraSpaceAfter: 10 });
  card(s, 8.8, 1.4, 3.9, 4.4, "New", null, { fill: "FBE3DC" });
  s.addText(bullets(["AlphaFold homogenises orthologs (30–40% more similar than crystals)", "Per-model 3Di noise inflates terminal branches", "AFDB clusters split ortholog groups, mostly by length"]), { x: 9.0, y: 2.0, w: 3.5, h: 3.7, fontSize: 17, color: C.text1, isTextBox: true, margin: 0, valign: "top", paraSpaceAfter: 10 });
}
{
  const s = content("Take-home messages", "Four numbers to remember.");
  const t = [["≈4.6×", "sequence changes faster than 3Di structure (internal branches)", C.accent1], ["0.75–1.0", "curvature: structure roughly proportional to sequence, mild saturation at most", C.text2], ["30–40%", "AlphaFold makes orthologs look more alike than experiments", C.accent3], ["70%", "of ortholog groups are split across AFDB structural clusters", C.accent2]];
  t.forEach((x, i) => { const col = i % 2, row = Math.floor(i / 2); card(s, 0.6 + col * 6.2, 1.4 + row * 2.7, 5.9, 2.5, "", null, {}); stat(s, 0.9 + col * 6.2, 1.6 + row * 2.7, 5.4, x[0], x[1], x[2]); });
}

// ================= 7 DECISION MAKING & FUTURE =================
section(7, "Decision making and the future", "What transfers to the next AI-led project");
{
  const s = content("The decision loop, and where it changed course", "Amber callouts mark where a step changed the plan: pilot (D22), synthetic tests (D31/D32), dry run (missing sensitivities), and the post-hoc refutation (per-model noise).");
  img(s, AS("s_workflow_decisions.png"), 0.4, 1.3, 12.5, 5.6, "Decision workflow");
}
{
  const s = content("What worked, and what did not", "From the milestone retrospectives and LESSONS.md.");
  card(s, 0.6, 1.4, 5.9, 4.6, "Worked", null, { fill: "DDF1EE" });
  s.addText(bullets(["Parallel scouting before design", "Small real pilot before freezing", "Synthetic worlds with known answers", "Pre-registration + deviation log", "Independent re-derivation and refutation attempts", "Research tree as memory across a node move"]), { x: 0.85, y: 2.0, w: 5.4, h: 3.9, fontSize: 18, color: C.text1, isTextBox: true, margin: 0, valign: "top", paraSpaceAfter: 10 });
  card(s, 6.8, 1.4, 5.9, 4.6, "Did not", null, { fill: "FBE3DC" });
  s.addText(bullets(["Noise floor from a single control type", "Plan-specified sensitivities unwired until the dry run", "Exploratory work competing with the critical path", "Reviewers shared blind spots (same model family)", "Harness quirks: blocked file writes, self-matching wait loops"]), { x: 7.05, y: 2.0, w: 5.4, h: 3.9, fontSize: 18, color: C.text1, isTextBox: true, margin: 0, valign: "top", paraSpaceAfter: 10 });
}
{
  const s = content("Carrying the decision process to future projects", "Rules distilled from this project's LESSONS.md.");
  const r = [["1", "Scout data, tools and literature in parallel before designing"], ["2", "Pilot the whole pipeline on small real data before freezing"], ["3", "Test every estimator on synthetic truths, incl. nuisance noise"], ["4", "Validate noise floors with two independent controls"], ["5", "Freeze, confirm once, log every deviation"], ["6", "Ask a fresh agent to refute surprising conclusions"]];
  r.forEach((x, i) => { const col = i % 2, row = Math.floor(i / 2); badge(s, 0.7 + col * 6.2, 1.6 + row * 1.7, x[0], col ? C.accent2 : C.accent1, 0.7); s.addText(x[1], { x: 1.6 + col * 6.2, y: 1.55 + row * 1.7, w: 5.0, h: 0.8, fontSize: 17, color: C.text1, isTextBox: true, margin: 0, valign: "middle" }); });
}
{
  const s = pres.addSlide({ masterName: "TITLE_DARK", sectionTitle: sectionName });
  s.addText("Conclusion: AI-led research", { placeholder: "title" });
  s.addText([
    { text: "In ~23 compute hours an AI team designed, pre-registered, ran and verified a study of 1.4 M protein pairs, and corrected its own headline.", options: { breakLine: true } },
    { text: "Limits: shared blind spots across agents; no human verification yet; predicted structures only.", options: { breakLine: true, color: C.accent2 } },
    { text: "Next: human audit · AlphaFold runs with varied MSA depth · PDB same-protein floors · time-calibrated structural phylogenetics · agentic labs at scale", options: { color: C.background2 } },
  ], { placeholder: "body" });
  s.addNotes("Close on the two-sided message: the AI-led process was fast, transparent and self-correcting, but its conclusions still need human verification.");
}
{
  const s = content("Acknowledgements and open slots", "Fill these before presenting.");
  card(s, 0.6, 1.4, 5.9, 4.0, "To fill in", null, {});
  s.addText(bullets(["Hackathon dates and participants: [__]", "Human contributors / project lead: [__]", "Project repository link: [__]", "Verification status by humans: [pending]"]), { x: 0.85, y: 2.0, w: 5.4, h: 3.3, fontSize: 17, color: C.text1, isTextBox: true, margin: 0, valign: "top", paraSpaceAfter: 10 });
  card(s, 6.8, 1.4, 5.9, 4.0, "Resources used", null, { fill: "FDEBD9" });
  s.addText(bullets(["AlphaFold DB & AFDB clusters (Steinegger, Beltrao labs; CC-BY 4.0)", "Foldseek, FoldMason, Q.3Di.AF", "Ancestra (hackathon Phase 1)", "TimeTree (data not redistributed)", "UniProt, PDB/SIFTS, TED, Pfam"]), { x: 7.05, y: 2.0, w: 5.4, h: 3.3, fontSize: 17, color: C.text1, isTextBox: true, margin: 0, valign: "top", paraSpaceAfter: 10 });
}

(async () => {
  const out = "project_slides.pptx";
  await pres.writeFile({ fileName: out });
  let applyTheme = null;
  for (const p of [process.env.PPTX_SKILL_DIR && path.join(process.env.PPTX_SKILL_DIR, "scripts/apply_theme.js")].filter(Boolean)) {
    if (fs.existsSync(p)) applyTheme = require(p).applyTheme;
  }
  if (applyTheme) await applyTheme(out, THEME);
  console.log("wrote", out, applyTheme ? "(theme applied)" : "(theme NOT applied)");
})();
