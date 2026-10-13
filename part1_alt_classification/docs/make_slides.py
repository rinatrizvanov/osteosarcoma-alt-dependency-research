"""
Generate an editable project_slides.pptx from scratch using python-pptx.
Run: python3 docs/make_slides.py
Output: docs/project_slides.pptx
"""

from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from pptx.util import Inches, Pt
import copy
from lxml import etree

# ── Palette ────────────────────────────────────────────────────────────────
NAVY    = RGBColor(0x0F, 0x34, 0x60)
RED     = RGBColor(0xE9, 0x45, 0x60)
TEAL    = RGBColor(0xA8, 0xD8, 0xEA)
WHITE   = RGBColor(0xFF, 0xFF, 0xFF)
LIGHT   = RGBColor(0xF0, 0xF4, 0xF8)
DARK    = RGBColor(0x1A, 0x1A, 0x2E)
GREY    = RGBColor(0x88, 0x88, 0x88)
GREEN   = RGBColor(0x27, 0xAE, 0x60)
AMBER   = RGBColor(0xFF, 0xC1, 0x07)

W = Inches(13.33)   # widescreen 16:9
H = Inches(7.5)

prs = Presentation()
prs.slide_width  = W
prs.slide_height = H

BLANK = prs.slide_layouts[6]   # completely blank

# ── Helpers ────────────────────────────────────────────────────────────────

def add_rect(slide, x, y, w, h, fill=None, line=None, line_w=Pt(0)):
    shape = slide.shapes.add_shape(1, x, y, w, h)   # MSO_SHAPE_TYPE.RECTANGLE = 1
    shape.line.width = line_w
    if fill:
        shape.fill.solid()
        shape.fill.fore_color.rgb = fill
    else:
        shape.fill.background()
    if line:
        shape.line.color.rgb = line
    else:
        shape.line.fill.background()
    return shape


def add_text(slide, text, x, y, w, h,
             size=18, bold=False, color=DARK, align=PP_ALIGN.LEFT,
             wrap=True, italic=False):
    txb = slide.shapes.add_textbox(x, y, w, h)
    txb.word_wrap = wrap
    tf = txb.text_frame
    tf.word_wrap = wrap
    p = tf.paragraphs[0]
    p.alignment = align
    run = p.add_run()
    run.text = text
    run.font.size = Pt(size)
    run.font.bold = bold
    run.font.italic = italic
    run.font.color.rgb = color
    return txb


def add_para(tf, text, size=18, bold=False, color=DARK,
             align=PP_ALIGN.LEFT, italic=False, space_before=0, bullet=False):
    from pptx.util import Pt
    p = tf.add_paragraph()
    p.alignment = align
    p.space_before = Pt(space_before)
    if bullet:
        p.level = 1
    run = p.add_run()
    run.text = text
    run.font.size = Pt(size)
    run.font.bold = bold
    run.font.italic = italic
    run.font.color.rgb = color
    return p


def title_slide(prs, title, subtitle="", byline=""):
    slide = prs.slides.add_slide(BLANK)
    # Background
    add_rect(slide, 0, 0, W, H, fill=NAVY)
    # Red accent bar
    add_rect(slide, 0, H - Inches(0.08), W, Inches(0.08), fill=RED)
    # Title
    t = slide.shapes.add_textbox(Inches(1), Inches(1.8), Inches(11.3), Inches(2))
    t.word_wrap = True
    tf = t.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.alignment = PP_ALIGN.LEFT
    run = p.add_run()
    run.text = title
    run.font.size = Pt(40)
    run.font.bold = True
    run.font.color.rgb = WHITE
    # Subtitle
    if subtitle:
        add_para(tf, subtitle, size=22, color=TEAL, space_before=6)
    # Byline
    t2 = slide.shapes.add_textbox(Inches(1), Inches(5.8), Inches(11), Inches(1))
    tf2 = t2.text_frame
    p2 = tf2.paragraphs[0]
    run2 = p2.add_run()
    run2.text = byline
    run2.font.size = Pt(16)
    run2.font.color.rgb = GREY
    return slide


def section_slide(prs, heading, sub=""):
    slide = prs.slides.add_slide(BLANK)
    add_rect(slide, 0, 0, W, H, fill=NAVY)
    add_rect(slide, 0, H - Inches(0.08), W, Inches(0.08), fill=RED)
    t = slide.shapes.add_textbox(Inches(1.5), Inches(2.5), Inches(10.3), Inches(2))
    t.word_wrap = True
    tf = t.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.alignment = PP_ALIGN.CENTER
    run = p.add_run()
    run.text = heading
    run.font.size = Pt(44)
    run.font.bold = True
    run.font.color.rgb = WHITE
    if sub:
        add_para(tf, sub, size=20, color=TEAL, align=PP_ALIGN.CENTER, space_before=8)
    return slide


def content_slide(prs, title_text):
    slide = prs.slides.add_slide(BLANK)
    # Top accent
    add_rect(slide, 0, 0, W, Inches(0.08), fill=NAVY)
    # Left red bar
    add_rect(slide, Inches(0.4), Inches(0.18), Inches(0.06), Inches(0.7), fill=RED)
    # Title
    t = slide.shapes.add_textbox(Inches(0.6), Inches(0.12), Inches(11.8), Inches(0.85))
    t.word_wrap = True
    tf = t.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.alignment = PP_ALIGN.LEFT
    run = p.add_run()
    run.text = title_text
    run.font.size = Pt(28)
    run.font.bold = True
    run.font.color.rgb = NAVY
    # Divider line
    add_rect(slide, Inches(0.4), Inches(1.0), Inches(12.5), Inches(0.025), fill=LIGHT)
    # Page number placeholder (bottom right)
    return slide


def bullet_box(slide, items, x=Inches(0.5), y=Inches(1.15),
               w=Inches(12.3), h=Inches(5.9),
               size=18, color=DARK, indent="• "):
    txb = slide.shapes.add_textbox(x, y, w, h)
    txb.word_wrap = True
    tf = txb.text_frame
    tf.word_wrap = True
    first = True
    for item in items:
        if first:
            p = tf.paragraphs[0]
            first = False
        else:
            p = tf.add_paragraph()
        p.space_before = Pt(4)
        run = p.add_run()
        is_sub = item.startswith("  ")
        run.text = ("  " if is_sub else "") + indent + item.lstrip()
        run.font.size = Pt(size - 2 if is_sub else size)
        run.font.color.rgb = color
    return txb


def two_col(slide, left_items, right_items,
            left_title="", right_title="",
            size=17):
    mid = Inches(6.7)
    top = Inches(1.15)
    col_w = Inches(5.9)
    col_h = Inches(5.8)
    if left_title:
        add_text(slide, left_title, Inches(0.5), top, col_w, Inches(0.4),
                 size=16, bold=True, color=NAVY)
        top2 = top + Inches(0.45)
    else:
        top2 = top
    bullet_box(slide, left_items, x=Inches(0.5), y=top2, w=col_w, h=col_h, size=size)
    add_rect(slide, mid - Inches(0.02), Inches(1.1), Inches(0.025), Inches(6.0), fill=LIGHT)
    if right_title:
        add_text(slide, right_title, mid + Inches(0.2), top, col_w, Inches(0.4),
                 size=16, bold=True, color=NAVY)
        top3 = top + Inches(0.45)
    else:
        top3 = top
    bullet_box(slide, right_items, x=mid + Inches(0.2), y=top3, w=col_w, h=col_h, size=size)


def table_slide(slide, headers, rows,
                x=Inches(0.4), y=Inches(1.1),
                w=Inches(12.5), col_widths=None):
    n_cols = len(headers)
    n_rows = len(rows)
    tbl = slide.shapes.add_table(n_rows + 1, n_cols, x, y, w, Inches(0.38 * (n_rows + 1))).table
    if col_widths:
        for i, cw in enumerate(col_widths):
            tbl.columns[i].width = cw

    def cell_text(cell, text, bold=False, bg=None, color=WHITE, size=15, align=PP_ALIGN.LEFT):
        cell.text = text
        p = cell.text_frame.paragraphs[0]
        p.alignment = align
        run = p.runs[0] if p.runs else p.add_run()
        run.font.size = Pt(size)
        run.font.bold = bold
        run.font.color.rgb = color
        if bg:
            cell.fill.solid()
            cell.fill.fore_color.rgb = bg

    for ci, h in enumerate(headers):
        cell_text(tbl.cell(0, ci), h, bold=True, bg=NAVY, color=WHITE, size=15)
    for ri, row in enumerate(rows):
        bg = LIGHT if ri % 2 == 0 else WHITE
        for ci, val in enumerate(row):
            cell_text(tbl.cell(ri + 1, ci), str(val), bg=bg, color=DARK, size=14)
    return tbl


# ══════════════════════════════════════════════════════════════════════════════
# Slide 1 — Title
# ══════════════════════════════════════════════════════════════════════════════
title_slide(prs,
    "Multi-Source ALT Classification\nin Osteosarcoma",
    "A Computational Feasibility Study",
    "Rinat Rizvanov  ·  June 2026  ·  Independent computational analysis")

# ══════════════════════════════════════════════════════════════════════════════
# Slide 2 — Agenda
# ══════════════════════════════════════════════════════════════════════════════
sl = content_slide(prs, "Agenda")
bullet_box(sl, [
    "Scientific motivation — OS, ALT, TCAB1",
    "Research aims and go/no-go framework",
    "Pipeline architecture and data sources",
    "Labeling framework — four-tier hierarchy",
    "Feasibility audit — Pass A and Pass B",
    "RNA availability crisis and three recovery routes",
    "Aim 2 discovery attempt — a null result (Notebook 11)",
    "Sex confounding diagnosis (Notebook 12, D1 & D2)",
    "Hypothesis-driven curated gene score — D3",
    "Robustness analysis (Notebook 13)",
    "Limitations and future directions",
], size=19)

# ══════════════════════════════════════════════════════════════════════════════
# Slide 3 — Section: Scientific Background
# ══════════════════════════════════════════════════════════════════════════════
section_slide(prs, "Scientific Background", "Why osteosarcoma, ALT, and TCAB1?")

# ══════════════════════════════════════════════════════════════════════════════
# Slide 4 — Osteosarcoma: The Problem
# ══════════════════════════════════════════════════════════════════════════════
sl = content_slide(prs, "Osteosarcoma: The Problem")
bullet_box(sl, [
    "Most common primary bone malignancy in children and adolescents",
    "5-year survival: ~65–70% localized  /  <20% metastatic or relapsed",
    "Standard of care = MAP chemo (methotrexate, doxorubicin, cisplatin) — unchanged since 1980s",
    "Extreme genomic instability: near-universal TP53 disruption, complex SVs",
    "",
    "Defining molecular feature:",
    "  ~50% of osteosarcomas use Alternative Lengthening of Telomeres (ALT)",
    "  instead of telomerase — highest ALT prevalence of any tumor type studied",
    "",
    "ALT makes OS an ideal system to study ALT as a therapeutic vulnerability",
], size=18)

# ══════════════════════════════════════════════════════════════════════════════
# Slide 5 — What Is ALT?
# ══════════════════════════════════════════════════════════════════════════════
sl = content_slide(prs, "What Is ALT?")
two_col(sl,
    left_items=[
        "Homologous-recombination-based telomere maintenance",
        "Replaces telomerase in ~50% of osteosarcomas",
        "",
        "Molecular hallmarks:",
        "  ALT-associated PML bodies (APBs)",
        "  C-circle assay — extrachromosomal telomeric repeat DNA",
        "  Loss of ATRX/DAXX → aberrant H3.3 loading at telomeres",
        "  Heterogeneous, ultra-long telomeres",
    ],
    right_items=[
        "ALT cells depend on:",
        "  HR machinery for telomere synthesis",
        "  FANCM-mediated break-induced replication at telomeres",
        "",
        "These dependencies are NOT shared by telomerase+ cells",
        "→ Potential synthetic lethal targets",
        "",
        "Experimental assays for ALT:",
        "  C-circle assay, APB immunofluorescence, TRAP (telomerase activity)",
    ],
    left_title="Biology", right_title="Therapeutic relevance",
    size=17)

# ══════════════════════════════════════════════════════════════════════════════
# Slide 6 — TCAB1 / WRAP53
# ══════════════════════════════════════════════════════════════════════════════
sl = content_slide(prs, "TCAB1 / WRAP53 — The Original Question")
bullet_box(sl, [
    "TCAB1 (Telomerase Cajal Body protein 1), encoded by WRAP53 at 17p13.1",
    "Localizes telomerase to Cajal bodies — required for efficient telomere extension",
    "",
    "In an ALT-positive tumor:",
    "  TCAB1 function may be irrelevant or actively harmful",
    "  Evidence: TCAB1-disrupted tumors may depend more on ALT-mediated repair",
    "",
    "Two research questions:",
    "  1.  Is TCAB1 genomic disruption (17p13.1 SV breakpoints) associated",
    "      with a distinct transcriptional program in osteosarcoma?",
    "  2.  Can a transcriptional signature distinguish ALT+ from ALT− lines?",
], size=18)

# ══════════════════════════════════════════════════════════════════════════════
# Slide 7 — Section: Architecture & Data
# ══════════════════════════════════════════════════════════════════════════════
section_slide(prs, "Pipeline Architecture & Data", "What I built and where the data came from")

# ══════════════════════════════════════════════════════════════════════════════
# Slide 8 — Repository Architecture
# ══════════════════════════════════════════════════════════════════════════════
sl = content_slide(prs, "Repository Architecture")
two_col(sl,
    left_items=[
        "src/  — all reusable logic",
        "  manifest_builder.py  — cohort loaders",
        "  label_rules.py       — deterministic labels",
        "  sv_annotator.py      — 17p13.1 SV calls",
        "  audit.py             — Pass A / B logic",
        "  aim2_signature.py    — LOOCV training",
        "",
        "notebooks/  00–13  — thin orchestration layers",
        "tests/  — 26 unit tests, all passing",
        "data/manifests/  — master_manifest.csv (115×46)",
        "data/processed/  — reports, sha256 hashes",
    ],
    right_items=[
        "Design principles:",
        "",
        "Notebooks import from src/ —",
        "  one fix propagates everywhere",
        "",
        "Every output is sha256-verified —",
        "  byte-identical re-runs confirmed",
        "",
        "Analysis spec frozen in label_spec_v1.0.md",
        "  before any modeling code was written",
        "",
        "Prevents circular method refinement",
    ],
    left_title="Structure", right_title="Principles",
    size=16)

# ══════════════════════════════════════════════════════════════════════════════
# Slide 9 — Five Cohorts
# ══════════════════════════════════════════════════════════════════════════════
sl = content_slide(prs, "Five Cohorts, 115 Samples")
table_slide(sl,
    headers=["Cohort", "N", "Label type", "RNA-seq"],
    rows=[
        ["DepMap 26Q1", "35 cell lines", "ALTitude predictions", "✅ 37 with RNA"],
        ["PCGP 2014", "34 pediatric tumors", "SV + TP53 (peer-reviewed)", "✗ public release"],
        ["de Nonneville 2022", "22 bone tumors", "C-circle / APB (peer-reviewed)", "✗ microarray only"],
        ["Lawlor 2026 GCCRI", "17 PDX", "SV / IHC (preprint)", "✗ internal only"],
        ["Taglialatela 2026", "7 cell lines", "ALT list (preprint)", "✗"],
    ],
    col_widths=[Inches(3.2), Inches(2.0), Inches(4.3), Inches(2.8)],
    y=Inches(1.15), w=Inches(12.4)
)
add_text(sl,
    "Key tension: Labels and expression data are in completely disjoint datasets.\n"
    "The DepMap RNA matrix (1,719 cell lines × 19,177 genes) had no ALT labels.\n"
    "The labeled cohorts had no RNA-seq available in their public releases.",
    Inches(0.5), Inches(4.5), Inches(12.3), Inches(1.8),
    size=17, color=NAVY, italic=True)

# ══════════════════════════════════════════════════════════════════════════════
# Slide 10 — Section: Labeling Framework
# ══════════════════════════════════════════════════════════════════════════════
section_slide(prs, "Labeling Framework", "Four-tier hierarchy, frozen before any modeling")

# ══════════════════════════════════════════════════════════════════════════════
# Slide 11 — ALT Status Hierarchy
# ══════════════════════════════════════════════════════════════════════════════
sl = content_slide(prs, "ALT Status Hierarchy")
table_slide(sl,
    headers=["Tier", "Label", "Source", "Used in"],
    rows=[
        ["1", "alt_gold_positive", "C-circle / APB / TRAP, peer-reviewed", "Primary + Pragmatic"],
        ["1", "alt_negative_high_conf", "Same experimental evidence", "Primary + Pragmatic"],
        ["2", "alt_silver", "ALTitude prediction only (no wet-lab confirmation)", "Pragmatic only"],
        ["3", "alt_unknown", "No usable evidence", "Excluded from all modeling"],
    ],
    col_widths=[Inches(0.6), Inches(2.6), Inches(5.5), Inches(3.5)],
    y=Inches(1.15), w=Inches(12.4)
)
bullet_box(sl, [
    "Three parallel definitions: Strict (gold only) → primary   |   Pragmatic (gold+silver) → sensitivity   |   Broad → no new samples",
    "",
    "Critical: alt_call_direction (positive/negative/unknown) is the arm-assignment column",
    "  Silver labels include BOTH ALT+ (n=6) and ALT− (n=5) predictions",
    "  Using alt_status directly would assign all silver to the positive arm — wrong",
    "  alt_call_direction was added in Notebook 10 after recognizing this during the merge",
], size=17, y=Inches(4.0))

# ══════════════════════════════════════════════════════════════════════════════
# Slide 12 — Section: Feasibility Audit
# ══════════════════════════════════════════════════════════════════════════════
section_slide(prs, "Feasibility Audit", "Two checks before any modeling code was written")

# ══════════════════════════════════════════════════════════════════════════════
# Slide 13 — Pass A and Pass B
# ══════════════════════════════════════════════════════════════════════════════
sl = content_slide(prs, "Pass B (Confounding) + Pass A (Availability)")
two_col(sl,
    left_items=[
        "Pass B: Confounding Census",
        "",
        "Checks if confounders are imbalanced",
        "across TCAB1 and ALT arms",
        "",
        "Critical finding (Aim 1):",
        "  PCGP 2014 = 100% of TCAB1-intact arm",
        "  Any DE result would be completely",
        "  confounded with cohort of origin",
        "",
        "Documented: Halt override: True",
        "Audit runs to completion for records",
    ],
    right_items=[
        "Pass A: RNA Availability Census",
        "",
        "Checks how many labeled samples",
        "have file-backed RNA-seq on disk",
        "",
        "Initial result:",
        "  Aim 1 RNA-backed disrupted = 0",
        "  Aim 2 RNA-backed ALT-labeled = 0",
        "",
        "Stage 0 diagnostics (nb06, nb07)",
        "were designed to catch this before",
        "any analysis code was written",
        "→ not a failure, a designed check",
    ],
    left_title="", right_title="",
    size=17)

# ══════════════════════════════════════════════════════════════════════════════
# Slide 14 — Section: RNA Crisis & Recovery
# ══════════════════════════════════════════════════════════════════════════════
section_slide(prs, "RNA Crisis and Recovery", "Three routes — one succeeded")

# ══════════════════════════════════════════════════════════════════════════════
# Slide 15 — Three Recovery Routes
# ══════════════════════════════════════════════════════════════════════════════
sl = content_slide(prs, "Three Recovery Routes")
table_slide(sl,
    headers=["Route", "Description", "Outcome"],
    rows=[
        ["A", "de Nonneville → ArrayExpress E-MTAB-11743",
         "Affymetrix microarray (not RNA-seq) → external validation only"],
        ["B ✅", "ALTitude predictions → DepMap cell-line-name join",
         "Succeeded: 37 RNA-backed OS lines matched → Aim 2 rescued"],
        ["C", "Lawlor PDX internal RNA-seq (PI collaboration)",
         "Pending — critical dependency for Aim 1"],
    ],
    col_widths=[Inches(1.0), Inches(4.8), Inches(6.4)],
    y=Inches(1.15), w=Inches(12.4)
)
bullet_box(sl, [
    "The column-binding bug (Notebook 09):",
    "  Initial run: NOT ESTIMABLE — suspicious given 37 RNA-backed lines",
    "  Root cause: code read Curated.Marker (ATRX mutation feature) not ALT_pred_class",
    "  Caught by ground-truth check: U2OS returned ambiguous_or_missing → wrong column",
    "  Fixed; corrected inventory: 31 RNA-backed lines / 12 ALT+ / 7 ALT− gold",
    "",
    "After Notebook 10 direction-aware merge:",
    "  Primary (gold): N = 20  (12 ALT+ / 8 ALT−)     Pragmatic (gold+silver): N = 31  (18 ALT+ / 13 ALT−)",
], size=17, y=Inches(3.9))

# ══════════════════════════════════════════════════════════════════════════════
# Slide 16 — Section: Aim 2 Modeling
# ══════════════════════════════════════════════════════════════════════════════
section_slide(prs, "Aim 2 Modeling", "Notebook 11 → 12 → 13")

# ══════════════════════════════════════════════════════════════════════════════
# Slide 17 — Notebook 11: Null Result
# ══════════════════════════════════════════════════════════════════════════════
sl = content_slide(prs, "Notebook 11: Whole-Transcriptome Discovery — Null")
bullet_box(sl, [
    "Nested LOOCV, elastic-net logistic regression, 19,177 genes",
    "Every preprocessing step is fold-local (no leakage):",
    "  VarianceGeneFilter (top K=500) → StandardScaler → ElasticNet (inner 5-fold C selection)",
], size=18, h=Inches(1.5))
table_slide(sl,
    headers=["Definition", "N", "AUROC [95% CI]", "Permutation p"],
    rows=[
        ["Primary (gold)", "20", "0.4688  [0.20, 0.75]", "0.60  (null)"],
        ["Pragmatic (gold+silver)", "31", "0.6880  [0.45, 0.90]", "0.06  (borderline)"],
    ],
    col_widths=[Inches(3.0), Inches(0.8), Inches(4.5), Inches(4.0)],
    y=Inches(2.9), w=Inches(12.4)
)
bullet_box(sl, [
    "Primary = worse than random chance",
    "Top features exposed the problem:",
    "  Ranks 1–4: XAGE1B, MAGEB2, XAGE1A, MAGEA10  (X-linked cancer-testis antigens)",
    "  Ranks 17–18: USP9Y, DDX3Y  (Y-chromosome genes)",
    "  Zero ALT-gene-list genes in top 20 of either definition",
    "The model was learning sex, not ALT biology",
], size=17, y=Inches(4.55))

# ══════════════════════════════════════════════════════════════════════════════
# Slide 18 — Notebook 12: Sex Confounding
# ══════════════════════════════════════════════════════════════════════════════
sl = content_slide(prs, "Notebook 12: Sex Confounding Diagnosis")
table_slide(sl,
    headers=["Design", "Method", "Primary AUROC", "Perm p"],
    rows=[
        ["D1 — sex covariate", "Append Female/Male indicator after gene filter", "0.5000", "0.53"],
        ["D2 — no Y/XCT genes", "Remove 41 Y-chr + 4 X-CT antigen columns", "0.3958", "0.66"],
        ["D3 — curated ALT genes", "19-gene prespecified list, no variance filter", "0.9479", "< 0.001"],
    ],
    col_widths=[Inches(2.8), Inches(4.8), Inches(2.4), Inches(2.2)],
    y=Inches(1.15), w=Inches(12.4)
)
bullet_box(sl, [
    "D1 and D2 confirm: the borderline pragmatic signal in nb11 depended entirely on sex-linked features",
    "After sex correction: no ALT signal survives the unguided whole-transcriptome approach",
    "",
    "D3 is a different kind of test — a hypothesis test, not a discovery attempt",
    "  Feature set fixed before any model fitting using a pre-specified biological gene list",
    "  Bypasses the problem of finding a needle in 19,177-gene haystack at N=20",
], size=17, y=Inches(3.85))

# ══════════════════════════════════════════════════════════════════════════════
# Slide 19 — D3: Curated Gene Score Results
# ══════════════════════════════════════════════════════════════════════════════
sl = content_slide(prs, "D3: Hypothesis-Driven Curated Gene Score")
bullet_box(sl, [
    "21-gene ALT pathway list: TERT, TERC, ATRX, DAXX, SMARCAL1, TP53, WRN, BLM, MRE11A,",
    "  RAD50, NBN, FANCD2, RAD51AP1, WRAP53, RAD52, BRCA1, BRCA2, FANCM, SLX4, GEN1, MUS81",
    "19 of 21 found in DepMap matrix  (TERC = ncRNA;  MRE11A = symbol mismatch in DepMap)",
], size=17, h=Inches(1.3))
table_slide(sl,
    headers=["Definition", "N", "AUROC [95% CI]", "AUPRC", "Sens@90%Spec", "Perm p"],
    rows=[
        ["Primary (gold)", "20", "0.9479  [0.81, 1.00]", "0.9755", "91.7%", "< 0.001"],
        ["Pragmatic", "31", "0.8291  [0.67, 0.96]", "0.9166", "77.8%", "0.010"],
    ],
    col_widths=[Inches(2.2), Inches(0.7), Inches(2.9), Inches(1.5), Inches(2.3), Inches(2.5)],
    y=Inches(2.6), w=Inches(12.4)
)
bullet_box(sl, [
    "Primary nonzero coefficients (full-data refit, C = 0.3162):",
    "  TERT  −0.689   FANCM  +0.683   BRCA2  +0.197   WRN  −0.153   NBN  +0.105",
    "",
    "Important framing: this is a test of a pre-formed hypothesis, not a discovery",
    "  The result says these 19 genes carry information about ALT status in this cohort",
    "  It does not establish generalization to new cohorts — external validation is needed",
], size=17, y=Inches(4.55))

# ══════════════════════════════════════════════════════════════════════════════
# Slide 20 — Gene Interpretation
# ══════════════════════════════════════════════════════════════════════════════
sl = content_slide(prs, "Gene Interpretation: TERT and FANCM")
two_col(sl,
    left_items=[
        "TERT  (negative coefficient)",
        "",
        "Catalytic subunit of telomerase",
        "Telomerase and ALT are mutually exclusive",
        "→ Lower TERT expression = ALT-positive",
        "",
        "Most stable predictor across all analyses:",
        "  Bootstrap freq = 1.00 (both definitions)",
        "  LOOCV fold freq = 1.00 (all three cohorts)",
        "",
        "WRAP53  (zero coefficient)",
        "Central to the original question, but",
        "WRAP53 mRNA adds no information",
        "beyond TERT + FANCM at N = 20",
    ],
    right_items=[
        "FANCM  (positive coefficient)",
        "",
        "Fanconi anemia helicase required for",
        "break-induced replication at ALT telomeres",
        "",
        "Bhowmick et al. (2016): FANCM depletion",
        "selectively inhibits ALT activity",
        "→ Proposed synthetic lethal target in ALT+",
        "",
        "Higher FANCM expression = ALT-positive",
        "",
        "Bootstrap freq = 1.00 (both definitions)",
        "Second most stable predictor overall",
    ],
    left_title="", right_title="",
    size=17)

# ══════════════════════════════════════════════════════════════════════════════
# Slide 21 — Notebook 13: Robustness
# ══════════════════════════════════════════════════════════════════════════════
sl = content_slide(prs, "Notebook 13: Robustness Analysis")
table_slide(sl,
    headers=["Cohort", "N", "AUROC", "Perm p (1,000 shuffles)"],
    rows=[
        ["Primary", "20", "0.9479", "< 0.001  (0/1,000)"],
        ["Pragmatic", "31", "0.8291", "0.010  (10/1,000)"],
        ["Pragmatic, known-sex only", "27", "0.8765", "0.003  (3/1,000)"],
    ],
    col_widths=[Inches(3.3), Inches(0.8), Inches(2.0), Inches(6.1)],
    y=Inches(1.15), w=Inches(12.4)
)
bullet_box(sl, [
    "Extended permutation test (1,000 shuffles vs 100 in nb12):",
    "  Primary highly robust (0/1,000);  pragmatic nominally significant — interpret cautiously",
    "  Implementation: precomputed fold-local StandardScaler once (label-independent) → ~10× faster",
    "",
    "Known-sex-only pragmatic (drop 4 unknown-sex samples):",
    "  AUROC improves 0.8291 → 0.8765;  CI tightens;  p = 0.003",
    "  TP53 coefficient shrinks substantially (−0.337 → +0.101) when unknown-sex dropped",
    "  Suggests TP53's pragmatic contribution was entangled with label uncertainty",
], size=17, y=Inches(3.65))

# ══════════════════════════════════════════════════════════════════════════════
# Slide 22 — Bootstrap Stability
# ══════════════════════════════════════════════════════════════════════════════
sl = content_slide(prs, "Bootstrap Coefficient Stability (200 Resamples)")
table_slide(sl,
    headers=["Gene", "Primary freq", "Primary class", "Pragmatic freq", "Pragmatic class"],
    rows=[
        ["TERT",    "1.00", "stable",     "1.00", "stable"],
        ["FANCM",   "1.00", "stable",     "1.00", "stable"],
        ["BRCA2",   "0.81", "stable",     "0.655","borderline"],
        ["TP53",    "0.445","unstable",   "0.935","stable"],
        ["WRN",     "0.57", "borderline", "0.705","borderline"],
        ["NBN",     "0.625","borderline", "0.785","borderline"],
        ["ATRX",    "0.235","unstable",   "0.555","borderline"],
        ["12 others","< 0.30","unstable", "—",    "—"],
    ],
    col_widths=[Inches(1.8), Inches(2.0), Inches(2.2), Inches(2.2), Inches(4.0)],
    y=Inches(1.15), w=Inches(12.4)
)
add_text(sl,
    "TERT and FANCM: the only genes stable in both definitions. "
    "LOOCV fold frequency confirms both appear in every held-out fold (freq = 1.00). "
    "The result is not driven by any single sample.",
    Inches(0.5), Inches(6.5), Inches(12.3), Inches(0.9),
    size=17, color=NAVY, bold=True)

# ══════════════════════════════════════════════════════════════════════════════
# Slide 23 — Limitations
# ══════════════════════════════════════════════════════════════════════════════
sl = content_slide(prs, "Limitations")
two_col(sl,
    left_items=[
        "Sample size",
        "  Primary N = 20  (12+ / 8−)",
        "  AUROC CI = [0.81, 1.00] — wide",
        "  Expected AUROC in new samples",
        "  could be substantially lower",
        "",
        "Label quality",
        "  Gold labels: different groups,",
        "  different protocols, different studies",
        "  Silver (ALTitude): computationally derived,",
        "  never independently validated",
        "",
        "Single platform",
        "  DepMap 26Q1 log2(TPM+1) only",
        "  Cell lines ≠ primary tumor profiles",
    ],
    right_items=[
        "This is confirmatory, not discovery",
        "  D3 tests prior biological knowledge",
        "  High AUROC validates the prior —",
        "  it does not discover new biology",
        "  Framing it as novel signature discovery",
        "  would overstate the finding",
        "",
        "Permutation test conservatism",
        "  C fixed at median observed value",
        "  (fully nested null is intractable ~9h)",
        "  Conservative: true p may be smaller",
        "",
        "No probability calibration",
        "  LOOCV scores are not calibrated probabilities",
        "  This is a signature score, not a classifier",
    ],
    size=16)

# ══════════════════════════════════════════════════════════════════════════════
# Slide 24 — Future Directions
# ══════════════════════════════════════════════════════════════════════════════
sl = content_slide(prs, "Future Directions")
table_slide(sl,
    headers=["Priority", "What", "Dependency"],
    rows=[
        ["Near-term",
         "External validation on de Nonneville microarray (22 tumors, C-circle labels)\n"
         "Apply 19-gene curated score after within-platform normalization",
         "None — can pursue independently"],
        ["PI collaboration",
         "Aim 1: Lawlor PDX RNA-seq\n"
         "8 TCAB1-disrupted vs 32 intact → differential expression",
         "Internal RNA-seq from Lawlor 2026 GCCRI group"],
        ["DepMap only",
         "FANCM dependency in DepMap CRISPR (Chronos scores)\n"
         "Does FANCM dependency correlate with ALT status?",
         "None — data already on disk"],
        ["Alternative Aim 1",
         "WRAP53 expression-quantile contrast within DepMap\n"
         "Stratify 35 OS lines by WRAP53 expression level",
         "None — but answers a different question than Aim 1"],
    ],
    col_widths=[Inches(2.0), Inches(6.8), Inches(3.4)],
    y=Inches(1.15), w=Inches(12.4)
)

# ══════════════════════════════════════════════════════════════════════════════
# Slide 25 — Summary Table
# ══════════════════════════════════════════════════════════════════════════════
sl = content_slide(prs, "Summary")
table_slide(sl,
    headers=["What was done", "Result"],
    rows=[
        ["5-cohort master manifest (115 samples, 46 columns)", "Built and audited"],
        ["Four-tier ALT label hierarchy, frozen spec", "Documented in label_spec_v1.0.md"],
        ["Stage 0 RNA availability audit", "Both aims initially inestimable"],
        ["Three recovery routes evaluated", "Route B rescued Aim 2"],
        ["Column-binding bug caught and fixed", "U2OS sanity check"],
        ["Aim 2 whole-transcriptome LOOCV (nb11)", "NULL — sex confound identified"],
        ["Sex confounding diagnosis (nb12 D1, D2)", "Confirmed: no ALT signal after correction"],
        ["Curated gene score, primary N=20 (nb12 D3)", "AUROC = 0.948, p < 0.001"],
        ["1,000-shuffle permutation test (nb13)", "Primary 0/1,000; pragmatic 10/1,000"],
        ["Bootstrap + LOOCV fold stability (nb13)", "TERT and FANCM stable in all analyses"],
        ["Aim 1 (TCAB1 differential expression)", "NOT ESTIMABLE — pending Route C"],
    ],
    col_widths=[Inches(7.4), Inches(4.8)],
    y=Inches(1.1), w=Inches(12.4)
)

# ══════════════════════════════════════════════════════════════════════════════
# Slide 26 — Closing
# ══════════════════════════════════════════════════════════════════════════════
title_slide(prs,
    "Key Takeaway",
    "The ALT transcriptional signal in osteosarcoma is real but concentrated in\n"
    "a small number of functionally motivated genes.\n"
    "It cannot be found by unguided gene selection at N = 20,\n"
    "but it is recoverable when you know where to look.\n\n"
    "TERT and FANCM are robust across every analysis.\n"
    "WRAP53 mRNA does not contribute to the classification signal.",
    "Rinat Rizvanov  ·  github.com/rinatrizvanov")

# ── Save ────────────────────────────────────────────────────────────────────
import pathlib
out = pathlib.Path(__file__).parent / "project_slides.pptx"
prs.save(out)
print(f"Saved: {out}")
