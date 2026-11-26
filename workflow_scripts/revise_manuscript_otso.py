"""Surgical Otso-review revision of the existing JSDM Studio manuscript DOCX."""

import sys
import re
from copy import deepcopy
from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt


def replace_paragraph(paragraph, value):
    original_rpr = None
    for run in paragraph.runs:
        if run.text.strip() and run._r.rPr is not None:
            original_rpr = deepcopy(run._r.rPr)
            break
    paragraph.text = value
    if original_rpr is not None and paragraph.runs:
        paragraph.runs[0]._r.insert(0, original_rpr)


def remove_paragraph(paragraph):
    element = paragraph._element
    element.getparent().remove(element)


def remove_table_row(table, row_index):
    element = table.rows[row_index]._tr
    element.getparent().remove(element)


def main(source, target):
    doc = Document(source)
    paragraphs = list(doc.paragraphs)
    if not paragraphs[99].text.startswith("Figure 1. Conceptual links"):
        raise ValueError("Unexpected source manuscript; Figure 1 anchor changed")
    if not paragraphs[106].text.startswith("Figure 4 Statistical Universal Benchmark"):
        raise ValueError("Unexpected source manuscript; Figure 4 anchor changed")

    replacements = {
        5: "Joint species distribution models (JSDMs) provide complementary approaches to analysing multispecies responses, prediction and residual dependence, but their software implementations differ in inputs, model structure, diagnostics and inferential targets. These differences make multi-engine analyses difficult to reproduce and easy to compare inappropriately.",
        7: "We demonstrate the workflow with four provenance-checked empirical cases spanning bryophytes, vascular plants, deadwood fungi and breeding birds. Each case illustrates a distinct data structure or modelling question, rather than one universal community-assembly framework. Software behaviour was additionally checked with synthetic workflows and a Universal Benchmark whose data-generating structure is now explicit. Recovery statistics are interpreted conditionally on that structure and the fitted model's estimand.",
        8: "JSDM Studio is an application and reproducibility layer, not a new JSDM. It helps users record and repeat environmental-response, trait, phylogenetic, spatial, detection and prediction analyses while keeping engine-specific assumptions and outputs visible. Ecological process claims require independent study design and evidence beyond a successful model fit.",
        11: "Joint species distribution models (JSDMs) analyse multiple species jointly, often estimating species-specific environmental responses and some form of residual dependence. Depending on the formulation, they may also represent spatial or temporal structure, species traits, phylogenetic relationships and imperfect detection (Pollock et al., 2014; Warton et al., 2015; Ovaskainen et al., 2017; Tikhonov et al., 2020). The resulting quantities are not interchangeable: a latent factor, a species-association matrix and an occupancy probability answer different questions. Residual association is not direct evidence of a biotic interaction, and spatial pattern alone does not identify dispersal (Dormann et al., 2018; Blanchet et al., 2020). These interpretive limits make transparent model specification and diagnostics essential.",
        15: "Here, we introduce JSDM Studio, an open-source graphical workbench for reproducible multi-engine joint species distribution modelling. It integrates Hmsc, Hmsc-HPC, jSDM, GJAM, spOccupancy, sjSDM and boral through engine-specific workflows for data checking, model specification, fitting, diagnosis, prediction and export. The software records configurations, data transformations, scripts and structured outputs without treating different model parameters as a single common estimand. We demonstrate its use with four empirical cases and assess software execution with synthetic tests. A documented benchmark allows sensitivity to the chosen data-generating structure to be examined, but a single simulated truth is not used to rank engines.",
        18: paragraphs[18].text.replace("(Figure 2)", "(Figure 1)"),
        28: "A multi-engine simulation is only interpretable after its data-generating mechanism and target estimands have been specified (Williams et al., 2024). The original Universal Benchmark used one correlated latent probit structure, which aligns more closely with some fitted binary models than with others. JSDM Studio now offers two explicitly named Studio-generated mechanisms: correlated latent probit occurrence and logistic occupancy with a separate replicated-detection process. Each mechanism produces matched site, species and predictor identifiers and exports the generating code, parameters and truth tables. These are structural templates written by JSDM Studio, not exact native simulations from any of the seven packages; even a matching link is only a partial within-structure check.",
        29: "The benchmark records whether an engine's configured link or observation process is partly aligned with the generator, or whether the fit is cross-structure. Environmental coefficients are evaluated only on a compatible scale, and prediction probabilities are assessed against the simulated truth at the training sites; the current comparison is not an independent hold-out test. Species-association quantities retain their engine-specific definitions and cannot be ranked as if they were the same parameter (Warton et al., 2015; Wilkinson et al., 2019, 2021). Both generating mechanisms can be run sequentially with the same design settings, but a rigorous within-engine validation would additionally require native, engine-specific simulation, replicated seeds and convergence assessment.",
        30: "Release validation also includes synthetic workflows spanning the seven execution targets (Table S2). These tests check execution, input transformation, diagnostic generation and output integrity; they do not establish scientific adequacy or comparative statistical superiority (Barros et al., 2023). The two-structure benchmark is a sensitivity design rather than a neutral league table. Manuscript-level recovery claims require repeated simulations, justified sampler settings and model-specific evaluation of parameter and prediction targets. The tested environments and dependencies are documented in Tables S2 and S3.",
        32: "We retained four empirical cases with traceable inputs: Tibetan bryophytes, Siskiyou vascular plants, deadwood fungi and Finnish breeding birds (Table 1). They address different questions about environmental responses, traits, phylogeny, residual structure and projection, without assigning all cases to a metacommunity paradigm. A previously labelled freshwater example was excluded because its taxonomic provenance could not be reconciled with the claimed source. Table 2 states each retained case's question and inferential limit; details are provided in Supporting Information S1.",
        34: "The Tibetan bryophyte case combined community observations for 190 species at 30 sites with climatic, solar and snow-related predictors, functional-trait axes, phylogenetic information and spatial structure (Gu et al., 2026). We asked how measured environmental gradients related to occurrence and whether species traits or phylogenetic relatedness helped explain variation among responses. Spatial residuals were treated as unexplained structure, not a direct estimate of dispersal limitation.",
        37: "The Siskiyou vascular-plant example analysed 75 species across 52 sites distributed along Whittaker's topographic moisture gradient (Damschen et al., 2010; Miller et al., 2018, 2019). We asked how species responses varied with dryness, whether measured traits explained some of that variation, and whether related species showed similar fitted responses. The case has one observed gradient and a taxonomy-derived surrogate tree; it does not by itself identify an assembly mechanism.",
        40: "The deadwood-fungus case used communities of wood-inhabiting fungi from Norway spruce logs to examine fitted responses to wood-decay stage and sequencing effort and residual dependence after those predictors were included (Ovaskainen et al., 2013; Ovaskainen & Abrego, 2020). Logs are the sampling units. Read counts are not direct fungal biomass, and residual association cannot establish pairwise interaction or stochastic colonisation.",
        36: "3.2 Siskiyou vascular plants: moisture-gradient responses and traits",
        39: "3.3 Deadwood fungi: decay-stage responses and residual dependence",
        45: "3.4 Finnish birds: climate, habitat and spatial prediction",
        46: "The Finnish bird case links environmental-response modelling to geographic prediction using the breeding-bird dataset employed in the Hmsc worked analysis (Tikhonov et al., 2020; Ovaskainen & Abrego, 2020). The analysis included 50 species surveyed on 137 routes, with habitat, spring climate, functional traits, phylogenetic information and a national prediction grid. It asks how climate and habitat relate to occurrence and whether traits or phylogeny help explain variation among species responses; these associations do not by themselves establish an assembly process.",
        48: "3.5 Cross-case interpretation",
        49: "Across the four retained cases, species responded differently to measured environmental gradients. Trait and phylogenetic patterns were context dependent, while residual spatial and association structures were less directly interpretable. We describe these as fitted statistical patterns rather than assigning the cases to one metacommunity paradigm or interpreting residuals as proof of dispersal or interaction. Such process claims require additional design or independent evidence.",
        50: "The cases illustrate why JSDM Studio harmonises workflows rather than conclusions. Environmental effects, trait mediation, phylogenetic parameters, residual dependence, observation processes and predictions remain distinct model components. The existing case-study outputs span all seven engines across different cases, but they are not a same-dataset seven-engine comparison. We therefore do not infer qualitative agreement among all seven models from those results. A future cross-engine worked example should use one defensible response target, display output availability for every engine and mark models requiring unavailable repeat-survey data as not applicable.",
        53: "The empirical examples show the value of keeping ecological questions and statistical outputs separate. Environmental responses, traits, phylogeny, spatial terms and residual associations can help formulate hypotheses about community structure, but their causal interpretation depends on the sampling design, the measured predictors and the model assumptions. JSDM Studio makes these choices and limitations visible across several modelling frameworks rather than prescribing one community-assembly theory for all applications. In particular, residual association need not indicate interaction, and spatial variation need not indicate dispersal.",
        54: "The Universal Benchmark illustrates another limit of cross-model interpretation. A single probit generator is closer to some fitted binary-link models than to a detection-aware logistic occupancy model; comparing raw recovery or runtime across those fits can therefore be misleading. The revised software records its generating structure and offers a second logistic occupancy-detection design, with results labelled as partly structure-aligned or cross-structure. This exposes sensitivity to model misspecification, but it does not yet constitute full native within-engine validation. Further work should simulate from each package's supported model structures where possible, repeat the design over seeds and sizes, and compare only matched prediction and parameter targets (Wilkinson et al., 2019, 2021; Williams et al., 2024).",
        91: "Table 1. Four empirical case-study datasets retained for the JSDM Studio worked examples. A previously labelled freshwater example was excluded pending source verification. B, binary occurrence; C, count or frequency; G, continuous or log-transformed conditional response.",
        95: "Table 2. Ecological questions, analytical focus and inferential limits across four JSDM Studio case studies.",
        117: "Ecological case studies and interpretation in JSDM Studio",
        119: "Four empirical case studies with traceable inputs are retained. They show how JSDM Studio accommodates contrasting data structures and ecological questions. The previously labelled freshwater example is withheld pending data-source verification. Interpretation is limited to the evidence available in each dataset and the estimands of its fitted models.",
        122: "S2 Ecological questions and analytical focus across case studies",
        123: "The retained cases differ in environmental gradients, response measurements, trait and phylogenetic information, and spatial support. Table S1 identifies the analytical focus and principal inferential limit for each case rather than assigning one metacommunity or assembly framework to all four.",
        124: "Table S1. Ecological questions, analytical focus and interpretive limits across four retained empirical cases.",
        277: "The Finnish bird case asked whether habitat and spring climate were associated with community composition and spatial species richness, whether migration strategy or body size helped explain differences among species responses, and whether those responses showed phylogenetic structure. The fitted Hmsc model indicated strong climate-related responses and additional habitat variation, whereas trait mediation was limited and phylogenetic structuring was negligible. Spatial projection showed a south-to-north richness gradient across Finland, but the reported prediction metrics are primarily in-sample and do not replace spatially blocked validation.",
        333: "Across the four retained cases, measured environmental gradients were associated with species-specific differences: climate, solar and snow variables in Tibetan bryophytes; topographic moisture in Siskiyou vascular plants; wood decay in deadwood fungi; and spring climate and habitat in Finnish birds. These fitted relationships are illustrative, not a quantitative cross-system comparison or direct evidence for a particular community-assembly process.",
        334: "Trait and phylogenetic results were more context dependent. Bryophyte trait-PC relationships and fitted phylogenetic structure were uncertain; leaf C:N modified the vascular-plant moisture response; and only limited trait or phylogenetic signal was supported for Finnish birds. The deadwood case does not include trait or phylogenetic inputs. None of these model summaries substitutes for independent trait or evolutionary-process evidence.",
        335: "The bryophyte Omega matrix and the deadwood Hmsc, GJAM and boral association summaries describe engine-specific residual dependence after fitted covariates. These matrices are not equivalent estimands and do not identify pairwise ecological interactions. Spatial random effects in the bryophyte and bird analyses likewise describe unexplained structure, not measured dispersal limitation.",
        336: "Projection is separate from estimation and validation. Bryophyte and vascular-plant examples project along supplied environmental gradients; the Finnish bird example projects to a national grid (Figures S4, S9 and S24-S25). These are model-based expectations at new scenarios, not independent field observations. In-sample bird fit metrics do not replace spatially blocked validation (Figure S26).",
        102: paragraphs[102].text.replace("Figure 2.", "Figure 1.", 1),
        104: "Figure 2. Hmsc worked-example outputs generated through JSDM Studio. (A) Species-specific environmental effects and variance partitioning. (B) Hmsc residual species-association matrix (Omega). (C) Posterior phylogenetic-signal parameter (rho). (D) Trait-related responses and trait-environment coefficients (Gamma). (E) Finnish bird spatial projections. This figure demonstrates the Hmsc workflow; it does not show equivalent outputs from all seven engines. Residual association is conditional on the fitted model and is not evidence of direct interaction.",
    }
    for index, value in replacements.items():
        replace_paragraph(paragraphs[index], value)

    for table_to_trim in doc.tables[:3]:
        remove_table_row(table_to_trim, 4)

    table = doc.tables[1]
    table.cell(0, 2).text = "Analytical focus and inferential limit"
    focus = [
        "Climate, snow and trait responses; phylogenetic structure and spatial residuals. Spatial pattern alone does not identify dispersal limitation.",
        "Species response to one topographic moisture gradient and trait mediation. The taxonomy-derived tree is a surrogate, not an observed molecular phylogeny.",
        "Decay-stage and sequencing-effort effects on fungal detections. Read counts are not direct biomass measurements; residual association is not interaction proof.",
        "Habitat and climate responses, trait variation and Finland-wide projection. Prediction grid cells are unsurveyed projections, not independent validation sites.",
    ]
    for row, value in zip(table.rows[1:], focus):
        row.cells[2].text = value
    table.autofit = False
    for column, width in zip(table.columns, (2.0, 3.4, 3.6)):
        column.width = Inches(width)
    for row in table.rows:
        for cell, width in zip(row.cells, (2.0, 3.4, 3.6)):
            cell.width = Inches(width)
            tc_mar = cell._tc.get_or_add_tcPr().first_child_found_in("w:tcMar")
            if tc_mar is None:
                tc_mar = OxmlElement("w:tcMar")
                cell._tc.get_or_add_tcPr().append(tc_mar)
            for side in ("top", "left", "bottom", "right"):
                element = tc_mar.find(qn(f"w:{side}"))
                if element is None:
                    element = OxmlElement(f"w:{side}")
                    tc_mar.append(element)
                element.set(qn("w:w"), "90")
                element.set(qn("w:type"), "dxa")
            for para in cell.paragraphs:
                para.alignment = WD_ALIGN_PARAGRAPH.LEFT
                para.paragraph_format.space_after = Pt(0)
                para.paragraph_format.line_spacing = 1.05
                for run in para.runs:
                    run.font.name = "Times New Roman"
                    run.font.size = Pt(9)

    supplemental = doc.tables[2]
    supplemental.cell(0, 3).text = "Analytical focus and limit"
    subset = [
        "30 sites; 190 species; 15 predictors; five trait PCs; source tree; ten projection scenarios",
        "52 sites; 75 species; one moisture gradient; three traits; surrogate tree",
        "Worked log subset; decay stage and sequencing effort; no spatial coordinates",
        "137 routes; 50 species; habitat and climate; traits; tree; Finland grid",
    ]
    question = [
        "How do climate, snow, traits and phylogeny relate to bryophyte occurrence?",
        "How do species and traits respond to the moisture gradient?",
        "How do fungal detections change with decay stage and read depth?",
        "Can habitat and climate predict bird occurrence across Finland?",
    ]
    limit = [
        "Spatial residuals do not identify dispersal.",
        "Tree is taxonomy-derived, not molecular.",
        "Read counts are not direct biomass; associations are not interactions.",
        "Grid cells are projections, not validation surveys.",
    ]
    result = [
        "Environmental responses and trait/phylogenetic summaries; see figures for uncertainty.",
        "Species responses vary along the measured moisture gradient.",
        "Decay-stage responses vary among taxa; residual signals need caution.",
        "Fitted occurrence surfaces project richness and species probabilities.",
    ]
    supplemental.autofit = False
    widths = (1.45, 1.75, 1.8, 1.75, 2.3, 0.55)
    for column, width in zip(supplemental.columns, widths):
        column.width = Inches(width)
    for row_index, row in enumerate(supplemental.rows):
        row.height = None
        row.height_rule = None
        if row_index:
            for col_index, value in zip((1, 2, 3, 4),
                                        (subset[row_index - 1], question[row_index - 1],
                                         limit[row_index - 1], result[row_index - 1])):
                row.cells[col_index].text = value
        for cell, width in zip(row.cells, widths):
            cell.width = Inches(width)
            tc_mar = cell._tc.get_or_add_tcPr().first_child_found_in("w:tcMar")
            if tc_mar is None:
                tc_mar = OxmlElement("w:tcMar")
                cell._tc.get_or_add_tcPr().append(tc_mar)
            for side in ("top", "left", "bottom", "right"):
                element = tc_mar.find(qn(f"w:{side}"))
                if element is None:
                    element = OxmlElement(f"w:{side}")
                    tc_mar.append(element)
                element.set(qn("w:w"), "70")
                element.set(qn("w:type"), "dxa")
            for para in cell.paragraphs:
                para.alignment = WD_ALIGN_PARAGRAPH.LEFT
                para.paragraph_format.space_after = Pt(0)
                para.paragraph_format.line_spacing = 1.05
                for run in para.runs:
                    run.font.name = "Times New Roman"
                    run.font.size = Pt(8.5)

    # Remove the theory illustration and the old one-generator benchmark plot.
    for index in list(range(239, 275)) + [44, 43, 42, 106, 105, 99, 98, 10]:
        remove_paragraph(paragraphs[index])
    for paragraph in doc.paragraphs:
        value = paragraph.text
        if value.startswith("S7 Finnish birds"):
            value = value.replace("S7 Finnish birds", "S6 Finnish birds", 1)
        elif value.startswith("S8 Cross-case synthesis"):
            value = value.replace("S8 Cross-case synthesis", "S7 Cross-case synthesis", 1)
        value = re.sub(r"S(2[0-6])\b", lambda match: f"S{int(match.group(1)) - 4}", value)
        if value != paragraph.text:
            replace_paragraph(paragraph, value)
    for table in doc.tables:
        for row in table.rows:
            for cell in row.cells:
                for paragraph in cell.paragraphs:
                    value = re.sub(r"S(2[0-6])\b", lambda match: f"S{int(match.group(1)) - 4}",
                                   paragraph.text)
                    if value != paragraph.text:
                        replace_paragraph(paragraph, value)
    for paragraph in list(doc.paragraphs):
        if paragraph.text.startswith(("Jeliazkov,", "Leibold,", "Mellado-Diaz,", "Mellado-D")):
            remove_paragraph(paragraph)
    doc.save(target)
    print(target)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("Usage: revise_manuscript_otso.py input.docx output.docx")
    main(sys.argv[1], sys.argv[2])
