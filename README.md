# XRD Reference Search & Phase ID

Automated phase identification for powder X-ray diffraction (XRD), built on [pymatgen](https://pymatgen.org).

Give it a powder XRD scan and the elements in your sample, and it returns every crystalline phase present (minor phases near 1 wt% included), how confident each assignment is, and the lattice parameters of each phase. It finds its own reference structures: first in a local library of CIF files, then in the [Crystallography Open Database (COD)](https://www.crystallography.net/cod/) and [Materials Project](https://materialsproject.org) only when the library has no convincing match. MOF collections ([CoRE MOF 2024](https://zenodo.org/records/15055758), [QMOF](https://doi.org/10.6084/m9.figshare.13147324)) and your own CIFs can be added. Everything it downloads joins the library, so each new sample starts with more.

The same code handles inorganic mixtures, battery cathode materials (NMC) and broad-peak MOFs.

Current version: **2026-10-23** (check with `print(X.__version__)`).

![Search and identification flow](docs/search_flow.png)

## Features

- **Library-first reference search** in stages (retry failed queries → local library → COD exact elements → Materials Project → related smaller systems → COD superset), stopping at the first stage whose answer passes three tests: pattern score, share of peak intensity explained, and cell shift. For a mixture, all identified phases are judged together.
- **Several structure sources with on/off switches** (`DATABASES`): COD and Materials Project online; CoRE MOF and QMOF imported once and used offline; your own CIFs from the CCDC (`CSD`), papers' supporting information (`Published`) or models (`local`). The ranking labels each candidate experimental, DFT, model or local.
- **Iterative phase identification** (Cullity's method): main phase first, then side phases ranked on the leftover peaks. Each phase gets a confidence label (clear / ambiguous / tentative), and a phase is never assigned from a single peak.
- **Consistency checks**: candidates are penalised for strong predicted lines missing from the scan, and masking is judged peak by peak so a weak predicted line cannot hide a strong measured peak.
- **Lab, synchrotron and broad-peak data**: Cu Kα1 or the Kα1/Kα2 doublet; broad patterns (MOFs, nanoparticles) are matched at their measured peak width.
- **Noise-based peak detection** that finds weak, broad humps and minor phases near 1 wt% in good data.
- **Lattice parameters for every identified phase** (cubic, tetragonal, hexagonal/trigonal, orthorhombic), using only each phase's own peaks with a shared 2θ zero offset.
- **Traceable results**: the CIFs of the identified phases and the best candidates are saved with `sources.csv` (database links and DOIs); plots draw a guide line through every measured peak, with unexplained peaks in red.
- **Fast repeat runs**: a disk cache of parsed structures and simulated patterns, parallel scoring on several CPU cores and vectorised matching make repeat runs about 5× faster, with identical rankings.
- **Built-in structure models** for R-3m layered oxides (LiCoO2, NMC111/532/622/811, LiNiO2) with optional Li/Ni cation mixing, plus M(OH)2 hydroxides and NiO.

## Repository contents

| Path | What it is |
|---|---|
| `xrd_refs.py` | All the code: database search, reference library, simulation, scoring, phase identification, lattice fits, plotting |
| `XRD_phase_ID.ipynb` | The notebook: settings, collection import, reference search, phase identification, lattice parameters, plots |
| `reference_library/` | Starter library of ~1,700 CIFs (830 COD, 894 Materials Project) and saved searches; see [its README](reference_library/README.md) |
| `examples/cpd-1a.txt` | Example scan: IUCr CPD round-robin sample 1a (fluorite + zincite + corundum); see [its README](examples/README.md) |
| `requirements.txt` | Python dependencies |

The bundled library covers the Al–Ca–F–Zn–O systems of the example and Li-ion cathode chemistry, so the example runs offline.

## Installation

```bash
git clone https://github.com/kumuditharathnayake/xrd-phase-id.git
cd xrd-phase-id
pip install -r requirements.txt
```

A free [Materials Project API key](https://materialsproject.org/api) is needed only for online Materials Project searches. Store it as an environment variable rather than in the notebook, so it is never shared by accident:

```bash
# Windows (then restart Jupyter)
setx MP_API_KEY yourkey
# macOS / Linux
export MP_API_KEY=yourkey
```

## Quick start

1. Open `XRD_phase_ID.ipynb` in Jupyter.
2. Run all cells. The preset analyses `examples/cpd-1a.txt` and writes tables, plots and reference CIFs to `results/cpd-1a/`.
3. For your own sample, set `EXPERIMENTAL_FILE`, `OUT`, `ELEMENTS`, `SYSTEMS`, `WAVELENGTH` and `DATABASES` in the settings cell.

Scan files need two numeric columns (2θ, intensity), separated by commas, tabs, spaces or semicolons; a header row is fine. Put every element of the sample in `ELEMENTS`, H included: a reference is kept only if all its elements are on the list. Restart the kernel after editing `xrd_refs.py`.

**MOF collections:** in section 0, set `IMPORT_NOW = {"CoRE_MOF": True, "QMOF": True}`, run the cell once, then set it back to `False`. The import downloads about 436 MB (32,800 CIFs) and takes several minutes; afterwards it works offline.

## Settings for different sample types

| Setting | CPD-1a (inorganic mixture) | NMC (battery cathode) | Cu-HHTP (MOF) |
|---|---|---|---|
| `ELEMENTS` | `["Al","Ca","F","Zn","O"]` | `["Li","Ni","Mn","Co","O"]` | `["Cu","C","H","O"]` |
| `SYSTEMS` | `[["Al","O"],["Ca","F"],["Zn","O"]]` | `None` | `[["Cu","C","H","O"]]` |
| `WAVELENGTH` | `"CuKa"` | `"CuKa1"` | `"CuKa"` |
| `SIM_RANGE` | `(5, 150)` | `(10, 80)` | `(2, 60)` |
| `GOOD_SCORE`, `GOOD_EXPLAINED_PCT`, `MAX_CELL_SHIFT_PCT` | 0.85, 80, 1.0 | 0.85, 80, 1.0 | 0.60, 70, 1.5 |
| `DATABASES` on | COD, MP | COD, MP | COD, CoRE_MOF, QMOF, Published |
| `USE_GENERATED_MODELS` | `False` | `True` | `False` |
| `IMPURITY_LIBRARY` | `[]` | `"Li-ion layered oxide (NMC/NCA/LCO)"` | `[["Cu","O"],["Cu"],["Cu","O","H"]]` |

Optional: `X.WORKERS` (CPU cores; default all but one, at most 12), `X.set_cache_dir(path)`, `X.clear_cache()`, `X.PEAK_WIDTH` (force the matching peak width in degrees).

## How it works

Background is removed with a rolling minimum. Each reference is simulated with pymatgen, broadened to the measured peak width and compared with the scan by cosine similarity of square-root intensities, so weak peaks still count, while its cell is stretched by up to ±1.5 %. Near-ties go to the reference needing the smallest stretch.

Phases are then found one at a time: the best match is the main phase, its peaks are masked, and every candidate is ranked again on the leftover peaks. A leftover phase is accepted only if it explains at least half of the leftover intensity and at least two leftover peaks, with no more expected-but-absent lines than peaks it explains. The loop ends when every peak is explained, nothing fits, or four phases are found. A hydroxide CIF that lists no H atoms is treated as the same phase as its H-containing entry, not as a separate oxide.

## Outputs

| File (in `OUT`) | Contents |
|---|---|
| `search_log.csv` | Each search stage, its best reference, score and why it stopped |
| `ranking_vs_experiment.csv` | Every candidate: source, kind, score, cell stretch, % explained, unexplained and missing peaks |
| `reference_summary.csv` | Cell, ordering, duplicate group and strongest peaks of the 60 best references |
| `phase_identification.csv` | One row per phase: role, confidence, peaks still unexplained |
| `lattice_parameters.csv` | Fitted cell, zero shift, number of peaks, rms |
| `phases.png`, `phases_zoom.png` | Stacked plot with peak guide lines; zoom on weak peaks with reflection ticks per phase |
| `reference_cifs/` | CIFs of the identified phases and the 10 best candidates, with `sources.csv` (links, DOIs) |

## Validation

**IUCr CPD round-robin sample 1a** (weighed 94.81 wt% fluorite, 4.04 wt% zincite, 1.15 wt% corundum): all three phases found.

| Step | Phase | Reference | Cell stretch | Fitted cell (literature) |
|---|---|---|---|---|
| 1 | CaF₂ fluorite | COD 9007060 | 0.00 % | a = 5.4625 Å (5.4626) |
| 2 | ZnO zincite | COD 2300112 | 0.00 % | a = 3.2485, c = 5.2065 Å (3.2498, 5.2066) |
| 3 | Al₂O₃ corundum (1.15 wt%) | MP mp-aaaaabrz | −0.75 % | a = 4.7616, c = 12.9856 Å (4.7589, 12.991) |

Corundum is still found in 10 out of 10 runs when every predicted line intensity is randomly changed by up to ±25 %.

**Known-answer tests on synthetic NMC111 patterns:**

| Test | Expected | Result |
|---|---|---|
| NMC111 + 8 % NiO + 6 % Ni(OH)₂ | Both side phases found | Both found, "clear" |
| Pure NMC111 | Nothing extra | Main phase only |
| Lab-type NMC111 + 4 % NiO | No phase from lines hidden under the main phase | NiO not assigned |

**Speed** (CPD-1a, 305-structure library, 2-core machine): first run 64 s → 46 s; repeat run 64 s → 8–11 s.

## Limitations

- Identifies phases and fits cells but does not yet quantify phase fractions (wt%).
- Only phases with a crystal structure in a switched-on source can be found.
- DFT cells (Materials Project, QMOF) are often 1–3 % off; experimental entries are preferred.
- A minor phase whose peaks all sit under the main phase cannot be found by the leftover method.
- Preferred orientation and stacking disorder are not modelled, which limits scores for layered materials and 2D MOFs.
- Monoclinic and triclinic cells are not fitted.

## Roadmap

- [ ] Semi-quantitative wt% (RIR), checked against the CPD-1a weighed values
- [ ] Comparison with DIFFRAC.EVA on the same scans
- [ ] Le Bail / Rietveld refinement of identified phases through GSAS-II scripting
- [ ] Preferred-orientation correction in simulated patterns
- [ ] Size/strain (Williamson–Hall) analysis from fitted peak widths

## Author

**Kumuditha Rathnayake**, Ph.D. candidate in Materials Chemistry, University of Illinois Chicago
[Website](https://kumuditharathnayake.github.io/) · [LinkedIn](https://www.linkedin.com/in/kumuditha-rathnayake-126602103) · [ORCID](https://orcid.org/0000-0002-7591-5630)
