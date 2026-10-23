# Reference library

Crystal structures (CIF) used as phase-identification references. The notebook loads this folder first and only queries the online databases when nothing here matches well; anything it downloads is added here. Each source can be switched on or off with `DATABASES` in the notebook.

**Included in this repository**

| Folder | Contents | Source and terms |
|---|---|---|
| `COD/` | `COD_<id>.cif` plus saved searches (`_cod_search_*.json`) | [Crystallography Open Database](https://www.crystallography.net/cod/), public domain (CC0). Cite: Gražulis et al., *Nucleic Acids Res.* 40, D420–D427 (2012). |
| `MP/` | `<mp-id>.cif` plus saved searches with E above hull (`_mp_search_*.json`) | [Materials Project](https://materialsproject.org), CC BY 4.0. Cite: Jain et al., *APL Mater.* 1, 011002 (2013). |
| `local/` | Your own CIFs and models | Found without any database search |

Chemical systems covered: Al–Ca–F–Zn–O (all subsystems, for the CPD example) and Li-ion cathode chemistry: Li–Ni–Mn–Co–O, with Li/Ni/Co/Mn oxides, hydroxides and carbonates and Li–Al–O and Li–S–O side phases.

**Added on your own computer (not in the repository)**

| Folder | How to get it |
|---|---|
| `CoRE_MOF/` | `X.import_collection("CoRE_MOF", LIBRARY)`: CoRE MOF 2024 from [Zenodo](https://zenodo.org/records/15055758), about 8,300 experimental MOF structures |
| `QMOF/` | `X.import_collection("QMOF", LIBRARY)`: QMOF database from [figshare](https://doi.org/10.6084/m9.figshare.13147324), about 24,500 DFT-optimised MOFs |
| `CSD/` | CIFs you download from the [CCDC](https://www.ccdc.cam.ac.uk/structures/); their licence does not allow redistribution |
| `Published/` | CIFs from papers' supporting information or data repositories |
| `_cache/` | Created automatically: parsed structures and simulated patterns, so repeat runs are fast |
