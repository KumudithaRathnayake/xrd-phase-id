# Reference library

Crystal structures (CIF) used as phase-identification references. The notebook loads this folder first and only queries the online databases when nothing here matches well; anything it downloads is added here.

| Folder | Contents | Source and terms |
|---|---|---|
| `COD/` | `COD_<id>.cif` plus saved searches (`_cod_search_*.json`) | [Crystallography Open Database](https://www.crystallography.net/cod/), public domain (CC0). Cite: Gražulis et al., *Nucleic Acids Res.* 40, D420–D427 (2012). |
| `MP/` | `<mp-id>.cif` plus saved searches with E above hull (`_mp_search_*.json`) | [Materials Project](https://materialsproject.org), CC BY 4.0. Cite: Jain et al., *APL Mater.* 1, 011002 (2013). |
| `local/` | Your own CIFs (literature, DFT) | Add your own; they are found without any database search |

Chemical systems covered: Al–Ca–F–Zn–O (all subsystems, for the CPD example) and Li-ion cathode chemistry: Li–Ni–Mn–Co–O, with Li/Ni/Co/Mn oxides, hydroxides and carbonates and Li–Al–O and Li–S–O side phases.
