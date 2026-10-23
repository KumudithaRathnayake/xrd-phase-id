"""
xrd_refs.py - find, build and simulate reference structures for XRD phase identification.

Three ways to get reference structures, all returning the same kind of record (a dict):

    search_cod(elements)          Crystallography Open Database, searched by ELEMENTS
                                  (not by exact formula, so fractional compositions such
                                  as LiNi0.33Mn0.33Co0.33O2 or LiNi0.8Mn0.1Co0.1O2 are found)
    search_mp(elements)           Materials Project (needs a free API key)
    load_local_cifs(folder)       your own CIF files
    build_layered_oxide(...)      generate an R-3m layered oxide (LiCoO2 / NMC / NCA type)
                                  for ANY composition, with optional Li/Ni cation mixing

Each record looks like:
    {"source": "COD", "ref_id": "cod-1234567", "label": "...", "formula": "...",
     "spacegroup": "R-3m", "structure": <pymatgen Structure or None>,
     "cif_path": "...", "error": None, ...extra metadata...}

Then:
    summary_table(refs)           one row per reference: lattice, symmetry, ordered/disordered,
                                  duplicate groups, strongest peaks
    simulate_peaks(structure)     stick pattern as a DataFrame (2theta, d, intensity, hkl)
    layered_oxide_metrics(s)      I(003)/I(104), (006)/(102) and (108)/(110) splittings
    plot_stack(refs)              offset plot of broadened simulated patterns

Tested with pymatgen 2026.9 (pymatgen-core). Works with or without the mp-api package.
"""
from __future__ import annotations

__version__ = "2026-10-23 (peak guide lines in plots; save_reference_cifs; measured peak width; cache + multi-core)"

import os
import re
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from pymatgen.analysis.diffraction.xrd import WAVELENGTHS, XRDCalculator
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import Lattice, Structure
from pymatgen.io.cif import CifParser, CifWriter
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

COD_SEARCH_URL = "https://www.crystallography.net/cod/result"
COD_CIF_URL = "https://www.crystallography.net/cod/{}.cif"


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
def pretty_formula(composition, disordered=True):
    """Display formula per formula unit, also for partial occupancies, e.g.
    Li3Mn0.3Co0.3Ni2.4O6 -> LiMn0.1Co0.1Ni0.8O2,  Li0.5Ni0.167Co0.167Mn0.167O -> LiNi0.333Co0.333Mn0.333O2."""
    from math import gcd

    comp = composition.element_composition.reduced_composition  # drop oxidation states (Mn3.54+ -> Mn)
    if not disordered:  # an ordered compound keeps its true reduced formula
        return comp.reduced_formula
    amounts = {str(el): v for el, v in comp.items()}
    # layered Li-TM oxides (LiMO2 family, incl. Li-rich and partly delithiated): write per O2
    if "O" in amounts and "Li" in amounts:
        cations = sum(v for el, v in amounts.items() if el != "O")
        if 0.76 <= cations / amounts["O"] <= 1.15 and 0.3 <= amounts["Li"] / cations <= 0.65:
            f = 2 / amounts["O"]
            return "".join(f"{el}{'' if abs(v * f - 1) < 0.01 else f'{v * f:.3g}'}" for el, v in amounts.items())
    total = sum(amounts.values())
    major = [el for el, v in amounts.items() if v / total > 0.2]  # e.g. Li and O in NMC
    # 1) scale up until the major elements are (nearly) whole numbers
    for n in range(1, 5):  # small multipliers only; otherwise keep the formula as given
        if all(abs(n * amounts[el] - round(n * amounts[el])) < 0.05 and round(n * amounts[el]) > 0 for el in major):
            break
    else:
        n = 1
    amounts = {el: v * n for el, v in amounts.items()}
    # 2) divide by the common factor of the major elements (Li3...O6 -> Li...O2)
    g = 0
    for el in major:
        g = gcd(g, int(round(amounts[el])))
    minor_sum = sum(v for el, v in amounts.items() if el not in major)
    # largest divisor that keeps at least one whole minor-element unit (no (Ni0.4..)O(H) nonsense)
    k = max([q for q in range(1, g + 1) if g % q == 0 and (minor_sum == 0 or minor_sum / q >= 0.99)], default=1)
    if k > 1:
        amounts = {el: v / k for el, v in amounts.items()}
    return "".join(f"{el}{'' if abs(v - 1) < 0.01 else f'{v:.3g}'}" for el, v in amounts.items())


def _record(source, ref_id, label, structure=None, cif_path=None, error=None, **meta):
    rec = {
        "source": source,
        "ref_id": ref_id,
        "label": label,
        "formula": pretty_formula(structure.composition, not structure.is_ordered) if structure is not None else meta.pop("formula", ""),
        "spacegroup": None,
        "structure": structure,
        "cif_path": str(cif_path) if cif_path else None,
        "error": error,
    }
    meta.pop("formula", None)
    sg = meta.pop("spacegroup", None)
    if structure is not None:
        rec["spacegroup"] = sg or _spacegroup(structure)
    rec.update(meta)
    return rec


def _spacegroup(structure, symprec=0.01):
    try:
        return SpacegroupAnalyzer(structure, symprec=symprec).get_space_group_symbol()
    except Exception:
        return None


def parse_cif_text(text, occupancy_tolerance=1.05):
    """Parse CIF text into a Structure (keeps partial occupancies)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        parser = CifParser.from_str(text, occupancy_tolerance=occupancy_tolerance)
        return parser.parse_structures(primitive=False)[0]


def write_cif(structure, path, symprec=0.01):
    """Write a CIF, with symmetry if it can be found, otherwise in P1."""
    path = Path(path)
    try:
        CifWriter(structure, symprec=symprec).write_file(path)
    except Exception:
        CifWriter(structure).write_file(path)
    return path


# ---------------------------------------------------------------------------
# Source 1: Crystallography Open Database (experimental structures)
# ---------------------------------------------------------------------------
def _http_get(url, params=None, timeout=60, attempts=5, wait=3.0):
    """GET with retries. COD in particular drops connections when it gets several requests
    in a row (ConnectionResetError 10054 on Windows); waiting a few seconds and retrying fixes it."""
    last = None
    for i in range(attempts):
        try:
            r = requests.get(url, params=params, timeout=timeout)
            if getattr(r, "status_code", 200) in (429, 500, 502, 503, 504):
                raise requests.HTTPError(f"HTTP {r.status_code}", response=r)
            r.raise_for_status()
            return r
        except (requests.ConnectionError, requests.Timeout, requests.HTTPError) as ex:
            if isinstance(ex, requests.HTTPError) and ex.response is not None and ex.response.status_code < 429:
                raise  # 4xx: a real error, retrying will not help
            last = ex
            if i < attempts - 1:
                time.sleep(wait * (i + 1))
    raise last


def _cod_saved_path(elements, out_dir, exact=True, exclude=()):
    tag = "-".join(elements) + ("_exact" if exact else "_atleast") + ("_no" + "-".join(exclude) if exclude else "")
    return Path(out_dir) / f"_cod_search_{tag}.json"


def search_cod(elements, out_dir="refs_cod", exact=True, exclude=(), max_entries=200, timeout=60,
               refresh=False, pause=0.3):
    """Search COD by elements and download every hit as a CIF.

    elements : e.g. ["Li", "Ni", "Mn", "Co", "O"]
    exact    : True  -> only entries containing exactly these elements
               False -> entries containing at least these elements
    exclude  : elements that must NOT be present, e.g. ["H", "C"] (max 4, COD limit)

    refresh  : False -> reuse the saved search result and CIFs from an earlier run (no network
               needed); True -> ask COD again.
    pause    : seconds between CIF downloads (be gentle - COD drops rapid connections)

    pymatgen's built-in COD().get_structure_by_formula() needs the exact formula,
    which almost never matches for solid solutions like NMC - hence the element search.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    params = {f"el{i + 1}": el for i, el in enumerate(elements)}
    params.update({f"nel{i + 1}": el for i, el in enumerate(exclude)})
    if exact:
        params["strictmin"] = len(elements)
        params["strictmax"] = len(elements)
    params["format"] = "json"

    import json

    saved = _cod_saved_path(elements, out_dir, exact, exclude)
    if saved.exists() and not refresh:
        entries = json.loads(saved.read_text(encoding="utf-8"))
        print(f"  (using saved COD search {saved.name}; refresh=True to query COD again)")
    else:
        entries = _http_get(COD_SEARCH_URL, params=params, timeout=timeout).json()
        saved.write_text(json.dumps(entries), encoding="utf-8")
    if len(entries) > max_entries:
        print(f"COD returned {len(entries)} entries; keeping the first {max_entries} (raise max_entries to get all).")
        entries = entries[:max_entries]

    refs = []
    for e in entries:
        cod_id = str(e["file"])
        cif_path = out_dir / f"COD_{cod_id}.cif"
        meta = {
            "formula": str(e.get("formula", "")).strip("- ").strip(),
            "title": e.get("title"),
            "year": e.get("year"),
            "doi": e.get("doi"),
            "cod_spacegroup": e.get("sg"),
        }
        try:
            if not cif_path.exists():
                r = _http_get(COD_CIF_URL.format(cod_id), timeout=timeout)
                cif_path.write_text(r.text, encoding="utf-8")
                time.sleep(pause)
            s = parse_cif_text(cif_path.read_text(encoding="utf-8"))
            refs.append(_record("COD", f"cod-{cod_id}", f"COD {cod_id}", s, cif_path, **meta))
        except Exception as ex:  # keep going; report the failure in the table
            refs.append(_record("COD", f"cod-{cod_id}", f"COD {cod_id}", None, cif_path, error=str(ex)[:200], **meta))
    return refs


# ---------------------------------------------------------------------------
# Source 2: Materials Project (computed structures)
# ---------------------------------------------------------------------------
def search_mp(elements, out_dir="refs_mp", api_key=None, max_e_above_hull=0.10, conventional=True,
              refresh=False):
    """Search Materials Project for the exact chemical system and save CIFs.

    api_key : your key, or set the environment variable MP_API_KEY (or PMG_MAPI_KEY).
    max_e_above_hull : eV/atom; 0 = only predicted-stable phases.

    Notes:
    - MP structures are DFT-relaxed: lattice parameters are typically off by ~1-2 %, so peaks
      sit slightly shifted from experiment. Use them for identification, not for exact positions.
    - MP structures are ORDERED; a solid solution such as NMC appears as an ordered supercell
      in a lower-symmetry space group. For Ni/Mn/Co the extra superstructure peaks are very weak
      (similar scattering power), but for Li/TM ordering (e.g. Li2MnO3) they are clearly visible.
    - refresh=False reuses the result of an earlier identical search saved in out_dir (no key or
      internet needed); refresh=True asks Materials Project again.
    """
    import json

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    chemsys = "-".join(sorted(elements))
    saved = out_dir / f"_mp_search_{chemsys}.json"
    if saved.exists() and not refresh:
        info = json.loads(saved.read_text(encoding="utf-8"))
        if info.get("max_e_above_hull", 0) >= max_e_above_hull:
            refs = []
            for row in info["entries"]:
                if row["e_above_hull"] is not None and row["e_above_hull"] > max_e_above_hull:
                    continue
                cif = out_dir / f"{row['material_id']}.cif"
                if cif.exists():
                    st = parse_cif_text(cif.read_text(encoding="utf-8"))
                    refs.append(_record("MP", row["material_id"], f"{row['material_id']} {row['formula_pretty']}", st,
                                        cif, e_above_hull=row["e_above_hull"], theoretical=row.get("theoretical")))
            print(f"  (using saved MP search {saved.name}; refresh=True to query Materials Project again)")
            return refs

    api_key = api_key or os.environ.get("MP_API_KEY") or os.environ.get("PMG_MAPI_KEY")
    if not api_key:
        raise ValueError("No Materials Project API key. Pass api_key=... or set MP_API_KEY.")
    fields = ["material_id", "formula_pretty", "structure", "energy_above_hull", "symmetry", "theoretical"]

    try:  # full client, if installed
        from mp_api.client import MPRester

        with MPRester(api_key) as mpr:
            docs = mpr.materials.summary.search(chemsys=chemsys, fields=fields)
        rows = [{f: getattr(d, f, None) for f in fields} for d in docs]
    except ImportError:  # lightweight client that ships with pymatgen
        from pymatgen.ext.matproj import MPRester

        with MPRester(api_key) as mpr:
            rows = mpr.materials.summary.search(chemsys=chemsys, _fields=fields)

    refs = []
    for row in rows:
        ehull = row.get("energy_above_hull")
        if ehull is not None and ehull > max_e_above_hull:
            continue
        s = row["structure"]
        if isinstance(s, dict):
            s = Structure.from_dict(s)
        if conventional:
            try:
                s = SpacegroupAnalyzer(s, symprec=0.01).get_conventional_standard_structure()
            except Exception:
                pass
        mpid = str(row["material_id"])
        cif_path = write_cif(s, out_dir / f"{mpid}.cif")
        refs.append(
            _record(
                "MP", mpid, f"{mpid} {row.get('formula_pretty', '')}", s, cif_path,
                e_above_hull=ehull, theoretical=row.get("theoretical"),
            )
        )
    refs.sort(key=lambda r: (r.get("e_above_hull") is None, r.get("e_above_hull")))
    saved.write_text(json.dumps({
        "max_e_above_hull": max_e_above_hull,
        "entries": [{"material_id": r["ref_id"], "formula_pretty": r["label"].split(" ", 1)[-1],
                     "e_above_hull": r.get("e_above_hull"), "theoretical": r.get("theoretical")} for r in refs],
    }), encoding="utf-8")
    return refs


# ---------------------------------------------------------------------------
# Source 3: your own CIFs
# ---------------------------------------------------------------------------
def load_local_cifs(folder, pattern="*.cif"):
    refs = []
    for p in sorted(Path(folder).glob(pattern)):
        try:
            s = parse_cif_text(p.read_text(encoding="utf-8", errors="replace"))
            refs.append(_record("local", p.stem, p.stem, s, p))
        except Exception as ex:
            refs.append(_record("local", p.stem, p.stem, None, p, error=str(ex)[:200]))
    return refs


# ---------------------------------------------------------------------------
# Source 4: generate layered oxides (LiCoO2 / NMC / NCA, R-3m, alpha-NaFeO2 type)
# ---------------------------------------------------------------------------
# Approximate literature lattice parameters (Angstrom) for pristine materials.
# Starting values only - refine a and c against your own data.
LAYERED_PRESETS = {
    "LiCoO2": ({"Co": 1.0}, 2.816, 14.05),
    "NMC111": ({"Ni": 1 / 3, "Mn": 1 / 3, "Co": 1 / 3}, 2.860, 14.22),
    "NMC532": ({"Ni": 0.5, "Mn": 0.3, "Co": 0.2}, 2.868, 14.23),
    "NMC622": ({"Ni": 0.6, "Mn": 0.2, "Co": 0.2}, 2.869, 14.22),
    "NMC811": ({"Ni": 0.8, "Mn": 0.1, "Co": 0.1}, 2.872, 14.20),
    "LiNiO2": ({"Ni": 1.0}, 2.878, 14.19),
}


def build_layered_oxide(tm, a, c, li_ni_mixing=0.0, z_O=0.2416, label=None):
    """Disordered R-3m layered oxide Li(TM)O2 in the hexagonal setting.

    tm            : transition-metal fractions on the TM site, e.g. {"Ni": .8, "Mn": .1, "Co": .1}
    a, c          : hexagonal lattice parameters (Angstrom)
    li_ni_mixing  : fraction x of Ni on the Li site (and x Li on the TM site), e.g. 0.05
    z_O           : oxygen z (6c). 0.2416 gives TM-O ~1.97 A, Li-O ~2.10 A for NMC111.

    Sites: Li 3a (0,0,0), TM 3b (0,0,1/2), O 6c (0,0,z).
    """
    tm = {k: float(v) for k, v in tm.items()}
    if abs(sum(tm.values()) - 1) > 1e-3:
        raise ValueError(f"TM fractions must sum to 1, got {sum(tm.values()):.3f}")
    x = float(li_ni_mixing)
    if x > 0 and tm.get("Ni", 0) < x:
        raise ValueError("li_ni_mixing cannot exceed the Ni fraction.")

    li_site = {"Li": 1 - x}
    tm_site = dict(tm)
    if x > 0:
        li_site["Ni"] = x
        tm_site["Ni"] = tm_site["Ni"] - x
        tm_site["Li"] = x
    li_site = {k: v for k, v in li_site.items() if v > 1e-9}
    tm_site = {k: v for k, v in tm_site.items() if v > 1e-9}

    s = Structure.from_spacegroup(
        "R-3m", Lattice.hexagonal(a, c), [li_site, tm_site, "O"], [[0, 0, 0], [0, 0, 0.5], [0, 0, z_O]]
    )
    if label is None:
        label = "Li(" + "".join(f"{el}{v:.2g}" for el, v in tm.items()) + ")O2"
        if x:
            label += f" ({x:.0%} Li/Ni mixing)"
    return _record("built", label.replace(" ", "_"), label, s, None, a=a, c=c, li_ni_mixing=x)


def build_layered_preset(name, li_ni_mixing=0.0, a=None, c=None):
    tm, a0, c0 = LAYERED_PRESETS[name]
    label = name + (f" ({li_ni_mixing:.0%} Li/Ni mixing)" if li_ni_mixing else "")
    return build_layered_oxide(tm, a or a0, c or c0, li_ni_mixing=li_ni_mixing, label=label)


def ordered_sqrt3_model(name="NMC111"):
    """An ORDERED sqrt3 x sqrt3 Ni/Mn/Co arrangement of NMC111 - the kind of model a computed
    database stores (lower symmetry, 3x larger cell). Compare it with the disordered R-3m model:
    the superstructure peaks near 20-24 deg (Cu Ka) are only ~0.1-0.2 % because Ni, Mn and Co
    scatter almost equally, so XRD alone cannot tell ordered from disordered NMC."""
    tm, a, c = LAYERED_PRESETS[name]
    s = build_layered_oxide({"Co": 1.0}, a, c)["structure"].copy()
    s.make_supercell([[2, 1, 0], [-1, 1, 0], [0, 0, 1]])
    tm_idx = [i for i, site in enumerate(s) if site.specie.symbol == "Co"]
    els = [el for el, frac in tm.items() for _ in range(round(frac * 3))]
    layers = {}
    for i in tm_idx:
        layers.setdefault(round(s[i].frac_coords[2], 3), []).append(i)
    for n, (_, idx) in enumerate(sorted(layers.items())):
        idx.sort(key=lambda i: tuple(np.round(s[i].frac_coords[:2] % 1, 3)))
        for j, i in enumerate(idx):
            s.replace(i, els[(j + n) % 3])
    return _record("built", f"{name}_ordered", f"{name} ordered (computed-style)", s, None)


def build_tm_hydroxide(tm=None, a=3.126, c=4.605, label=None):
    """beta-Ni(OH)2 / brucite-type M(OH)2, P-3m1 - the usual NMC precursor phase.

    tm : metal fractions, e.g. {"Ni": .8, "Mn": .1, "Co": .1} (default pure Ni).
    Default a, c are for beta-Ni(OH)2; Co(OH)2 ~ (3.18, 4.65), Mn(OH)2 ~ (3.32, 4.73).
    Real precursors are often turbostratic: broad, asymmetric hk0 bands - the simulated
    peaks will be much sharper than the measured ones.
    """
    tm = tm or {"Ni": 1.0}
    s = Structure.from_spacegroup(
        "P-3m1", Lattice.hexagonal(a, c), [dict(tm), "O", "H"],
        [[0, 0, 0], [1 / 3, 2 / 3, 0.2208], [1 / 3, 2 / 3, 0.4300]],
    )
    label = label or "M(OH)2 " + "".join(f"{el}{v:.2g}" for el, v in tm.items())
    return _record("built", label.replace(" ", "_"), label, s, None, a=a, c=c)


# ---------------------------------------------------------------------------
# Simulation and comparison
# ---------------------------------------------------------------------------
def _hkl3(hkl):
    """pymatgen gives 4-index (h,k,i,l) for hexagonal cells; reduce to (h,k,l)."""
    hkl = tuple(int(v) for v in hkl)
    return (hkl[0], hkl[1], hkl[3]) if len(hkl) == 4 else hkl


def simulate_peaks(structure, two_theta=(10, 80), wavelength="CuKa"):
    """Stick pattern as a DataFrame: two_theta, d, intensity (max=100), hkl (list of (h,k,l))."""
    calc = XRDCalculator(wavelength=wavelength)
    p = calc.get_pattern(structure, two_theta_range=two_theta)
    return pd.DataFrame(
        {
            "two_theta": p.x,
            "d": p.d_hkls,
            "intensity": p.y,
            "hkl": [[_hkl3(h["hkl"]) for h in hk] for hk in p.hkls],
        }
    )


def _hex_key(hkl):
    """Family key for hexagonal (h,k,l): in-plane (h,k,i=-h-k) up to permutation/sign, and |l|.
    So (1,1,0) and (2,-1,0) - both written by different programs for the same 110 peak - match."""
    h, k, l = hkl
    return tuple(sorted((abs(h), abs(k), abs(h + k)))), abs(l)


def _find(peaks, hkl):
    """Row of the peak whose hkl family contains `hkl` (hexagonal indexing), or None."""
    target = _hex_key(hkl)
    for _, row in peaks.iterrows():
        if any(_hex_key(h) == target for h in row["hkl"]):
            return row
    return None


def layered_oxide_metrics(structure, wavelength="CuKa", peaks=None):
    """Classic quality indicators for R-3m layered oxides (hexagonal indexing).

    I003_I104   : higher = better-ordered layers / less Li/Ni mixing (typically > 1.2 is 'good')
    split_006_102, split_108_110 : 2theta separation (deg); clear splitting = well-developed layering
    peaks : an already simulated stick pattern to reuse (simulated here if missing reflections)
    """
    need = [(0, 0, 3), (1, 0, 4), (0, 0, 6), (1, 0, 2), (1, 0, 8), (1, 1, 0)]
    rows = {h: _find(peaks, h) for h in need} if peaks is not None else {}
    if peaks is None or any(v is None for v in rows.values()):
        pk = simulate_peaks(structure, (10, 80), wavelength)
        rows = {h: _find(pk, h) for h in need}
    out = {}
    if rows[(0, 0, 3)] is not None and rows[(1, 0, 4)] is not None:
        out["I003_I104"] = rows[(0, 0, 3)]["intensity"] / rows[(1, 0, 4)]["intensity"]
        out["2theta_003"] = rows[(0, 0, 3)]["two_theta"]
        out["2theta_104"] = rows[(1, 0, 4)]["two_theta"]
    for a, b in [((0, 0, 6), (1, 0, 2)), ((1, 0, 8), (1, 1, 0))]:
        if rows[a] is not None and rows[b] is not None:
            out[f"split_{''.join(map(str, a))}_{''.join(map(str, b))}"] = abs(rows[a]["two_theta"] - rows[b]["two_theta"])
    return out


def _peak_label(hkls, hexagonal=False):
    """'003', '104', '110' ... For hexagonal cells, use the usual positive form (2,-1,0) -> 110."""
    out = []
    for h in hkls[:2]:
        if hexagonal and len(h) == 3:
            (a, b, _), l = _hex_key(h)
            h = (b, a, l)
        out.append("".join(str(v) for v in h))
    return "/".join(out)


def _any_cached_peaks(ref, wavelength):
    """Widest stick pattern already simulated for this reference at this wavelength, or None."""
    best = None
    for (lo, hi, wl), pk in ref.get("_peaks_cache", {}).items():
        if wl == str(wavelength) and (best is None or hi - lo > best[0]):
            best = (hi - lo, pk)
    return None if best is None else best[1]


def summary_table(refs, two_theta=(10, 80), wavelength="CuKa", group_duplicates=True, n_peaks=3,
                  reuse_simulations=True):
    """One row per reference. Duplicates (same structure within tolerance) share a 'group'.

    reuse_simulations : use patterns already simulated (e.g. by the ranking) instead of simulating
                        again - much faster; the strongest peaks are then those inside the range that
                        was simulated (your scan range).
    """
    rows = []
    good = [r for r in refs if r["structure"] is not None]

    groups = {}
    if group_duplicates and good:
        # duplicates must have the same composition and symmetry: compare only within those buckets
        matcher = StructureMatcher(ltol=0.2, stol=0.3, angle_tol=5, primitive_cell=True, scale=False)
        buckets = {}
        for r in good:
            key = (r["structure"].composition.element_composition.reduced_formula, r.get("spacegroup"))
            buckets.setdefault(key, []).append(r)
        gi = 0
        for members in buckets.values():
            if len(members) == 1:
                groups[members[0]["ref_id"]] = gi
                gi += 1
                continue
            try:
                grouped = matcher.group_structures([m["structure"] for m in members])
            except Exception:
                grouped = [[m["structure"]] for m in members]
            for grp in grouped:
                for st in grp:
                    for m in members:
                        if m["structure"] is st:
                            groups[m["ref_id"]] = gi
                gi += 1

    for r in refs:
        s = r["structure"]
        row = {k: r.get(k) for k in ("source", "ref_id", "label", "formula", "spacegroup")}
        row["kind"] = r.get("kind") or {"COD": "experimental", "MP": "DFT", "built": "model"}.get(r.get("source"), "")
        for k in ("e_above_hull", "year", "title", "li_ni_mixing", "error"):
            if r.get(k) is not None:
                row[k] = r[k]
        if s is not None:
            lat = s.lattice
            row.update(
                a=round(lat.a, 4), b=round(lat.b, 4), c=round(lat.c, 4),
                alpha=round(lat.alpha, 2), beta=round(lat.beta, 2), gamma=round(lat.gamma, 2),
                vol_per_atom=round(lat.volume / s.composition.num_atoms, 3),
                ordered=s.is_ordered, n_sites=len(s), group=groups.get(r["ref_id"]),
            )
            try:
                pk_all = _any_cached_peaks(r, wavelength) if reuse_simulations else None
                if pk_all is None:
                    pk_all = _cached_peaks(r, two_theta, wavelength)
                pk = pk_all.nlargest(n_peaks, "intensity").sort_values("two_theta")
                hexa = s.lattice.is_hexagonal()
                row["strongest_peaks"] = ", ".join(
                    f"{t:.2f}({_peak_label(h, hexa)})" for t, h in zip(pk["two_theta"], pk["hkl"])
                )
            except Exception as ex:
                row["error"] = f"XRD sim failed: {ex}"[:200]
            if r.get("spacegroup") == "R-3m":  # hexagonal-setting indices only mean 003/104 here
                pk_all = _any_cached_peaks(r, wavelength)
                row.update({k: round(v, 3) for k, v in layered_oxide_metrics(s, wavelength, pk_all).items()})
        rows.append(row)
    df = pd.DataFrame(rows)
    if "group" in df:
        df["group"] = df["group"].astype("Int64")
    return df


def pseudo_voigt_profile(grid, positions, intensities, fwhm=0.15, eta=0.5):
    """Broaden sticks with a pseudo-Voigt (eta=0 Gaussian, 1 Lorentzian). Normalised to max 100."""
    grid = np.asarray(grid, float)
    y = np.zeros_like(grid)
    sig = fwhm / (2 * np.sqrt(2 * np.log(2)))
    gam = fwhm / 2
    for x0, amp in zip(positions, intensities):
        dx = grid - x0
        g = np.exp(-0.5 * (dx / sig) ** 2)
        lor = 1 / (1 + (dx / gam) ** 2)
        y += amp * (eta * lor + (1 - eta) * g)
    return y / y.max() * 100 if y.max() > 0 else y


def plot_stack(refs, two_theta=(10, 80), wavelength="CuKa", fwhm=0.15, eta=0.5, offset=135,
               label_hkl_above=15, ax=None, experimental=None, peak_lines=True, peak_min_rel=0.0,
               label_peaks_above=3.0, unexplained=None, tol_deg=None):
    """Offset plot of broadened simulated patterns (first ref at the bottom).

    experimental : optional (two_theta_array, intensity_array) plotted at the very bottom.
    label_hkl_above : label sticks stronger than this (relative %) with their hkl.
    peak_lines   : dashed vertical line through every measured peak (needs experimental), from the
                   scan up through all reference patterns, so you can see which reference line
                   sits under which measured peak. peak_min_rel: only peaks >= this % of the
                   strongest; label_peaks_above: write the 2theta of peaks >= this % at the top.
    unexplained  : 2theta values of peaks no identified phase explains (e.g. from
                   unexplained_peaks); their lines are drawn red and always labelled.
    """
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots(figsize=(11, 1.6 + 1.3 * (len(refs) + (experimental is not None))))
    grid = np.arange(two_theta[0], two_theta[1], 0.01)
    k = 0
    if experimental is not None:
        x, y = map(np.asarray, experimental)
        y = (y - y.min()) / (y.max() - y.min()) * 100
        ax.plot(x, y, color="black", lw=1.2)
        ax.text(two_theta[1], 0.45 * offset, "experimental  ", ha="right", fontsize=9)
        k = 1
    for r in refs:
        if r["structure"] is None:
            continue
        lines = _emission_lines(wavelength)            # Ka1 + Ka2 for a lab doublet ('CuKa')
        pk = simulate_peaks(r["structure"], two_theta, lines[0][0] if len(lines) > 1 else wavelength)
        prof = sum(w * pseudo_voigt_profile(grid, p_["two_theta"], p_["intensity"], fwhm, eta) / 100
                   for (lam, w), p_ in ((ln, pk if i == 0 else simulate_peaks(r["structure"], two_theta, ln[0]))
                                       for i, ln in enumerate(lines)))
        prof = prof / prof.max() * 100 if prof.max() > 0 else prof
        base = k * offset
        line, = ax.plot(grid, prof + base, lw=1.2)
        text = r.get("plot_label") or f"{r['label']}  [{r['spacegroup']}]"
        ax.text(two_theta[1], base + 0.45 * offset, text + "  ", ha="right",
                fontsize=9, color=line.get_color())
        hexa = r["structure"].lattice.is_hexagonal()
        labelled = []
        for _, p in pk[pk["intensity"] >= label_hkl_above].sort_values("intensity", ascending=False).iterrows():
            if any(abs(p["two_theta"] - t) < 0.6 for t in labelled):
                continue                      # overlapping reflections (low symmetry): label the strongest only
            labelled.append(p["two_theta"])
            ax.text(p["two_theta"], base + min(p["intensity"], 100) + 3, _peak_label(p["hkl"], hexa),
                    rotation=90, ha="center", va="bottom", fontsize=7, color=line.get_color())
        k += 1
    ax.set_xlim(*two_theta)
    top = (k - 1) * offset + 125
    if experimental is not None and peak_lines:
        xe, ye = np.asarray(experimental[0], float), np.asarray(experimental[1], float)
        pk_e = find_experimental_peaks(xe, subtract_background(xe, ye)[0])
        pk_e = pk_e[(pk_e.rel_intensity >= peak_min_rel) & (pk_e.two_theta >= two_theta[0])
                    & (pk_e.two_theta <= two_theta[1])]
        bad = np.asarray(list(unexplained) if unexplained is not None else [], float)
        tol = tol_deg if tol_deg is not None else max(0.05, 0.5 * float(np.median(pk_e.fwhm)) if len(pk_e) else 0.1)
        band = 0.32 * offset                       # a strip above the top pattern for the 2theta labels
        top += band
        min_sep = 0.011 * (two_theta[1] - two_theta[0])
        labelled = []
        order = sorted(zip(pk_e.two_theta, pk_e.rel_intensity),
                       key=lambda p_: (not (len(bad) and np.min(np.abs(bad - p_[0])) <= tol), -p_[1]))
        for t_, h_ in order:                       # unexplained first, then strongest first
            red = bool(len(bad) and np.min(np.abs(bad - t_)) <= tol)
            ax.axvline(t_, ls=(0, (4, 3)), lw=0.9 if red else 0.6, color="red" if red else "0.45",
                       alpha=0.9 if red else 0.55, zorder=0)
            if (red or h_ >= label_peaks_above) and all(abs(t_ - s) >= min_sep for s in labelled):
                labelled.append(t_)
                ax.text(t_, top - 3, f"{t_:.2f}", rotation=90, ha="center", va="top", fontsize=6.5,
                        color="red" if red else "0.3", zorder=3,
                        bbox=dict(boxstyle="square,pad=0.1", fc="white", ec="none", alpha=0.8))
        if len(bad):
            ax.plot([], [], ls=(0, (4, 3)), color="red", label="measured peak not explained")
            ax.plot([], [], ls=(0, (4, 3)), color="0.45", label="measured peak")
            ax.legend(loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=2, fontsize=7, frameon=False,
                      borderaxespad=0.2)
    ax.set_ylim(-5, top)
    ax.set_yticks([])
    ax.set_xlabel(r"2$\theta$ / degrees" + (f"  ({wavelength})" if isinstance(wavelength, str) else ""))
    ax.set_ylabel("Intensity + offset (a.u.)")
    return ax


# ---------------------------------------------------------------------------
# Experimental data and automatic ranking
# ---------------------------------------------------------------------------
def load_xy(path, x_col=0, y_col=1):
    """Read a powder scan (.txt/.xy/.csv; comma, tab, space or semicolon; header or not).
    x_col / y_col count only the numeric columns."""
    df = pd.read_csv(path, sep=r"[\s,;]+", engine="python", header=None, comment="#", skip_blank_lines=True)
    df = df.apply(pd.to_numeric, errors="coerce").dropna(axis=1, how="all")
    df = df.dropna(subset=[df.columns[x_col], df.columns[y_col]])
    x = df.iloc[:, x_col].to_numpy(float)
    y = df.iloc[:, y_col].to_numpy(float)
    order = np.argsort(x)
    return x[order], y[order]


def _wavelength_value(wavelength):
    return WAVELENGTHS[wavelength] if isinstance(wavelength, str) else float(wavelength)


def _emission_lines(wavelength, ka2_ratio=0.5):
    """[(lambda, relative weight), ...] actually present in the data.
    'CuKa' (also Co, Mo, Fe, Cr, Ag) = the Ka1/Ka2 doublet of a lab tube without a Ka1
    monochromator, so each reflection is matched as two peaks (Ka2 at half intensity).
    'CuKa1' or a number = a single line."""
    if isinstance(wavelength, str) and wavelength.endswith("Ka") and f"{wavelength}1" in WAVELENGTHS:
        return [(WAVELENGTHS[f"{wavelength}1"], 1.0), (WAVELENGTHS[f"{wavelength}2"], ka2_ratio)]
    return [(_wavelength_value(wavelength), 1.0)]


def subtract_background(x, y, window_deg=2.0):
    """Rolling-minimum background (fine for crystalline samples; amorphous humps wider than
    `window_deg` are treated as background). Returns (y - background >= 0, background)."""
    from scipy.ndimage import minimum_filter1d, uniform_filter1d

    step = np.median(np.diff(x))
    n = max(3, int(window_deg / step)) | 1
    bg = uniform_filter1d(minimum_filter1d(y, n, mode="nearest"), n, mode="nearest")
    bg = np.minimum(bg, y)
    return y - bg, bg


def find_experimental_peaks(x, y_bg, min_rel=None, snr=6.0, floor_rel=0.001, smooth_pts=2.0):
    """Peaks in a background-subtracted scan: 2theta, relative intensity (%), FWHM (deg).

    Default (min_rel=None): a peak must stand `snr` x the noise above its surroundings
    (noise estimated from the data), and at least `floor_rel` of the strongest peak. This finds
    weak broad humps of minor phases in good data without picking up noise in poor data.
    min_rel=0.02 instead uses a fixed 2 %-of-maximum threshold.
    """
    from scipy.ndimage import gaussian_filter1d
    from scipy.signal import find_peaks, peak_widths

    step = np.median(np.diff(x))
    ys = gaussian_filter1d(y_bg, smooth_pts)
    if min_rel is None:
        hf = ys - gaussian_filter1d(ys, 10 * smooth_pts)
        noise = 1.4826 * np.median(np.abs(hf - np.median(hf)))
        threshold = max(snr * noise, floor_rel * ys.max())
    else:
        threshold = min_rel * ys.max()
    idx, _ = find_peaks(ys, prominence=threshold, width=1)
    fwhm = peak_widths(ys, idx, 0.5)[0] * step
    out = pd.DataFrame({"two_theta": x[idx], "rel_intensity": ys[idx] / ys.max() * 100, "fwhm": fwhm})
    out.attrs["detection_limit_rel"] = 100 * threshold / ys.max() if ys.max() > 0 else 0.0
    return out


# ---------------------------------------------------------------------------
# Speed: a disk cache in the library + several CPU cores
# ---------------------------------------------------------------------------
# CACHE_DIR : where parsed structures and simulated patterns are kept between sessions.
#             None = set automatically to <library>/_cache the first time a library is loaded;
#             False = no disk cache.
# WORKERS   : CPU cores used for simulating and scoring. None = all cores but one (at most 12);
#             1 = everything on one core (as before).
CACHE_DIR = None
WORKERS = None
PARALLEL_MIN = 12           # fewer jobs than this run on one core (starting workers costs ~2-5 s)
_CACHE_VERSION = 1
_POOL = None
_POOL_N = 0
_PARALLEL_BROKEN = False


def set_cache_dir(path):
    """Use this folder for the disk cache (False / None to switch it off for this session)."""
    global CACHE_DIR
    CACHE_DIR = Path(path) if path else False


def _cache_root():
    if not CACHE_DIR:
        return None
    root = Path(CACHE_DIR)
    try:
        root.mkdir(parents=True, exist_ok=True)
    except OSError:
        return None
    return root


def _auto_cache_dir(library_dir):
    """First library loaded in a session decides the cache folder, unless the user chose one."""
    global CACHE_DIR
    if CACHE_DIR is None:
        CACHE_DIR = Path(library_dir) / "_cache"


def clear_cache(what="all"):
    """Delete the disk cache ('peaks', 'structures' or 'all'); it is rebuilt as needed."""
    import shutil

    root = _cache_root()
    if root is None:
        return
    for sub in (("peaks", "structures") if what == "all" else (what,)):
        shutil.rmtree(root / sub, ignore_errors=True)


def _cache_read(path):
    import pickle

    try:
        with open(path, "rb") as fh:
            obj = pickle.load(fh)
        return obj if isinstance(obj, dict) and obj.get("v") == _CACHE_VERSION else None
    except Exception:
        return None


def _cache_write(path, obj):
    """Atomic write (a half-written file is never read back)."""
    import pickle

    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_name(path.name + f".{os.getpid()}.tmp")
        with open(tmp, "wb") as fh:
            pickle.dump({"v": _CACHE_VERSION, **obj}, fh, protocol=pickle.HIGHEST_PROTOCOL)
        os.replace(tmp, path)
    except Exception:
        pass


def _structure_key(structure):
    """Fingerprint of a structure (lattice, species/occupancies, positions) for the pattern cache."""
    import hashlib

    h = hashlib.md5()
    h.update(np.round(structure.lattice.matrix, 5).tobytes())
    h.update(np.round(np.mod(structure.frac_coords, 1.0), 5).tobytes())
    h.update("|".join(site.species.formula for site in structure).encode())
    return h.hexdigest()


def load_cif_cached(path):
    """(structure, spacegroup) for a CIF file - parsed once, then read back from the disk cache
    (keyed by the file's content, so a copied or moved file is still a cache hit)."""
    import hashlib

    path = Path(path)
    data = path.read_bytes()
    root = _cache_root()
    key = hashlib.md5(data).hexdigest()
    if root is not None:
        hit = _cache_read(root / "structures" / f"{key}.pkl")
        if hit is not None:
            return hit["structure"], hit["spacegroup"]
    st = parse_cif_text(data.decode("utf-8", errors="replace"))
    sg = _spacegroup(st)
    if root is not None:
        _cache_write(root / "structures" / f"{key}.pkl", {"structure": st, "spacegroup": sg})
    return st, sg


def _n_workers():
    if WORKERS is not None:
        return max(1, int(WORKERS))
    return max(1, min(12, (os.cpu_count() or 2) - 1))


def _pool(n):
    """One pool of worker processes for the session (Windows/Jupyter safe: the work functions
    live in this module, which the workers import)."""
    global _POOL, _POOL_N
    from concurrent.futures import ProcessPoolExecutor

    if _POOL is None or _POOL_N != n:
        if _POOL is not None:
            _POOL.shutdown(wait=False, cancel_futures=True)
        _POOL, _POOL_N = ProcessPoolExecutor(max_workers=n), n
    return _POOL


def shutdown_workers():
    """Stop the worker processes (they also stop when the kernel restarts)."""
    global _POOL, _POOL_N
    if _POOL is not None:
        _POOL.shutdown(wait=False, cancel_futures=True)
    _POOL, _POOL_N = None, 0


def _parallel_map(func, jobs, min_jobs=None):
    """func over jobs on several cores; falls back to one core if workers cannot start."""
    global _PARALLEL_BROKEN
    n = _n_workers()
    if n <= 1 or len(jobs) < (PARALLEL_MIN if min_jobs is None else min_jobs) or _PARALLEL_BROKEN:
        return [func(j) for j in jobs]
    try:
        chunk = max(1, len(jobs) // (n * 4))
        return list(_pool(n).map(func, jobs, chunksize=chunk))
    except Exception as ex:                       # e.g. a blocked process start - keep working on one core
        _PARALLEL_BROKEN = True
        shutdown_workers()
        print(f"  (several cores not available here - continuing on one core: {str(ex)[:80]})")
        return [func(j) for j in jobs]


def _simulate_job(job):
    structure, lo, hi, wavelength = job
    try:
        return simulate_peaks(structure, (lo, hi), wavelength)
    except Exception as ex:
        return ex


def _peaks_disk_path(ref, wavelength):
    root = _cache_root()
    if root is None or ref.get("structure") is None:
        return None
    if "_key" not in ref:
        ref["_key"] = _structure_key(ref["structure"])
    wl = re.sub(r"[^A-Za-z0-9.]", "_", str(wavelength))
    return root / "peaks" / f"{ref['_key']}_{wl}.pkl"


def _memory_peaks(ref, lo, hi, wavelength):
    for (clo, chi, wl), pk in ref.get("_peaks_cache", {}).items():
        if wl == str(wavelength) and clo <= lo + 1e-9 and chi >= hi - 1e-9:
            return pk[(pk.two_theta >= lo) & (pk.two_theta <= hi)]
    return None


def _disk_peaks(ref, lo, hi, wavelength):
    p = _peaks_disk_path(ref, wavelength)
    if p is None or not p.exists():
        return None
    hit = _cache_read(p)
    if hit is None or not (hit["lo"] <= lo + 1e-9 and hit["hi"] >= hi - 1e-9):
        return None
    ref.setdefault("_peaks_cache", {})[(hit["lo"], hit["hi"], str(wavelength))] = hit["df"]
    return hit["df"][(hit["df"].two_theta >= lo) & (hit["df"].two_theta <= hi)]


def _store_peaks(ref, lo, hi, wavelength, pk):
    ref.setdefault("_peaks_cache", {})[(lo, hi, str(wavelength))] = pk
    p = _peaks_disk_path(ref, wavelength)
    if p is not None:
        old = _cache_read(p) if p.exists() else None
        if old is None or (hi - lo) >= (old["hi"] - old["lo"]):     # keep the widest range on disk
            _cache_write(p, {"lo": lo, "hi": hi, "df": pk})


def precompute_peaks(refs, two_theta, wavelength="CuKa"):
    """Make sure every reference has its pattern for this 2theta range: memory, then disk cache,
    then simulate what is missing - on several cores when there is enough to do."""
    lo, hi = float(two_theta[0]), float(two_theta[1])
    todo = []
    for r in refs:
        if r.get("structure") is None:
            continue
        if _memory_peaks(r, lo, hi, wavelength) is not None or _disk_peaks(r, lo, hi, wavelength) is not None:
            continue
        todo.append(r)
    if not todo:
        return 0
    if len(todo) >= PARALLEL_MIN and _n_workers() > 1:
        print(f"  simulating {len(todo)} patterns on {_n_workers()} cores ...")
    res = _parallel_map(_simulate_job, [(r["structure"], lo, hi, wavelength) for r in todo])
    for r, pk in zip(todo, res):
        if not isinstance(pk, Exception):
            _store_peaks(r, lo, hi, wavelength, pk)
    return len(todo)


# ---------------------------------------------------------------------------
# Peak width used for matching
# ---------------------------------------------------------------------------
# Sharp lab peaks (FWHM up to ~0.25 deg) are matched with a fixed 0.15 deg width - that works for
# CPD and NMC data. Broad peaks (MOFs, nanoparticles, poorly crystalline samples: FWHM >= 0.3 deg)
# are matched with their measured width, otherwise a reference with the right peaks but the
# "wrong" sharpness loses to big low-symmetry cells that have a line everywhere.
PEAK_WIDTH = None          # None = automatic; a number (deg) forces it, e.g. X.PEAK_WIDTH = 0.5
BROAD_PEAKS_FROM = 0.30    # measured FWHM (deg) above which the measured width is used


def peak_width(x, y):
    """Peak width (deg) used to match references to this scan (see PEAK_WIDTH)."""
    if PEAK_WIDTH:
        return float(PEAK_WIDTH)
    x, y = np.asarray(x, float), np.asarray(y, float)
    pk = find_experimental_peaks(x, subtract_background(x, y)[0])
    if pk.empty:
        return 0.15
    w = float(pk.nlargest(5, "rel_intensity").fwhm.median())
    return float(min(w, 1.0)) if w >= BROAD_PEAKS_FROM else 0.15


def _widths(x, y, width_deg=None, tol_deg=None, mask_deg=None):
    """(blur width, peak-matching tolerance, masking half-width) for this scan."""
    w = float(width_deg) if width_deg is not None else peak_width(x, y)
    tol = float(tol_deg) if tol_deg is not None else max(0.15, 0.5 * w)
    mask = float(mask_deg) if mask_deg is not None else max(0.25, 0.5 * w)
    return w, tol, mask


def _cached_peaks(ref, two_theta, wavelength):
    """simulate_peaks with a per-reference cache (memory, then the disk cache); reuses any cached
    wider range."""
    lo, hi = float(two_theta[0]), float(two_theta[1])
    pk = _memory_peaks(ref, lo, hi, wavelength)
    if pk is None:
        pk = _disk_peaks(ref, lo, hi, wavelength)
    if pk is None:
        pk = simulate_peaks(ref["structure"], (lo, hi), wavelength)
        _store_peaks(ref, lo, hi, wavelength, pk)
    return pk


class _Matcher:
    """Shared machinery for scoring references against one scan."""

    def __init__(self, x, y, wavelength, two_theta, width_deg, max_shift_pct, step_pct, min_rel):
        from scipy.ndimage import gaussian_filter1d

        lo, hi = two_theta if two_theta else (x.min(), x.max())
        self.lo, self.hi = max(lo, x.min()), min(hi, x.max())
        self.lines = _emission_lines(wavelength)
        self.lam = self.lines[0][0]
        self.wavelength = wavelength
        self.width = width_deg
        self.dx = width_deg / 4
        self.edges = np.arange(self.lo, self.hi + self.dx, self.dx)
        self.grid = 0.5 * (self.edges[1:] + self.edges[:-1])
        sig = width_deg / 2.3548 / self.dx
        k = np.arange(-int(4 * sig) - 1, int(4 * sig) + 2)
        self.kernel = np.exp(-0.5 * (k / sig) ** 2)

        y_bg, _ = subtract_background(x, y)
        self.x, self.y_bg = x, y_bg
        step = np.median(np.diff(x))
        y_s = gaussian_filter1d(y_bg, width_deg / 2.3548 / step)  # same blur as the simulation
        # remove the noise floor, otherwise sqrt-weighting turns noise into "signal"
        floor = np.median(y_s) + 3 * 1.4826 * np.median(np.abs(y_s - np.median(y_s)))
        self.noise_floor = floor
        self.exp_lin = np.clip(np.interp(self.grid, x, y_s) - floor, 0, None)
        self.peaks = find_experimental_peaks(x, y_bg, min_rel)
        self.detection_limit = float(self.peaks.attrs.get("detection_limit_rel", 0.0))   # % of max
        self.peaks = self.peaks[(self.peaks.two_theta >= self.lo) & (self.peaks.two_theta <= self.hi)]
        self.scales = 1 + np.arange(-max_shift_pct, max_shift_pct + 1e-9, step_pct) / 100
        self.sim_range = (max(self.lo - 3, 1.0), min(self.hi + 3, 170.0))

    def sticks(self, ref):
        pk = _cached_peaks(ref, self.sim_range, self.wavelength)
        return pk["d"].to_numpy(float), pk["intensity"].to_numpy(float)

    def line_sticks(self, d, inten, scale):
        """Peak positions for every emission line (Ka1 and Ka2 for a lab doublet).
        Returns (2theta, weighted intensity, index into d, wavelength) arrays."""
        tts, Is, idx, lams = [], [], [], []
        for lam, w in self.lines:
            arg = lam / (2 * d * scale)
            ok = arg < 1
            tts.append(np.degrees(2 * np.arcsin(arg[ok])))
            Is.append(inten[ok] * w)
            idx.append(np.nonzero(ok)[0])
            lams.append(np.full(ok.sum(), lam))
        return np.concatenate(tts), np.concatenate(Is), np.concatenate(idx), np.concatenate(lams)

    def profile(self, d, inten, scale):
        return self.profiles(d, inten, np.array([scale]))[0]

    def profiles(self, d, inten, scales):
        """Blurred stick patterns for several cell scales at once (rows = scales)."""
        from scipy.ndimage import convolve1d

        nb = len(self.edges) - 1
        out = np.zeros((len(scales), nb))
        for lam, w in self.lines:
            arg = lam / (2 * d[None, :] * np.asarray(scales)[:, None])
            ok = arg < 1
            tt = np.degrees(2 * np.arcsin(np.where(ok, arg, 0.0)))
            k = np.floor((tt - self.edges[0]) / self.dx).astype(np.int64)
            ok &= (k >= 0) & (k < nb)
            rows = np.broadcast_to(np.arange(len(scales))[:, None], k.shape)
            flat = (rows * nb + k)[ok]
            out += np.bincount(flat, weights=np.broadcast_to(inten * w, k.shape)[ok],
                               minlength=len(scales) * nb).reshape(len(scales), nb)
        return convolve1d(out, self.kernel, axis=1, mode="constant", cval=0.0)

    def scale_scores(self, d, inten, target_lin, weight=None):
        """Cosine score for every trial cell scale (-inf where undefined)."""
        w = np.ones_like(target_lin) if weight is None else weight
        t = np.sqrt(target_lin) * w
        tn = np.linalg.norm(t)
        if tn == 0:
            return np.full(len(self.scales), -np.inf)
        P = np.sqrt(np.clip(self.profiles(d, inten, self.scales), 0, None)) * w
        pn = np.linalg.norm(P, axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(pn > 0, (P @ t) / (tn * pn), -np.inf)

    def best_scale(self, d, inten, target_lin, weight=None):
        c = self.scale_scores(d, inten, target_lin, weight)
        i = int(np.argmax(c))
        if not np.isfinite(c[i]):
            return (-1.0, 1.0)
        return (float(c[i]), float(self.scales[i]))

    def best_scale_on_peaks(self, d, inten, target_lin, weight, peaks, tol_deg=0.15, near_best=0.8):
        """For the LEFTOVER signal: a minor phase's score hardly changes with the cell scale, so noise
        can pick a stretch that misses every leftover peak. Among the scales scoring within
        near_best of the best, take the one whose visible lines land on the most leftover peaks
        (then the highest score)."""
        c = self.scale_scores(d, inten, target_lin, weight)
        i_best = int(np.argmax(c))
        if not np.isfinite(c[i_best]) or not len(peaks):
            return self.best_scale(d, inten, target_lin, weight)
        cand = np.flatnonzero(c >= near_best * c[i_best]) if c[i_best] > 0 else np.array([i_best])
        if len(cand) == 1:
            return (float(c[cand[0]]), float(self.scales[cand[0]]))
        pk = np.sort(peaks.two_theta.to_numpy())
        any_hit = np.zeros((len(cand), len(pk)), bool)      # scale x leftover peak: a line lands on it
        for lam, w in self.lines:                            # all candidate scales at once
            dv = d[inten * w >= 0.5]
            if not len(dv):
                continue
            arg = lam / (2 * dv[None, :] * self.scales[cand][:, None])
            tt = np.where(arg < 1, np.degrees(2 * np.arcsin(np.clip(arg, 0, 1))), np.inf)
            any_hit |= (np.abs(tt[:, :, None] - pk[None, None, :]) <= tol_deg).any(axis=1)
        n = any_hit.sum(axis=1)
        order = sorted(range(len(cand)), key=lambda j: (n[j], c[cand[j]]), reverse=True)
        i = cand[order[0]]
        return (float(c[i]), float(self.scales[i]))


def _height_ratio(tt, I, peaks, tol_deg=0.15):
    """Measured peak height (% of max) per unit of predicted line intensity, as the median over the
    measured peaks this phase has a line under. None if no peak matches."""
    if not len(peaks) or not len(tt):
        return None
    order = np.argsort(tt)
    ts, Is = tt[order], I[order]
    ratios = []
    for t0, h in zip(peaks.two_theta.to_numpy(), peaks.rel_intensity.to_numpy()):
        lo_i, hi_i = np.searchsorted(ts, t0 - tol_deg), np.searchsorted(ts, t0 + tol_deg, side="right")
        if hi_i > lo_i and Is[lo_i:hi_i].max() > 0:
            ratios.append(h / Is[lo_i:hi_i].max())
    return float(np.median(ratios)) if ratios else None


def _mask_weight(m, phases, target, mask_deg=0.25, min_rel=0.05, max_reflections=150, visible_frac=0.5,
                 max_underprediction=30.0):
    """1 everywhere except within mask_deg of the reflections of the given phases.

    A line (>= min_rel %) is masked when the phase could plausibly produce what is measured there:
      - it would be visible on its own (the phase put on the measured scale by the median height
        ratio of its peaks reaches visible_frac x the detection limit), or
      - the measured peak under it is at most max_underprediction x the line's scaled height -
        weak lines are often under-predicted (fluorite 200 in CPD-1a: measured ~14x the simulation),
        but a 0.05 % line of a big MOF cell cannot account for a peak at 35 % of the maximum.
    At most the max_reflections strongest lines (stable order). Each phase is aligned (isotropic
    cell shift) on what the previous phases left over."""
    weight = np.ones_like(target)
    pk_t = m.peaks.two_theta.to_numpy()
    pk_h = m.peaks.rel_intensity.to_numpy()
    for ph in phases:
        d, inten = m.sticks(ph)
        _, sc = m.best_scale(d, inten, target, weight)
        tt, I, _, _ = m.line_sticks(d, inten, sc)
        keep = I >= min_rel
        k = _height_ratio(tt, I, m.peaks, max(0.15, 0.5 * m.width))
        if k is not None and len(pk_t):
            scaled = k * I
            visible = scaled >= visible_frac * m.detection_limit
            # strongest measured peak within the masked window of each line
            order = np.argsort(pk_t)
            ts, hs = pk_t[order], pk_h[order]
            lo_i = np.searchsorted(ts, tt - mask_deg)
            hi_i = np.searchsorted(ts, tt + mask_deg, side="right")
            hmax = np.array([hs[a:b].max() if b > a else 0.0 for a, b in zip(lo_i, hi_i)])
            keep &= visible | (hmax <= max_underprediction * scaled)
        tt, I = tt[keep], I[keep]
        for t in tt[np.argsort(-I, kind="stable")[:max_reflections]]:
            weight[np.abs(m.grid - t) <= mask_deg] = 0.0
    return weight


def score_against_experiment(refs, x, y, wavelength="CuKa", two_theta=None, width_deg=None,
                             max_shift_pct=1.5, step_pct=0.05, min_rel=None, tol_deg=None,
                             subtract=None, mask_deg=None, progress=True):
    """Rank references by how well they reproduce a measured scan.

    Each reference is simulated, blurred to `width_deg`, and compared with the
    background-subtracted scan (cosine similarity on sqrt intensities, so weak peaks count).
    The cell is allowed to expand/shrink isotropically by up to `max_shift_pct` %, so small
    lattice differences (composition, DFT error, temperature) do not sink a correct phase.

    subtract : list of refs (e.g. [best_phase]). Regions within `mask_deg` of their peaks are
               blanked out first, so the ranking tells you what explains the LEFTOVER peaks
               (Cullity's iterative method for mixtures). Minor phases whose peaks all overlap
               the main phase cannot be found this way.

    Returns a DataFrame sorted by score (1 = perfect) with:
      cell_scale_pct      how much larger (+) / smaller (-) the measured cell is than the reference
      explained_pct       share of measured peak intensity that has a reference peak within tol_deg
      unexplained_peaks   measured peaks (2theta, %) the reference does NOT account for
      missing_strong      reference peaks >= 10 % that are absent from the scan
      missing_expected    reference lines that SHOULD be visible but are not: the reference is
                          scaled to the measured heights of the peaks it explains, and a line
                          predicted at >= 2x the detection limit with no measured peak counts.
                          This works for minor phases too (their weak lines are below the noise
                          and are not counted) - a wrong phase that fits a few peaks by chance
                          predicts many strong lines that are not there.
    Near-ties (score within 0.005) are ordered by the smallest cell stretch, so an experimental
    entry with the exact cell wins over a DFT one that needs 1 % stretching.
    """
    x, y = np.asarray(x, float), np.asarray(y, float)
    width_deg, tol_deg, mask_deg = _widths(x, y, width_deg, tol_deg, mask_deg)
    m = _Matcher(x, y, wavelength, two_theta, width_deg, max_shift_pct, step_pct, min_rel)
    target = m.exp_lin.copy()

    note, weight = "", None
    if subtract:
        # Cullity's iterative method: blank out everything the known phase(s) explain,
        # then rank the references on what is left.
        weight = _mask_weight(m, subtract, target, mask_deg)
        note = " ".join(r["ref_id"] for r in subtract)
        keep = np.interp(m.peaks.two_theta, m.grid, weight) > 0.5
        peaks = m.peaks[keep]
    else:
        peaks = m.peaks

    rows = []
    good = [r for r in refs if r.get("structure") is not None]
    precompute_peaks(good, m.sim_range, wavelength)          # cache / several cores
    sticks = []
    for r in good:
        try:
            sticks.append(m.sticks(r))
        except Exception as ex:
            sticks.append(None)
            rows.append({"ref_id": r["ref_id"], "source": r["source"], "error": str(ex)[:120]})
    idx = [i for i, s_ in enumerate(sticks) if s_ is not None]
    n_cores = _n_workers()
    if n_cores > 1 and len(idx) >= PARALLEL_MIN and not _PARALLEL_BROKEN:
        chunks = [list(c) for c in np.array_split(idx, min(len(idx), 3 * n_cores)) if len(c)]
        jobs = [(m, target, weight, peaks, tol_deg, [sticks[i] for i in c]) for c in chunks]
        res = [row for part in _parallel_map(_score_job, jobs, min_jobs=2) for row in part]
    else:
        res = []
        for k, i in enumerate(idx):
            if progress and k and k % 100 == 0:
                print(f"  scored {k}/{len(idx)}")
            res.append(_score_one(m, target, weight, peaks, tol_deg, *sticks[i]))
    for i, sc_row in zip(idx, res):
        r = good[i]
        rows.append({
            "ref_id": r["ref_id"], "source": r["source"],
            "kind": r.get("kind") or {"COD": "experimental", "MP": "DFT", "built": "model"}.get(r["source"], ""),
            "label": r["label"], "formula": r["formula"], "spacegroup": r["spacegroup"], **sc_row})
    if not rows:                                   # no candidates: an empty table with the usual columns
        return pd.DataFrame(columns=["ref_id", "source", "kind", "label", "formula", "spacegroup", "score",
                                     "cell_scale_pct", "explained_pct", "n_unexplained", "unexplained_peaks",
                                     "missing_strong", "missing_expected"])
    out = pd.DataFrame(rows).sort_values("score", ascending=False, na_position="last").reset_index(drop=True)
    if not subtract:
        out = _prefer_exact_cell(out)
    if note:
        out.attrs["subtracted"] = note.strip()
    return out


def _near(a, b, tol):
    """For each value in a: is there a value in b within tol? (vectorised)"""
    a = np.asarray(a, float)
    b = np.sort(np.asarray(b, float))
    if not len(b) or not len(a):
        return np.zeros(len(a), bool)
    i = np.clip(np.searchsorted(b, a), 1, len(b) - 1)
    return (np.abs(b[i - 1] - a) <= tol) | (np.abs(b[i] - a) <= tol) | (np.abs(b[0] - a) <= tol)


def _score_one(m, target, weight, peaks, tol_deg, d, inten):
    """Score one reference (its d-spacings and intensities) against the scan held by matcher m."""
    if weight is None:
        score, sc = m.best_scale(d, inten, target, weight)
    else:                                   # leftover signal: prefer the scale that hits leftover peaks
        score, sc = m.best_scale_on_peaks(d, inten, target, weight, peaks, tol_deg)
    tt, I, _, lam = m.line_sticks(d, inten, sc)
    inr = (tt >= m.lo) & (tt <= m.hi)
    tt, I, lam = tt[inr], I[inr], lam[inr]
    visible = tt[I >= 0.5]
    expl = _near(peaks.two_theta.to_numpy(), visible, tol_deg) if len(peaks) else np.array([], bool)
    unexpl = peaks[~expl] if len(peaks) else peaks
    tot = peaks.rel_intensity.sum() if len(peaks) else 0
    # "missing" = predicted at least as strong as this phase's weakest DETECTED line (and >= 10 %),
    # yet absent: for a minor phase, lines below the noise are not counted against it
    all_pk = m.peaks.two_theta.to_numpy()
    own_pk = peaks.two_theta.to_numpy() if len(peaks) else np.array([])
    seen = _near(tt, all_pk, tol_deg)
    own = _near(tt, own_pk, tol_deg)                       # lines this phase explains
    floor_I = max(10.0, I[own & (I >= 0.5)].min()) if np.any(own & (I >= 0.5)) else 10.0
    missing = [t for t, i_, s_ in zip(tt, I, seen) if i_ >= floor_I and not s_]
    return {
        "score": round(score, 4),
        "cell_scale_pct": round((sc - 1) * 100, 2),
        "explained_pct": round(100 * peaks.rel_intensity[expl].sum() / tot, 1) if tot else np.nan,
        "n_unexplained": int(len(unexpl)),
        "unexplained_peaks": ", ".join(f"{t:.2f}({i:.0f}%)" for t, i in
                                       unexpl.sort_values("rel_intensity", ascending=False)
                                       [["two_theta", "rel_intensity"]].head(6).to_numpy()),
        "missing_strong": len(missing),
        "missing_expected": _missing_expected(tt, I, lam, own, seen, peaks, m.lam, m.detection_limit, tol_deg),
    }


def _score_job(job):
    """Worker: score a chunk of references (several cores)."""
    m, target, weight, peaks, tol_deg, items = job
    return [_score_one(m, target, weight, peaks, tol_deg, d, inten) for d, inten in items]


def _missing_expected(tt, I, lam, own, seen, peaks, lam1, det_limit, tol_deg, factor=2.0, b_iso=0.5):
    """Count Ka1 lines that should be above the detection limit but have no measured peak.
    The reference is put on the measured scale with the median height ratio of the lines it
    explains; a mild Debye-Waller damping (B = b_iso A^2) keeps high-angle lines realistic."""
    if not len(peaks) or not np.any(own):
        return 0
    s2 = (np.sin(np.radians(tt / 2)) / lam) ** 2
    Id = I * np.exp(-2 * b_iso * s2)
    ratios = []
    for t, h in zip(peaks.two_theta.to_numpy(), peaks.rel_intensity.to_numpy()):
        near = np.abs(tt - t) <= tol_deg
        if np.any(near) and Id[near].max() > 0:
            ratios.append(h / Id[near].max())
    if not ratios:
        return 0
    k = float(np.median(ratios))
    thr = factor * max(det_limit, 1e-6)
    ka1 = np.abs(lam - lam1) < 1e-6
    return int(np.sum(ka1 & (k * Id >= thr) & ~seen))


def _prefer_exact_cell(rank, tol=0.005):
    """Among the near-tied best entries (score within tol of the top), put the one needing the
    smallest cell stretch first: score differences that small mean nothing."""
    if rank.empty or "score" not in rank or rank.score.isna().all():
        return rank
    top = rank.score.iloc[0]
    tied = rank[(rank.score >= top - tol) & rank.cell_scale_pct.notna()]
    if len(tied) < 2:
        return rank
    best = tied.cell_scale_pct.abs().idxmin()
    order = [best] + [i for i in rank.index if i != best]
    return rank.loc[order].reset_index(drop=True)


def pick_for_plot(ranking, n_top=4, per_source=1, sources=None, table=None):
    """ref_ids to plot: the n_top best overall plus the best `per_source` of each source,
    skipping duplicates (same 'group' in the summary table, if given)."""
    seen_groups, picked = set(), []

    def take(row):
        g = None
        if table is not None and "group" in table:
            gg = table.loc[table.ref_id == row.ref_id, "group"]
            g = None if gg.empty or pd.isna(gg.iloc[0]) else int(gg.iloc[0])
        if row.ref_id in picked or (g is not None and g in seen_groups):
            return False
        picked.append(row.ref_id)
        if g is not None:
            seen_groups.add(g)
        return True

    rk = ranking.dropna(subset=["score"])
    if sources is None:                      # every source present in the ranking, best first
        sources = list(dict.fromkeys(rk.source))
    for row in rk.itertuples():
        if len(picked) >= n_top:
            break
        take(row)
    for src in sources:
        k = 0
        for row in rk[rk.source == src].itertuples():
            if k >= per_source:
                break
            if row.ref_id in picked or take(row):
                k += 1
    return picked


def _cell_system(lattice, tol=1e-3):
    a, b, c = lattice.abc
    al, be, ga = lattice.angles
    eq = lambda u, v: abs(u - v) <= tol * max(u, v)
    right = all(abs(t - 90) < 0.01 for t in (al, be, ga))
    if right and eq(a, b) and eq(b, c):
        return "cubic"
    if right and eq(a, b):
        return "tetragonal"
    if abs(al - 90) < 0.01 and abs(be - 90) < 0.01 and abs(ga - 120) < 0.01 and eq(a, b):
        return "hexagonal"
    if right:
        return "orthorhombic"
    return None


def _cell_row(system, h, k, l):
    """Design-matrix row: 1/d^2 = row . coefficients (coefficients 1/a^2, 1/c^2, ...)."""
    if system == "cubic":
        return (h * h + k * k + l * l,)
    if system == "tetragonal":
        return (h * h + k * k, l * l)
    if system == "hexagonal":
        return (4 / 3 * (h * h + h * k + k * k), l * l)
    return (h * h, k * k, l * l)


def fit_cell(x, y, ref, wavelength="CuKa", other_phases=(), tol_deg=None, min_rel=None, refine_zero=True,
             zero_shift=None):
    """Lattice parameters of one phase from YOUR peak positions, indexed with `ref`.

    Works for cubic, tetragonal, hexagonal/trigonal (hexagonal setting, incl. R-3m) and
    orthorhombic references; monoclinic/triclinic are refused.
    In a mixture, pass the other identified phases as `other_phases`: measured peaks that one of
    them could also explain are left out, so only this phase's own peaks are used.
    For a Ka1/Ka2 doublet ('CuKa'), resolved Ka1 and Ka2 peaks are each fitted with their own
    wavelength, blended low-angle peaks with the intensity-weighted mean.
    zero_shift : fix the 2theta zero offset (deg) instead of fitting it - use the main phase's value
                 for minor phases, which have too few peaks to fit cell and offset together.
    Returns (result dict, DataFrame of the peaks used).
    """
    x, y = np.asarray(x, float), np.asarray(y, float)
    system = _cell_system(ref["structure"].lattice)
    if system is None:
        raise ValueError(f"{ref['ref_id']}: monoclinic/triclinic cells are not supported by this fit")
    w_fit, tol_deg, _ = _widths(np.asarray(x, float), np.asarray(y, float), None, tol_deg, None)
    m = _Matcher(x, y, wavelength, None, w_fit, 1.5, 0.05, min_rel)

    phases = [ref] + [o for o in other_phases if o["ref_id"] != ref["ref_id"]]

    def sticks_of(r):
        # align each phase with the OTHER phases blanked out, as identify_phases does - on the
        # whole pattern a minor phase would be pulled onto the main phase's peaks
        d, inten = m.sticks(r)
        rest = [p for p in phases if p["ref_id"] != r["ref_id"]]
        w = _mask_weight(m, rest, m.exp_lin) if rest else None
        _, sc = m.best_scale(d, inten, m.exp_lin, w)
        return d, inten, sc

    d_ref, I_ref, sc = sticks_of(ref)
    pk = _cached_peaks(ref, (max(m.lo - 3, 1.0), min(m.hi + 3, 170.0)), wavelength)
    hkls = list(pk["hkl"])
    tt, I, idx, lam = m.line_sticks(d_ref, I_ref, sc)
    others = []
    for o in other_phases:
        if o["ref_id"] == ref["ref_id"]:
            continue
        od, oi, osc = sticks_of(o)
        ott, oI, _, _ = m.line_sticks(od, oi, osc)
        others.append(ott[oI >= 1.0])
    others = np.concatenate(others) if others else np.array([])
    lam_avg = sum(l * w for l, w in m.lines) / sum(w for _, w in m.lines)

    rows = []
    for _, p in m.peaks.iterrows():
        near = np.abs(tt - p.two_theta) <= tol_deg
        near &= I >= 0.5
        if not near.any():
            continue
        if others.size and np.any(np.abs(others - p.two_theta) <= tol_deg):
            continue                                   # shared with another phase: ambiguous, skip
        j = np.flatnonzero(near)[np.argmin(np.abs(tt[near] - p.two_theta))]
        fam = {tuple(round(v, 9) for v in _cell_row(system, *hh)) for hh in hkls[idx[j]]}  # signed: h*k matters
        if len(fam) != 1:
            continue                                   # different reflections at the same angle
        same = np.flatnonzero(idx == idx[j])           # this reflection's Ka1 and Ka2 sticks
        blended = len(same) > 1 and abs(tt[same].max() - tt[same].min()) < max(p.fwhm, 0.1)
        win = (x > p.two_theta - max(p.fwhm, 0.02)) & (x < p.two_theta + max(p.fwhm, 0.02))
        yw = m.y_bg[win]
        top_ = yw >= 0.5 * yw.max()
        centre = float(np.sum(x[win][top_] * yw[top_]) / np.sum(yw[top_]))
        (row,) = fam
        rows.append({"hkl": _peak_label([hkls[idx[j]][0]], hexagonal=(system == "hexagonal")),
                     "line": "Ka1+Ka2" if blended else ("Ka2" if lam[j] > m.lines[0][0] + 1e-6 else "Ka1"),
                     "lambda": lam_avg if blended else lam[j], "two_theta_obs": centre,
                     "rel_intensity": p.rel_intensity, "row": row})
    used = pd.DataFrame(rows)
    if used.empty:
        raise ValueError("No peaks could be indexed for this phase.")
    used = used.drop_duplicates(subset=["two_theta_obs"])
    A = np.array([list(r) for r in used.row])
    if len(used) < A.shape[1] + (0 if zero_shift is not None else 1) or np.linalg.matrix_rank(A) < A.shape[1]:
        raise ValueError(f"Too few independent indexed peaks ({len(used)}) for a {system} cell.")

    def solve(zero):
        th = np.radians((used.two_theta_obs.to_numpy() - zero) / 2)
        inv_d2 = (2 * np.sin(th) / used["lambda"].to_numpy()) ** 2
        coef, *_ = np.linalg.lstsq(A, inv_d2, rcond=None)
        calc = np.degrees(2 * np.arcsin(np.clip(used["lambda"].to_numpy() / 2 * np.sqrt(np.clip(A @ coef, 0, None)), 0, 1))) + zero
        return coef, calc

    if zero_shift is not None:
        zeros = [float(zero_shift)]
    else:
        zeros = np.arange(-0.1, 0.1001, 0.0005) if refine_zero and len(used) > A.shape[1] + 1 else [0.0]
    best = min(zeros, key=lambda z: np.sum((solve(z)[1] - used.two_theta_obs.to_numpy()) ** 2))
    coef, calc = solve(best)
    if np.any(coef <= 0):
        raise ValueError("Fit gave an unphysical cell - peaks were probably mis-indexed.")
    p = 1 / np.sqrt(coef)
    if system == "cubic":
        a_ = b_ = c_ = p[0]
    elif system in ("tetragonal", "hexagonal"):
        a_, b_, c_ = p[0], p[0], p[1]
    else:
        a_, b_, c_ = p[0], p[1], p[2]
    vol = a_ * b_ * c_ * (np.sqrt(3) / 2 if system == "hexagonal" else 1)
    used["two_theta_calc"] = calc
    used["diff_deg"] = used.two_theta_obs - calc
    res = {"phase": ref["ref_id"], "formula": ref.get("formula"), "system": system,
           "a": a_, "b": b_, "c": c_, "volume": vol,
           "zero_shift_deg": best, "zero_fixed": zero_shift is not None, "rms_deg": float(np.sqrt(np.mean(used.diff_deg ** 2))),
           "n_peaks": len(used)}
    if system == "hexagonal":
        res["c_over_a"] = c_ / a_
    res["indexed_with"] = ref["ref_id"]
    if res["rms_deg"] > 0.02:
        res["warning"] = ("Poor fit (rms > 0.02 deg): peaks were probably mis-indexed, or the "
                          "phase is a minor one with few clean peaks.")
    return res, used.drop(columns=["row"])


def fit_hexagonal_cell(x, y, ref, wavelength="CuKa", tol_deg=0.15, min_rel=0.02, refine_zero=True):
    """Hexagonal / R-3m lattice parameters (kept for older notebooks; see fit_cell)."""
    return fit_cell(x, y, ref, wavelength, tol_deg=tol_deg, min_rel=min_rel, refine_zero=refine_zero)

# ---------------------------------------------------------------------------
# Impurity / secondary-phase search
# ---------------------------------------------------------------------------
# COD element searches (exact element sets) for the side phases usually met in a family of
# materials. Add your own lists; each entry is one COD query.
IMPURITY_LIBRARIES = {
    "Li-ion layered oxide (NMC/NCA/LCO)": [
        ["Ni", "O", "H"], ["Co", "O", "H"], ["Mn", "O", "H"],   # M(OH)2, MOOH precursors / hydrates
        ["Ni", "O"], ["Co", "O"], ["Mn", "O"],                  # NiO rock salt, Co3O4, Mn3O4, MnO2 ...
        ["Li", "O", "H"], ["Li", "C", "O"], ["Li", "O"],        # LiOH, Li2CO3, Li2O (surface residues)
        ["Ni", "C", "O"], ["Co", "C", "O"], ["Mn", "C", "O"],   # carbonate precursors
        ["Li", "S", "O"],                                        # Li2SO4 (sulfate route)
        ["Li", "Mn", "O"], ["Li", "Ni", "O"], ["Li", "Co", "O"],  # spinel, Li2MnO3, Li2NiO2 ...
        ["Al", "O"], ["Li", "Al", "O"],                          # Al2O3 / LiAlO2 (NCA, coatings)
    ],
}


def built_impurity_models():
    """Offline fallbacks (no database needed): brucite-type hydroxides and NiO rock salt.
    Approximate literature cells."""
    nio = Structure.from_spacegroup("Fm-3m", Lattice.cubic(4.177), ["Ni", "O"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    return [
        build_tm_hydroxide({"Ni": 1.0}, 3.126, 4.605, label="Ni(OH)2"),
        build_tm_hydroxide({"Co": 1.0}, 3.183, 4.652, label="Co(OH)2"),
        build_tm_hydroxide({"Mn": 1.0}, 3.322, 4.734, label="Mn(OH)2"),
        _record("built", "NiO", "NiO (rock salt)", nio, None),
    ]


def search_impurities(library="Li-ion layered oxide (NMC/NCA/LCO)", out_dir="refs_impurities",
                      max_per_search=80, include_built=None, refresh=False, pause_between=2.0):
    """Fetch candidate side phases from COD for every query in IMPURITY_LIBRARIES[library]
    (or a list of element lists you pass instead). Failed queries are reported and skipped.
    Returns a list of reference records tagged with r["role"] = "impurity"."""
    queries = IMPURITY_LIBRARIES[library] if isinstance(library, str) else library
    refs, failed = [], []
    for q in queries:
        online = refresh or not _cod_saved_path(q, Path(out_dir) / "COD", True).exists()
        try:
            got = search_cod(q, out_dir=Path(out_dir) / "COD", exact=True, max_entries=max_per_search,
                             refresh=refresh)
            ok = [r for r in got if r["structure"] is not None]
            print(f"  impurities {'-'.join(q):10s}: {len(ok)} structures")
            refs += ok
        except Exception as ex:
            failed.append("-".join(q))
            print(f"  impurities {'-'.join(q):10s}: FAILED ({str(ex)[:80]})")
        if online:
            time.sleep(pause_between)      # COD drops rapid back-to-back queries (no wait when reused)
    if include_built is None:      # the built Li-ion side-phase models only belong to the named Li-ion library
        include_built = isinstance(library, str)
    if include_built:
        refs += built_impurity_models()
    for r in refs:
        r["role"] = "impurity"
    if failed:
        print("Re-run later to fill in:", ", ".join(failed), "(finished queries are reused)")
    return list({r["ref_id"]: r for r in refs}.values())


def unexplained_peaks(phases, x, y, wavelength="CuKa", width_deg=None, min_rel=None, mask_deg=None):
    """Measured peaks that none of `phases` accounts for (DataFrame: two_theta, rel_intensity, fwhm)."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    width_deg, _, mask_deg = _widths(x, y, width_deg, None, mask_deg)
    m = _Matcher(x, y, wavelength, None, width_deg, 1.5, 0.05, min_rel)
    w = _mask_weight(m, phases, m.exp_lin.copy(), mask_deg)
    keep = np.interp(m.peaks.two_theta, m.grid, w) > 0.5
    return m.peaks[keep].reset_index(drop=True)


def _has_h(ref):
    return any(el.symbol == "H" for el in ref["structure"].composition.elements)


def _same_phase(r1, r2, vol_tol=0.04):
    """Same crystalline phase as far as X-rays can tell: same space group, same composition
    ignoring hydrogen (many older CIFs omit H atoms - e.g. Ni(OH)2 stored as 'NiO2'), and the
    same volume per non-H atom within vol_tol. A different cell (real NiO2, a ~ 2.8 A vs
    Ni(OH)2, a ~ 3.1 A) therefore stays a different phase."""
    from pymatgen.core import Composition

    if r1 is r2 or r1["ref_id"] == r2["ref_id"]:
        return True
    if r1.get("spacegroup") != r2.get("spacegroup"):
        return False

    def no_h(ref):
        comp = ref["structure"].composition.element_composition
        return Composition({el: v for el, v in comp.items() if el.symbol != "H"})

    c1, c2 = no_h(r1), no_h(r2)
    if not c1.fractional_composition.almost_equals(c2.fractional_composition, rtol=0.05, atol=0.02):
        return False
    v1 = r1["structure"].volume / c1.num_atoms
    v2 = r2["structure"].volume / c2.num_atoms
    return abs(v1 - v2) / max(v1, v2) <= vol_tol


def identify_phases(candidates, x, y, wavelength="CuKa", impurity_candidates=(), max_phases=4,
                    min_score=0.30, min_explained_pct=50.0, min_peaks=2, min_rel=None,
                    ambiguity_margin=0.03, minor_phase_peaks=3, verbose=True):
    """Automatic phase identification by Cullity's iterative procedure:

      1. the best-scoring candidate is the main phase;
      2. its peaks are blanked out; if measured peaks remain, ALL candidates (main + impurity
         pool) are ranked on the leftover signal;
      3. the winner is accepted if it accounts for >= min_explained_pct % of the leftover peak
         intensity AND for at least `min_peaks` leftover peaks (a single peak can be matched by
         almost anything), and its leftover score is >= min_score - or, for a minor phase, it
         explains >= minor_phase_peaks leftover peaks; then repeat from 2.
         Otherwise the leftover peaks are reported with the most likely candidates.
    Stops when nothing is left unexplained, nothing qualifies, or max_phases is reached.

    Returns (phases, report): the accepted reference records, and a DataFrame with one row
    per step (score, what it explained, what was still unexplained, strong peaks it predicts
    that are NOT in the scan -> a warning sign for a wrong assignment).
    """
    pool = list({r["ref_id"]: r for r in list(candidates) + list(impurity_candidates)}.values())
    pool = [r for r in pool if r.get("structure") is not None]
    rank = score_against_experiment(list(candidates), x, y, wavelength, min_rel=min_rel, progress=False)
    rank = rank.dropna(subset=["score"])
    if rank.empty:
        raise ValueError("identify_phases: no candidate structures to compare with the scan - check ELEMENTS, "
                         "DATABASES and that the library holds structures with these elements")
    by_id = {r["ref_id"]: r for r in pool}
    first = rank.iloc[0]
    phases = [by_id[first.ref_id]]
    second = rank.iloc[1] if len(rank) > 1 else None
    rows = [{"step": 1, "ref_id": first.ref_id, "formula": first.formula, "spacegroup": first.spacegroup,
             "score": first.score, "explained_pct": first.explained_pct,
             "missing_strong": first.missing_strong, "missing_expected": first.missing_expected,
             "role": "main phase",
             "confidence": "" if second is None else f"next best {second.ref_id} {second.formula} ({second.score:.2f})"}]

    while len(phases) < max_phases:
        left = unexplained_peaks(phases, x, y, wavelength, min_rel=min_rel)
        rows[-1]["unexplained_after"] = ", ".join(f"{t:.2f}({i:.0f}%)" for t, i in
                                                  left[["two_theta", "rel_intensity"]].to_numpy()) or "none"
        if left.empty:
            break
        rest = [r for r in pool if not any(_same_phase(r, p) for p in phases)]
        rk = score_against_experiment(rest, x, y, wavelength, subtract=phases, min_rel=min_rel, progress=False)
        rk = rk.dropna(subset=["score"]).reset_index(drop=True)
        if rk.empty:
            rows[-1]["not_assigned"] = (f"{len(left)} leftover peak(s) not assigned: no other candidate "
                                        "structures (switch more DATABASES on or add CIFs)")
            if verbose:
                print("  ->", rows[-1]["not_assigned"])
            break
        # on leftovers, a phase must explain MANY of the leftover peaks, not just fit one of them well:
        # rank by pattern score x fraction of leftover peak intensity it accounts for
        # small penalty for needing a big cell stretch (an exact-cell match is more convincing)
        # ... and must not predict clearly visible lines that are absent (consistency): a big
        # low-symmetry cell can place lines on any 4 leftover peaks, but then predicts dozens of
        # strong lines that are not in the scan
        n_hit = len(left) - rk.n_unexplained.fillna(len(left))
        rk["consistency"] = (n_hit / (n_hit + rk.missing_expected.fillna(0))).where(n_hit > 0, 0.0)
        rk["leftover_rank"] = (rk.score * rk.explained_pct.fillna(0) / 100 * rk.consistency
                               * (1 - 0.05 * rk.cell_scale_pct.abs() / 1.5))
        rk = rk.sort_values("leftover_rank", ascending=False).reset_index(drop=True)
        same_as_top = [_same_phase(by_id[rid], by_id[rk.ref_id.iloc[0]]) for rid in rk.ref_id]
        rk["same_phase_as_top"] = np.array(same_as_top, dtype=bool)
        note = ""
        top = rk.iloc[0]
        if not _has_h(by_id[top.ref_id]):
            # an equivalent entry WITH hydrogen atoms tells the chemistry better (Ni(OH)2, not "NiO2")
            alt = rk[rk.same_phase_as_top & (rk.leftover_rank >= top.leftover_rank - ambiguity_margin)
                     & np.array([_has_h(by_id[rid]) for rid in rk.ref_id], dtype=bool)]
            if len(alt):
                note = (f"{top.ref_id} lists no H atoms ({top.formula}); same structure as "
                        f"{alt.iloc[0].ref_id} {alt.iloc[0].formula} - X-rays barely see H")
                top = alt.iloc[0]
        n_expl = len(left) - int(top.n_unexplained)
        # a minor phase shows only its strongest lines above the noise, so its pattern score on the
        # leftover signal stays low; three or more leftover peaks fitting one cell are evidence enough
        score_ok = top.score >= min_score or (n_expl >= minor_phase_peaks and top.score >= 0.05)
        consistent = top.missing_expected <= max(n_expl, 1)
        ok = (score_ok and (top.explained_pct or 0) >= min_explained_pct and n_expl >= min_peaks
              and consistent)
        # how sure is the pick? near-ties and one-peak matches are only suggestions
        # a rival explains (mostly) the SAME leftover peaks about as well; a candidate explaining
        # OTHER leftover peaks is not a rival - it is the next phase
        def _explained(row):
            miss = {round(float(t.split("(")[0]), 2) for t in str(row.unexplained_peaks).split(", ") if "(" in t}
            return {round(t, 2) for t in left.two_theta} - miss
        top_set = _explained(top)
        rivals = rk[(rk.leftover_rank >= top.leftover_rank - ambiguity_margin)
                    & ~rk.same_phase_as_top]          # other entries of the same phase are not rivals
        rivals = rivals[np.array([len(_explained(r) & top_set) >= 0.5 * max(len(top_set), 1)
                                  for r in rivals.itertuples()], dtype=bool)]
        if n_expl < 2:
            confidence = "tentative: rests on a single peak"
        elif len(rivals):
            confidence = "ambiguous: as good -> " + ", ".join(
                f"{r.ref_id} {r.formula} ({r.leftover_rank:.2f})" for r in rivals.head(3).itertuples())
        else:
            other = rk[~rk.same_phase_as_top]
            confidence = (f"clear (next best other phase {other.iloc[0].ref_id} {other.iloc[0].formula}: "
                          f"{other.iloc[0].leftover_rank:.2f} vs {top.leftover_rank:.2f})") if len(other) else "clear"
        if verbose:
            print(f"step {len(phases) + 1}: {len(left)} leftover peaks; best candidate {top.ref_id} "
                  f"({top.formula}, {top.spacegroup}, cell {top.cell_scale_pct:+.2f} %) leftover score {top.score:.2f}, explains "
                  f"{top.explained_pct:.0f}% ({n_expl} peaks), {int(top.missing_expected)} expected lines absent "
                  f"-> {'ACCEPTED' if ok else 'rejected'}; {confidence}")
        if not ok:
            if n_expl < min_peaks:
                why = f"explains only {n_expl} leftover peak(s) - one peak cannot identify a phase"
            elif not consistent:
                why = (f"best candidate predicts {int(top.missing_expected)} clearly visible lines that are "
                       f"absent - the right reference is probably not in the candidate list")
            else:
                why = f"too weak (score {top.score:.2f}, explains {top.explained_pct:.0f}%)"
            rows[-1]["not_assigned"] = (f"{len(left)} leftover peak(s) not assigned: {why}. Possible: "
                                        + "; ".join(f"{r.formula} [{r.spacegroup}] {r.ref_id}"
                                                    for r in rk.head(5).itertuples()))
            if verbose:
                print("  ->", rows[-1]["not_assigned"])
            break
        phases.append(by_id[top.ref_id])
        rows.append({"step": len(phases), "ref_id": top.ref_id, "formula": top.formula,
                     "spacegroup": top.spacegroup, "score": top.score, "explained_pct": top.explained_pct,
                     "missing_strong": top.missing_strong, "missing_expected": top.missing_expected,
                     "role": "secondary phase" if by_id[top.ref_id].get("role") != "impurity" else "impurity",
                     "cell_scale_pct": top.cell_scale_pct, "peaks_explained": n_expl, "confidence": confidence,
                     "note": note})
    else:
        left = unexplained_peaks(phases, x, y, wavelength, min_rel=min_rel)
        rows[-1]["unexplained_after"] = ", ".join(f"{t:.2f}({i:.0f}%)" for t, i in
                                                  left[["two_theta", "rel_intensity"]].to_numpy()) or "none"
    return phases, pd.DataFrame(rows)



# ---------------------------------------------------------------------------
# Offline structure collections: MOF databases and CIFs you download yourself
# ---------------------------------------------------------------------------
# Each collection lives in its own library subfolder (library/<name>/) with an index file
# (_index.json: elements, formula, cell of every CIF). The index is what makes them usable:
# filtering 20 000 MOFs by elements takes milliseconds, and only the matching CIFs are parsed.
COLLECTIONS = {
    "CoRE_MOF": {
        "title": "CoRE MOF 2024 - experimental MOF structures (computation-ready: solvent removed; "
                 "plus not-computation-ready entries)",
        "kind": "experimental",
        "download": "https://zenodo.org/records/15055758/files/CoREMOF2024DB_SI_20250204.zip?download=1",
        "page": "https://zenodo.org/records/15055758",
        "note": "about 44 MB, 8 300 structures, CC BY 4.0",
    },
    "QMOF": {
        "title": "QMOF database - DFT-optimised MOF structures",
        "kind": "DFT",
        "figshare_article": 13147324,
        "page": "https://doi.org/10.6084/m9.figshare.13147324",
        "note": "large download; DFT cells can differ from experiment by 1-3 %",
    },
    "CSD": {
        "title": "CIFs you downloaded from the CCDC (Access Structures / your CSD licence)",
        "kind": "experimental",
        "manual": True,
        "page": "https://www.ccdc.cam.ac.uk/structures/",
    },
    "Published": {
        "title": "CIFs from papers' supporting information and data repositories (Zenodo, Mendeley ...)",
        "kind": "experimental",
        "manual": True,
    },
}

_INDEX = "_index.json"


def _element_symbol(token):
    """'Cu1' -> 'Cu', 'CU1' -> 'Cu', 'O2A' -> 'O', 'C10' -> 'C', 'Zn2+' -> 'Zn'; None if not an element."""
    from pymatgen.core.periodic_table import Element

    m = re.match(r"([A-Za-z]{1,2})", token.strip().strip("'\""))
    if not m:
        return None
    s = m.group(1)
    for cand in ((s[0].upper() + s[1:].lower()) if len(s) == 2 else None, s[0].upper()):
        if cand and Element.is_valid_symbol(cand):
            return cand
    return None


def _cif_quick_info(text):
    """Elements, formula and cell of a CIF without building the structure (fast, for indexing).
    Elements come from the atom-site loop (what pymatgen will actually build); the formula from
    _chemical_formula_sum when present."""
    info = {"elements": [], "formula": "", "a": None, "b": None, "c": None, "spacegroup": "", "n_sites": 0}

    def num(key):
        m = re.search(rf"^\s*{key}\s+([-+0-9.]+)", text, re.M)
        return float(m.group(1)) if m else None

    info["a"], info["b"], info["c"] = num("_cell_length_a"), num("_cell_length_b"), num("_cell_length_c")
    m = re.search(r"^\s*_(?:symmetry_space_group_name_H-M|space_group_name_H-M_alt)\s+(['\"]?)(.+?)\1\s*$", text, re.M)
    if m:
        info["spacegroup"] = m.group(2).strip()
    m = re.search(r"^\s*_chemical_formula_sum\s+(['\"]?)(.+?)\1\s*$", text, re.M)
    if m:
        info["formula"] = m.group(2).strip()

    # atom-site loop: header lines starting with _atom_site_, then rows until a blank line / loop_ / _tag
    lines = text.splitlines()
    els = set()
    i = 0
    while i < len(lines):
        if lines[i].strip().lower() == "loop_":
            j, heads = i + 1, []
            while j < len(lines) and lines[j].strip().startswith("_"):
                heads.append(lines[j].strip().split()[0].lower())
                j += 1
            if any(h.startswith("_atom_site_") and not h.startswith("_atom_site_aniso") for h in heads):
                col = (heads.index("_atom_site_type_symbol") if "_atom_site_type_symbol" in heads
                       else heads.index("_atom_site_label") if "_atom_site_label" in heads else None)
                n = 0
                while j < len(lines):
                    ln = lines[j].strip()
                    if not ln or ln.lower() == "loop_" or ln.startswith("_") or ln.startswith("#"):
                        break
                    parts = ln.split()
                    if col is not None and len(parts) > col:
                        el = _element_symbol(parts[col])
                        if el and el not in ("D", "T"):
                            els.add(el)
                        elif el in ("D", "T"):
                            els.add("H")
                    n += 1
                    j += 1
                info["n_sites"] += n
            i = j
            continue
        i += 1
    if not els and info["formula"]:
        els = {e for e in (_element_symbol(t) for t in re.findall(r"[A-Za-z]{1,2}", info["formula"])) if e}
    info["elements"] = sorted(els)
    return info


def _safe_name(rel):
    """A flat, unique file name from a path inside an archive: 'CR/ASR/ABAVIJ.cif' -> 'CR__ASR__ABAVIJ.cif'."""
    parts = [p for p in re.split(r"[\\/]+", rel) if p not in ("", ".", "..")]
    name = "__".join(parts)
    name = re.sub(r"[^A-Za-z0-9._+-]", "_", name)
    if len(name) > 150:                       # Windows path limits
        stem, ext = os.path.splitext(name)
        import hashlib
        name = stem[:110] + "_" + hashlib.md5(name.encode()).hexdigest()[:8] + ext
    return name


def _extract_cifs(archive, dest, prefix="", verbose=True):
    """Copy every .cif out of a zip / tar(.gz) archive into dest (nested archives included).
    Returns the list of written file names."""
    import io
    import tarfile
    import zipfile

    written = []

    def handle(name, data):
        low = name.lower()
        if low.endswith(".cif"):
            out = _safe_name(prefix + name)
            (dest / out).write_bytes(data)
            written.append(out)
        elif low.endswith(".zip") or low.endswith((".tar", ".tar.gz", ".tgz")):
            written.extend(_extract_cifs(io.BytesIO(data), dest, prefix + os.path.splitext(name)[0] + "/",
                                         verbose=False))

    src = archive
    if isinstance(src, (str, Path)):
        src = Path(src)
        is_zip = zipfile.is_zipfile(src)
    else:
        is_zip = zipfile.is_zipfile(src)
        src.seek(0)
    if is_zip:
        with zipfile.ZipFile(src) as z:
            members = [m for m in z.infolist() if not m.is_dir()]
            wanted = [m for m in members if m.filename.lower().endswith((".cif", ".zip", ".tar", ".tar.gz", ".tgz"))]
            # one top folder around everything (e.g. 'CoREMOF2024DB/') adds nothing to the names - drop it
            tops = {m.filename.replace("\\", "/").split("/")[0] for m in wanted}
            strip = (next(iter(tops)) + "/") if len(tops) == 1 and all("/" in m.filename.replace("\\", "/") for m in wanted) else ""
            for k, mem in enumerate(wanted):
                fn = mem.filename.replace("\\", "/")
                handle(fn[len(strip):] if strip and fn.startswith(strip) else fn, z.read(mem))
                if verbose and k and k % 5000 == 0:
                    print(f"  ... {k}/{len(members)} archive members read, {len(written)} CIFs so far")
    else:
        mode = "r:*"
        with (tarfile.open(src, mode) if isinstance(src, Path) else tarfile.open(fileobj=src, mode=mode)) as t:
            for mem in t:
                if mem.isfile() and mem.name.lower().endswith((".cif", ".zip", ".tar", ".tar.gz", ".tgz")):
                    handle(mem.name, t.extractfile(mem).read())
    return written


def _download(url, dest_file, verbose=True, chunk=1 << 20):
    """Stream a (large) file to disk with progress; skips the download if it is already there."""
    dest_file = Path(dest_file)
    if dest_file.exists() and dest_file.stat().st_size > 0:
        if verbose:
            print(f"  using the file already downloaded: {dest_file}")
        return dest_file
    dest_file.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest_file.with_suffix(dest_file.suffix + ".part")
    with requests.get(url, stream=True, timeout=120) as r:
        r.raise_for_status()
        total = int(r.headers.get("content-length", 0))
        done, next_report = 0, 0.1
        with open(tmp, "wb") as fh:
            for block in r.iter_content(chunk):
                fh.write(block)
                done += len(block)
                if verbose and total and done / total >= next_report:
                    print(f"  downloaded {done / 1e6:.0f} / {total / 1e6:.0f} MB")
                    next_report += 0.1
    tmp.replace(dest_file)
    return dest_file


def _figshare_zip_url(article_id, prefer="qmof_database"):
    """Download URL of the main zip in a figshare record (the one whose name contains `prefer`,
    else the largest zip)."""
    r = _http_get(f"https://api.figshare.com/v2/articles/{article_id}/files", timeout=60)
    files = [f for f in r.json() if str(f.get("name", "")).lower().endswith(".zip")]
    if not files:
        raise RuntimeError(f"figshare record {article_id} lists no zip file")
    pick = [f for f in files if prefer in f["name"].lower()] or files
    best = max(pick, key=lambda f: f.get("size", 0))
    return best["download_url"], best["name"], best.get("size", 0)


def index_collection(name, library_dir, verbose=True):
    """(Re)build library/<name>/_index.json from the CIFs in that folder. Run it after copying
    CIFs into a manual collection (CSD, Published) - load_collection also does it automatically
    when it sees new or removed files."""
    import json

    folder = Path(library_dir) / name
    folder.mkdir(parents=True, exist_ok=True)
    meta = COLLECTIONS.get(name, {})
    old = {}
    if (folder / _INDEX).exists():
        try:
            old = {e["file"]: e for e in json.loads((folder / _INDEX).read_text(encoding="utf-8"))["entries"]}
        except Exception:
            old = {}
    entries, bad = [], 0
    cifs = sorted(folder.glob("*.cif"))
    for k, cif in enumerate(cifs):
        if cif.name in old and old[cif.name].get("mtime") == int(cif.stat().st_mtime):
            entries.append(old[cif.name])
            continue
        try:
            info = _cif_quick_info(cif.read_text(encoding="utf-8", errors="replace"))
        except Exception:
            bad += 1
            continue
        parts = cif.stem.split("__")
        entries.append({"file": cif.name, "id": parts[-1], "subset": "/".join(parts[:-1]),
                        "mtime": int(cif.stat().st_mtime), **info})
        if verbose and k and k % 2000 == 0:
            print(f"  indexed {k}/{len(cifs)}")
    (folder / _INDEX).write_text(json.dumps({"collection": name, "kind": meta.get("kind", "experimental"),
                                             "indexed": time.strftime("%Y-%m-%d %H:%M"),
                                             "entries": entries}), encoding="utf-8")
    if verbose:
        print(f"{name}: {len(entries)} CIFs indexed" + (f", {bad} unreadable" if bad else ""))
    return entries


def import_collection(name, library_dir, source=None, download=True, keep_archive=True, verbose=True):
    """Put a structure collection into the library, once. Afterwards it works offline.

      name    : "CoRE_MOF", "QMOF", "CSD", "Published", or any name of your own
      source  : None  -> download it (CoRE_MOF from Zenodo, QMOF from figshare)
                a .zip / .tar.gz you downloaded yourself (nested archives are fine)
                a folder of CIFs (e.g. files saved from the CCDC or a paper's SI) - copied in
                a URL of a zip
    Every CIF is copied to library/<name>/ (sub-folders become part of the name, e.g.
    'CR__ASR__ABAVIJ.cif', so the subset stays visible) and indexed.
    """
    import shutil

    lib = Path(library_dir)
    dest = lib / name
    dest.mkdir(parents=True, exist_ok=True)
    meta = COLLECTIONS.get(name, {})

    if source is None:
        if meta.get("manual"):
            if verbose:
                print(f"{name}: copy your CIFs into {dest} (or pass source=<folder or zip>), then it is indexed.")
            return index_collection(name, lib, verbose)
        if not download:
            raise ValueError(f"{name}: no source given and download=False")
        if "download" in meta:
            url, fname = meta["download"], meta["download"].split("/")[-1].split("?")[0]
        elif "figshare_article" in meta:
            url, fname, size = _figshare_zip_url(meta["figshare_article"])
            if verbose and size:
                print(f"{name}: {fname} is {size / 1e6:.0f} MB")
        else:
            raise ValueError(f"{name}: unknown collection - pass source=<zip, folder or URL>")
        if verbose:
            print(f"{name}: downloading {url}")
        source = _download(url, lib / "_downloads" / fname, verbose)
    elif isinstance(source, str) and source.lower().startswith(("http://", "https://")):
        fname = source.split("/")[-1].split("?")[0] or f"{name}.zip"
        source = _download(source, lib / "_downloads" / fname, verbose)

    source = Path(source)
    if source.is_dir():
        n = 0
        for cif in source.rglob("*.cif"):
            shutil.copy2(cif, dest / _safe_name(str(cif.relative_to(source))))
            n += 1
        if verbose:
            print(f"{name}: copied {n} CIFs from {source}")
    else:
        if verbose:
            print(f"{name}: extracting CIFs from {source.name} ...")
        got = _extract_cifs(source, dest, verbose=verbose)
        if verbose:
            print(f"{name}: {len(got)} CIFs extracted to {dest}")
        if not keep_archive and source.parent == lib / "_downloads":
            source.unlink()
    return index_collection(name, lib, verbose)


def collection_status(library_dir):
    """Which collections are in the library, how many structures each, experimental or DFT."""
    import json

    rows = []
    lib = Path(library_dir)
    for name in sorted(set(COLLECTIONS) | {p.parent.name for p in lib.glob(f"*/{_INDEX}")}):
        f = lib / name / _INDEX
        n, when = 0, ""
        if f.exists():
            try:
                d = json.loads(f.read_text(encoding="utf-8"))
                n, when = len(d["entries"]), d.get("indexed", "")
            except Exception:
                pass
        rows.append({"collection": name, "structures": n, "indexed": when,
                     "kind": COLLECTIONS.get(name, {}).get("kind", "experimental"),
                     "what": COLLECTIONS.get(name, {}).get("title", "your own collection")})
    return pd.DataFrame(rows)


def load_collection(name, library_dir, allowed_elements=None, max_sites=None, require=None, verbose=True):
    """Structures of one collection whose elements are all in allowed_elements (and include every
    element in `require`, if given). Only the matching CIFs are parsed. max_sites skips huge
    cells (number of atom-site rows in the CIF) if simulation gets slow."""
    import json

    lib = Path(library_dir)
    _auto_cache_dir(lib)
    folder = lib / name
    if not folder.exists():
        return []
    f = folder / _INDEX
    n_cif = sum(1 for _ in folder.glob("*.cif"))
    entries = None
    if f.exists():
        try:
            entries = json.loads(f.read_text(encoding="utf-8"))["entries"]
        except Exception:
            entries = None
    if entries is None or len(entries) != n_cif:
        entries = index_collection(name, lib, verbose)
    allowed = set(allowed_elements) if allowed_elements else None
    need = set(require) if require else set()
    hits = [e for e in entries if e["elements"]
            and (allowed is None or set(e["elements"]) <= allowed)
            and need <= set(e["elements"])
            and (max_sites is None or (e.get("n_sites") or 0) <= max_sites)]
    kind = COLLECTIONS.get(name, {}).get("kind", "experimental")
    refs, bad = [], 0
    for e in hits:
        path = folder / e["file"]
        try:
            st, sg = load_cif_cached(path)
        except Exception:
            bad += 1
            continue
        rid = f"{name}:{Path(e['file']).stem}"
        refs.append(_record(name, rid, f"{name} {e.get('id', '')} {e.get('subset', '')}".strip(), st, path,
                            kind=kind, subset=e.get("subset", ""), spacegroup=sg))
    if verbose:
        print(f"  {name}: {len(hits)} of {len(entries)} structures match the elements"
              + (f" ({bad} could not be parsed)" if bad else ""))
    return refs


# ---------------------------------------------------------------------------
# Reference library: use what you already have, search databases only when needed
# ---------------------------------------------------------------------------
def _db_on(databases, name, default=True):
    """Is database `name` switched on? databases = {"COD": True, "MP": False, ...}; names not
    listed use `default`."""
    if not databases:
        return default
    return bool(databases.get(name, default))


def load_library(library_dir, allowed_elements=None, verbose=True, databases=None):
    """Load every CIF already in a reference library folder (from earlier COD / Materials
    Project searches, any project) - no internet needed.

    databases : switches, e.g. {"COD": True, "MP": False, "local": True, "CoRE_MOF": True,
                "QMOF": False, "CSD": True, "Published": True}; anything not listed is on.
                Collections (library/<name>/ with an _index.json) are loaded through their
                index, so only structures with your elements are parsed.

    library_dir/COD/COD_<id>.cif   -> ref_id "cod-<id>"
    library_dir/MP/<mp-id>.cif     -> ref_id "<mp-id>" (E above hull restored from the saved searches)
    library_dir/local/*.cif        -> your own CIFs (optional)
    allowed_elements : keep only structures whose elements are all in this set
                       (e.g. the elements of your sample), so unrelated chemistry is skipped.
    """
    import json

    lib = Path(library_dir)
    _auto_cache_dir(lib)
    allowed = set(allowed_elements) if allowed_elements else None
    mp_meta = {}
    for f in (lib / "MP").glob("_mp_search_*.json"):
        try:
            for row in json.loads(f.read_text(encoding="utf-8"))["entries"]:
                mp_meta[row["material_id"]] = row
        except Exception:
            pass

    refs, skipped, bad = [], 0, 0
    for sub, source in (("COD", "COD"), ("MP", "MP"), ("local", "local")):
        if not _db_on(databases, sub):
            continue
        for cif in sorted((lib / sub).glob("*.cif")):
            try:
                st, sg = load_cif_cached(cif)
            except Exception:
                bad += 1
                continue
            if allowed is not None and not {e.symbol for e in st.composition.elements} <= allowed:
                skipped += 1
                continue
            if source == "COD":
                rid = "cod-" + cif.stem.split("_", 1)[-1]
                refs.append(_record("COD", rid, f"COD {rid[4:]}", st, cif, kind="experimental", spacegroup=sg))
            elif source == "MP":
                meta = mp_meta.get(cif.stem, {})
                refs.append(_record("MP", cif.stem, f"{cif.stem} {meta.get('formula_pretty', '')}".strip(), st, cif,
                                    e_above_hull=meta.get("e_above_hull"), theoretical=meta.get("theoretical"),
                                    kind="DFT", spacegroup=sg))
            else:
                refs.append(_record("local", cif.stem, cif.stem, st, cif, kind="local", spacegroup=sg))
    if verbose:
        print(f"Library {lib}: {len(refs)} structures loaded"
              + (f", {skipped} skipped (other elements)" if skipped else "")
              + (f", {bad} unreadable" if bad else ""))
    # offline collections (CoRE MOF, QMOF, CSD downloads, published CIFs ...)
    names = sorted(set(n for n in COLLECTIONS if (lib / n).is_dir()) | {p.parent.name for p in lib.glob(f"*/{_INDEX}")})
    for name in names:
        if not _db_on(databases, name):
            if verbose:
                print(f"  {name}: switched off")
            continue
        refs += load_collection(name, lib, allowed_elements, verbose=verbose)
    return refs


def good_match(ranking, min_score=0.85, min_explained_pct=80.0, max_cell_shift_pct=1.0):
    """Is the best reference a convincing main phase?
    - pattern score >= min_score, and it explains >= min_explained_pct % of the peak intensity
    - its cell needed <= max_cell_shift_pct % stretching: right structure type but clearly different
      composition (e.g. LiCoO2 for an NMC, +1 %) is NOT good enough - go and look for a closer one."""
    if ranking is None or ranking.empty or "score" not in ranking:
        return False
    top = ranking.iloc[0]
    return bool(top.score >= min_score and (top.explained_pct or 0) >= min_explained_pct
                and abs(top.cell_scale_pct) <= max_cell_shift_pct)


def _related_systems(elements, size):
    """Smaller chemical systems worth searching for a sample: keep the anion (O, S, ...) and,
    if the sample contains a metal, at least one metal - so for a MOF (Cu-C-H-O) we search
    Cu-C-O, Cu-H-O, Cu-O ... but NOT C-H-O, which is every organic molecule in COD."""
    from itertools import combinations

    from pymatgen.core import Element

    els = list(elements)
    anions = {"O", "S", "Se", "F", "Cl", "N", "P"} & set(els)
    metals = {e for e in els if Element(e).is_metal or Element(e).is_metalloid}
    out = []
    for c in combinations(els, size):
        c = set(c)
        if anions and not anions & c:
            continue
        if metals and not metals & c:
            continue
        out.append([e for e in els if e in c])
    return out


def auto_reference_search(elements, x, y, library_dir, wavelength="CuKa", mp_api_key=None,
                          extra_refs=(), min_score=0.85, min_explained_pct=80.0, max_cell_shift_pct=1.0,
                          related=True, max_e_above_hull=0.10, systems=None, mixture=True,
                          min_mixture_main_score=0.5, databases=None, verbose=True):
    """Find references for a scan, going online only when the library is not good enough.

      stage 1  library only (offline): every CIF already downloaded, filtered to your elements
      stage 2  COD, exact element set, asked again (refresh) to catch new entries
      stage 3  Materials Project for the exact system, then COD + MP for related, smaller systems
               level by level (4-, 3-, 2-element subsets that keep the anion: for Li-Ni-Mn-Co-O
               that is Li-Ni-Mn-O ..., then end members Li-Ni-O, Li-Co-O ..., then Ni-O ...)
      stage 4  COD "at least these elements" (also structures with extra elements, e.g. dopants)

    It stops at the first stage whose best match has score >= min_score, explains
    >= min_explained_pct % of your peak intensity and needs <= max_cell_shift_pct % cell stretch. Everything downloaded goes into the library,
    so the next sample starts with more. extra_refs (e.g. your generated models) take part in
    every stage.

    Mixtures (several phases):
      systems  : the chemical systems you expect, e.g. [["Al","O"], ["Ca","F"], ["Zn","O"]] for
                 corundum + fluorite + zincite. Any system with no structure in the library is
                 fetched before stage 1, and stage 2/3 search these systems instead of relying
                 on the whole element list (COD has nothing containing Al-Ca-F-Zn-O together).
      mixture  : when no single phase is good enough, identify the phases (as identify_phases)
                 and accept the stage if together they explain >= min_explained_pct % of the peak
                 intensity, each within max_cell_shift_pct, the main one scoring
                 >= min_mixture_main_score.

    databases : switch sources on/off, e.g. {"COD": True, "MP": False, "CoRE_MOF": True,
                "QMOF": False, "CSD": True, "Published": True, "local": True}. COD and MP off
                means neither their library folders nor the online search are used; the
                collections are offline (import them once with import_collection).

    Returns (refs, ranking, log): all references considered, their ranking vs the scan, and a
    small table saying what each stage found.
    """
    import json
    from itertools import combinations

    lib = Path(library_dir)
    cod_dir, mp_dir = lib / "COD", lib / "MP"
    allowed = set(elements)
    if systems and all(isinstance(e, str) for e in systems):
        # ["Cu", "C", "H", "O"] means ONE system; a list of systems is [["Cu","C","H","O"], ["Cu","O"]]
        systems = [list(systems)]
    for sysm in systems or []:
        bad = set(sysm) - allowed
        if bad:
            raise ValueError(f"SYSTEMS entry {sysm} has elements not in ELEMENTS: {sorted(bad)}")
    log = []

    def evaluate(stage, refs):
        uniq = list({r["ref_id"]: r for r in list(refs) + list(extra_refs)}.values())
        uniq = [r for r in uniq if r.get("structure") is not None]
        if not uniq:
            log.append({"stage": stage, "n_refs": 0, "good_enough": False})
            if verbose:
                print(f"[{stage}] no references yet")
            return uniq, pd.DataFrame(), False
        rk = score_against_experiment(uniq, x, y, wavelength, progress=False)
        ok = good_match(rk, min_score, min_explained_pct, max_cell_shift_pct)
        top = rk.iloc[0]
        entry = {"stage": stage, "n_refs": len(uniq), "best": top.ref_id, "formula": top.formula,
                 "spacegroup": top.spacegroup, "score": top.score, "cell_scale_pct": top.cell_scale_pct,
                 "explained_pct": top.explained_pct, "good_enough": ok}
        msg = (f"score {top.score:.3f}, cell {top.cell_scale_pct:+.2f} %, explains {top.explained_pct:.0f}%")
        if not ok and mixture and top.score >= min_mixture_main_score:
            # several phases together?
            phases, rep = identify_phases(uniq, x, y, wavelength, verbose=False)
            if len(phases) > 1:
                m = _Matcher(np.asarray(x, float), np.asarray(y, float), wavelength, None, 0.15, 1.5, 0.05, None)
                left = unexplained_peaks(phases, x, y, wavelength)
                total = m.peaks.rel_intensity.sum()
                tot_expl = 100 * (1 - left.rel_intensity.sum() / total) if total else 0.0
                shifts = [abs(top.cell_scale_pct)] + [abs(v) for v in rep.cell_scale_pct.dropna()]
                ok = tot_expl >= min_explained_pct and max(shifts) <= max_cell_shift_pct
                entry.update(phases=", ".join(f"{p['formula']} ({p['ref_id']})" for p in phases),
                             mixture_explained_pct=round(tot_expl, 1), good_enough=ok)
                msg += f"; as a mixture of {len(phases)} phases explains {tot_expl:.0f}%"
        log.append(entry)
        if verbose:
            print(f"[{stage}] {len(uniq)} refs -> best {top.ref_id} {top.formula} [{top.spacegroup}] "
                  f"{msg} -> {'GOOD ENOUGH, stop' if ok else 'not good enough'}")
        return uniq, rk, ok

    # queries that failed (server down, connection reset ...) are remembered in the library and
    # tried again at the start of the next run - a failed search must not silently mean "nothing there"
    failed_file = lib / "_failed_queries.json"
    try:
        failed = json.loads(failed_file.read_text(encoding="utf-8")) if failed_file.exists() else {}
    except Exception:
        failed = {}
    failed_now = []

    def _save_failed():
        lib.mkdir(parents=True, exist_ok=True)
        failed_file.write_text(json.dumps(failed, indent=1), encoding="utf-8")

    def _note(db, els, exact, ex=None):
        key = f"{db}:{'-'.join(els)}:{'exact' if exact else 'atleast'}"
        if ex is None:
            if failed.pop(key, None) is not None:
                _save_failed()
            return
        failed[key] = {"db": db, "elements": list(els), "exact": exact, "error": str(ex)[:200],
                       "when": time.strftime("%Y-%m-%d %H:%M")}
        failed_now.append(key)
        _save_failed()

    def cod(els, exact=True, refresh=False):
        if not _db_on(databases, "COD"):
            return []
        try:
            out = [r for r in search_cod(els, out_dir=cod_dir, exact=exact, refresh=refresh) if r["structure"] is not None]
        except Exception as ex:
            print(f"  COD {'-'.join(els)} FAILED (will be retried next run): {str(ex)[:80]}")
            _note("COD", els, exact, ex)
            return []
        _note("COD", els, exact)
        return out

    def mp(els):
        if not _db_on(databases, "MP"):
            return []
        if not (mp_api_key or os.environ.get("MP_API_KEY") or os.environ.get("PMG_MAPI_KEY")):
            return []
        try:
            out = search_mp(els, out_dir=mp_dir, api_key=mp_api_key, max_e_above_hull=max_e_above_hull)
        except Exception as ex:
            print(f"  MP {'-'.join(els)} FAILED (will be retried next run): {str(ex)[:80]}")
            _note("MP", els, True, ex)
            return []
        _note("MP", els, True)
        return out

    def _finish(refs, rk):
        out_log = pd.DataFrame(log)
        if failed_now:
            msg = ("WARNING - these database queries failed, so the result may be incomplete "
                   "(a phase can be missed because its reference was never downloaded): "
                   + ", ".join(sorted(set(failed_now))) + ". Run again later; they are retried first.")
            print(msg)
            out_log.attrs["failed_queries"] = sorted(set(failed_now))
        return refs, rk, out_log

    # stage 1: what we already have (expected systems missing from the library are fetched first)
    if verbose and databases:
        print("Databases: " + ", ".join(f"{k} {'on' if v else 'OFF'}" for k, v in databases.items()))
    if verbose:
        _w = peak_width(x, y)
        if _w > 0.15 + 1e-9:
            print(f"Broad peaks: matching with the measured peak width {_w:.2f} deg (X.PEAK_WIDTH to override)")
    refs = load_library(lib, allowed_elements=allowed, verbose=verbose, databases=databases)
    retry = [v for v in failed.values() if set(v["elements"]) <= allowed and _db_on(databases, v["db"])]
    if retry and verbose:
        print(f"  retrying {len(retry)} query(ies) that failed last time: "
              + ", ".join(f"{v['db']} {'-'.join(v['elements'])}" for v in retry))
    for v in retry:
        refs += cod(v["elements"], exact=v["exact"]) if v["db"] == "COD" else mp(v["elements"])
    for sysm in systems or []:
        have = any({e.symbol for e in r["structure"].composition.elements} == set(sysm) for r in refs)
        if not have:
            if not (_db_on(databases, "COD") or _db_on(databases, "MP")):
                if verbose:
                    print(f"  library has no {'-'.join(sysm)} structure (COD and MP are off, so nothing is fetched)")
                continue
            if verbose:
                print(f"  library has no {'-'.join(sysm)} structure yet - fetching it")
            refs += cod(list(sysm))
            refs += mp(list(sysm))
    refs, rk, ok = evaluate("1 library (offline)", refs)
    if ok:
        return _finish(refs, rk)

    # stage 2: ask COD again for exactly these elements (and for each expected system)
    refs += cod(elements, exact=True, refresh=True)
    for sysm in systems or []:
        refs += cod(list(sysm), exact=True, refresh=True)
    refs, rk, ok = evaluate("2 COD exact, refreshed", refs)
    if ok:
        return _finish(refs, rk)

    # stage 3: Materials Project for the exact system, then related (smaller) systems level by
    # level - 4-element, 3-element (end members such as LiNiO2, LiCoO2), 2-element - each keeping
    # the anion; stop as soon as something fits
    if related:
        refs += mp(elements)
        for sysm in systems or []:
            refs += mp(list(sysm))
        refs, rk, ok = evaluate("3a Materials Project, exact system", refs)
        if ok:
            return _finish(refs, rk)
        for size in range(len(elements) - 1, 1, -1):
            subs = _related_systems(elements, size)
            if not subs:
                continue
            for sub in subs:
                refs += cod(sub)
                refs += mp(sub)
                time.sleep(1.0)
            refs, rk, ok = evaluate(f"3b related {size}-element systems ({len(subs)})", refs)
            if ok:
                return _finish(refs, rk)

    # stage 4: structures that contain these elements plus others (dopants, mixed anions ...)
    refs += cod(elements, exact=False)
    refs, rk, ok = evaluate("4 COD 'at least these elements'", refs)
    if not ok and verbose:
        print("No convincing match in the databases. The phase may be new or not in COD/MP: "
              "build a model (ASE / X.build_...), add your own CIFs to <library>/local, or check "
              "the elements and the wavelength setting.")
    return _finish(refs, rk)



_ACS_NEW = {"ja": "jacs", "am": "acsami", "ic": "acs.inorgchem", "cm": "acs.chemmater", "cg": "acs.cgd",
            "jo": "acs.joc", "om": "acs.organomet", "ma": "acs.macromol", "nl": "acs.nanolett"}


def doi_from_name(name):
    """Probable DOI hidden in a CIF file name from a paper's supporting information (as used by
    CoRE MOF and many SI downloads): 'ja4015666_si_004' -> 10.1021/ja4015666,
    'am8b02015_si_009' -> 10.1021/acsami.8b02015, 'c6ce00407e1' -> 10.1039/c6ce00407e,
    'S0885715619000228sup002' -> 10.1017/S0885715619000228. '' if no pattern fits."""
    s = str(name)
    for part in re.split(r"__|/|\\", s):
        m = re.match(r"^(ja|cm|cg|ic|jo|om|ma|nl)(\d{6,7}[a-z]?)(?![a-z0-9])", part)
        if m:
            return f"10.1021/{m.group(1)}{m.group(2)}"
        m = re.match(r"^(ja|am|ic|cm|cg|jo|om|ma|nl)(\d{1,2}b\d{5})", part)
        if m and m.group(1) in _ACS_NEW:
            return f"10.1021/{_ACS_NEW[m.group(1)]}.{m.group(2)}"
        m = re.match(r"^([b-d]\d[a-z]{2}\d{5}[a-z])", part)
        if m:
            return f"10.1039/{m.group(1)}"
        m = re.match(r"^(S\d{16})", part)
        if m:
            return f"10.1017/{m.group(1)}"
    return ""


def _source_link(r):
    src, rid = r.get("source", ""), str(r.get("ref_id", ""))
    if src == "COD":
        return f"https://www.crystallography.net/cod/{rid.split('-', 1)[-1]}.html"
    if src == "MP":
        return f"https://next-gen.materialsproject.org/materials/{rid}"
    if src == "CoRE_MOF":
        return COLLECTIONS["CoRE_MOF"]["page"]
    if src == "QMOF":
        return COLLECTIONS["QMOF"]["page"]
    return ""


def save_reference_cifs(ranking, refs, out_dir, n_top=10, phases=(), folder="reference_cifs"):
    """Copy the CIFs of the identified phases and of the n_top best-ranked references next to the
    results (out_dir/reference_cifs/), with sources.csv saying where each one comes from:
    database, its page (COD / Materials Project), the original file in the library, the CSD
    refcode or QMOF id where the name holds one, and a probable DOI read from SI file names
    (check it - it is a guess from the name). Returns the sources table."""
    import shutil

    out = Path(out_dir) / folder
    out.mkdir(parents=True, exist_ok=True)
    by_id = {r["ref_id"]: r for r in refs}
    for p in phases:
        by_id.setdefault(p["ref_id"], p)
    rank_of = {rid: i + 1 for i, rid in enumerate(ranking.ref_id)} if ranking is not None and len(ranking) else {}
    row_of = {r.ref_id: r for r in ranking.itertuples()} if ranking is not None and len(ranking) else {}
    chosen = [p["ref_id"] for p in phases]
    for rid in (list(ranking.ref_id) if ranking is not None and len(ranking) else []):
        if len(chosen) >= len(phases) + n_top:
            break
        if rid not in chosen:
            chosen.append(rid)
    phase_ids = {p["ref_id"]: i + 1 for i, p in enumerate(phases)}
    rows = []
    for k, rid in enumerate(chosen, 1):
        r = by_id.get(rid)
        if r is None or r.get("structure") is None:
            continue
        fname = f"{k:02d}_" + re.sub(r"[^A-Za-z0-9._+-]", "_", rid) + ".cif"
        src_file = r.get("cif_path")
        try:
            if src_file and Path(src_file).exists():
                shutil.copy2(src_file, out / fname)
            else:
                write_cif(r["structure"], out / fname)
        except Exception as ex:
            fname = f"(not written: {str(ex)[:60]})"
        stem = Path(src_file).stem if src_file else rid
        refcode = ""
        m = re.search(r"core_([A-Z]{6}(?:\d{2})?)", stem)
        if m:
            refcode = m.group(1)
        qmof = re.search(r"(qmof-[0-9a-f]{7})", stem)
        doi_cif = ""
        try:
            txt = Path(src_file).read_text(encoding="utf-8", errors="replace") if src_file and Path(src_file).exists() else ""
            m = re.search(r"^\s*_(?:journal_paper_doi|citation_doi|journal_DOI)\s+['\"]?([^\s'\"]+)", txt, re.M | re.I)
            if m and m.group(1) not in ("?", "."):
                doi_cif = m.group(1)
        except Exception:
            pass
        rk = row_of.get(rid)
        rows.append({
            "file": fname,
            "role": (f"identified phase {phase_ids[rid]}" if rid in phase_ids else "candidate"),
            "rank": rank_of.get(rid), "ref_id": rid, "source": r.get("source"),
            "kind": r.get("kind") or {"COD": "experimental", "MP": "DFT", "built": "model"}.get(r.get("source"), ""),
            "formula": r.get("formula"), "spacegroup": r.get("spacegroup"),
            "score": getattr(rk, "score", None), "cell_scale_pct": getattr(rk, "cell_scale_pct", None),
            "explained_pct": getattr(rk, "explained_pct", None),
            "link": _source_link(r),
            "doi_in_cif": doi_cif,
            "probable_doi": "" if doi_cif else doi_from_name(stem),
            "csd_refcode": refcode, "qmof_id": qmof.group(1) if qmof else "",
            "original_file": str(src_file or ""),
        })
    table = pd.DataFrame(rows)
    table.to_csv(out / "sources.csv", index=False)
    print(f"{len(table)} CIFs saved to {out} (see sources.csv)")
    return table


def add_to_library(sources, library_dir):
    """Copy CIFs and saved searches from earlier result folders into one shared library, once.
    Finds (recursively) COD_*.cif and _cod_search_*.json -> library/COD,
    mp-*.cif and _mp_search_*.json -> library/MP. Existing files are not overwritten."""
    import shutil

    lib = Path(library_dir)
    (lib / "COD").mkdir(parents=True, exist_ok=True)
    (lib / "MP").mkdir(parents=True, exist_ok=True)
    n = {"COD": 0, "MP": 0}
    for src in ([sources] if isinstance(sources, (str, Path)) else sources):
        src = Path(src)
        if not src.exists():
            print(f"  (skipped {src}: not found)")
            continue
        for f in src.rglob("*"):
            if not f.is_file():
                continue
            if f.name.startswith("COD_") and f.suffix == ".cif" or f.name.startswith("_cod_search_"):
                dest = lib / "COD" / f.name
            elif f.name.startswith("mp-") and f.suffix == ".cif" or f.name.startswith("_mp_search_"):
                dest = lib / "MP" / f.name
            else:
                continue
            if lib in dest.parents and dest.resolve() == f.resolve():
                continue
            if not dest.exists():
                shutil.copy2(f, dest)
                n[dest.parent.name] += 1
    print(f"Library {lib}: added {n['COD']} COD files and {n['MP']} MP files")
