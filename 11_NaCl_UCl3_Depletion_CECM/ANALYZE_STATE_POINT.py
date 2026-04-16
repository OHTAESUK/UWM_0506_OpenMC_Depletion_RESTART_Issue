#!/usr/bin/env python3
# ============================================================
# OpenMC SINGLE Statepoint Analyzer (Spyder-friendly)
#
# HOW TO USE (Spyder):
#   1) CHANGE N_STATEPOINT VALUE
#   2) Run (F5)
#
# Example:
#   N_STATEPOINT = 10
#   -> reads openmc_simulation_n10.h5 ONLY
#
# Outputs:
#   - keff
#   - Lambda_eff, beta_eff
#   - beta_i (all delayed neutron groups)
#   - neutron spectrum (energy-filtered flux tally)
#
# Author: Taesuk Oh
# ============================================================

import numpy as np
import openmc
import matplotlib.pyplot as plt

# ============================================================
# USER INPUT (<<< ONLY CHANGE THIS SECTION >>>)
# ============================================================
N_STATEPOINT = 5
PLOT_SPECTRUM = True   

# ============================================================
# Load statepoint
# ============================================================
sp_file = f"openmc_simulation_n{N_STATEPOINT}.h5"

print("\n========================================================")
print(" OpenMC Single Statepoint Analyzer (Spyder)")
print("========================================================")
print(f" Using statepoint: {sp_file}")

sp = openmc.StatePoint(sp_file)

# ============================================================
# 1) k-effective
# ============================================================
print("\n[1] k-effective")
print("--------------------------------------------------------")

try:
    k_mean = sp.keff.nominal_value
    k_std  = sp.keff.std_dev
except AttributeError:
    k_mean, k_std = sp.keff

print(f"  k-eff = {k_mean:.6f} ± {k_std:.3e}")

# ============================================================
# 2) Adjoint-weighted IFP kinetics
# ============================================================
print("\n[2] Adjoint-weighted kinetics (IFP)")
print("--------------------------------------------------------")

def get_tally_safe(sp, name):
    try:
        return sp.get_tally(name=name)
    except Exception:
        return None

t_core = get_tally_safe(sp, "ifp_kinetics")
t_grp  = get_tally_safe(sp, "ifp_beta_by_group")

if t_core is None or t_grp is None:
    print("  IFP tallies not available.")
else:
    df = t_core.get_pandas_dataframe()

    def pull(score):
        r = df[df["score"] == score]
        return float(r["mean"]), float(r["std. dev."])

    S_t, S_t_s = pull("ifp-time-numerator")
    S_b, S_b_s = pull("ifp-beta-numerator")
    S_d, S_d_s = pull("ifp-denominator")

    eps = 1e-300

    beta_eff = S_b / S_d
    beta_eff_s = beta_eff * np.sqrt(
        (S_b_s / max(S_b, eps))**2 +
        (S_d_s / max(S_d, eps))**2
    )

    Lambda = (S_t / S_d) * k_mean
    Lambda_s = Lambda * np.sqrt(
        (S_t_s / max(S_t, eps))**2 +
        (S_d_s / max(S_d, eps))**2 +
        (k_std  / max(k_mean, eps))**2
    )

    print(f"  Lambda_eff = {Lambda:.5e} ± {Lambda_s:.2e} s")
    print(f"  beta_eff  = {beta_eff:.5e} ± {beta_eff_s:.2e}")

    # --- delayed neutron group breakdown
    df_g = t_grp.get_pandas_dataframe()
    gcol = [c for c in df_g.columns if "group" in str(c).lower()][0]
    df_g = df_g[df_g["score"] == "ifp-beta-numerator"].sort_values(gcol)

    num   = df_g["mean"].to_numpy()
    num_s = df_g["std. dev."].to_numpy()

    beta_i = num / S_d
    beta_i_s = beta_i * np.sqrt(
        (num_s / np.maximum(num, eps))**2 +
        (S_d_s / max(S_d, eps))**2
    )

    print("  beta_i by delayed group:")
    for i, (b, s) in enumerate(zip(beta_i, beta_i_s), start=1):
        print(f"    beta_{i} = {b:.5e} ± {s:.2e}")

# ============================================================
# 3) Neutron spectrum
# ============================================================

print("\n[3] Neutron spectrum (energy-filtered flux tally)")
print("--------------------------------------------------------")

def find_energy_flux_tally(sp):
    fallback = None
    for t in sp.tallies.values():
        if not any(isinstance(f, openmc.EnergyFilter) for f in t.filters):
            continue
        if "flux" in [s.lower() for s in t.scores]:
            return t
        fallback = t
    return fallback

t_spec = find_energy_flux_tally(sp)

if t_spec is None:
    print("  No energy-filtered tally found.")
else:
    ef = next(f for f in t_spec.filters if isinstance(f, openmc.EnergyFilter))
    E_bins = np.array(ef.bins)

    E_lo, E_hi = E_bins[:, 0], E_bins[:, 1]
    flux = t_spec.mean.flatten()

    mask = (E_lo > 0.0) & (E_hi > E_lo)
    E_lo = E_lo[mask]
    E_hi = E_hi[mask]
    flux = flux[mask]

    E_mid = np.sqrt(E_lo * E_hi)
    dlnE  = np.log(E_hi / E_lo)
    spec  = flux / dlnE   # E*Phi(E)

    print("  Group |   E_low[eV]     E_high[eV]      E_mid[eV]     Flux        E*Phi(E)")
    print("  -----------------------------------------------------------------------------")
    for i in range(len(E_mid)):
        print(f"  {i+1:3d} | {E_lo[i]:11.3e}  {E_hi[i]:11.3e}  "
              f"{E_mid[i]:11.3e}  {flux[i]:11.3e}  {spec[i]:11.3e}")

    if PLOT_SPECTRUM:
        plt.figure(figsize=(7, 5))
        plt.loglog(E_mid, spec / np.max(spec), lw=1.8)
        plt.xlabel("Energy [eV]")
        plt.ylabel(r"Normalized $E\Phi(E)$")
        plt.grid(True, which="both", ls="--", alpha=0.6)
        plt.tight_layout()
        plt.show()

'''
# ============================================================
# 4) Isotopes — FULL depletion inventory (ZAID order, XS-only)
#     XS-only = nuclides listed in <library materials="...">
# ============================================================

print("\n[4] Isotope inventory (FULL DUMP, ZAID-sorted, XS-ONLY)")
print("--------------------------------------------------------")

import os
import re
import xml.etree.ElementTree as ET
import openmc.deplete
import openmc.data

# ------------------------------------------------------------
# Load depletion results
# ------------------------------------------------------------
res_dep = openmc.deplete.Results("depletion_results.h5")

step = N_STATEPOINT
sr = res_dep[step]

FUEL_MAT_KEY = "6"

if FUEL_MAT_KEY not in sr.index_mat:
    raise ValueError(
        f"Material key '{FUEL_MAT_KEY}' not found.\n"
        f"Available keys: {list(sr.index_mat.keys())}"
    )

# ------------------------------------------------------------
# Parse XS nuclides from cross_sections.xml (ENDF/B-VIII style)
# ------------------------------------------------------------
xs_xml = os.environ.get("OPENMC_CROSS_SECTIONS")
if xs_xml is None:
    raise RuntimeError("OPENMC_CROSS_SECTIONS is not set.")

tree = ET.parse(xs_xml)
root = tree.getroot()

xs_nuclides = set()

for lib in root.findall("library"):
    mats = lib.get("materials")
    if mats:
        for name in mats.split():
            xs_nuclides.add(name)

if not xs_nuclides:
    raise RuntimeError(
        "Failed to extract XS nuclides from <library materials='...'>.\n"
        "cross_sections.xml format not recognized."
    )

# ------------------------------------------------------------
# ZAID sorting key: (Z, A, metastable)
# ------------------------------------------------------------
def zaid_key(nuc):
    m = re.match(r"([A-Za-z]+)(\d+)(?:_m(\d+))?$", nuc)
    if not m:
        return (999, 999, 0)
    sym, A, mstate = m.groups()
    Z = openmc.data.ATOMIC_NUMBER.get(sym, 999)
    A = int(A)
    mstate = int(mstate) if mstate else 0
    return (Z, A, mstate)

# ------------------------------------------------------------
# Filter depletion nuclides by XS availability
# ------------------------------------------------------------
nuclide_names = [
    n for n in sr.index_nuc.keys()
    if n in xs_nuclides
]
nuclide_names.sort(key=zaid_key)

# ------------------------------------------------------------
# Print BUMAT-style table
# ------------------------------------------------------------
print(f"Depletion step : {step}")
print(f"Material key   : {FUEL_MAT_KEY}")
print(f"XS source      : {xs_xml}")
print("")
print("  Nuclide        Number density [atoms/b-cm]")
print("  ------------------------------------------")

total_nd = 0.0
for nuc in nuclide_names:
    nd = sr[FUEL_MAT_KEY, nuc]
    total_nd += nd
    print(f"  {nuc:<12s}  {nd:14.6e}")

print("  ------------------------------------------")
print(f"  TOTAL        {total_nd:14.6e}")

print(f"\n✔ XS-only ZAID-sorted isotope dump completed "
      f"({len(nuclide_names)} nuclides)\n")

# ============================================================
# 5) Element-summed inventory (for chemistry/materials people)
#     - Sum over isotopes (e.g., U234 + U235 + ... → U)
# ============================================================

print("\n[5] Element inventory (ISOTOPES SUMMED)")
print("--------------------------------------------------------")

from collections import defaultdict

# ------------------------------------------------------------
# Extract element symbol from nuclide name
# e.g., U235 -> U
#       Am242_m1 -> Am
# ------------------------------------------------------------
def element_from_nuclide(nuc):
    m = re.match(r"([A-Za-z]+)", nuc)
    if m:
        return m.group(1)
    return None

element_nd = defaultdict(float)

for nuc in nuclide_names:
    nd = sr[FUEL_MAT_KEY, nuc]
    elem = element_from_nuclide(nuc)
    if elem is not None:
        element_nd[elem] += nd

# Sort elements by atomic number
def element_key(elem):
    return openmc.data.ATOMIC_NUMBER.get(elem, 999)

elements_sorted = sorted(element_nd.keys(), key=element_key)

print("")
print("  Element        Number density [atoms/b-cm]")
print("  ------------------------------------------")

total_elem_nd = 0.0
for elem in elements_sorted:
    nd = element_nd[elem]
    total_elem_nd += nd
    print(f"  {elem:<10s}  {nd:14.6e}")

print("  ------------------------------------------")
print(f"  TOTAL        {total_elem_nd:14.6e}")

print(f"\n✔ Element-summed inventory completed "
      f"({len(elements_sorted)} elements)\n")

# ============================================================
# 6) Reaction rates (Serpent-style micro detector comparison)
# ============================================================

print("\n[6] Reaction rates by nuclide (Serpent-style)")
print("--------------------------------------------------------")

# ------------------------------------------------------------
# Read tallies
# ------------------------------------------------------------
def get_tally_safe(name):
    try:
        return sp.get_tally(name=name)
    except Exception:
        return None

t_flux = get_tally_safe("flux_total")
t_salt = get_tally_safe("rr_salt_MT1_MT102")
t_fuel = get_tally_safe("rr_fuel_MT1_MT18_MT102")

# ------------------------------------------------------------
# Load material number densities
# ------------------------------------------------------------
materials = openmc.Materials.from_xml("materials.xml")

fuel_mat = None
for m in materials:
    if m.name == "NaCl_UCl3":
        fuel_mat = m
        break

if fuel_mat is None:
    raise RuntimeError("Material 'NaCl_UCl3' not found in materials.xml")

# nuclide atom densities [atoms/b-cm]
atom_dens = fuel_mat.get_nuclide_atom_densities()

# ------------------------------------------------------------
# Extract global flux
# ------------------------------------------------------------
if t_flux is None:
    print("  flux_total tally not found.")
else:

    df_flux = t_flux.get_pandas_dataframe()

    flux = float(df_flux["mean"].iloc[0])
    flux_sd = float(df_flux["std. dev."].iloc[0])

    print(f"\n  Global flux = {flux:.6e} ± {flux_sd:.2e}\n")

    # --------------------------------------------------------
    # Helper to process RR tallies
    # --------------------------------------------------------
    def process_rr_tally(tally):

        df = tally.get_pandas_dataframe()

        rows = []

        for _, r in df.iterrows():

            nuc = r["nuclide"]
            score = r["score"]

            rr = float(r["mean"])
            rr_sd = float(r["std. dev."])

            # ------------------------------------------------
            # Atom density for this nuclide
            # ------------------------------------------------
            if nuc not in atom_dens:
                continue

            N = atom_dens[nuc]

            # ------------------------------------------------
            # Derived quantities
            # ------------------------------------------------

            # Effective macroscopic XS
            xs_macro = rr / flux

            # Density-removed response (microscopic-like RR)
            rr_div_N = rr / N

            # Effective microscopic XS
            xs_micro = rr / (flux * N)

            rows.append((nuc, score, rr, rr_sd, xs_macro, rr_div_N, xs_micro))

        return rows

    rows = []

    if t_salt is not None:
        rows.extend(process_rr_tally(t_salt))

    if t_fuel is not None:
        rows.extend(process_rr_tally(t_fuel))

    # --------------------------------------------------------
    # Desired print order (Serpent detector style)
    # --------------------------------------------------------
    desired_order = [
        ("Na23",  "total"),
        ("Na23",  "(n,gamma)"),
        ("Cl35",  "total"),
        ("Cl35",  "(n,gamma)"),
        ("Cl37",  "total"),
        ("Cl37",  "(n,gamma)"),
        ("U235",  "total"),
        ("U235",  "(n,gamma)"),
        ("U235",  "fission"),
        ("U238",  "total"),
        ("U238",  "(n,gamma)"),
        ("U238",  "fission"),
    ]

    row_map = {
        (nuc, score): (rr, rr_sd, xs_macro, rr_div_N, xs_micro)
        for nuc, score, rr, rr_sd, xs_macro, rr_div_N, xs_micro in rows
    }

    # --------------------------------------------------------
    # Print table
    # --------------------------------------------------------
    print("  Nuclide   Reaction      RR              RR_sd           RR/flux        RR/N           RR/(flux*N)")
    print("  ---------------------------------------------------------------------------------------------------------")

    for key in desired_order:

        if key not in row_map:
            continue

        nuc, score = key

        rr, rr_sd, xs_macro, rr_div_N, xs_micro = row_map[key]

        print(
            f"  {nuc:<8s}  {score:<10s}  "
            f"{rr:12.5e}   {rr_sd:12.3e}   "
            f"{xs_macro:12.5e}   {rr_div_N:12.5e}   {xs_micro:12.5e}"
        )

    print("  ---------------------------------------------------------------------------------------------------------")

    # --------------------------------------------------------
    # Explanation of columns
    # --------------------------------------------------------
    print("\n  Column meaning:")
    print("    RR            : OpenMC reaction rate tally  ~ ∫ φ N σ dE dV")
    print("    RR/flux       : Effective macroscopic XS    (≈ Σ_eff)")
    print("    RR/N          : Density-removed response")
    print("    RR/(flux*N)   : Effective microscopic XS    (≈ σ_eff)")

    print("\n  Note:")
    print("    RR/(flux*N) is usually the closest quantity to Serpent")
    print("    positive-MT detector responses for homogeneous systems.")
'''

print("\n✔ Analysis complete.\n")
