#!/usr/bin/env python3
# ============================================================
# OpenMC depletion post-processing (NO PLOTS)
#   - EFPD
#   - Burnup [GWd/MTU]  (user-defined, explicit; U-only basis)
#   - keff ± sigma
#   - beta_eff (IFP total) ± sigma
# ============================================================

import numpy as np
import glob
import re
import xml.etree.ElementTree as ET
import openmc
import openmc.deplete
from openmc.data import atomic_mass

# ============================================================
# User-fixed constants (MATCH INPUT)
# ============================================================

POWER_W = 20.0  # total power (W), constant

# Identify the depleted fuel material in materials.xml
FUEL_MAT_NAME = "NaCl_UCl3"   # change if your materials.xml uses a different name
# Alternatively, use id:
FUEL_MAT_ID = None            # e.g., 1

# ============================================================
# 1) Initial uranium mass (MTU) from materials.xml (U-only)
# ============================================================

def uranium_mass_MTU_from_materials_xml(
    materials_xml="materials.xml",
    material_name=None,
    material_id=None
) -> float:
    tree = ET.parse(materials_xml)
    root = tree.getroot()

    target = None
    for m in root.findall("material"):
        if material_name is not None and m.get("name") == material_name:
            target = m
            break
        if material_id is not None and m.get("id") == str(material_id):
            target = m
            break

    if target is None:
        raise ValueError("Target material not found in materials.xml")

    # volume is required for depletion normalization (and for MTU)
    vol_attr = target.get("volume")
    if vol_attr is None:
        raise ValueError("Target material has no 'volume' attribute in materials.xml")
    volume_cm3 = float(vol_attr)

    dens = target.find("density")
    if dens is None:
        raise ValueError("Target material has no <density> in materials.xml")

    rho_val = float(dens.get("value"))
    rho_units = dens.get("units", None)

    # OpenMC normally writes g/cm3 when you set_density("g/cm3", ...)
    if rho_units not in (None, "g/cm3"):
        raise ValueError(f"Unsupported density units in materials.xml: {rho_units}")

    total_mass_g = rho_val * volume_cm3

    # Convert ao/wo to a consistent mass fraction using atomic masses
    total_mass_units = 0.0
    uranium_mass_units = 0.0

    for nuc in target.findall("nuclide"):
        name = nuc.get("name")
        ao = nuc.get("ao")
        wo = nuc.get("wo")

        if ao is None and wo is None:
            raise ValueError(f"Nuclide {name} has neither 'ao' nor 'wo' in materials.xml")

        if ao is not None:
            frac = float(ao)
            m_unit = frac * atomic_mass(name)
        else:
            frac = float(wo)
            m_unit = frac  # already weight-like

        total_mass_units += m_unit
        if name.startswith("U"):
            uranium_mass_units += m_unit

    if total_mass_units <= 0.0:
        raise ValueError("Computed total_mass_units <= 0 from materials.xml")

    u_mass_frac = uranium_mass_units / total_mass_units
    u_mass_g = total_mass_g * u_mass_frac
    u_mass_MTU = u_mass_g / 1.0e6  # g -> metric tons uranium

    return u_mass_MTU

u_mass_MTU = uranium_mass_MTU_from_materials_xml(
    materials_xml="materials.xml",
    material_name=FUEL_MAT_NAME,
    material_id=FUEL_MAT_ID
)

# print(f"[INFO] Initial U mass = {u_mass_MTU:.6e} MTU (U-only, from materials.xml)")

# ============================================================
# 2) Load depletion results
# ============================================================

res = openmc.deplete.Results("depletion_results.h5")

time_sec, keff_pair = res.get_keff()  # times [s] by default
time_sec = np.asarray(time_sec)
keff_pair = np.asarray(keff_pair)

keff   = keff_pair[:, 0]
keff_s = keff_pair[:, 1]

efpd = time_sec / 86400.0

# cumulative energy -> burnup (GWd/MTU)
time_day = time_sec / 86400.0
burnup = (POWER_W * time_day) / 1.0e9 / u_mass_MTU

# ============================================================
# 3) Statepoint files (numeric order)
# ============================================================

def idx(f):
    m = re.search(r"n(\d+)", f)
    return int(m.group(1)) if m else -1

sp_files = sorted(glob.glob("openmc_simulation_n*.h5"), key=idx)

# ============================================================
# 4) beta_eff (IFP total)
# ============================================================

beta_eff   = []
beta_eff_s = []
efpd_ifp   = []
burn_ifp   = []
keff_ifp   = []
keff_s_ifp = []

N = len(efpd)  # only use steps that exist in depletion_results

for i, f in enumerate(sp_files[:N]):
    with openmc.StatePoint(f) as sp:
        t  = sp.get_tally(name="ifp_kinetics")
        df = t.get_pandas_dataframe()

    Sb = df[df["score"] == "ifp-beta-numerator"]
    Sd = df[df["score"] == "ifp-denominator"]

    Sb_m = float(Sb["mean"].values[0])
    Sb_s = float(Sb["std. dev."].values[0])
    Sd_m = float(Sd["mean"].values[0])
    Sd_s = float(Sd["std. dev."].values[0])

    beta = Sb_m / Sd_m
    beta_s = beta * np.sqrt(
        (Sb_s / Sb_m)**2 +
        (Sd_s / Sd_m)**2
    )

    beta_eff.append(beta)
    beta_eff_s.append(beta_s)

    efpd_ifp.append(efpd[i])
    burn_ifp.append(burnup[i])
    keff_ifp.append(keff[i])
    keff_s_ifp.append(keff_s[i])

# ============================================================
# 5) Print table
# ============================================================

print("\nEFPD       keff      dkeff      beta_eff     dbeta_eff")
print("-------------------------------------------------------------------------")

for t, k, ks, be, bes in zip(
    efpd_ifp, keff_ifp, keff_s_ifp, beta_eff, beta_eff_s
):
    print(f"{t:6.1f}   {k:.6f}  {ks:.6f}  {be:.6e}  {bes:.6e}")