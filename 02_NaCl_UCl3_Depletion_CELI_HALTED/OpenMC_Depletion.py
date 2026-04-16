#!/usr/bin/env python3
# ============================================================
# OpenMC Depletion Script
#
# CONTROL MODE (DOCUMENTED, VERIFIED):
#   - Burnup-driven depletion using timestep_units = "MWd/kg"
#     (OpenMC CECMIntegrator official API)
#   - Power density specified directly via power_density [W/gHM]
#
# SERPENT CORRESPONDENCE (EXPLICIT):
#   dep butot   -> cumulative burnup [GWd/MTU]
#   set powdens -> power_density [W/gHM]
#
# IMPORTANT (DOCUMENTED BEHAVIOR):
#   - OpenMC timesteps are INCREMENTAL, NOT cumulative
#   - GWd/MTU == MWd/kg (numerically identical)
#
# PHYSICS-IMPACTING CHANGES (LABELED):
#   - Burnup control via MWd/kg
#   - Power normalization via power_density
# ============================================================

import os
# ============================================================
# 0) Cross section path (explicit, no assumptions)
# ============================================================
os.environ["OPENMC_CROSS_SECTIONS"] = "/nuclear-data/endfb71/cross_sections.xml"

import openmc
import openmc.deplete

# ============================================================
# 1) Load steady-state model (user-provided XMLs)
# ============================================================

materials = openmc.Materials.from_xml("materials.xml")
geometry  = openmc.Geometry.from_xml("geometry.xml")
settings  = openmc.Settings.from_xml("settings.xml")
tallies   = openmc.Tallies.from_xml("tallies.xml")

model = openmc.Model(
    materials=materials,
    geometry=geometry,
    settings=settings,
    tallies=tallies
)

# ============================================================
# 2) Serpent dep butot (CUMULATIVE) [GWd/MTU]
# ============================================================

burnup_cumulative = [
      5.0,  10.0,  20.0,  30.0,  40.0,
     50.0,  60.0,  70.0,  80.0,  90.0,
    100.0, 110.0, 120.0, 130.0, 140.0,
    150.0, 160.0, 170.0, 180.0, 190.0,
    200.0, 210.0, 220.0, 230.0, 240.0,
    250.0
]

# ============================================================
# 3) Convert cumulative -> INCREMENTAL (REQUIRED BY OpenMC)
# ============================================================
# VERIFIED FROM DOC:
#   "timesteps values are NOT cumulative"

burnup_steps = [burnup_cumulative[0]] + [
    burnup_cumulative[i] - burnup_cumulative[i - 1]
    for i in range(1, len(burnup_cumulative))
]

# Units:
#   GWd/MTU == MWd/kg  (definition-based identity)

# ============================================================
# 4) Power density (Serpent: set powdens)
# ============================================================
# DOCUMENTED API:
#   power_density is in [W / gHM]

POWER_DENSITY_W_PER_GHM = 20.0   # [PHYSICS-IMPACTING CHANGE]

# ============================================================
# 5) Depletion operator
# ============================================================

CHAIN_FILE = "/home/toh8/OpenMC/nuclear-data/chains/chain_endfb71_sfr_Lithium6_v2.xml"

operator = openmc.deplete.CoupledOperator(
    model,
    chain_file=CHAIN_FILE,
    normalization_mode="fission-q"
)

# ============================================================
# 6) Integrator (burnup-controlled, DOCUMENTED MODE)
# ============================================================

integrator = openmc.deplete.CELIIntegrator(
    operator=operator,
    timesteps=burnup_steps,        # INCREMENTAL
    timestep_units="MWd/kg",       # burnup unit (DOCUMENTED)
    power_density=POWER_DENSITY_W_PER_GHM,
    solver="cram48"
)

# ============================================================
# 7) Run depletion
# ============================================================

print("==============================================")
print("[INFO] OpenMC burnup-driven depletion")
print("[INFO] timestep_units = MWd/kg (burnup control)")
print(f"[INFO] power_density = {POWER_DENSITY_W_PER_GHM} W/gHM")
print(f"[INFO] number of steps = {len(burnup_steps)}")
print("==============================================")

integrator.integrate(
    path="depletion_results.h5",
    write_rates=False
)

print("==============================================")
print("[INFO] Depletion complete")
print("[INFO] Output: depletion_results.h5")
print("==============================================")
