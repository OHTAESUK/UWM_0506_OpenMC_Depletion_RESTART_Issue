#!/usr/bin/env python3

# ============================================================
# OpenMC Depletion Restart Script (SAFE DAY-BASED CONTINUATION)
#
# Strategy:
#   - Read EXACT previous elapsed times from depletion_results.h5
#   - Reconstruct previous incremental day steps
#   - Append additional day steps manually
#   - Use continue_timesteps=True
#   - Use CECM integrator
#
# Notes:
#   - This is a day-based continuation strategy.
#   - It avoids trying to recover burnup steps from Results.
#   - It reuses the exact previous time history stored in results.
# ============================================================

import os
import warnings

os.environ["OPENMC_CROSS_SECTIONS"] = "/nuclear-data/endfb71/cross_sections.xml"

import openmc
import openmc.deplete

warnings.filterwarnings("ignore")

# -------------------------------
# 1) Load model
# -------------------------------
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

# -------------------------------
# 2) Load previous depletion results
# -------------------------------
prev_results = openmc.deplete.Results("depletion_results.h5")

print("==============================================")
print("[INFO] Restart mode activated")
print("==============================================")

# -------------------------------
# 3) Recover EXACT previous elapsed times [days]
# -------------------------------
prev_times = prev_results.get_times("d")   # cumulative times including t=0

if len(prev_times) < 2:
    raise RuntimeError(
        "depletion_results.h5 does not contain enough depletion steps to restart."
    )

# Remove initial t=0
prev_times = prev_times[1:]

# Convert cumulative -> incremental
prev_steps = [prev_times[0]] + [
    prev_times[i] - prev_times[i - 1]
    for i in range(1, len(prev_times))
]

n_prev = len(prev_steps)

print(f"[INFO] Completed steps = {n_prev}")
print(f"[INFO] Previous final time = {prev_times[-1]} d")

# -------------------------------
# 4) Append additional day steps
# -------------------------------
# Choose how many additional depletion steps to run
N_EXTRA = 11

# Reuse the last previous dt
last_dt = prev_steps[-1]
extra_steps = [last_dt] * N_EXTRA

time_steps = prev_steps + extra_steps

print(f"[INFO] Last previous dt = {last_dt} d")
print(f"[INFO] Extra steps added = {N_EXTRA}")
print(f"[INFO] Total schedule length = {len(time_steps)}")

# -------------------------------
# 5) Power (must match original run)
# -------------------------------
power_density = 20.0  # W/gHM

hm_mass = sum(mat.fissionable_mass for mat in materials)
TOTAL_POWER = power_density * hm_mass

print(f"[INFO] Power density = {power_density} W/gHM")
print(f"[INFO] Total power   = {TOTAL_POWER} W")

# -------------------------------
# 6) Operator
# -------------------------------
CHAIN_FILE = "/home/toh8/OpenMC/nuclear-data/chains/chain_endfb71_sfr_Lithium6_v2.xml"

operator = openmc.deplete.CoupledOperator(
    model,
    chain_file=CHAIN_FILE,
    normalization_mode="fission-q",
    prev_results=prev_results
)

# -------------------------------
# 7) Integrator (CECM)
# -------------------------------
integrator = openmc.deplete.CECMIntegrator(
    operator,
    time_steps,
    power=TOTAL_POWER,
    timestep_units="d",
    solver="cram48",
    continue_timesteps=True
)

# -------------------------------
# 8) Run
# -------------------------------
print("==============================================")
print("[INFO] Starting restart with CECMIntegrator...")
print("==============================================")

integrator.integrate(
    final_step=True,
    output=True,
    path="depletion_results_restart.h5"
)

print("==============================================")
print("[INFO] Restart completed")
print("[INFO] Output: depletion_results_restart.h5")
print("==============================================")