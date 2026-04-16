#!/usr/bin/env python3

# ============================================================
# OpenMC Depletion Restart Script (ACTUAL FINAL FIX)
#
# Strategy:
#   - Use EXACT previous timesteps (no reconstruction)
#   - Append additional steps manually
#   - Avoid floating-point mismatch
# ============================================================

import os
import sys
import warnings

os.environ["OPENMC_CROSS_SECTIONS"] = "/nuclear-data/endfb71/cross_sections.xml"

import openmc
import openmc.deplete

warnings.filterwarnings("ignore")

# -------------------------------
# Load model
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
# Load previous results
# -------------------------------
prev_results = openmc.deplete.Results("depletion_results.h5")

print("[INFO] Restart mode activated")

# -------------------------------
# 🔥 EXACT previous timesteps
# -------------------------------
prev_times = prev_results.get_times("d")  # cumulative
prev_times = prev_times[1:]               # remove t=0

# convert to incremental
prev_steps = [prev_times[0]] + [
    prev_times[i] - prev_times[i-1]
    for i in range(1, len(prev_times))
]

n_prev = len(prev_steps)

print(f"[INFO] Completed steps = {n_prev}")

# -------------------------------
# 🔥 ADD NEW STEPS HERE
# -------------------------------
# 원하는 만큼 추가
N_EXTRA = 10

last_dt = prev_steps[-1]
extra_steps = [last_dt] * N_EXTRA

time_steps = prev_steps + extra_steps

print(f"[INFO] Total steps = {len(time_steps)}")

# -------------------------------
# Power (must match original)
# -------------------------------
power_density = 20.0  # W/gHM

hm_mass = sum(mat.fissionable_mass for mat in materials)
TOTAL_POWER = power_density * hm_mass

# -------------------------------
# Operator
# -------------------------------
CHAIN_FILE = "/home/toh8/OpenMC/nuclear-data/chains/chain_endfb71_sfr_Lithium6_v2.xml"

operator = openmc.deplete.CoupledOperator(
    model,
    chain_file=CHAIN_FILE,
    normalization_mode="fission-q",
    prev_results=prev_results
)

# -------------------------------
# Integrator
# -------------------------------
integrator = openmc.deplete.CELIIntegrator(
    operator,
    time_steps,
    power=TOTAL_POWER,
    timestep_units="d",
    solver="cram48",
    continue_timesteps=True
)

# -------------------------------
# Run
# -------------------------------
print("[INFO] Starting restart...")

integrator.integrate(
    final_step=True,
    output=True,
    path="depletion_results_restart.h5"
)

print("[INFO] Restart completed")