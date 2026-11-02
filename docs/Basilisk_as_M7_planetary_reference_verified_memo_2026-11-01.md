# Basilisk as M7 planetary reference: verified claims memo

Date: 2026-11-01, 21:15 CST
Status: the three-researcher deep-research run was stopped before any notes were written, so this is not the full report. It is a direct verification pass over the claims your agent's synthesis made (SPICE dual path, Basilisk modules, Vizard transport, real-time mode, coupling precedents, budget criteria), checked against primary documentation and papers tonight. Every row below carries its source. Rows marked "not confirmed" should not be banked into NeuralCompose memory yet.

## 1. SPICE dual path

| Claim | Verdict | Evidence |
|---|---|---|
| Basilisk `spiceInterface` reads JPL SPICE kernels and emits one `SpicePlanetStateMsgPayload` per registered body, plus `SpiceTimeMsgPayload`, spacecraft states, and attitude or translation reference messages | Confirmed | [spiceInterface module page](https://avslab.github.io/basilisk/Documentation/simulation/environment/spiceInterface/spiceInterface.html): outputs `spiceTimeOutMsg`, `planetStateOutMsgs`, `scStateOutMsgs`, `attRefStateOutMsgs`, `transRefStateOutMsgs`; bodies registered with `addPlanetNames` and `addSpacecraftNames`; kernels via `loadSpiceKernel` or `addKernelPath(s)` |
| Basilisk planet states are in SI meters, not SPICE kilometers | Confirmed | [SpicePlanetStateMsgPayload](https://hanspeterschaub.info/basilisk/Documentation/architecture/msgPayloadDefC/SpicePlanetStateMsgPayload.html): `PositionVector[3]` documented as "m", `VelocityVector[3]` as "m/s", `J2000Current` as "s" (TDB seconds past J2000), `J20002Pfix[3][3]` and `J20002Pfix_dot[3][3]` give the planet-fixed orientation and its rate |
| Default planet-fixed frames are the IAU frames | Confirmed | spiceInterface page: default planet frames are "IAU_" + planet name; optional per-channel SPICE query reconstruction with 60 s knot spacing and cubic Hermite interpolation |
| SPICE.jl is a JuliaAstro wrapper of the NAIF toolkit with CSPICE function names | Confirmed for naming, not for units | [SPICE.jl README](https://github.com/JuliaAstro/SPICE.jl): `furnsh`, `utc2et`, `spkpos` used as in CSPICE; tests derived from SpiceyPy's suite; install is `Pkg.add("SPICE")`. The README does not state units. CSPICE itself returns km and km/s, and nothing in the README claims a conversion, so treat SPICE.jl output as km until you check `spkpos` on a known body yourself |
| Consequence for M7 | Design rule, not a sourced fact | If WorldScenario mirrors Basilisk's message contract, convert SPICE.jl positions by 1e3 to meters and velocities by 1e3 to m/s at the boundary, and carry the IAU frame DCM separately. Do the conversion in exactly one function with a unit test against a Basilisk-generated state for the same epoch |

## 2. Vizard and transport

| Claim | Verdict | Evidence |
|---|---|---|
| Vizard is a Unity application | Confirmed | [Vizard overview](https://hanspeterschaub.info/basilisk/Vizard/Vizard.html): "The Vizard Unity-based Basilisk visualization is able to display in a three-dimensional view the Basilisk simulation data." License listed as Freeware |
| Transport is protobuf over ZeroMQ | Confirmed | [vizInterface module source](https://hanspeterschaub.info/bskOlderDocs/bsk_2_6_0/_sources/Documentation/simulation/vizard/vizInterface/vizInterface.rst.txt): "the serialized message protobuffer" is sent each update; BSK `zmq_bind`s the broadcast socket and `zmq_connect`s the two-way socket; handshake strings are `PING`, `SIM_UPDATE`, `OK`, `REQUEST_INPUT`, `VIZARD_INPUT`, `REQUEST_IMAGE_` |
| Three live modes with default ports | Confirmed | [Live communication page](https://hanspeterschaub.info/basilisk/Vizard/vizardAdvanced/vizardLiveComm.html): `liveStream=True` (lockstep, waits for OK each step, supports user input and opNav cameras) on `tcp://localhost:5556`; `broadcastStream=True` (one-way publish/subscribe, multiple Vizards, BSK does not wait) on `tcp://localhost:5570`; `noDisplay=True` (headless opNav, cannot combine with the other two). Legacy `opNavMode` 1 and 2 map to `liveStream` and `noDisplay` |
| Vizard can also play back saved files | Confirmed in general | vizInterface page: data can be recorded to a file for later playback or streamed live; the overview page shows a timeline slider |

Design consequence for M7: the broadcast socket on 5570 is the right tap for a read-only Julia or NeuralCompose HUD subscriber, because it never blocks the simulation. The two-way socket on 5556 is for the case where WorldScenario must feed events back. Either way the Julia side needs Basilisk's `.proto` definitions compiled for Julia (ProtoBuf.jl) and a ZMQ.jl SUB socket; nothing in the docs suggests a JSON side channel.

## 3. Real-time mode

| Claim | Verdict | Evidence |
|---|---|---|
| Basilisk has a module that slaves the simulation to wall clock with an acceleration factor | Confirmed | [simSynch module](https://www.hanspeterschaub.info/basilisk/_sources/Documentation/simulation/simSynch/simSynch.rst.txt): "The clock synchronization module is used to slave the simulation to realtime"; acceleration factor can be changed mid-run if the time-initialized flag is also reset |
| Variable names `accelFactor`, `displayTime`, `timeInitialized` | Not confirmed from the page | The documentation page defers to `simSynch.h`; I did not read the header. The older module was named `clock_synch` (bsk 1.x docs). Check the header before quoting names |
| Example scenario streaming live to Vizard | Confirmed to exist | [scenarioBasicOrbitStream](https://avslab.github.io/basilisk/examples/scenarioBasicOrbitStream.html) is the canonical live-stream example |

## 4. Coupling precedents (mission or geometry models to subsurface habitability)

| Claim | Verdict | Evidence |
|---|---|---|
| PlanetProfile couples SPICE geometry to interior structure and magnetic induction | Confirmed | [PlanetProfile on PyPI](https://pypi.org/project/PlanetProfile/): 1D self-consistent interior models for ocean worlds and rocky dwarf planets; install instructions pull `spiceypy` (conda-forge) and `MoonMag` (pip); citation Styczinski, Vance and Melwani Daswani 2023, Earth and Space Science 10(8), doi 10.1029/2022ea002748. This is the clearest published SPICE-to-interior coupling and the closest architectural analogue to what WorldScenario wants |
| Enceladus: Affholder et al. 2021, Nature Astronomy, Bayesian comparison of abiotic vs biotic CH4 against Cassini plume data | Confirmed (via Astrobites summary; the HAL PDF was bot-blocked tonight) | [Astrobites summary](https://astrobites.org/2021/06/15/bayesian-enceladus-plume/) and [HAL record](https://hal.sorbonne-universite.fr/hal-03384808): hydrothermal vent chemistry model with and without methanogenesis; 50,000 prior draws, about 32 percent habitable, split into two subsets of 8,008; observables are H2 and CH4 escape rates and the H2:CH4 ratio; abiotic-only serpentinization cannot match the observed CH4 escape rate; biotic posterior dominates once the prior probability of methanogenesis given habitability exceeds about 0.35. Caveats: only serpentinization modeled as abiotic source; primordial CH4 not modeled |
| Mars: Sauterey et al. 2022, Nature Astronomy, H2-based methanogens and early Mars habitability with global cooling feedback | Confirmed title and venue only | [Nature Astronomy s41550-022-01786-w](https://www.nature.com/articles/s41550-022-01786-w.pdf), [arXiv 2210.04948](https://arxiv.org/abs/2210.04948v1). I did not read the body tonight. Your agent's synthesis should cite this for Mars, not for Enceladus |
| Europa precedent | Not checked tonight | No source fetched. Do not bank a Europa coupling precedent until one is read |

## 5. Budget closure: what separates a habitat model from a visualization

These are the criteria the verified sources actually implement. Where your agent listed seven, I could not see its text after context compaction, so I cannot tell you which of its seven match; compare against this list.

1. Energy supply is computed from a specified catabolic reaction and the local chemistry, not assumed. NutMEG needs "only an overall catabolic reaction and some environmental information" to start ([Higgins and Cockell 2020, J. R. Soc. Interface 17:20200588](https://pmc.ncbi.nlm.nih.gov/articles/PMC7653372)).
2. Maintenance power is an explicit, temperature-dependent sink per cell. NutMEG's P_M (W per cell) rises roughly exponentially with temperature, calibrated against methanogen growth data (phymet2) and the Tijhuis et al. cell-size prediction at 300 K, with Lever et al. 2015 as the thermodynamic floor (protein replacement after 2 to 10 percent amino acid racemization). The fitted range is 288 to 371 K; below 280 K the lower bound goes to zero, which the authors flag as non-physical. So any M7 planetary run colder than 280 K is outside the calibrated range and must say so.
3. Growth is what remains after maintenance, with a maintenance efficiency parameter. The paper's results show the mechanism: energy- and nutrient-saturated methanogens from 280 to 330 K give final biomass falling from about 1e9 to 1e6 cells per litre, total CH4 from about 62 uM to 3 uM, peak growth rate rising from 0.007 to 0.14 per hour (doubling 100 h to 5 h). A visualization has no such dependence; a habitat model must reproduce it.
4. Nutrient limitation is enforced as a separate gate (CO2, H2, phosphorus in NutMEG). At 330 K, 90 percent of optimal CO2 yields little to no growth. Enceladus phosphorus may be as low as about 1e-10 M, which is exactly the kind of number WorldScenario needs to carry.
5. ATP yield is treated as an uncertainty range (0.5 to 1.5 mol ATP per mol CO2), not a constant.
6. Material closure is tested against an observable flux. Affholder et al. close the model against Cassini H2 and CH4 escape rates and their ratio, then score hypotheses by Bayes factor. For Mars the analogue would be atmospheric CH4 or H2 constraints; for Europa there is no equivalent plume flux yet.
7. Interior state is self-consistent with geometry. PlanetProfile derives ocean depth, conductivity and induction response from the same body parameters that SPICE positions, so the habitat volume is not a free parameter.

A visualization satisfies none of these; Vizard draws states it is handed. Basilisk itself satisfies the orbital and attitude dynamics side and knows nothing about chemistry. The M7 claim worth making is that WorldScenario sits between them: SPICE geometry in, NutMEG-style power and nutrient gating in the middle, and a flux observable out that can be scored.

## 6. Recommendation on banking memory

Bank only the Confirmed rows above, with their URLs, into NeuralCompose memory. Hold the SPICE.jl unit claim until you have run `spkpos` once and read the magnitude. Hold the `accelFactor` name until `simSynch.h` is read. Do not bank a Europa precedent or the "7 criteria" wording from the unverified synthesis; use the list in section 5, which is traceable.

## 7. If you want the full report later

The deep-research version (Basilisk capabilities, precedents, budget closure, each as a sourced section, then a report writer) needs about 45 to 60 minutes of agent time and would read the primary Affholder PDF, the Sauterey body, a Europa source, and the Basilisk headers. Say so and I will relaunch it; nothing from tonight's cancelled run is reusable because no notes were written.

## Sources

- [Basilisk spiceInterface module](https://avslab.github.io/basilisk/Documentation/simulation/environment/spiceInterface/spiceInterface.html)
- [SpicePlanetStateMsgPayload](https://hanspeterschaub.info/basilisk/Documentation/architecture/msgPayloadDefC/SpicePlanetStateMsgPayload.html)
- [Basilisk vizInterface module source (bsk 2.6.0)](https://hanspeterschaub.info/bskOlderDocs/bsk_2_6_0/_sources/Documentation/simulation/vizard/vizInterface/vizInterface.rst.txt)
- [Vizard live communication](https://hanspeterschaub.info/basilisk/Vizard/vizardAdvanced/vizardLiveComm.html)
- [Vizard overview](https://hanspeterschaub.info/basilisk/Vizard/Vizard.html)
- [Basilisk simSynch module](https://www.hanspeterschaub.info/basilisk/_sources/Documentation/simulation/simSynch/simSynch.rst.txt)
- [scenarioBasicOrbitStream](https://avslab.github.io/basilisk/examples/scenarioBasicOrbitStream.html)
- [SPICE.jl (JuliaAstro)](https://github.com/JuliaAstro/SPICE.jl)
- [PlanetProfile on PyPI](https://pypi.org/project/PlanetProfile/)
- [Higgins and Cockell 2020, NutMEG, J. R. Soc. Interface](https://pmc.ncbi.nlm.nih.gov/articles/PMC7653372)
- [Affholder et al. 2021, HAL record](https://hal.sorbonne-universite.fr/hal-03384808) and [Astrobites summary](https://astrobites.org/2021/06/15/bayesian-enceladus-plume/)
- [Sauterey et al. 2022, Nature Astronomy](https://www.nature.com/articles/s41550-022-01786-w.pdf) and [arXiv 2210.04948](https://arxiv.org/abs/2210.04948v1)
