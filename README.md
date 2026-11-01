# Tashkent traffic simulation

A macroscopic morning-traffic model for Tashkent with an animated map viewer.

It builds trip demand per mahalla (585 zones) from population, age structure, jobs, school and university enrollment and 2GIS place data. A gravity model distributes the trips. Every driver takes the fastest path (Dijkstra) over the OpenStreetMap road graph. Roads slow down as they fill, and the loop repeats until the times and flows agree.

**This is a model, not a measurement.** Car shares, trip timing, lane counts and road capacities are assumptions, listed in [`docs/ASSUMPTIONS.md`](docs/ASSUMPTIONS.md). Nothing is calibrated against observed traffic.

## See the animation (no setup)

Open `viewer/index.html` in a browser (it needs internet only for the fonts). It shows cars moving along the roads from 06:00 to 10:00. The model is solved every 15 minutes; the "Update" dropdown sets how often the picture refreshes (15, 5 or 1 minute). The 5 and 1-minute options blend between the 15-minute results, so they look smoother but add no new model information. Zoom with the wheel, drag to pan, drag the slider to jump in time.

## Run the whole pipeline

```bash
pip install -r requirements.txt
# put the raw inputs in data/inputs/ (see below), then:
bash run_pipeline.sh          # about 6 minutes with ITERS=1, about 25 minutes with the default iterations
```

Each step is a standalone script in `pipeline/` that prints its usage in the file header. `run_pipeline.sh` shows the order and arguments.

| Step | Script | Output |
|---|---|---|
| 1 | `build_attraction.py` | classify 2GIS places, attraction weights per category |
| 2 | `build_mahalla_zones.py` | one row per mahalla (population, jobs, attraction) |
| 3 | `build_od.py` | first OD matrices for work, health, retail (straight-line proxy times) |
| 4 | `build_student_production.py` | school-age and student population per mahalla |
| 5 | `build_university_attraction.py` | university students per mahalla |
| 6 | `build_student_od.py` | school, college and university OD, `pa_vectors.npz` |
| 7 | `build_graph.py` | drivable road graph from the OSM `.pbf` |
| 8 | `build_travel_times.py` | free-flow zone-to-zone times (Dijkstra) |
| 9 | `congestion_loop.py` | one-hour equilibrium (`am0809` or `am0708`) |
| 10 | `congestion_slices.py` | sixteen 15-minute slices, 06:00 to 10:00 |
| 11-12 | `viewer/prep_viewer_data.py`, `viewer/assemble_viewer.py` | the viewer page |

`pipeline/build_lane_checklist.py` writes a spreadsheet of the 100 busiest streets so lane counts can be checked by hand. Pass the filled file back with `LANES_FILE=...` to steps 9 and 10.

## Raw inputs (not in the repo)

Place these in `data/inputs/`:

| File | What it is |
|---|---|
| `combined_places_2gis.geojson` | 2GIS places for Tashkent (scraped, so not redistributed here) |
| `tashkent_districts_master.geojson` | 12 district polygons, `d_id` 2601-2612 |
| `tashkent_mahallas_master.geojson` | 585 mahalla polygons with population and employment |
| `sdmx/sdmx_data_{248,589,600,607,617,626}.csv` | Statistics Agency tables: total population (248), ages 6-7, 8-15, 16-17, 18-19, 20-24 (589, 600, 607, 617, 626) |
| `OTMdagi_talabalar_soni_hududlar_kesimida.xlsx` | students per university by home region |
| `*.osm.pbf` | Uzbekistan extract from Geofabrik (OpenStreetMap) |

`data/derived/` holds small outputs so the viewer and later steps work without the raw data, including `road_edges.csv.gz` (the road graph from OSM, ODbL licence).

## Parameters

Every assumption can be overridden with an environment variable. The main ones:

| Variable | Default | Meaning |
|---|---|---|
| `CAR_WORK` | 0.45 | share of commuters driving |
| `OCC_WORK` | 1.15 | people per car for commuters |
| `BETA_WORK`, `BETA_SCHOOL`, `BETA_COLLEGE`, `BETA_UNIVERSITY` | 0.08, 0.30, 0.10, 0.03 | distance decay per minute |
| `ITERS` | 15 (loop), 4 (slices) | equilibrium iterations |
| `LANES_FILE` | none | filled lane checklist |
| `UNIV_TOTAL` | 354000 | university students making a daily trip |
| `TARGET_SCHOOL`, `TARGET_COLLEGE`, `TARGET_UNIVERSITY` | 10, 20, 30 | used by `build_student_od.py` only (mean trip time that sets decay) |

## Data sources and licences

- Road network: OpenStreetMap contributors, ODbL.
- Population by district and age: Statistics Agency of Uzbekistan (SDMX tables).
- School pupils by district (1 Jan 2023) and university students (2022/23): stat.uz.
- Places: 2GIS. Check their terms before sharing the raw or place-level data.

## What it is not

The viewer is built from the macro model. It is not a microsimulation. Each 15-minute step is solved separately with no queues carried over, and the dots follow modelled flows and speeds, not individual vehicles. Buses, trucks, taxis, the metro, and cars from outside the city are not modelled.

## Trip layer (origin -> destination)

`pipeline/build_trips.py` draws individual cars from the OD matrices of every 15-minute slice (1 drawn car = `SAMPLE` real cars, default 80), gives each a departure time and its congested shortest path, and writes `trips.json`. Pass that file as the 4th argument of `viewer/assemble_viewer.py`; the viewer then offers "Trips origin -> destination": each car moves along its own route, a green ring marks where it left, a white dot where it arrived, and the district selector shows only cars leaving or arriving in one district. Run it with the same `CAR_WORK` / `OCC_WORK` / `BETA_*` / `LANES_FILE` as `congestion_slices.py`, otherwise it stops or the trips will not match the slices. Intrazonal trips (inside one mahalla) are not drawn because they never load the road graph.

## Zone graph (districts and mahallas)

`pipeline/build_od_graph.py` + `viewer/od_graph_template.html` give a second page, `viewer/od_graph.html`. Every district (or mahalla, switchable) is a circle placed at its centre; the circle area is the number of people currently in it, so it shrinks as people leave for work or school and grows where they arrive. Arcs show who leaves in the current 15 minutes. People counted: everyone with a modelled morning trip (work, school, college, university; all modes), starting at home at 06:00. Arrival uses the congested car time for everybody (an assumption: walkers and bus riders are really slower). Run it with the same `CAR_WORK` / `BETA_*` / `WORK_PROFILE` as the slices.

## Dashboard (three windows on one page)

`viewer/dashboard_template.html` puts the car-flow simulation, the district graph and the mahalla graph on one page with a shared clock (06:00-11:00; the road model itself ends at 10:00, the graphs show arrivals until 11:00). Drag the bars between windows to resize them, use the maximize button on any window, or pick another layout (roads + two graphs, side by side, stacked). Clicking a district or mahalla circle filters the car trips to that district. Build it with:

    python3 viewer/assemble_viewer.py viewer/dashboard_template.html data/work/viewer_data.json viewer/dashboard.html data/work/trips.json data/work/od_graph.json

## Traffic-signal delay (optional, off by default)

`pipeline/build_signals.py` extracts the signalised junctions and signalised pedestrian crossings from the OSM extract (about 1,900 junction nodes and 2,000 crossings in the study area). Run the slices with `SIGNALS_FILE=data/work/signals.csv` to add a fixed delay to every road section that enters such a node: `SIGNAL_DELAY` seconds at a junction (default 20) and `CROSSING_DELAY` at a crossing (default 5). Both are ASSUMPTIONS. Without the file the model has no intersection delay and an empty network averages about 58 km/h; with 40 s / 10 s it averages about 38 km/h and the busiest slice about 32 km/h. The delay is a fixed cost, so it does not tell you how much of a real slowdown comes from queues; calibrate it against measured off-peak and peak speeds.
