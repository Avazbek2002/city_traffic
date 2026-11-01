#!/usr/bin/env bash
# Full pipeline, step by step. Put the raw inputs in data/inputs/ (see README), then: bash run_pipeline.sh
# Every step can be run on its own; each script prints its usage at the top of the file.
set -euo pipefail
IN=${IN:-data/inputs}; W=${W:-data/work}; mkdir -p "$W"
PBF=${PBF:-$(ls "$IN"/*.pbf | head -1)}
PLACES=$IN/combined_places_2gis.geojson
DIST=$IN/tashkent_districts_master.geojson
MAH=$IN/tashkent_mahallas_master.geojson
SDMX=$IN/sdmx                                   # sdmx_data_{248,589,600,607,617,626}.csv
OTM=$IN/OTMdagi_talabalar_soni_hududlar_kesimida.xlsx

python3 pipeline/build_attraction.py "$PLACES" "$DIST" "$W/attr"                                        # 1 classify 2GIS places, weights per category
python3 pipeline/build_mahalla_zones.py "$W/attr/places_assigned.csv" "$MAH" "$W/mahalla_zones.csv"    # 2 one row per mahalla (OD zone)
python3 pipeline/build_od.py "$W/mahalla_zones.csv" "$W/od"                                             # 3 first OD matrices (work, health, retail; proxy times)
python3 pipeline/build_student_production.py "$W/mahalla_zones.csv" "$SDMX" "$W/mahalla_student_production.csv"   # 4 age-structured production
python3 pipeline/build_university_attraction.py "$OTM" "$W/mahalla_zones_places.csv" "$W"               # 5 universities -> mahallas
python3 pipeline/build_student_od.py "$W/mahalla_zones.csv" "$W/mahalla_student_production.csv" "$W/mahalla_zones_places.csv" "$W/mahalla_university_students.csv" "$W/students"   # 6 school/college/university OD + pa_vectors.npz
python3 pipeline/build_graph.py "$PBF" "$MAH" "$W/road_edges.csv.gz"                                    # 7 drivable road graph from the OSM extract
python3 pipeline/build_travel_times.py "$W/road_edges.csv.gz" "$W/mahalla_zones.csv" "$W"               # 8 free-flow zone-to-zone times (Dijkstra)
python3 pipeline/congestion_loop.py "$W/road_edges.csv.gz" "$W/mahalla_zones.csv" "$W/students/pa_vectors.npz" "$W/loop" am0809    # 9 one-hour equilibrium (am0809 or am0708)
python3 pipeline/congestion_slices.py "$W/road_edges.csv.gz" "$W/mahalla_zones.csv" "$W/students/pa_vectors.npz" "$W/slices"        # 10 sixteen 15-minute slices (~15 min)
python3 viewer/prep_viewer_data.py "$W/slices/slices_flow.npz" "$W/slices/slice_stats.csv" "$DIST" "$W/viewer_data.json"             # 11 pack for the viewer
python3 pipeline/build_trips.py "$W/road_edges.csv.gz" "$W/mahalla_zones.csv" "$W/students/pa_vectors.npz" "$W/slices/slices_flow.npz" "$W/trips.json"   # 12 sample individual cars (origin -> destination + route); same env vars as step 10
python3 viewer/assemble_viewer.py viewer/viewer_template.html "$W/viewer_data.json" viewer/index.html "$W/trips.json"                    # 13 build the page
python3 pipeline/build_od_graph.py "$W/road_edges.csv.gz" "$W/mahalla_zones.csv" "$W/students/pa_vectors.npz" "$W/slices/slices_flow.npz" "$DIST" "$W/od_graph.json"   # 14 zone graph data (people per district / mahalla over time)
python3 viewer/assemble_viewer.py viewer/od_graph_template.html "$W/od_graph.json" viewer/od_graph.html                                                          # 15 build the graph page
echo "Done. Open viewer/index.html (roads) and viewer/od_graph.html (district / mahalla graph) in a browser."
