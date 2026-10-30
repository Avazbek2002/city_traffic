# Assumptions

Everything below is a placeholder to be replaced by data. "Where" gives the file and variable to change.

| Area | Assumption | Where |
|---|---|---|
| Work demand | one commuter per job (`employment` column, meaning unconfirmed: jobs in the mahalla or employed residents) | `build_od.py`, `congestion_*.py` |
| Commuter mode | 35% drive, 1.15 people per car | `CAR_WORK`, `OCC_WORK` |
| Student mode | school 10% car (drop-off, both legs), college 2%, university 2%; the rest walk or use public transport | `build_student_od.py` `MODE`, `congestion_*.py` `CAR` |
| School production | ages 6-15 plus half of ages 16-17 | `build_student_od.py` `SHARE_1617_SCHOOL` |
| College production | half of ages 16-17 plus 15% of ages 18-19 | `SHARE_1617_SCHOOL`, `COLLEGE_1819` |
| School attraction | district pupil totals (2023) spread over that district's 2GIS schools | `PUPILS_D` |
| College attraction | 2GIS college listings (only about 15) | `build_student_od.py` |
| University | 354k daily students; 23% live in Tashkent, the rest assumed half near campus and half spread like residents | `UNIV_TOTAL`, `UNIV_RES_SHARE`, `NONRES_DORM` |
| University matching | OTM names matched to 2GIS by hand; 16 universities have no location and are dropped | `build_university_attraction.py` `MATCH`, `UNLOC` |
| Attraction weights | per 2GIS category (government 50, school 1, bazaar 100, ...) | `build_attraction.py` |
| Age bands | inferred from the numbers: 6-7, 8-15, 16-17, 18-19, 20-24 | `build_student_production.py` |
| District mapping | polygon `d_id` to official district matched by population | `build_student_production.py` `DID` |
| Trip timing | share of each purpose's trips per 15-minute slice | `congestion_slices.py` `PROFILE` |
| Distance decay | exponential, per-minute beta by purpose | `BETA_*` |
| Lane counts | 15% of road length tagged in OSM; otherwise trunk 3, primary 2, secondary 2, others 1 per direction | `build_graph.py`, `LANES_FILE` |
| Speeds | OSM `maxspeed` on 44% of length; otherwise 80/60/50/45/40/30/30 km/h by class | `build_graph.py` `SPEED` |
| Capacity | veh/h per lane by class (motorway 1800 ... residential 500), links x0.8 | `congestion_*.py` `CAP_LANE` |
| Delay | BPR, a=0.15, b=4; no signals or turn delays | `BPR_A`, `BPR_B` |
| Route choice | all drivers take the least-time path, averaged over iterations (MSA) | `congestion_*.py` |
| Intrazonal trips | not loaded on the network | `congestion_*.py` |
