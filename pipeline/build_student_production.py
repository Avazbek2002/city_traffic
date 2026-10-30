"""Age-structured production per mahalla.
usage: python3 build_student_production.py mahalla_zones.csv DEMO_DIR OUT.csv
DEMO_DIR holds sdmx_data_{248,589,600,607,617,626}.csv (SDMX Tashkent-city district rows 1726262..1726294).
Age bands INFERRED by cohort tracking (confirm!): 589=6-7, 600=8-15, 607=16-17, 617=18-19, 626=20-24.
d_id -> official district matched by population (2026)."""
import sys, pandas as pd
mz, D, out = sys.argv[1:4]
BANDS={'589':'age6_7','600':'age8_15','607':'age16_17','617':'age18_19','626':'age20_24'}
DID={2601:'Mirabad',2602:'Mirzo Ulugbek',2603:'Yunusabad',2604:'Yakkasaray',2605:'Shaykhantakhur',
     2606:'Chilanzar',2607:'Sergeli',2608:'Yashnabad',2609:'Almazar',2610:'Uchtepa',2611:'Bektemir',2612:'Yangikhayot'}
def row(k):
    d=pd.read_csv(f'{D}/sdmx_data_{k}.csv').set_index('Klassifikator_en')
    return d['2026']
tot=row('248')*1000
sh=pd.DataFrame({v:row(k) for k,v in BANDS.items()})
sh=sh.div(tot,axis=0)
rows={}
for did,name in DID.items():
    rows[did]=sh.loc[name+' district']
sh=pd.DataFrame(rows).T; sh.index.name='district_id'
m=pd.read_csv(mz).merge(sh.reset_index(),on='district_id')
for b in BANDS.values(): m[b]=m.population*m[b]
m['age6_24']=m[list(BANDS.values())].sum(axis=1)
m['school_age_6_17']=m.age6_7+m.age8_15+m.age16_17
m['student_prod_share']=m.age6_24/m.age6_24.sum()
m['plain_prod_share']=m.population/m.population.sum()
m['prod_ratio_vs_plain']=m.student_prod_share/m.plain_prod_share
keep=['mahalla_id','mahalla_name','district_id','population']+list(BANDS.values())+['age6_24','school_age_6_17','student_prod_share','plain_prod_share','prod_ratio_vs_plain']
m[keep].round(5).to_csv(out,index=False)
print(m[keep[4:9]+['age6_24']].sum().round(0)); print(m.prod_ratio_vs_plain.describe())
print((sh.sum(axis=1)).round(3))
