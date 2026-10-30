"""University attraction per mahalla: students per university (OTM xlsx) matched to 2GIS higher_ed listings.
usage: python3 build_university_attraction.py OTM.xlsx mahalla_zones_places.csv OUT_DIR
MATCH maps OTM code -> regex on 2GIS brand. TASHKENT = OTM codes judged to be located in Tashkent city (manual; review!)."""
import sys, re, pandas as pd
xl, pk, out = sys.argv[1:4]
MATCH={302:'узбекского языка и литературы',306:'Международная исламская',314:'Национальный университет Республики Узбекистан',
315:'технический университет им. Ислама',324:'государственный экономический',325:'хореографии',326:'мировых языков',
327:'востоковедения',328:'архитектурно-строительный',329:'Бехзода',331:'текстильной',334:'государственный транспортный|^Ташкентский транспортный',
337:'химико-технологический',339:'Европейский Медицинский',350:'мировой экономики и дипломатии',355:'педагогический университет Узбекистана',
369:'Ташкентский фармацевтический',377:'искусств и культуры',380:'Аль-Хорезми',384:'юридический университет',395:'журналистики',
400:'правоохранительных',407:'TGFU',414:'Юнуса Раджабий',415:'международный университет Кимё',416:'Central Asian university',
417:'Губкина',418:'Ломоносова',419:'Плеханова',420:'Туринский',421:'Вестминстерский',422:'Сингапурский',423:'Инха',424:'Bucheon',
425:'Мифи',426:'Аджу',427:'Amity',429:'Webster',430:'МЭИ',431:'Менделеева',433:'экономики и технологий, ООО Toshkent',
435:'социально-инновационный',436:'IT Park',437:'University of Education',439:'фармацевтического образования',441:'Profi',
442:'Oriental',443:'эстрадного',444:'Diplomat',448:'Pisa',449:'менеджмента и экономики',457:'Alfraganus',461:'ТМС',463:'Vatan',
465:'Герасимова',470:'Yangi Asr',472:'USAT',475:'Renessans',476:'Perfect',485:'точных и социальных',489:'Mehnat',490:'Stars',
491:'British',493:'Nordic',496:'Турецкий',498:'Japan digital',499:'Afrasiab',500:'UDEA',501:'Toshkent Amaliy',506:'PDP',507:'^Team',
509:'International school of finance',511:'University of Business and Science',513:'Millat',514:'New Uzbekistan',
516:'финансового управления',517:'^Tashkent International University$',520:'исламский институт',522:'Банковско-финансовая',
524:'Oxus',530:'Tashkent metropolitan',533:'Сарбонский',535:'American university',547:'государственный медицинский',
555:'Академия государственного управления',562:'Высшая налоговая',}
UNLOC=[332,376,378,403,434,471,478,481,525,529,545,554,556,566,567,387,554]  # Tashkent-based by my reading, no 2GIS match
TASHKENT=sorted(set(MATCH)|set(UNLOC))
d=pd.read_excel(xl); d=d[d['OTM kodi'].astype(str).str.isdigit()].copy(); d['code']=d['OTM kodi'].astype(int)
d=d[d.code.isin(TASHKENT)]
h=(pd.read_csv(pk) if str(pk).endswith('.csv') else pd.read_pickle(pk)); h=h[h.subtype=='higher_ed'].copy()
h['b']=h.brand.fillna(h.name).str.replace(r',\s*(ООО|университет|институт|офис|филиал.*)$','',regex=True).str.strip()
rows=[];alloc=[]
for _,r in d.iterrows():
    rx=MATCH.get(r.code); g=h[h.b.str.contains(rx,regex=True,na=False)] if rx else h.iloc[0:0]
    g=g[g.mahalla_id.notna()]
    rows.append(dict(code=r.code,name=r.OTM,students=r.Jami,from_tashkent_city=r['Toshkent shahri'],n_listings=len(g),brands=';'.join(sorted(g.b.unique()))[:120]))
    for _,x in g.iterrows(): alloc.append(dict(code=r.code,mahalla_id=x.mahalla_id,students=r.Jami/len(g)))
u=pd.DataFrame(rows); u.to_csv(f'{out}/university_match.csv',index=False)
a=pd.DataFrame(alloc).groupby('mahalla_id').students.sum().rename('univ_students').reset_index()
a.to_csv(f'{out}/mahalla_university_students.csv',index=False)
print(len(u),'universities; students',int(u.students.sum()),'; located',int(u[u.n_listings>0].students.sum()),'; unlocated',int(u[u.n_listings==0].students.sum()))
print('students from Tashkent city in these:',int(u.from_tashkent_city.sum()))
print(u.sort_values('students',ascending=False).head(12)[['name','students','n_listings']].to_string())
print(u[u.n_listings==0].sort_values('students',ascending=False)[['code','name','students']].head(15).to_string())
print(len(a),'mahallas get university students; top:',a.sort_values('univ_students').tail(5).to_dict('records'))
