// Registered population and diabetes register size by ICB, from the OHID Fingertips API
// (https://fingertips.phe.org.uk/api). Run in a browser console on any fingertips.phe.org.uk page
// (same-origin fetch). Returns the CSV text of data/icb_populations.csv and its SHA-256.
// Extracted 1 October 2026; data/icb_populations.csv has SHA-256
// 5dfe6d1f875c903201d89c802ad4f7f87f3eb45c1a337a3b49b83f3c70c6bb78 (values can change if OHID revises the data).
//   Area type 221 = "ICBs, former STPs" (the 42 ICBs, April 2023 codes); rows for Persons only.
//   popYYYY: GP-registered population, the denominator of indicator 93468 (largest value across the age
//     groups of that year). Values for 2021 and 2022 differ from later years by up to 9% in some ICBs, so
//     the analysis uses the mean of 2023 to 2025.
//   dmYYYY_YY: number of patients on the QOF diabetes register (indicator 241, aged 17+), by financial year.
// The populations are used only to describe dispensing volume per registered patient (code/need_summary_numbers.py).
function parseCSV(t){const rows=[];let row=[],f='',q=false;for(let i=0;i<t.length;i++){const c=t[i];if(q){if(c=='"'){if(t[i+1]=='"'){f+='"';i++}else q=false}else f+=c}else{if(c=='"')q=true;else if(c==','){row.push(f);f=''}else if(c=='\n'){row.push(f);rows.push(row);row=[];f=''}else if(c!='\r')f+=c}}if(f.length||row.length){row.push(f);rows.push(row)}return rows}
async function get(id){for(let k=0;k<4;k++){const r=await fetch('/api/all_data/csv/by_indicator_id?indicator_ids='+id+'&child_area_type_id=221&parent_area_type_id=15');const rows=parseCSV(await r.text());const h=rows[0];const ix=n=>h.indexOf(n);if(r.status==200&&ix('Area Type')>=0){const out=rows.slice(1).filter(x=>x.length>5&&x[ix('Area Type')]=='ICBs'&&x[ix('Sex')]=='Persons'&&x[ix('Category Type')]=='').map(x=>({ons:x[ix('Area Code')],ods:(x[ix('Area Name')].match(/- ([A-Z0-9]{3})$/)||[])[1],age:x[ix('Age')],tp:x[ix('Time period')],cnt:x[ix('Count')],den:x[ix('Denominator')]}));if(out.length)return out;}await new Promise(s=>setTimeout(s,1500));}throw new Error('indicator '+id+' failed');}
const pop={}, dm={}, ons={};
for(const x of await get(93468)){ if(['2021','2022','2023','2024','2025'].includes(x.tp) && x.den!==''){ pop[x.ods]=pop[x.ods]||{}; const d=Math.round(+x.den); pop[x.ods][x.tp]=Math.max(pop[x.ods][x.tp]||0,d); ons[x.ods]=x.ons; } }
for(const x of await get(241)){ if(['2021/22','2022/23','2023/24','2024/25'].includes(x.tp) && x.cnt!==''){ dm[x.ods]=dm[x.ods]||{}; dm[x.ods][x.tp]=Math.round(+x.cnt); } }
const Y=['2021','2022','2023','2024','2025'], F=['2021/22','2022/23','2023/24','2024/25'];
let csv='ods,ons,'+Y.map(y=>'pop'+y).join(',')+','+F.map(f=>'dm'+f.replace('/','_')).join(',')+'\n';
for(const o of Object.keys(pop).sort()){ csv+=o+','+ons[o]+','+Y.map(y=>pop[o][y]??'').join(',')+','+F.map(f=>(dm[o]||{})[f]??'').join(',')+'\n'; }
const buf=await crypto.subtle.digest('SHA-256',new TextEncoder().encode(csv));
({sha256:[...new Uint8Array(buf)].map(b=>b.toString(16).padStart(2,'0')).join(''), csv});
