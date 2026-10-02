// ICB-level need covariates from the OHID Fingertips API (https://fingertips.phe.org.uk/api).
// Run in a browser console on any fingertips.phe.org.uk page (same-origin fetch). It returns the
// CSV text of data/icb_covariates.csv and its SHA-256. Extracted 1 October 2026; data/icb_covariates.csv has
// SHA-256 fe78a076dbf3c5ad440502b32f64265fee96bbaf206cb4f7caf9c548fca2baca (values can change if
// OHID revises the source data).
//   Area type 221 = "ICBs, former STPs" (the 42 ICBs, April 2023 codes); rows for Persons only.
//   Age structure: indicator 93468 (proportion of GP-registered population by age group), mean of
//     2023, 2024 and 2025 (the years available for all 42 ICBs; 15+ only 2023 and 2025).
//   Deprivation: 94240 (IMD 2025 score, all 42 ICBs) and 93553 (IMD 2019 score, 29 ICBs only).
//   QOF prevalence: mean of 2021/22 to 2024/25 (depression 3 years, obesity 2023/24 only).
function parseCSV(t){const rows=[];let row=[],f='',q=false;for(let i=0;i<t.length;i++){const c=t[i];if(q){if(c=='"'){if(t[i+1]=='"'){f+='"';i++}else q=false}else f+=c}else{if(c=='"')q=true;else if(c==','){row.push(f);f=''}else if(c=='\n'){row.push(f);rows.push(row);row=[];f=''}else if(c!='\r')f+=c}}if(f.length||row.length){row.push(f);rows.push(row)}return rows}
const ages={'0-4 yrs':'age_0_4','5-14 yrs':'age_5_14','<18 yrs':'age_u18','15+ yrs':'age_15plus','65+ yrs':'age_65plus','75+ yrs':'age_75plus','85+ yrs':'age_85plus'};
const qof=[[241,'qof_dm'],[219,'qof_hyp'],[90933,'qof_ast'],[253,'qof_copd'],[273,'qof_chd'],[262,'qof_hf'],[280,'qof_af'],[258,'qof_ckd'],[848,'qof_dep'],[212,'qof_stroke'],[247,'qof_dem'],[90581,'qof_smi'],[224,'qof_epi'],[276,'qof_cancer'],[90443,'qof_osteo'],[91269,'qof_ra'],[94136,'qof_obes'],[92590,'qof_pad'],[200,'qof_ld'],[93797,'qof_ndh'],[91280,'qof_smok']];
const QP=['2021/22','2022/23','2023/24','2024/25']; const AP=['2023','2024','2025'];
const data={}; const ons={};
function add(ods,col,v){data[ods]=data[ods]||{}; data[ods][col]=data[ods][col]||[]; data[ods][col].push(v);}
async function get(id){for(let k=0;k<4;k++){const r=await fetch('/api/all_data/csv/by_indicator_id?indicator_ids='+id+'&child_area_type_id=221&parent_area_type_id=15');const rows=parseCSV(await r.text());const h=rows[0];const ix=n=>h.indexOf(n);if(r.status==200&&ix('Area Type')>=0){const out=rows.slice(1).filter(x=>x.length>5&&x[ix('Area Type')]=='ICBs'&&x[ix('Sex')]=='Persons'&&x[ix('Category Type')]=='').map(x=>({ons:x[ix('Area Code')],ods:(x[ix('Area Name')].match(/- ([A-Z0-9]{3})$/)||[])[1],age:x[ix('Age')],tp:x[ix('Time period')],v:x[ix('Value')]}));if(out.length)return out;}await new Promise(s=>setTimeout(s,1500));}throw new Error('indicator '+id+' failed');}
for(const x of await get(93468)){if(ages[x.age]&&AP.includes(x.tp)&&x.v!==''){add(x.ods,ages[x.age],+x.v);ons[x.ods]=x.ons;}}
for(const x of await get(94240)){if(x.v!==''){add(x.ods,'imd2025',+x.v);ons[x.ods]=x.ons;}}
for(const x of await get(93553)){if(x.v!==''){add(x.ods,'imd2019',+x.v);}}
for(const [id,col] of qof){for(const x of await get(id)){if(QP.includes(x.tp)&&x.v!==''){add(x.ods,col,+x.v);}}}
const cols=[...Object.values(ages),'imd2025','imd2019',...qof.map(q=>q[1])];
const odsList=Object.keys(data).sort();
let csv='ods,ons,'+cols.join(',')+'\n';
for(const o of odsList){csv+=o+','+ons[o]+','+cols.map(c=>{const a=data[o][c];return a&&a.length?(a.reduce((s,v)=>s+v,0)/a.length).toFixed(4):'';}).join(',')+'\n';}
const buf=await crypto.subtle.digest('SHA-256',new TextEncoder().encode(csv));
const sha256=[...new Uint8Array(buf)].map(b=>b.toString(16).padStart(2,'0')).join('');
({csv, sha256});
