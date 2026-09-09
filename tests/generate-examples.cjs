/* Development helper only: node tests/generate-examples.cjs. No runtime dependency. */
'use strict';
const fs=require('node:fs');
const path=require('node:path');
const dir=path.join(__dirname,'..','examples');
fs.mkdirSync(dir,{recursive:true});
const tone=t=>Math.sin(2*Math.PI*2*t)+.5*Math.sin(2*Math.PI*10*t)+.2*Math.sin(2*Math.PI*25*t);
const n=4000,dt=.001;
function write(name,header,row,delimiter=','){
 const lines=[header.join(delimiter)];
 for(let i=0;i<n;i++) lines.push(row(i).map(x=>x.toPrecision(15)).join(delimiter));
 fs.writeFileSync(path.join(dir,name),lines.join('\n')+'\n');
}
write('three-tones.csv',['time','value'],i=>[i*dt,tone(i*dt)]);
write('dc-offset.csv',['time','value'],i=>[i*dt,3+tone(i*dt)]);
write('linear-trend.csv',['time','value'],i=>[i*dt,2+.8*i*dt+tone(i*dt)]);
write('nonuniform.csv',['time','value'],i=>{const t=i*dt+.00012*Math.sin(i*.7);return[t,tone(t)];});
write('spectral-leakage.csv',['time','value'],i=>[i*dt,Math.sin(2*Math.PI*2.13*i*dt)]);
write('chirp.csv',['time','value'],i=>{const t=i*dt;return[t,Math.sin(2*Math.PI*(5*t+.5*18.75*t*t))];});
write('multicolumn.tsv',['time','ax','ay','az'],i=>{const t=i*dt;return[t,Math.sin(2*Math.PI*2*t),.5*Math.sin(2*Math.PI*10*t),.2*Math.sin(2*Math.PI*25*t)];},'\t');
console.log('Generated 7 deterministic example files (4000 rows each).');
