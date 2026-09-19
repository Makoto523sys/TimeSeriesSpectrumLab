// Browser engine is the parity reference; no generated fixtures or dependencies are needed.
const fs=require('node:fs'),vm=require('node:vm'),path=require('node:path');
const html=fs.readFileSync(path.join(__dirname,'../../index.html'),'utf8');
const ctx={};vm.createContext(ctx);vm.runInContext(html.match(/<script id="numerical-engine">([\s\S]*?)<\/script>/)[1]+';this.E=createEngine();',ctx);
const cases=[
 {n:129,window:'hann',filter:'none'},
 {n:128,window:'rectangular',filter:'lowpass',highCut:18},
 {n:131,window:'blackman',filter:'highpass',lowCut:5,detrend:true},
 {n:160,window:'hamming',filter:'bandpass',lowCut:4,highCut:20,integrate:1},
 {n:128,window:'hann',filter:'none',integrate:2,removeMean:false},
 {n:129,window:'rectangular',filter:'none',resample:true,irregular:true}
];
for(let ci=0;ci<cases.length;ci++){
 const c=cases[ci],t=Array.from({length:c.n},(_,i)=>i*.01+(c.irregular?.0001*Math.sin(i):0)),x=t.map(v=>.4+Math.sin(2*Math.PI*3*v)+.3*Math.cos(2*Math.PI*12*v)+.2*v);
 const o={removeMean:true,integrate:0,padding:'double',segmentLength:33,overlap:.5,stft:true,stftLength:32,stftOverlap:.5,stftWindow:'hann',peakDistance:1,...c};
 if(!o.integrate)o.frs={min:1,max:25,count:5,spacing:'log',damping:[.005,.05,.2],steps:100,tail:2,unit:'gal'};
 const r=ctx.E.analyze({t,x,options:o});
 const data={t,x,window:o.window,filter:o.filter,low_cut:o.lowCut||1,high_cut:o.highCut||20,detrend:!!o.detrend,remove_mean:o.removeMean,integrate:o.integrate,resample:!!o.resample,processed:Array.from(r.processed),fft:Array.from(r.amplitude.values),frequency:Array.from(r.amplitude.frequency),psd:Array.from(r.psd.values),stft:r.stft.values.flatMap(v=>Array.from(v)),stft_times:Array.from(r.stft.times),peaks:r.peaks.map(p=>[p.frequency,p.period,p.amplitude,p.psd]),frs:r.frs?r.frs.curves.flatMap(c=>Array.from(r.frs.frequency,(_,j)=>[c.sa[j],c.psa[j],c.sd[j],c.sv[j],c.psv[j]])):[]};
 console.log(`[case${ci+1}]`);for(const [k,v]of Object.entries(data))console.log(`${k} = ${JSON.stringify(v)}`);
}
