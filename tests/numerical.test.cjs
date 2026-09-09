'use strict';
// Development-only independent verification; application needs no Node.js.
const fs=require('node:fs');
const path=require('node:path');
const vm=require('node:vm');
const assert=require('node:assert/strict');
const {performance}=require('node:perf_hooks');
const root=path.resolve(__dirname,'..');
let source,sourceName;
const htmlPath=path.join(root,'index.html');
if(fs.existsSync(htmlPath)){
 const html=fs.readFileSync(htmlPath,'utf8');
 const match=html.match(/<script\b[^>]*\bid=["']numerical-engine["'][^>]*>([\s\S]*?)<\/script>/i);
 if(match){source=match[1];sourceName='index.html#numerical-engine';}
}
if(!source){source=fs.readFileSync(path.join(root,'dev','engine.js'),'utf8');sourceName='dev/engine.js';}
const sandbox={console,Float64Array,Float32Array,Uint32Array,Int32Array,Math,Number,Array,Object,JSON,Date,Set,Map,Error,Infinity,NaN};
vm.createContext(sandbox);
vm.runInContext(source+'\n;globalThis.__engine = createEngine();',sandbox);
const E=sandbox.__engine;
const results=[];
function check(name,fn){const start=performance.now();try{const detail=fn();results.push({name,passed:true,ms:performance.now()-start,detail:detail||''});console.log('PASS '+name+(detail?' — '+detail:''));}catch(error){results.push({name,passed:false,ms:performance.now()-start,detail:error.stack});console.error('FAIL '+name+'\n'+error.stack);}}
function close(actual,expected,absolute=1e-10,relative=1e-10){assert.ok(Number.isFinite(actual),'nonfinite actual '+actual);assert.ok(Math.abs(actual-expected)<=absolute+relative*Math.abs(expected),`${actual} != ${expected}; error ${Math.abs(actual-expected)}`);}
function signal(n,fs,terms,dc=0){return Float64Array.from({length:n},(_,i)=>dc+terms.reduce((s,[a,f,p=0])=>s+a*Math.sin(2*Math.PI*f*i/fs+p),0));}
function directDFT(x){const re=[],im=[];for(let k=0;k<x.length;k++){let r=0,j=0;for(let q=0;q<x.length;q++){const a=-2*Math.PI*k*q/x.length;r+=x[q]*Math.cos(a);j+=x[q]*Math.sin(a);}re.push(r);im.push(j);}return{re,im};}
check('Independent DFT comparison: lengths 1–35, 63, 64, 65, 97',()=>{
 let error=0;for(const n of [...Array.from({length:35},(_,i)=>i+1),63,64,65,97]){
 const x=Float64Array.from({length:n},(_,i)=>Math.sin(i*.72)+.37*Math.cos(i*i*.19));const got=E.fft(x),ref=directDFT(x);
 for(let k=0;k<n;k++){error=Math.max(error,Math.abs(got.real[k]-ref.re[k]),Math.abs(got.imag[k]-ref.im[k]));close(got.real[k],ref.re[k],2e-10);close(got.imag[k],ref.im[k],2e-10);}}
 return `maximum absolute complex-component error ${error.toExponential(3)}`;
});
check('CSV, tab, whitespace, quoted headers, no header, selected columns',()=>{
 for(const delim of [',','\t',' ']){const p=E.parse('time'+delim+'ax'+delim+'ay\n0'+delim+'2'+delim+'3\n.1'+delim+'4'+delim+'5');assert.deepEqual(Array.from(p.headers),['time','ax','ay']);assert.deepEqual(Array.from(p.rows,r=>Number(r[2])),[3,5]);}
 const plain=E.parse('0 2 3\n.1 4 5');assert.equal(plain.rows.length,2);assert.equal(plain.headers[0],'Column 1');
 const quotes=E.parse('\uFEFF"time","pressure, outlet"\r\n0,2\r\n.1,3');assert.equal(quotes.headers[1],'pressure, outlet');
 assert.throws(()=>E.parse('t,x\n0,1,2'),/columns/);assert.throws(()=>E.parse('"time,x\n0,1'),/quote/);
 const p=E.parse('x,time,y\n'+Array.from({length:16},(_,i)=>`${i},${i*.01},${2*i}`).join('\n'));
 const out=E.analyze({t:p.rows.map(r=>Number(r[1])),x:p.rows.map(r=>Number(r[2])),options:{stft:false}});close(out.summary.fs,100);close(out.summary.max,30);
});
check('Sampling metadata and population statistics',()=>{
 const t=Float64Array.from({length:4000},(_,i)=>7+i*.001),q=E.inspect(t);close(q.dt,.001);close(q.fs,1000);close(q.nyquist,500);close(q.resolution,.25);close(q.duration,3.999);close(q.start,7);close(q.end,10.999);close(q.minDt,.001);close(q.maxDt,.001);
 const x=signal(4000,1000,[[1,2],[.5,10],[.2,25]],3),s=E.stats(x);close(s.mean,3);close(s.std,Math.sqrt(.645));close(s.rms,Math.sqrt(9.645));
});
check('Known three-tone amplitude and frequency: all four windows',()=>{
 const x=signal(4000,1000,[[1,2],[.5,10],[.2,25]]);for(const win of ['none','hann','hamming','blackman']){const a=E.spectrum(x,1000,win);for(const [f,v] of [[2,1],[10,.5],[25,.2]]){const k=Math.round(f/a.binWidth);close(a.frequency[k],f);close(a.values[k],v,1e-9);}}
});
check('DC and even Nyquist are not doubled; odd final bin is doubled',()=>{
 for(const n of [63,64,99,100,127,128]){const fs=n;const x=Float64Array.from({length:n},(_,i)=>2+1.7*Math.cos(2*Math.PI*Math.floor(n/2)*i/n));const a=E.spectrum(x,fs,'none');assert.equal(a.values.length,Math.floor(n/2)+1);close(a.values[0],2);close(a.values[Math.floor(n/2)],1.7,1e-9);close(a.frequency.at(-1),Math.floor(n/2));}
});
check('Mean removal and least-squares linear detrend',()=>{
 const x=Float64Array.from({length:101},(_,i)=>3+.2*i);const y=E.preprocess(x,{detrend:true});assert.ok(Math.max(...y.map(Math.abs))<1e-12);
 const tone=signal(4000,1000,[[1,2]],3);close(E.spectrum(E.preprocess(tone,{removeMean:true}),1000,'none').values[0],0);close(E.spectrum(tone,1000,'none').values[0],3);
 const mixed=Float64Array.from({length:101},(_,i)=>Math.sin(i*.2)+3+.2*i),d=E.preprocess(mixed,{detrend:true});close(d.reduce((a,b)=>a+b,0),0,1e-10);close(d.reduce((a,b,i)=>a+(i-50)*b,0),0,1e-9);
});
check('Periodic window coherent gains and correction',()=>{
 for(const n of [99,100,128])for(const [name,gain] of [['none',1],['hann',.5],['hamming',.54],['blackman',.42]]){const x=signal(n,n,[[2.3,7]]);const a=E.spectrum(x,n,name);close(a.coherentGain,gain);close(a.values[7],2.3,1e-9);}
});
check('Zero padding changes grid, preserves normalization and native resolution',()=>{
 const n=1000,fs=1000,x=Float64Array.from({length:n},(_,i)=>2*Math.cos(2*Math.PI*125*i/fs));
 for(const padding of ['none','next','double']){const a=E.spectrum(x,fs,'hann',padding);close(a.values[Math.round(125/a.binWidth)],2,1e-9);close(a.resolution,1);assert.equal(a.nfft,padding==='none'?1000:padding==='next'?1024:2048);close(a.binWidth,fs/a.nfft);}
});
check('PSD Parseval identity for odd/even lengths and every window',()=>{
 for(const n of [99,100,127,128,1000])for(const win of ['none','hann','hamming','blackman']){const x=Float64Array.from({length:n},(_,i)=>1+.8*Math.sin(i*.8)+.1*Math.cos(i*i*.03));const p=E.welch(x,1000,n,0,win);const w=E.windowValues(n,win);let weighted=0,energy=0;for(let i=0;i<n;i++){weighted+=(x[i]*w[i])**2;energy+=w[i]**2;}close(p.values.reduce((a,b)=>a+b,0)*p.binWidth,weighted/energy,1e-9);}
});
check('Welch segment averaging, overlap, tonal power, and PSD units normalization',()=>{
 const x=signal(4096,1024,[[2,64],[.5,128]],3);for(const win of ['none','hann']){const p=E.welch(x,1024,256,.5,win);assert.equal(p.segments,31);assert.equal(p.step,128);close(p.binWidth,4);close(p.values.reduce((a,b)=>a+b,0)*p.binWidth,11.125,1e-9);}
 const p=E.welch(x,1024,256,.5,'none');close(p.values[16],2/4);close(p.values[32],.125/4);close(p.values[0],9/4);
});
check('Welch arbitrary signal equals independently averaged window-weighted power',()=>{
 const n=1237,fs=500,x=Float64Array.from({length:n},(_,i)=>Math.sin(i*i*.071)+.3*Math.cos(i*.023)+.0002*i);
 for(const length of [99,128,3000])for(const overlap of [0,.5,.73]){const p=E.welch(x,fs,length,overlap,'hamming'),L=Math.min(n,length),hop=Math.max(1,Math.round(L*(1-overlap)));const w=E.windowValues(L,'hamming');let reference=0,count=0;const energy=w.reduce((a,b)=>a+b*b,0);
 for(let start=0;start+L<=n;start+=hop){let power=0;for(let i=0;i<L;i++)power+=(x[start+i]*w[i])**2;reference+=power/energy;count++;}
 close(p.values.reduce((a,b)=>a+b,0)*p.binWidth,reference/count,1e-9);assert.equal(p.segments,count);assert.equal(p.segmentLength,L);}
});
check('Finite-data and settings failures are explicit',()=>{
 const t=Float64Array.from({length:16},(_,i)=>i*.01),x=new Float64Array(16);
 for(const options of [{overlap:1},{overlap:-.1},{segmentLength:7},{padding:'invalid'},{tolerance:-1},{peakThreshold:2},{integrate:3},{filter:'bandpass',lowCut:30,highCut:20}])assert.throws(()=>E.analyze({t,x,options}));
 assert.throws(()=>E.analyze({t,x:new Float64Array(15)}),/equal length/);assert.throws(()=>E.inspect([0,1,2,3,4,5,6,Infinity]),/finite/);
});
check('Large representable PSD avoids raw FFT-square overflow; invalid fs is rejected',()=>{
 const n=10000,fs=1000,amplitude=1e152,x=signal(n,fs,[[amplitude,2]]);const p=E.welch(x,fs,n,0,'hann');
 close(p.values[20],amplitude*amplitude*n/(3*fs),0,1e-9);assert.ok(p.values.every(Number.isFinite));
 for(const bad of [0,-1,Infinity,NaN]){assert.throws(()=>E.spectrum(x,bad),/frequency/);assert.throws(()=>E.welch(x,bad),/frequency/);}
});
check('Local peaks, period conversion, DC exclusion and separation',()=>{
 const x=signal(4000,1000,[[1,2],[.5,10],[.2,25]],3),a=E.spectrum(E.preprocess(x,{removeMean:true}),1000,'hann'),p=E.welch(E.preprocess(x,{removeMean:true}),1000,1000,.5,'hann');const peaks=E.peaks(a,p,{peakThreshold:.1,peakCount:5,peakDistance:1});assert.deepEqual(Array.from(peaks,q=>q.frequency),[2,10,25]);for(let i=0;i<peaks.length;i++)close(peaks[i].period,1/[2,10,25][i]);
 const synthetic={frequency:Float64Array.from([0,1,2,3,4,5,6]),values:Float64Array.from([100,0,8,7,0,6,0])};assert.deepEqual(Array.from(E.peaks(synthetic,synthetic,{peakThreshold:0,peakDistance:4}),q=>q.frequency),[2]);
});
check('Nonuniform sampling blocks FFT; linear resampling preserves linear signal',()=>{
 const t=Float64Array.from({length:100},(_,i)=>i*.01+.0004*Math.sin(i*.9)),x=Float64Array.from(t,v=>2+3*v);assert.ok(E.inspect(t).jitter>.01);assert.throws(()=>E.analyze({t,x,options:{stft:false}}),/Non-uniform/);
 const a=E.analyze({t,x,options:{resample:true,removeMean:false,stft:false}});assert.ok(a.processedSummary.jitter<1e-10);assert.ok(a.warnings.some(s=>/interpolation/i.test(s)));for(let i=0;i<t.length;i++)close(a.x[i],2+3*a.t[i]);
 assert.throws(()=>E.inspect([0,.1,.2,.3,.3,.5,.6,.7]),/strictly/);assert.throws(()=>E.analyze({t:[0,1,2,3,4,5,6,7],x:[0,1,2,NaN,4,5,6,7]}),/finite/);
});
check('STFT chirp ridge tracks instantaneous frequency and segment centers',()=>{
 const fs=1000,n=4000,t=Float64Array.from({length:n},(_,i)=>i/fs),x=Float64Array.from(t,v=>Math.sin(2*Math.PI*(5*v+.5*18.75*v*v)));const a=E.analyze({t,x,options:{stft:true,stftLength:256,stftOverlap:.75}});assert.ok(a.stft);close(a.stft.times[0],.1275);assert.equal(a.stft.step,64);let maxError=0;
 a.stft.values.forEach((v,j)=>{let k=1;for(let i=2;i<v.length;i++)if(v[i]>v[k])k=i;const err=Math.abs(a.stft.frequency[k]-(5+18.75*a.stft.times[j]));maxError=Math.max(maxError,err);assert.ok(err<=fs/256,`ridge error ${err}`);});return `maximum chirp ridge error ${maxError.toFixed(4)} Hz (bin width 3.90625 Hz)`;
});
check('Centered FIR low/high/band filters reject stop bands with interior phase retained',()=>{
 const fs=1000,n=8192,x=signal(n,fs,[[1,40],[1,200]]);
 function projection(y,f){let s=0,c=0;for(let i=512;i<n-512;i++){s+=y[i]*Math.sin(2*Math.PI*f*i/fs);c+=y[i]*Math.cos(2*Math.PI*f*i/fs);}return {amp:2*Math.hypot(s,c)/(n-1024),phase:Math.atan2(c,s)};}
 for(const [filter,lowCut,highCut,pass,stop] of [['lowpass',0,100,40,200],['highpass',100,0,200,40],['bandpass',150,250,200,40]]){const y=E.filterSignal(x,fs,{filter,lowCut,highCut});const p=projection(y,pass),s=projection(y,stop);assert.ok(Math.abs(p.amp-1)<.01);assert.ok(s.amp<.01);assert.ok(Math.abs(p.phase)<.01);}
});
check('Trapezoidal integration: constant acceleration, zero initial conditions',()=>{
 const t=Float64Array.from({length:100},(_,i)=>i*.01),x=new Float64Array(100).fill(2);for(const integrate of [1,2]){const a=E.analyze({t,x,options:{removeMean:false,window:'none',integrate,stft:false}});for(let i=0;i<t.length;i++)close(a.processed[i],integrate===1?2*t[i]:t[i]*t[i]);assert.ok(a.warnings.some(s=>/drift/.test(s)));}
});
check('All shipped example files parse and analyze as intended',()=>{
 for(const file of fs.readdirSync(path.join(root,'examples'))){if(!/\.(csv|tsv)$/.test(file))continue;const p=E.parse(fs.readFileSync(path.join(root,'examples',file),'utf8'));assert.equal(p.rows.length,4000);const data={t:p.rows.map(r=>Number(r[0])),x:p.rows.map(r=>Number(r[1])),options:{resample:file==='nonuniform.csv',stft:file==='chirp.csv'}};const a=E.analyze(data);assert.equal(a.summary.n,4000);if(file==='three-tones.csv')assert.deepEqual(Array.from(a.peaks.slice(0,3),p=>p.frequency),[2,10,25]);}
});
const benchmarks=[];
if(!process.argv.includes('--no-benchmark')) for(const n of [10000,100000,300000]) check(`Performance ${n.toLocaleString()} samples: FFT + Welch + STFT`,()=>{
 const t=Float64Array.from({length:n},(_,i)=>i*.001),x=signal(n,1000,[[1,2],[.5,10],[.2,25]]);const start=performance.now();const a=E.analyze({t,x,options:{stft:true,stftLength:1024,stftOverlap:.5,segmentLength:1024,overlap:.5}});const ms=performance.now()-start;assert.ok(a.stft);assert.equal(a.amplitude.values.length,Math.floor(n/2)+1);benchmarks.push({n,ms,frames:a.stft.times.length});return `${ms.toFixed(1)} ms; ${a.stft.times.length} STFT frames`;
});
const failed=results.filter(r=>!r.passed);
const report=`# Numerical validation\n\nExecuted ${new Date().toISOString()} using ${process.version}, ${process.platform}/${process.arch}.\n\nSource tested: \`${sourceName}\`.\n\n${results.length-failed.length}/${results.length} checks passed. Tests use independent direct DFT references only for small vectors; production uses FFT.\n\n| Check | Result | Detail |\n|---|---|---|\n${results.map(r=>`| ${r.name} | ${r.passed?'PASS':'FAIL'} | ${String(r.detail).replace(/\n/g,' ').replace(/\|/g,' / ')} |`).join('\n')}\n\n## Performance\n\nWall-clock measurements are from this development container, not a browser guarantee. Each includes full-record arbitrary-length FFT, overlapping Welch PSD, peak detection, and STFT; timings exclude file parsing, Web Worker transfer, rendering, and user interaction.\n\n| Samples | Elapsed (ms) | STFT frames |\n|---:|---:|---:|\n${benchmarks.map(r=>`| ${r.n} | ${r.ms.toFixed(1)} | ${r.frames} |`).join('\n')}\n\n## Reproduce\n\nRun \`node tests/numerical.test.cjs\` from a checkout. Node is used only for developer verification; end users open \`index.html\` directly. Regenerate sample files using \`node tests/generate-examples.cjs\`.\n\nTests prefer the shipped \`index.html\` numerical-engine script and otherwise use \`dev/engine.js\`. GUI/browser verification is documented separately, if performed. Population standard deviation is tested; segment detrending is not assumed by Welch tests. PSD integrals use discrete bin sums multiplied by bin width (no half-weight endpoints).\n`;
fs.writeFileSync(path.join(__dirname,'VALIDATION.md'),report);
console.log(`\n${results.length-failed.length}/${results.length} checks passed; source ${sourceName}`);
process.exitCode=failed.length?1:0;
