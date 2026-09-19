'use strict';
const fs=require('node:fs'),path=require('node:path'),vm=require('node:vm'),assert=require('node:assert/strict');
const html=fs.readFileSync(path.join(__dirname,'../index.html'),'utf8'),ctx={};
vm.createContext(ctx);vm.runInContext(html.match(/<script id="numerical-engine">([\s\S]*?)<\/script>/)[1]+';this.E=createEngine();',ctx);
const E=ctx.E,results=[],tau=2*Math.PI;
const wave=(n,dt,f)=>Float64Array.from({length:n},(_,i)=>Math.sin(tau*f*i*dt));
function near(a,b,rtol=1e-10,atol=1e-12){assert.ok(Math.abs(a-b)<=atol+rtol*Math.abs(b),`${a} != ${b}`);}
function check(name,fn){try{const detail=fn()||'';results.push({name,pass:true,detail});console.log('PASS '+name+' '+detail);}catch(e){results.push({name,pass:false,detail:e.message});console.error('FAIL '+name,e);}}
check('Zero acceleration gives zero response for every measure and damping',()=>{const r=E.responseSpectrum(new Float64Array(50),.01,{max:40,count:7,damping:[0,.005,.05,.9]});for(const c of r.curves)for(const key of ['sa','psa','sd','sv','psv'])assert.ok(c[key].every(v=>v===0));});
check('Undamped resonant sine agrees with closed-form transient maxima',()=>{
 const dt=.0005,n=4001,w=tau*2,r=E.responseSpectrum(wave(n,dt,2),dt,{min:2,max:3,count:2,damping:[0],tail:0});
 // u = t*cos(w*t)/(2*w) - sin(w*t)/(2*w*w), v = -t*sin(w*t)/2.
 let sd=0,sv=0;for(let i=0;i<n;i++){const t=i*dt;sd=Math.max(sd,Math.abs(t*Math.cos(w*t)/(2*w)-Math.sin(w*t)/(2*w*w)));sv=Math.max(sv,Math.abs(-t*Math.sin(w*t)/2));}
 near(r.curves[0].sd[0],sd,2e-4);near(r.curves[0].sv[0],sv,2e-4);near(r.curves[0].sa[0],w*w*sd,2e-4);near(r.curves[0].sa[0],r.curves[0].psa[0]);return `Sa relative error ${(Math.abs(r.curves[0].sa[0]/(w*w*sd)-1)).toExponential(3)}`;
});
check('Damped step agrees with analytic overshoot including near-critical damping',()=>{
 for(const z of [.02,.05,.2,.9,.999]){const r=E.responseSpectrum(new Float64Array(20001).fill(1),.0005,{min:2,max:3,count:2,damping:[z],tail:0});const exact=1+Math.exp(-Math.PI*z/Math.sqrt(1-z*z));near(r.curves[0].psa[0],exact,2e-4);}
});
check('Independent RK4 reference for broadband piecewise-linear input and free tail',()=>{
 const dt=.01,n=501,input=Float64Array.from({length:n},(_,i)=>(Math.sin(tau*1.3*i*dt)+.3*Math.cos(tau*7.1*i*dt))*Math.sin(Math.PI*i/(n-1))**2);
 const r=E.responseSpectrum(input,dt,{min:1.5,max:9,count:3,spacing:'linear',damping:[.005,.05,.2],steps:200,tail:2});let maxError=0;
 // Independent integrator: very fine RK4 with linear forcing and identical endpoint convention.
 for(const curve of r.curves)for(let j=0;j<r.frequency.length;j++){
  const w=tau*r.frequency[j],c=2*curve.damping*w,k=w*w,sub=100,h=dt/sub;let u=0,v=0,sa=0,sd=0,sv=0;
  const end=(n-1)*dt,steps=Math.ceil((end+2/r.frequency[j])/h);
  function force(t){if(t>=end)return 0;const ix=Math.min(n-2,Math.floor(t/dt)),q=t/dt-ix;return input[ix]*(1-q)+input[ix+1]*q;}
  const acc=(t,u,v)=>-force(t)-c*v-k*u;
  for(let i=0;i<steps;i++){const t=i*h,k1u=v,k1v=acc(t,u,v),k2u=v+h*k1v/2,k2v=acc(t+h/2,u+h*k1u/2,v+h*k1v/2),k3u=v+h*k2v/2,k3v=acc(t+h/2,u+h*k2u/2,v+h*k2v/2),k4u=v+h*k3v,k4v=acc(t+h,u+h*k3u,v+h*k3v);u+=h*(k1u+2*k2u+2*k3u+k4u)/6;v+=h*(k1v+2*k2v+2*k3v+k4v)/6;sa=Math.max(sa,Math.abs(-c*v-k*u));sd=Math.max(sd,Math.abs(u));sv=Math.max(sv,Math.abs(v));}
  for(const [key,ref] of [['sa',sa],['sd',sd],['sv',sv]]){const err=Math.abs(curve[key][j]/ref-1);maxError=Math.max(maxError,err);near(curve[key][j],ref,.004);}
 }return `maximum relative error ${maxError.toExponential(3)}`;
});
check('Linearity, sign invariance and all supported input unit conversions',()=>{
 const x=wave(2001,.001,3),o={min:1,max:20,count:5,tail:1},ref=E.responseSpectrum(x,.001,o);
 for(const [unit,scale] of [['m/s2',1],['g',9.80665],['gal',.01],['mm/s2',.001]]){const r=E.responseSpectrum(Float64Array.from(x,v=>-3*v/scale),.001,{...o,unit});near(r.pfa,3*ref.pfa);for(let c=0;c<2;c++)for(const key of ['sa','psa','sd','sv','psv'])for(let j=0;j<5;j++)near(r.curves[c][key][j],3*ref.curves[c][key][j],1e-10);}
});
check('Absolute and pseudo acceleration differ at finite damping',()=>{const r=E.responseSpectrum(wave(20001,.001,2),.001,{min:2,max:3,count:2,damping:[.2],tail:0});assert.ok(Math.abs(r.curves[0].sa[0]/r.curves[0].psa[0]-1)>.02);});
check('Free vibration captures delayed peak; nonzero final sample resets equilibrium',()=>{
 const dt=.001,x=Float64Array.from({length:51},(_,i)=>Math.sin(Math.PI*i/50)),o={min:1,max:2,count:2,damping:[0],steps:100};const short=E.responseSpectrum(x,dt,{...o,tail:0}),tail=E.responseSpectrum(x,dt,{...o,tail:2});assert.ok(tail.curves[0].sa[0]>2*short.curves[0].sa[0]);
 // Rectangular pulse duration .05s: exact undamped residual amplitude = 2*sin(w*duration/2).
 const step=E.responseSpectrum(new Float64Array(51).fill(1),dt,{...o,tail:2});near(step.curves[0].sa[0],2*Math.sin(tau*.05/2),2e-4);
});
check('Internal time-step refinement converges on the same sampled input',()=>{
 const x=wave(501,.01,7),o={min:7,max:20,count:2,damping:[.02],tail:1};const a=E.responseSpectrum(x,.01,{...o,steps:20}),b=E.responseSpectrum(x,.01,{...o,steps:100}),c=E.responseSpectrum(x,.01,{...o,steps:500});const ea=Math.abs(a.curves[0].sa[0]-c.curves[0].sa[0]),eb=Math.abs(b.curves[0].sa[0]-c.curves[0].sa[0]);assert.ok(eb<ea/5);return `20-step error ${ea.toExponential(3)}; 100-step error ${eb.toExponential(3)}`;
});
check('FRS uses preprocessed acceleration without FFT window; integrated input rejected',()=>{
 const dt=.001,x=wave(1000,dt,3),t=Float64Array.from(x,(_,i)=>i*dt),frs={min:1,max:20,count:5,unit:'m/s2',tail:1};
 const a=E.analyze({t,x,options:{stft:false,frs,window:'hann'}}),b=E.analyze({t,x,options:{stft:false,frs,window:'rectangular'}});assert.deepEqual(a.frs,b.frs);
 const direct=E.responseSpectrum(a.processed,a.processedSummary.dt,frs);assert.deepEqual(a.frs,direct);
 assert.throws(()=>E.analyze({t,x,options:{stft:false,frs,integrate:1}}),/FRS/);
 const uneven=Float64Array.from(t,(v,i)=>v+.0002*Math.sin(i));assert.throws(()=>E.analyze({t:uneven,x,options:{stft:false,frs,tolerance:1}}),/Non-uniform/);
 assert.ok(E.analyze({t:uneven,x,options:{stft:false,frs,resample:true}}).frs);
});
check('Invalid ranges, damping, units, nonfinite data and work limits fail explicitly',()=>{
 const x=new Float64Array(50);for(const options of [{unit:''},{unit:'toString'},{min:0},{max:51},{min:5,max:4},{count:1},{count:1001},{spacing:'bad'},{damping:[]},{damping:[1]},{damping:[NaN]},{damping:[-.1]},{damping:Array(9).fill(.05)},{steps:19},{steps:Infinity},{tail:-1},{tail:NaN}])assert.throws(()=>E.responseSpectrum(x,.01,options),/FRS/);
 assert.throws(()=>E.responseSpectrum([0,NaN],.01),/FRS/);assert.throws(()=>E.responseSpectrum(x,0),/FRS/);
 assert.throws(()=>E.responseSpectrum(new Float64Array(100000),.01,{count:1000}),/1億/);
});
const passed=results.filter(r=>r.pass).length;
fs.writeFileSync(path.join(__dirname,'FRS_VALIDATION.md'),`# FRS numerical validation\n\n${passed}/${results.length} checks passed. Tested the engine embedded in index.html using ${process.version}.\n\n| Check | Result | Detail |\n|---|---|---|\n${results.map(r=>`| ${r.name} | ${r.pass?'PASS':'FAIL'} | ${r.detail.replaceAll('|','/')} |`).join('\n')}\n\nReproduce: \`node tests/frs.test.cjs\`. References are analytic transient solutions and an independent fine-step RK4 integrator, not another call to Newmark. Numerical verification is not validation of input records, seismic design compliance, or full browser compatibility.\n`);
console.log(`${passed}/${results.length} FRS checks passed`);process.exitCode=passed===results.length?0:1;
