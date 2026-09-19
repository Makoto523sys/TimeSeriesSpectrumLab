Base.@kwdef struct FRSOptions
    minimum::Float64=0.1
    maximum::Float64=50.0
    count::Int=200
    spacing::String="log"
    damping::Vector{Float64}=[.02,.05] # ratios, not percent
    steps_per_period::Int=100
    residual_cycles::Float64=5.0
    unit::String="m/s2"
end
struct FRSRow
    frequency::Float64
    period::Float64
    damping::Float64
    sa::Float64
    psa::Float64
    sd::Float64
    sv::Float64
    psv::Float64
    substeps::Int
end
function estimate_frs(n,dt,o::FRSOptions)
    require(Sys.WORD_SIZE==64,"Use 64-bit Julia for large-data analysis")
    require(n>=2 && finitepositive(dt),"FRS needs >=2 samples and positive dt")
    require(finitepositive(o.minimum) && isfinite(o.maximum) && o.minimum<o.maximum<=1/(2dt)*(1+1e-12),"FRS frequency must satisfy 0 < min < max <= Nyquist")
    require(o.count>=2 && o.spacing in ("log","linear"),"Invalid FRS frequency grid")
    require(!isempty(o.damping) && all(z->isfinite(z)&&0<=z<1,o.damping),"Damping ratios must satisfy 0 <= zeta < 1")
    require(o.steps_per_period>=20 && isfinite(o.residual_cycles) && o.residual_cycles>=0,"Require >=20 substeps/period and nonnegative residual cycles")
    require(o.unit in ("m/s2","g","gal","mm/s2"),"Explicit acceleration unit required")
    frequency=o.spacing=="log" ? exp.(range(log(o.minimum),log(o.maximum),length=o.count)) : collect(range(o.minimum,o.maximum,length=o.count))
    frequency[1]=o.minimum; frequency[end]=o.maximum
    subdivisions=Int[]; tails=Int[]; work=big(0)
    for f in frequency
        q=dt*f*o.steps_per_period
        require(isfinite(q) && q<typemax(Int),"Substep count exceeds machine integer range")
        sub=max(1,ceil(Int,q)); tail=o.residual_cycles/f/(dt/sub)
        require(isfinite(tail) && tail<typemax(Int),"Residual step count exceeds machine integer range")
        nt=ceil(Int,tail); push!(subdivisions,sub); push!(tails,nt)
        work+=big(n-1)*sub+nt
    end
    (;frequency,subdivisions,tails,work=work*length(o.damping))
end
@inline function newmark_step(u,v,a,force,h,c,k,denominator)
    pu=u+h*v+h*h*a/4; pv=v+h*a/2
    na=(-force-c*pv-k*pu)/denominator
    (pu+h*h*na/4,pv+h*na/2,na)
end
# O(1) working storage per oscillator: never allocate histories or interpolated input.
function oscillator(input,dt,f,z,sub,tail)
    w=TAU*f; k=w*w; c=2z*w; h=dt/sub; denominator=1+c*h/2+k*h*h/4
    require(all(isfinite,(w,k,c,h,denominator)) && h>0,"FRS coefficients exceed floating-point range")
    u=0.0; v=0.0; a=-input[1]; ma=0.0; mu=0.0; mv=0.0
    @inbounds for i in 2:length(input)
        previous=input[i-1]; delta=input[i]-previous
        for q in 1:sub
            u,v,a=newmark_step(u,v,a,previous+delta*(q/sub),h,c,k,denominator)
            ma=max(ma,abs(-c*v-k*u)); mu=max(mu,abs(u)); mv=max(mv,abs(v))
            q % 8192 == 0 && GC.safepoint()
        end
        i % 8192 == 0 && GC.safepoint()
    end
    a=-c*v-k*u # abrupt zero floor acceleration, continuous u and v
    for i in 1:tail
        u,v,a=newmark_step(u,v,a,0.0,h,c,k,denominator)
        ma=max(ma,abs(-c*v-k*u)); mu=max(mu,abs(u)); mv=max(mv,abs(v))
        i % 8192 == 0 && GC.safepoint()
    end
    require(all(isfinite,(ma,mu,mv,k*mu,w*mu)),"FRS response overflow; rescale input")
    FRSRow(f,1/f,z,ma,k*mu,mu,mv,w*mu,sub)
end
function atomic_toml(path,data)
    tmp=path*".tmp"
    open(tmp,"w") do io; TOML.print(io,data;sorted=true); end
    mv(tmp,path;force=true)
end
rowdict(r::FRSRow)=Dict{String,Any}(string(k)=>getfield(r,k) for k in fieldnames(FRSRow))
"FRS without a workload ceiling. Checkpoints commit complete oscillators; rerun to resume."
function response_spectrum(acceleration,dt; options=FRSOptions(),parallel=true,
                           checkpoint_directory=nothing,progress=(done,total)->nothing)
    o=options; grid=estimate_frs(length(acceleration),dt,o)
    scale=Dict("m/s2"=>1.0,"g"=>9.80665,"gal"=>.01,"mm/s2"=>.001)[o.unit]
    input=Float64.(acceleration).*scale
    require(all(isfinite,input),"FRS input contains nonfinite acceleration")
    pfa=maximum(abs,input)
    total=Base.checked_mul(o.count,length(o.damping)); rows=Vector{FRSRow}(undef,total)
    warnings=String[]
    o.maximum*dt>.05 && push!(warnings,"FRS input has <20 samples/period near the upper limit; interpolation cannot recover lost bandwidth.")
    (length(input)-1)*dt*o.minimum<1 && push!(warnings,"Longest oscillator period exceeds record duration; check baseline and record length.")
    o.residual_cycles==0 && push!(warnings,"Record-end free vibration is excluded.")
    o.residual_cycles>0 && abs(input[end])>.01pfa && push!(warnings,"Nonzero final acceleration is switched instantly to zero for the residual response.")
    fingerprint=bytes2hex(sha256(string("FRS-v1-Newmark;",dt,';',repr(grid.frequency),';',repr(o),';',bytes2hex(sha256(reinterpret(UInt8,input))))))
    lockdir=nothing
    if checkpoint_directory !== nothing
        mkpath(checkpoint_directory); lockdir=joinpath(checkpoint_directory,"run.lock")
        require(!ispath(lockdir),"Checkpoint is locked: $lockdir. If no process is using it, remove only this lock directory before restarting.")
        mkdir(lockdir)
    end
    try
        if checkpoint_directory !== nothing
            file=joinpath(checkpoint_directory,"metadata.toml")
            if isfile(file)
                require(get(TOML.parsefile(file),"fingerprint","")==fingerprint,"Checkpoint belongs to different input or FRS settings; choose a new directory")
            else
                require(all(name->name=="run.lock",readdir(checkpoint_directory)),"Nonempty checkpoint directory without metadata")
                atomic_toml(file,Dict("fingerprint"=>fingerprint,"jobs"=>total,"method"=>"Newmark beta=1/4 gamma=1/2"))
            end
        end
        pending=Int[]
        for job in 1:total
            j=(job-1)%o.count+1; d=(job-1)÷o.count+1
            file=checkpoint_directory===nothing ? "" : joinpath(checkpoint_directory,"$job.toml")
            if isfile(file)
                data=TOML.parsefile(file)
                require(get(data,"fingerprint","")==fingerprint,"Checkpoint fingerprint mismatch in job $job")
                r=FRSRow((data[string(k)] for k in fieldnames(FRSRow))...)
                require(r.frequency==grid.frequency[j] && r.damping==o.damping[d] && r.substeps==grid.subdivisions[j] && r.period==1/r.frequency && all(v->isfinite(v)&&v>=0,(r.sa,r.psa,r.sd,r.sv,r.psv)),"Corrupt checkpoint job $job")
                rows[job]=r
            else
                push!(pending,job)
            end
        end
        resumed=total-length(pending); done=Ref(resumed); progress_lock=ReentrantLock()
        progress(resumed,total)
        function runjob(job)
            j=(job-1)%o.count+1; d=(job-1)÷o.count+1
            r=oscillator(input,Float64(dt),grid.frequency[j],o.damping[d],grid.subdivisions[j],grid.tails[j])
            rows[job]=r
            if checkpoint_directory !== nothing
                data=rowdict(r); data["fingerprint"]=fingerprint
                atomic_toml(joinpath(checkpoint_directory,"$job.toml"),data)
            end
            lock(progress_lock) do
                done[]+=1; progress(done[],total)
            end
        end
        if parallel && Threads.nthreads()>1
            Threads.@threads :dynamic for job in pending; runjob(job); end
        else
            for job in pending; runjob(job); end
        end
        (;frequency=grid.frequency,rows,pfa,dt,options=o,work=grid.work,warnings,resumed,fingerprint)
    finally
        lockdir!==nothing && rm(lockdir) # only our empty lock; never deletes checkpoints
    end
end
