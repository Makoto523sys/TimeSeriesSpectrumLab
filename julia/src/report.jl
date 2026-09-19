escape_html(s)=replace(string(s),'&'=>"&amp;",'<'=>"&lt;",'>'=>"&gt;",'"'=>"&quot;",'\''=>"&#39;")
# Display-only min/max envelopes preserve extrema; numerical CSV/binary files stay complete.
function envelope_indices(y; buckets=700)
    n=length(y); n<=2buckets && return collect(1:n)
    indices=Int[]
    for b in 0:buckets-1
        lo=b*n÷buckets+1; hi=(b+1)*n÷buckets
        imin=lo; imax=lo
        @inbounds for i in lo:hi
            y[i]<y[imin] && (imin=i); y[i]>y[imax] && (imax=i)
        end
        append!(indices,sort!(unique([lo,imin,imax,hi])))
    end
    indices
end
const COLORS=["#12619c","#bc5623","#24803c","#8b4aa0","#9b7300","#b23e62","#377d83","#575d6a"]
function svg_lines(series,xlabel,ylabel;logx=false)
    io=IOBuffer(); width=1000; height=360; left=95; top=35; pw=865; ph=260
    prepared=[(x=s.x[envelope_indices(s.y)],y=s.y[envelope_indices(s.y)],name=s.name) for s in series]
    transform(v)=logx ? log10(v) : v
    xmin=minimum(transform(v) for s in prepared for v in s.x); xmax=maximum(transform(v) for s in prepared for v in s.x)
    ymin=minimum(v for s in prepared for v in s.y); ymax=maximum(v for s in prepared for v in s.y)
    xmin==xmax && (xmax=xmin+1)
    dy=ymax-ymin; dy==0 && (dy=max(1,abs(ymax)*.1)); ymin-=.06dy; ymax+=.06dy
    println(io,"<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 $width $height' role='img' aria-label='",escape_html(ylabel),"'><rect width='100%' height='100%' fill='white'/><g font-family='sans-serif' font-size='12'>")
    for i in 0:5
        x=left+pw*i/5; y=top+ph-ph*i/5; xv=xmin+(xmax-xmin)*i/5; yv=ymin+(ymax-ymin)*i/5
        println(io,"<path d='M $x $top V $(top+ph) M $left $y H $(left+pw)' stroke='#e3e8ed' fill='none'/>")
        println(io,"<text x='$x' y='$(top+ph+20)' text-anchor='middle'>",@sprintf("%.4g",logx ? 10.0^xv : xv),"</text><text x='$(left-8)' y='$(y+4)' text-anchor='end'>",@sprintf("%.4g",yv),"</text>")
    end
    for (i,s) in enumerate(prepared)
        print(io,"<polyline fill='none' stroke='",COLORS[mod1(i,length(COLORS))],"' stroke-width='1.3' points='")
        for j in eachindex(s.x)
            px=left+pw*(transform(s.x[j])-xmin)/(xmax-xmin); py=top+ph-ph*(s.y[j]-ymin)/(ymax-ymin)
            @printf(io,"%.2f,%.2f ",px,py)
        end
        println(io,"'><title>",escape_html(s.name),"</title></polyline>")
    end
    println(io,"<text x='$(left+pw/2)' y='345' text-anchor='middle'>",escape_html(xlabel),logx ? " (log)" : "","</text><text transform='translate(18 165) rotate(-90)' text-anchor='middle'>",escape_html(ylabel),"</text></g></svg>")
    String(take!(io))
end
function write_report(path,t,x,y,fft,psd,oldfft,oldpsd,peaks,stft,frs,meta)
    unit=meta["config"]["input"]["unit"]; integration=meta["config"]["preprocess"]["integrate"]
    processed_unit=integration==0 ? unit : "($unit)·s^$integration"
    atomic_file(path) do io
        println(io,"<!doctype html><html lang='ja'><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'><title>Time Series Spectrum Lab — Julia</title><style>body{font:15px system-ui;color:#183149;background:#eef2f6;max-width:1150px;margin:25px auto;padding:0 18px}section{background:white;border:1px solid #d4dfe8;padding:18px;margin:18px 0;border-radius:6px}svg{width:100%;height:auto}table{border-collapse:collapse}td,th{padding:7px 12px;border-bottom:1px solid #ddd;text-align:right}p{line-height:1.7}.legend{display:flex;gap:18px;flex-wrap:wrap}summary{cursor:pointer;font-weight:bold}</style><h1>Time Series Spectrum Lab — Julia</h1><p>数値結果は全点を保存。グラフは表示用の最小・最大包絡、STFTプレビューは区画内最大値です。CSVの再現条件と単位は run.toml に記録しています。</p>")
        println(io,"<section><h2>入力・処理後の統計量</h2><table><tr><th>量</th><th>入力</th><th>処理後</th></tr>")
        for key in ("n","dt","duration","fs","nyquist","mean","min","max","std","rms")
            println(io,"<tr><td>$key</td><td>",escape_html(meta["input_summary"][key]),"</td><td>",escape_html(meta["processed_summary"][key]),"</td></tr>")
        end
        println(io,"</table></section>")
        function chart(title,name,series,xlabel,ylabel;logx=false)
            svg=svg_lines(series,xlabel,ylabel;logx)
            atomic_file(joinpath(dirname(path),name*".svg")) do f; write(f,svg); end
            push!(meta["output_files"],name*".svg")
            println(io,"<section><h2>",escape_html(title),"</h2><a href='$name.svg'>SVG保存</a>",svg,"<div class='legend'>")
            for (i,s) in enumerate(series)
                println(io,"<span style='color:",COLORS[mod1(i,length(COLORS))],"'>",escape_html(s.name),"</span>")
            end
            println(io,"</div></section>")
        end
        chart("入力時刻歴","history-input",[(x=t,y=x,name="入力（必要時再サンプリング）")],"Time [s]","Input [$unit]")
        chart("処理後時刻歴","history-processed",[(x=t,y=y,name="平均除去・detrend・フィルタ・積分後")],"Time [s]","Processed [$processed_unit]")
        for (s,old,label,name) in [(fft,oldfft,"Amplitude","fft"),(psd,oldpsd,"PSD","welch")]
            s===nothing && continue
            units=label=="PSD" ? "($processed_unit)^2/Hz" : processed_unit
            chart(label,name,[(x=s.frequency,y=s.values,name="処理後")],"Frequency [Hz]","$label [$units]")
            if old!==nothing
                oldunits=label=="PSD" ? "($unit)^2/Hz" : unit
                chart(label*"（フィルタ・積分前）",name*"-before",[(x=old.frequency,y=old.values,name="平均除去・detrend後、フィルタ・積分前")],"Frequency [Hz]","$label [$oldunits]")
            end
        end
        if !isempty(peaks)
            println(io,"<section><h2>卓越周波数</h2><table><tr><th>Hz</th><th>周期 [s]</th><th>Amplitude</th><th>Welch PSD（補間）</th></tr>")
            for p in peaks; println(io,"<tr>",join(["<td>"*@sprintf("%.7g",v)*"</td>" for v in values(p)]),"</tr>"); end
            println(io,"</table></section>")
        end
        if stft!==nothing
            a=stft.preview; ny,nx=size(a); hi=maximum(a); lo=hi-80
            println(io,"<section><h2>STFT preview</h2><p>時間は左から右、周波数は下から上。色: PSD [dB re 1 ($processed_unit)²/Hz]、範囲 $lo ～ $(hi)。区画内最大値表示。</p><svg viewBox='0 0 1000 340'><rect width='100%' height='100%' fill='white'/>")
            for j in 1:nx,i in 1:ny
                q=clamp((a[i,j]-lo)/80,0,1); r=round(Int,20+235q); g=round(Int,40+160q); b=round(Int,120-90q)
                println(io,"<rect x='$(80+(j-1)*860/nx)' y='$(20+(ny-i)*260/ny)' width='$(860/nx+.1)' height='$(260/ny+.1)' fill='rgb($r,$g,$b)'/>")
            end
            println(io,"<g font-family='sans-serif' font-size='13'><text x='80' y='305'>",@sprintf("%.5g",stft.start)," s</text><text x='850' y='305'>",@sprintf("%.5g",stft.stop)," s</text><text x='5' y='28'>",@sprintf("%.4g",stft.frequency[end])," Hz</text><text x='25' y='280'>0 Hz</text></g></svg></section>")
        end
        if frs!==nothing
            println(io,"<section><h2>FRS</h2><p>PFA = ",frs.pfa," m/s² · ",frs.work," integration steps · 完了済み再利用 ",frs.resumed,"件。SaとPSaは異なる定義です。FRSの出力単位はSIです。</p></section>")
            for (key,label,units) in [(:sa,"Sa absolute","m/s²"),(:psa,"PSa pseudo","m/s²"),(:sd,"Sd relative","m"),(:sv,"Sv relative","m/s"),(:psv,"PSv pseudo","m/s")]
                println(io,key==:sa ? "<details open>" : "<details>","<summary>",label,"</summary>")
                series=[(x=frs.frequency,y=[getfield(r,key) for r in frs.rows[(d-1)*length(frs.frequency)+1:d*length(frs.frequency)]],name="ζ=$(100z)%") for (d,z) in enumerate(frs.options.damping)]
                key in (:sa,:psa) && push!(series,(x=frs.frequency,y=fill(frs.pfa,length(frs.frequency)),name="PFA"))
                chart(label,"frs-"*string(key),series,"Natural frequency [Hz]","$label [$units]";logx=true)
                periods=[(x=reverse(1.0./s.x),y=reverse(s.y),name=s.name) for s in series]
                chart(label*" / period","frs-"*string(key)*"-period",periods,"Natural period [s]","$label [$units]";logx=true)
                println(io,"</details>")
            end
        end
        println(io,"<section><h2>条件・出力</h2><ul>")
        for warning in meta["warnings"]; println(io,"<li>",escape_html(warning),"</li>"); end
        println(io,"</ul><p><a href='run.toml'>run.toml（全設定・単位・入力SHA-256）</a></p><ul>")
        for file in meta["output_files"]; println(io,"<li><a href='",escape_html(file),"'>",escape_html(file),"</a></li>"); end
        println(io,"</ul></section></html>")
    end
end
