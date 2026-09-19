# Julia版 — 大規模時刻歴・スペクトル解析

ブラウザ版と同じ数値処理を、64-bit Juliaで実行するCLI／ライブラリです。**FRSの1億ステップ上限、524,288点の入力上限、Welch/STFTのブラウザ用上限は設けていません。** 計算量に応じてCPU時間・RAM・ディスク容量は必要です。無制限のデータが任意のPCで動くという意味ではありません。

入力UIが未確定でも使えるよう、数値計算とCSV読み込み、設定、出力を分離しています。後からGUIやWeb画面を追加しても `src/` の計算関数を再利用できます。既存の `index.html` は変更しません。

## 始め方（Windows・WSL・Linux）

64-bit Julia 1.10以上を利用します。以下はリポジトリ直下で実行します。初回だけFFTWなどの依存を取得するため通信が必要です。入力データを外部へ送信する処理はありません。

```sh
julia --project=julia -e 'using Pkg; Pkg.instantiate()'
julia --project=julia --threads=auto julia/run.jl julia/examples/full-analysis.toml
```

Windowsのコマンドプロンプトでは初回コマンドを次に置き換えてください。

```bat
julia --project=julia -e "using Pkg; Pkg.instantiate()"
```

初回はコンパイル時間がかかります。結果は `julia/results/example/` に出力されます。`report.html` を開くと、時刻歴、FFT、Welch PSD、卓越周波数、STFT、各減衰比のFRSを確認できます。グラフはSVGでも保存されます。HTMLレポートは折り畳み式の確認用表示で、元のWebアプリの編集・ズーム操作やPNGダウンロードUIを移植したものではありません。

CPUを使いすぎたくない場合は `--threads=4` などを指定します。FRSでは各固有振動数・減衰比の組合せを独立に並列計算します。時刻方向の逐次積分は保持します。FFTW内部のスレッド数は変更していません。

### 大きい記録でFRSだけ必要な場合

`julia/examples/frs-only.toml` をコピーして入力パス・列・加速度単位を編集します。

```sh
julia --project=julia --threads=auto julia/run.jl julia/examples/frs-only.toml
```

この設定ではFFT／Welch／STFTと全時刻歴CSV出力を無効にし、FRSに必要な計算だけを行います。解析の種類を個別に有効・無効化できるため、FRSのためだけに大きなFFTを計算する必要はありません。

## 入力

設定はTOMLです。**TOML内のパスはその設定ファイルからの相対パス**です。絶対パスも使用できます。Windowsのパスは `C:/analysis/acc.csv` のように `/` を使うか、TOMLの単一引用符で `'C:\analysis\acc.csv'` と記述してください。

```toml
[input]
path = "my-floor-acceleration.csv"
delimiter = "auto"     # auto / comma / tab / space
header = "auto"        # auto / yes / no
skip_lines = 0         # 先頭のメタ情報行を読み飛ばす場合のみ明示
# 列番号はブラウザの内部番号と異なり、1始まり
time_column = 1
value_column = 2
time_scale = 1.0       # 入力時刻を秒に換算: msなら0.001
quantity = "Acceleration"
unit = "Gal"          # 一般解析の表示ラベル。数値換算はしない
```

ヘッダー・区切りの自動判定、引用符付きCSV、空白・タブ区切り、複数列に対応します。ファイルを1行ずつ読み、選択した時間列と値列だけを数値配列として保持します。複数物理行にまたがる引用符付きフィールドは未対応です。欠損・非数値・無限大、列数の不一致、重複・逆順時刻はエラーにし、勝手に削除・並べ替えしません。

時刻列がなく、サンプリング間隔が分かっている場合は次のように指定できます。

```toml
[input]
path = "acc.csv"
time_column = 0
value_column = 1
dt = 0.01             # 常に秒。100 Hzなら0.01 s
time_scale = 1.0       # time_column=0では使わない
header = "yes"
skip_lines = 0         # 実ファイルのメタ情報行数に合わせる
```

`time_column=0` は最初の値の時刻を0として生成します。気象庁形式などに専用の自動判定は設けていません。NS/EW/UDのどの列か、メタ情報の行数、単位とサンプリング周波数を実ファイルで確認してください。FRSを床応答として評価するには、床の**絶対加速度時刻歴**が必要です。地盤入力なら地盤入力の応答スペクトルになります。

## 設定とブラウザ版との対応

| 処理 | Julia設定／出力 | 数値上の扱い |
|---|---|---|
| サンプリング検査・統計 | `run.toml` | 母標準偏差、RMS、dt・Nyquist・jitter等 |
| 非等間隔処理 | `preprocess.resample` | 明示選択時のみ、端点・点数を保つ線形補間 |
| 平均除去・detrend | `remove_mean`, `detrend` | 記録全体に適用 |
| FIR | `filter`, `low_cut`, `high_cut` | 同じ最大257タップ、反射境界、中心合わせ |
| 数値積分 | `integrate = 0/1/2` | 初期値ゼロの台形則 |
| FFT | `[fft]`, `fft.csv` | 任意長、片側片振幅、周期窓、coherent gain補正 |
| Zero padding | `padding = "none"/"next"/"double"` | 窓を掛けた元の点数で正規化 |
| Welch PSD | `[welch]`, `welch.csv` | 同じ片側・窓エネルギー正規化、端数区間を除外 |
| 卓越周波数 | `[peaks]`, `peaks.csv` | 局所ピーク・DC除外・最小間隔、PSDは補間参考値 |
| STFT | `[stft]` | 1フレームずつFFT、全行列をメモリに持たない |
| 処理前後比較 | `fft-before.csv`, `welch-before.csv`, HTML | 平均除去・detrend後、フィルタ／積分前と比較 |
| FRS | `[frs]`, `frs.csv` | Sa・PSa・Sd・Sv・PSv、Newmark、残留自由振動 |
| グラフ | `report.html`, `*.svg` | 表示だけ最小・最大包絡、数値出力は全点 |

窓は `rectangular` / `hann` / `hamming` / `blackman`。Welch/STFTの `overlap` は百分率ではなく0～1未満の比率です。例えば50%は `0.5`。区間長は8以上、入力長を超える場合は入力長に合わせます。FFTとWelchで異なる窓を指定することもできますが、ブラウザ版との比較時は同じ窓にしてください。ピーク検出にはFFTとWelchの両方が必要です。

FRS使用時の例:

```toml
[frs]
enabled = true
unit = "gal"              # m/s2 / g / gal / mm/s2。SIへ数値換算
minimum = 0.1              # Hz
maximum = 50.0             # Hz、入力のNyquist以下
count = 200
spacing = "log"            # log / linear
damping_percent = [0.5, 2.0, 5.0]
steps_per_period = 100
residual_cycles = 5.0
parallel = true
checkpoint_directory = "checkpoints/case01"
```

ブラウザ版のFRSと同じ運動方程式・初期条件・終端条件を使用します。

\[
\ddot u+2\zeta\omega\dot u+\omega^2u=-a_f(t),\qquad
S_a=\max|-2\zeta\omega\dot u-\omega^2u|,\quad
PS_a=\omega^2\max|u|.
\]

入力を線形補間し、Newmark β=1/4、γ=1/2で内部細分化します。初期相対変位・速度はゼロ。記録終了後は床加速度をゼロへ切り替え、変位・速度を連続にして、指定周期数の自由振動を追います。FRSにFFT窓を掛けません。台形積分とFRSの同時使用はエラーです。FFTの前処理と同様、平均除去／detrend／FIRはFRSにも影響します。

**刻みを細かくしただけでは入力のaliasingは直りません。** 低減衰ピークは固有振動数格子への依存も確認してください。SaとPSa、SvとPSvは別の量です。ピーク拡幅・包絡化や規格適合の設計スペクトル化は実施しません。

## 大規模FRSのメモリ・途中再開

FRSは入力の数値配列を全ジョブで共有し、各振動子では現在の変位・速度・加速度と最大値だけを保持します。**全振動子の応答時刻歴や、細分化済み入力配列を生成しません。** 作業量の総計はBigIntで数え、1億ステップを超えても拒否しません。個々のループ数・配列サイズは64-bit整数と実メモリの範囲に制約されます。

- FRS自体のメモリは概ね入力点数に比例する配列と、振動数×減衰比の結果表です。
- CLIの前処理は時間・入力・処理後等のFloat64配列を持つため、ファイル全体をディスクだけで処理するout-of-coreソルバではありません。
- 全記録FFTは全データとFFT作業領域のRAMを使います。RAMに収まらない全記録FFTは本実装の対象外です。
- Welch/STFTはセグメント用バッファとFFT計画を再利用します。STFT全結果は逐次書き出します。

チェックポイントは1振動数・1減衰比の計算完了ごとに一時ファイルから置き換えて保存します。中断すると、その時点で完了した振動子までを再利用できます。計算中の1振動子の時刻途中からは再開しません。

同じ設定をもう一度実行すると、入力加速度内容・時間刻み・FRS設定のSHA-256識別子が一致するチェックポイントを読み込み、残りだけ計算します。入力や設定が違う場合はエラーにします。出力先を再使用する場合は `[output] overwrite = true` を明示してください。チェックポイントだけを再使用し、出力先を別名にすることもできます。

```toml
[output]
directory = "results/case01"
overwrite = true
history = false
report = true
```

`Ctrl+C`で中止できます。スレッド処理の中断が反映されるまで時間がかかることがあります。プロセス強制終了・電源断後に `analysis.lock` やチェックポイント内の `run.lock` が残る場合は、**同じ解析が動いていないことを確認して、その空のロックディレクトリだけを削除**して再実行してください。`.toml` の完了済みチェックポイントは残します。同じチェックポイント／出力先を複数プロセスで同時に使わないでください。

計算の完了状態は `run.toml` の `status = "complete"` で確認します。`running`／`incomplete` は全出力が完成したことを意味しません。個々のCSV／STFTファイルも `.partial` から完成時に置き換えます。`overwrite=true` は明示した出力先の生成物を上書きします。今回生成されたファイルは `run.toml` の `output_files` が正であり、以前の別設定で生成したファイルが残っていても今回の結果に含めないでください。

## STFTの保存形式

既定の `format = "binary"` は `stft.f64` にFloat64のPSDを保存します。周波数binが連続し、その次に次フレームが続きます。サイズは `8 × bin数 × frame数` bytes。byte order、shape、窓・hop・スケールは `stft-metadata.toml` に保存します。周波数と時刻はそれぞれ別CSVです。

同じbyte orderの環境ではJuliaで以下のように読み込めます（この読み方は全STFTをRAMに展開します。大規模時は必要部分だけread/seekしてください）。

```julia
using TOML
m = TOML.parsefile("results/stft-metadata.toml")
A = Matrix{Float64}(undef, m["bins"], m["frames"])
open("results/stft.f64", "r") do io
    read!(io, A)
end
```

`format = "csv"` では `time_s,frequency_Hz,PSD` のlong形式で全セルを書き出します。大きなSTFTではテキスト出力が計算より遅くなり、ディスク使用量も増えます。HTMLには最大160×120区画のプレビューだけを載せます。各区画の最大dB値であり、全数値結果は変更しません。

## Juliaから直接呼ぶ

```julia
using TimeSeriesSpectrumLab

t = collect(0.0:0.001:60.0)
a = sin.(2pi .* 3 .* t)                  # m/s²
opts = FRSOptions(minimum=0.1, maximum=50.0, count=400,
                  damping=[0.005, 0.02, 0.05], # APIは比率、TOMLは百分率
                  unit="m/s2", steps_per_period=100)
r = response_spectrum(a, 0.001; options=opts, parallel=true,
                      checkpoint_directory="checkpoints/example")
r.rows[1].sa
```

`read_signal`, `inspect_sampling`, `signal_stats`, `preprocess`, `filter_signal`, `integrate_signal`, `amplitude_spectrum`, `welch_psd`, `stft_each`, `find_peaks` も独立して呼べます。`stft_each` のコールバックへ渡すPSD配列はフレームごとに再利用するため、後で保持する場合は `copy` してください。`run_analysis(config)` は辞書による呼び出しにも対応します（辞書の相対パスはその時の作業ディレクトリ基準）。

## 検証

```sh
julia --project=julia --threads=4 julia/test/runtests.jl
# 以下の比較試験だけ、開発者用にNode.jsも必要
julia --project=julia --threads=4 julia/test/browser_parity.jl
# 非ゼロの100万1点入力、1億ステップ超を本当に積分
julia --project=julia --threads=4 julia/scripts/benchmark_frs.jl
julia --project=julia --threads=4 julia/scripts/benchmark_spectra.jl
```

閉形式解、単位換算、フィルタ・積分、PSDエネルギー整合、並列／逐次一致、再開と異なる入力の拒否、CSV／STFT／HTML出力を検証します。ブラウザ比較は実際の `index.html` 内の計算コードをNodeで実行します。検証結果・大規模実行の計測値は [VALIDATION.md](VALIDATION.md) を参照してください。

## 依存と参考

数値定義はリポジトリの [README](../README.md) と共通です。FFTには [FFTW.jl](https://juliamath.github.io/FFTW.jl/stable/) を使います。スレッド数の指定は [Julia公式ドキュメント](https://docs.julialang.org/en/v1/manual/multi-threading/)、設定書式は [TOML公式API](https://docs.julialang.org/en/v1/stdlib/TOML/) を参照してください。
