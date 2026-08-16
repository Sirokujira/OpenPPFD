# OpenPPFD

A radiative-transfer solver for **plant factory lighting** in C: photosynthetic
photon flux density (PPFD), spectral design, and canopy light distribution.
Companion project to [OpenFDTD](https://github.com/Sirokujira/OpenFDTD) — same
build conventions, same portability rules, no external libraries.

植物工場 (人工光型植物工場) の照明設計ソルバー (C 言語)。
[OpenFDTD](https://github.com/Sirokujira/OpenFDTD) の姉妹プロジェクトで、
ビルド規約・移植性規則を共有します。外部ライブラリ (HDF5/LAPACK/BLAS)
には依存しません。

## なぜ FDTD/RCWA と別なのか

植物工場の照明は**波動光学の問題ではありません**。扱うのは

- **光量子束密度 (PPFD)** — µmol/m²/s。エネルギー (W) ではなく光子の数
- **スペクトル設計** — 赤 660nm / 青 450nm / 遠赤 730nm の配分、R:FR 比
- **輻射伝達** — 反射壁での相互反射と、群落 (葉群) の Beer-Lambert 減衰

であって、波長スケールの構造による回折・干渉ではありません。したがって
FDTD (`OpenFDTD`) や RCWA (`OpenRCWA`) ではなく、光線・輻射伝達の枠組みで
解きます。

さらに、**一般照明の測光量がそのまま使えません**。一般照明は
lm / lx / lm/W (標準比視感度 V(λ) 基準、緑 555nm がピーク) ですが、
光合成は 400-700nm の光子を数え (PPFD)、量子収率は McCree 作用曲線
(赤 600-680nm がピーク) に従います。赤主体の園芸用 LED は
**lm/W では低く、µmol/J では高い**ので、一般照明の指標で評価すると
設計を誤ります。OpenPPFD は両方を必ず並べて出力します。

```
photon efficacy   =      2.69344 umol/J     <- 園芸用の指標
luminous efficacy =      63.8918 lm/W       <- 一般照明の指標 (同じ器具)
```

## Theory / 理論概要

1. **直接光** — 光源は配光 `I(θ) = Φ(m+1)/(2π) cos^m θ` (m=1 で Lambert、
   `iso` で等方) を持つ点光源。面光源は nu × nv 個の点光源へ分割します。
   実器具のカタログ配光 (IESNA LM-63 / EULUMDAT) も読めます
   (`ies` / `ldt` オプション)。ファイルから使うのは**配光の形だけ**で、
   放射束は入力ファイルの W 欄 (または `ppf` / `lumens`) から取ります
   — cd や lm は V(λ) 基準の測光量なので、SPD 抜きに W や µmol/s へは
   換算できないためです。配光は `∫Î dΩ = 1` に厳密に正規化するので、
   実測配光を使ってもエネルギー収支の検証はそのまま成り立ちます。
2. **相互反射** — チャンバ 6 面を矩形パッチに分割し、波長ビンごとに
   ラジオシティ方程式 `B_i = ρ_i (E_i^dir + Σ_j F_ij τ_ij B_j)` を
   Jacobi 反復で解きます。形態係数は「微小面 → 平面多角形」の閉形式を
   受光パッチ上で Gauss 求積したもの (放射側は厳密)。
3. **遮蔽物** — 棚板・トレイ・バッフルを軸並行直方体 (厚さ 0 の板も可) で
   置けます。表面はパッチに分割されてラジオシティに参加し、直接光には
   2 値の影を落とします。多段ラックのように棚板が断面いっぱいの配置では
   パッチ対の可視率が 0 か 1 にしかならず、形態係数は遮蔽があっても
   厳密なままです (段どうしの分離が厳密に成り立つ)。
4. **群落** — 葉面積密度 a = LAI/厚さ、葉群投影係数 G の一様スラブとして、
   経路長 s に対し `τ = exp(-G a s sqrt(1-ω_λ))`。ω_λ = 葉の反射率+透過率
   (Goudriaan の散乱補正) で、ω=0 のとき厳密な Beer-Lambert 則になります。
   赤・青は強く吸収され、緑と遠赤は群落深部まで届きます。
   `leafscatter = on` にすると葉が散乱した光をキャビティへ戻します
   (下記)。
5. **測定面** — 水平・上向きの量子センサ相当。各セルで直接光 + 壁からの
   間接光を合成し、PPFD / YPFD / ePAR / lx / DLI / R:FR を出します。

### 現在の制限

- 幾何は**直方体チャンバ 1 個 + 軸並行の遮蔽物**のみ。斜めの反射板や
  曲面は置けません。
- 反射は **Lambert 拡散 + 軸並行壁の鏡面反射** (`specular`)。一般の BRDF や
  斜めの反射板は未対応で、鏡面壁は群落・遮蔽物と併用できません (下記)。
- 群落の散乱光は `leafscatter = on` で戻せます (既定は off = 従来動作)。
  拡散場は水平一様な 1 次元カラムとして解くので、群落の高さ範囲に埋没した
  側壁パッチとの横方向のやりとりは表現できません (エネルギーは保存します)。
  群落スラブは 1 個だけなので、多段ラックの各段に別々の群落は置けません。
- パッチ間の群落透過率はパッチ重心どうしを結ぶ線分で評価する近似です。
- 部分的にしか遮られないパッチ対の可視率は標本化 (最大 4×4 × 4×4) で
  求める近似です。標本化に落ちた対の数は `ppfd.log` に出ます。
- 組み込みの McCree 作用曲線は公表曲線の**近似値**です
  (`actionspectrum = file` で上書きしてください。下記参照)。

## Build / ビルド

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j"$(nproc)"

# 検証 (すべて解析解・定義値との比較)
sh data/sample/ppfd_check.sh bin/oppfd /tmp/ppfd-check
```

- Linux / macOS / Windows (MSVC + Ninja) の 3 OS を CI で検証しています。
- OpenMP は任意です (見つかれば形態係数・直接光・測定面を並列化。
  無くてもビルド・実行可能)。

## Usage / 実行

```bash
oppfd [-n <threads>] input.ppfd
```

出力 (カレントディレクトリ):

| ファイル | 内容 |
|---|---|
| `ppfd.log` | 実行ログ。器具の総量・エネルギー収支・測定面ごとの表。`=== normal end ===` で正常終了 |
| `ppfd_summary.csv` | 測定面ごとの 1 行要約 (平均/最小/最大 PPFD、均斉度、DLI、YPFD、lx、R:FR) |
| `ppfd_map_<target>.csv` | 測定面のセルごとの分布 (GUI での等高線表示用) |
| `ppfd_spectrum.csv` | 測定面の面平均スペクトル + 作用曲線 + V(λ) |

## Input format / 入力形式

テキスト形式。1 行目は `OpenPPFD 1 0`、最終行は `end`。`#` はコメント行。
書式は `キー = 値...`。単位は長さ [m]、波長 [nm]、放射束 [W]。
**未知のキーは無視されます** (前方互換)。**省略時の既定値は、キーを
書かない場合と挙動が一致する**よう選んであります (後方互換)。

### 解析条件

| キー | 書式 | 既定 |
|---|---|---|
| `title` | `title = 任意文字列` | 空 |
| `wavelength` | `wavelength = λ開始 λ終了 Δλ` | `380 780 5` |
| `chamber` | `chamber = Lx Ly Lz` | `1 1 1` |
| `patchdiv` | `patchdiv = nu nv` (1 面あたりの分割数) | `6 6` |
| `photoperiod` | `photoperiod = 時間/日` (DLI 用) | `16` |
| `solver` | `solver = 最大反復数 収束判定` | `200 1e-6` |
| `quadrature` | `quadrature = msub` (パッチ直接光の複合分割数、0=自動) | `0` |
| `specbounce` | `specbounce = N` (鏡面反射の鏡像段数) | `2` |
| `actionspectrum` | `actionspectrum = file <path>` | 組み込み McCree |

### 材料 (Lambert 拡散面 / 葉)

```
material = <名前> gray  <ρ> [τ]
material = <名前> table <λ:ρ[:τ]> <λ:ρ[:τ]> ...
material = <名前> file  <path>            # "λ ρ [τ]" のテキスト
```

τ (透過率) は群落の葉材料でのみ使い、散乱係数 ω = ρ + τ になります。
テーブルは線形補間、範囲外は端の値でクランプします。

### スペクトル

```
spectrum = <名前> mono      <λ>                        # 単色
spectrum = <名前> gauss     <中心λ> <半値全幅>
spectrum = <名前> sum       <中心λ> <半値全幅> <重み> ...  # 3 つ組の繰り返し (白色 LED 等)
spectrum = <名前> blackbody <T[K]>
spectrum = <名前> table     <λ:値> ...                  # 分光放射束密度 [/nm 相対]
spectrum = <名前> file      <path>                      # "λ 値" のテキスト (実測 SPD)
```

内部では波長ビンごとの重み (総和 1) として保持します。ピークの 1e-10 未満の
裾は切り落として再正規化するので、放射束は厳密に保存されます。

### 光源

```
led   = <x> <y> <z> <放射束W> <スペクトル名> [オプション]
array = <x0> <y0> <z> <nx> <ny> <ピッチx> <ピッチy> <1灯あたり放射束W> <スペクトル名> [オプション]
```

| オプション | 意味 | 既定 |
|---|---|---|
| `beam <m>` | 配光 `I ∝ cos^m θ` | `1` (Lambert) |
| `iso` | 等方 `I = Φ/4π` | — |
| `ies <path>` | IESNA LM-63 の実測配光ファイル | — |
| `ldt <path>` | EULUMDAT の実測配光ファイル | — |
| `rot <deg>` | 配光を配光軸まわりに回す (`ies`/`ldt` 用) | `0` |
| `dir <dx> <dy> <dz>` | 配光軸 | `0 0 -1` (真下) |
| `size <w> <h> [nu nv]` | 面光源の寸法と分割数 | 点光源 |
| `ppf <µmol/s>` | 放射束をカタログの PPF で与える (W 欄より優先) | — |
| `lumens <lm>` | 放射束をカタログの光束で与える (W 欄より優先) | — |
| `input <W>` | 消費電力 (µmol/J・lm/W の算出に使う) | 未指定 |

#### 実測配光ファイル (IES / EULUMDAT)

```
led = 0.6 0.3 0.4 0 spd_red ies bar660.ies ppf 250 input 90
```

- 使うのは**配光の形だけ**です。ファイル中の cd / lm は放射束の決定には
  使いません (SPD 抜きに lm → W → µmol/s は決まらないため)。放射束は
  W 欄か `ppf` / `lumens` オプションで与えてください。ファイル記載の
  定格光束は参考値として `ppfd.log` に出ます。
- 角度の定義は IES 光度分布型 C / EULUMDAT と同じで、γ = 0 が配光軸
  (`dir`、既定は真下)。C 面は既定の `dir` のとき **C = 0 が +x、
  C = 90 が +y** で、`rot <deg>` で器具を配光軸まわりに回せます。
- IES は水平角の対称性 (0-90 / 0-180 / 0-360)、EULUMDAT は
  `Isym` = 0/1/2/3/4 を展開します。表の γ 範囲外の光度は 0 です。
- 配光は `∫Î dΩ = 1` へ正規化します。この積分は双線形補間の**厳密な**
  積分 (γ 方向は `∫(a+bγ)sinγ dγ` の閉形式) で評価するので、台形則の
  O(Δγ²) 誤差が放射束に漏れません。実測配光でも `closure error` は
  1e-6 のままです。
- 光度分布型 B / A (IES) は未対応で、型 C として扱い警告を出します。

### 壁・遮蔽物・群落・測定面・帯域

```
wall     = <面> <材料名>        # 面 = xmin/xmax/ymin/ymax/zmin(floor)/zmax(ceiling)/side/all
occluder = <名前> <x0> <x1> <y0> <y1> <z0> <z1> <材料名> [div <nu> <nv>]
canopy   = <ztop> <zbot> <LAI> <G> <葉材料名>
leafscatter = on | off          # 群落の散乱光をキャビティへ戻す (既定 off)
leaflayers  = <N>               # 二流の層数 (既定 8)
target   = <名前> <z> <nx> <ny> [<x0> <x1> <y0> <y1>]
band     = <名前> <λ1> <λ2>     # 光量子束の割合を出す帯域 (既定 : UV-A/青/緑/赤/遠赤)
```

`wall` を書かない面は完全吸収 (ρ=0) です。`target` の範囲を省略すると
チャンバの床面全体になります。

#### 鏡面反射する壁 (`specular`)

```
material = mirrorfilm gray 0.05 specular 0.90
```

材料の末尾に `specular <ρs>` を付けると、その面は `ρs` を**鏡面反射**します
(拡散反射 `ρd` とは別枠。吸収は `1 − ρd − ρs`)。チャンバは軸並行の直方体
なので、鏡面壁は**光源の鏡像**と等価です — 直接光についてはこれで厳密で、
近似は鏡像の段数の打ち切りだけです。段数は `specbounce = <N>` (既定 2)、
打ち切りで捨てた放射束は `ppfd.log` に出ます。鏡面が 1 面だけなら鏡像は
自分自身の面で折り返さないので、`specbounce = 1` で厳密になります
(ログの `truncated` が 0)。

**拡散壁との混在**も正しく扱います。拡散壁からの放射が鏡面壁で正反射して
戻る経路は**拡張形態係数** (放射側パッチの鏡像) で運びます。軸並行の
直方体では鏡像パッチへの直線が必ず鏡面の矩形内で平面を横切るので、
この拡張は厳密です。検証は「完全鏡面の天井 = 鏡映した幾何 (高さ 2 倍 +
光源の鏡像)」の等価性で縛っており、両者はビット単位で一致します。
`ρd + ρs > 1` の非物理な材料は入力段でエラーにします。

`ρs` は**波長非依存のスカラ**です。反射フィルムや研磨アルミはほぼ中性
なので、分光にせず 1 個の値で持ちます。

**制限** : 鏡像からの直線は「折り返した経路」を表すので、群落
(`canopy`) や遮蔽物 (`occluder`) があると経路長も遮蔽判定も合いません。
黙って誤差を出さないよう、この組み合わせは入力段でエラーにします。

#### 遮蔽物 (棚板・トレイ・バッフル)

`occluder` は軸並行の直方体です。**1 辺の長さを 0 にすると厚さ 0 の板**に
なり、表裏 2 面だけがパッチになります (棚板はこちらが確実)。表面は
`wall` と同じくラジオシティに参加するので、エネルギー収支は閉じたまま
です (`absorbed by walls` は遮蔽物の吸収を含みます)。`div` は 1 面あたりの
分割数で、省略時は `patchdiv` に従います。

置き方の注意 :

- **棚板の面が壁パッチの境界に乗るよう `patchdiv` を選ぶこと。** 例えば
  `chamber = 1.2 0.6 0.9` で棚板を z = 0.45 に置くなら `patchdiv = 8 8`
  (0.9 × 4/8 = 0.45)。パッチが棚板をまたぐと、そのパッチ 1 枚を通して
  上下段の間に光が漏れます。
- **遮蔽物の面を壁と同一平面で接触させないこと** (床の上に厚さのある箱を
  直に置く等)。同一平面のパッチどうしは互いに見えないものとして扱うので、
  その隙間の光が行き場を失います。床から浮かせるか厚さ 0 の板にします。
- 部分的にしか遮られないパッチ対の可視率は標本化で求めます。棚板が断面
  いっぱいの配置ではどの対も「全部見える / 全部見えない」に決まるので
  形態係数は厳密なままです。ログの
  `partially shadowed patch pairs` が 0 ならそうなっています。
- 影の境界ではパッチ求積の被積分関数が不連続になるため、境界にかかる
  パッチは自動的に分割を上げます (実行数はログに出ます)。それでも
  `closure error` は 1e-4 台が下限で、これは形態係数ではなく直接光の
  求積の限界です。

### 作用曲線の差し替え

組み込みの McCree (1972) 相対量子収率は**公表曲線をディジタイズした
近似値** (誤差の目安 ±0.03、625nm で 1.0 に正規化) です。厳密な値が
必要な場合は

```
actionspectrum = file mccree.txt      # 各行 "λ[nm] 相対量子収率"
```

で上書きしてください。YPFD の検証は「625nm 単色なら YPFD = PPFD」という
曲線の値に依存しない構造的性質で行っているので、差し替えても検証は通ります。

## Validation / 検証 (`data/sample/`)

`ppfd_check.sh` が実行する判定。すべて厳密解 (SI 定義または閉形式) が
期待値で、誤差の出どころは形態係数の求積と面光源の分割だけです。

| ケース | 解析解 | 実測誤差 | 許容 |
|---|---|---|---|
| `unit.ppfd` PPF | 1 W @555nm → 555×8.359346e-3 = 4.639437 µmol/s | +0.0001% | 0.01% |
| `unit.ppfd` 光束 | 1 W @555nm → 683 lm (lm の定義) | 0.0000% | 0.01% |
| `unit.ppfd` 逆二乗+余弦 | E(r) = (Φ/π)h²/(h²+r²)² を 2 点で | +0.0001% | 0.2% |
| `unit.ppfd` DLI | PPFD × 時間 × 3600 × 1e-6 | −0.0003% | 0.001% |
| `cavity.ppfd` 収支 (ρ=0.9) | 壁の吸収 = 光源の放射束 (形状によらず厳密) | −0.0009% | 0.5% |
| `cavity.ppfd` 収支 (ρ=0) | 同上 (直接光ソルバー単体) | 0.0000% | 0.5% |
| `cavity.ppfd` 平均壁面照度 | Φ/(A_tot(1−ρ)) = 3.086420 W/m² | −0.001% | 0.5% |
| `canopy.ppfd` 透過率 | exp(−G·LAI) = exp(−1), exp(−2) | −0.0001% | 1% |
| `panel.ppfd` 面光源軸上 | E = ΦF/A_s (矩形形態係数の閉形式) | +0.047% | 1% |
| `unit625.ppfd` YPFD | 正規化点 625nm で YPFD = PPFD | 0.0000% | 0.5% |
| `photometry.ppfd` PPF | 配光ファイルによらず PPF = 4.639437 µmol/s | +0.0001% | 0.01% |
| `photometry.ppfd` 等方 IES | `iso.ies` → E = Φ/(4πh²) (= `iso` と一致) | +0.0001% | 0.2% |
| `photometry.ppfd` Lambert IES | `lambert.ies` → E = Φ/(πh²) (= `beam 1` と一致) | +0.0026% | 0.2% |
| `photometry.ppfd` 収支 | 実測配光でも壁の吸収 = 放射束 (正規化が厳密) | 5e-7 | 1e-6 |
| `photometry.ppfd` C 面 | `asym.ies` の +x / +y の比 = (1+0.6s)/(1−0.6s) | −0.029% | 0.2% |
| `photometry.ppfd` IES=LDT | 同一配光を IES と EULUMDAT で書いて一致 | 0 | 1e-6 |
| `photometry.ppfd` rot 90 | `rot 90` で床面分布が転置される | 0 | 1e-6 |
| `unit.ppfd` ppf/lumens | `ppf 10` → PPF = 10、`lumens 683` → 1 W (555nm) | 0.0000% | 1e-6 |
| `shelf.ppfd` 収支 | 遮蔽物込みでも壁の吸収 = 放射束 | 0.0000% | 0.01% |
| `shelf.ppfd` 平均照度 | Φ/(A(1−ρ)) = 0.625 W/m² (A は棚板の表裏込み) | 0.0000% | 0.01% |
| `shelf.ppfd` 遮断 | 仕切られた下室の PPFD = 0 (厳密) | 0 | 0 |
| `shelf.ppfd` 可視率 | 断面いっぱいの棚板なら標本化に落ちる対が 0 | 0 | 0 |
| `shadow.ppfd` 本影 | 本影内の PPFD = 0 (厳密) | 0 | 0 |
| `shadow.ppfd` 影の外 | E = (Φ/π)h²/(h²+r²)² | +0.0001% | 0.2% |
| `shadow.ppfd` 収支 | 影の分は板が吸収するので収支は閉じる | −2.5e-4 | 1e-3 |
| `rack2.ppfd` 段の独立 | 下段の LED を消しても上段の分布が変わらない | 0 | 1e-6 |
| `rack2.ppfd` 段の遮断 | 上段の LED だけなら下段の PPFD = 0 | 0 | 0 |
| `mirror.ppfd` 鏡像 | E = Φ/(4π)(1/d₁² + ρs/d₂²) = 1.772134 | −0.0003% | 0.2% |
| `mirror.ppfd` 打ち切り | 鏡面 1 面なら捨てる放射束 = 0 | 0 | 0 |
| `mirror.ppfd` 収支 | 鏡面成分は鏡像が運ぶので壁の吸収 = 放射束 | −1.4e-11 | 1e-6 |
| `mirror.ppfd` ρs=0 | 鏡像が消えて直接光だけの解析解に戻る | +0.0001% | 0.2% |
| `mirror2.ppfd` 収支 | 鏡面天井 + 拡散床でも収支が閉じる (修正前 −7.7%) | 5e-7 | 1e-5 |
| `mirror2.ppfd` 等価 | 完全鏡面の折り返し = 鏡映した幾何 (参照解とビット一致) | 0 | 1e-6 |
| `mirror2.ppfd` 材料 | ρd + ρs > 1 は入力段で拒否 | — | — |
| 層カーネル (`--selftest`) | `R+T+A=1` / ω=0 で `T=e^{-t}`・`R=0` / ω=1 で `A=0` / t→∞ で `(1-u)/(1+u)` | 4e-16 | 1e-13 |
| `canopy.ppfd` 散乱 ω=0 | `leafscatter = on` でも出力がビット一致 | 0 | 0 |
| `canopy.ppfd` 散乱 収支 | 群落の吸収を実計算しても収支が閉じる | 4.6e-9 | 1e-5 |
| `canopy.ppfd` 散乱 ω=1 | 無吸収葉なら群落の吸収が厳密に 0 | 0 | 0 |
| `rack.ppfd` | 群落を通ると PPFD と R:FR が下がる (単調性) | — | — |

閉キャビティのエネルギー収支は特に効く検証です。反射率にも形状にも
よらず厳密に成り立つので、形態係数の求積誤差 (現状 1e-5) がそのまま
`closure error` に出ます。

## Example / 実用例

`data/sample/rack.ppfd` は栽培面 1.2 × 0.6 m のラック 1 段に、
赤/青/白/遠赤の LED バー 4 本 × 12 灯を配置し、LAI = 3 のレタス相当群落を
置いたものです。

```
--- source totals ---
radiant flux      =       39.504 W
PPF (400-700nm)   =      192.958 umol/s
luminous flux     =      4577.21 lm
input power       =        71.64 W
photon efficacy   =      2.69344 umol/J
luminous efficacy =      63.8918 lm/W

--- target "canopytop" (z = 0.25 m) ---
PPFD [umol/m2/s]             260.4
DLI [mol/m2/day]            15.001
uniformity min/avg          0.73445
R:FR (660/730)               15.311

--- target "canopybase" (z = 0.06 m) ---
PPFD [umol/m2/s]             32.65
R:FR (660/730)               5.124
```

群落を通ると PPFD が 1/8 に落ちる一方、R:FR は 15.3 → 5.1 まで下がります
(赤は葉に吸収され遠赤は透過するため)。これは徒長・光形態形成に直結する
量で、一般照明の lx では見えません。

### 多段ラック

`data/sample/rack2.ppfd` は同じ栽培面を 2 段にし、z = 0.45 m の棚板
(トレイ、ρ = 0.4) で仕切ったものです。棚板が無いと下段の LED の光が
上段へ回り込んでしまうので、多段ラックの設計には遮蔽判定が要ります。

```
--- target "tier2" (z = 0.5 m) ---      上段 (棚板の 5cm 上)
PPFD [umol/m2/s]             85.74
uniformity min/avg          0.81968

--- target "tier1" (z = 0.25 m) ---     下段 (群落上面)
PPFD [umol/m2/s]              63.46
uniformity min/avg           0.5899
```

同じ器具・同じ間隔でも、下段は群落の上面が反射の弱い葉であるぶん
均斉度が落ちます。上段の分布は下段の LED を消しても 1 ビットも
変わりません (棚板が断面いっぱいなので段どうしは厳密に独立)。

## Sibling projects / 姉妹リポジトリ

| リポジトリ | 手法 | バイナリ |
|---|---|---|
| [OpenFDTD](https://github.com/Sirokujira/OpenFDTD) | 3 次元 FDTD | `ofd` |
| [OpenRCWA](https://github.com/Sirokujira/OpenRCWA) | 周期構造 RCWA | `orcwa` |
| [OpenBPM](https://github.com/Sirokujira/OpenBPM) | 導波路 BPM | `obpm` |
| [OpenPEEC](https://github.com/Sirokujira/OpenPEEC) | 準静的 PEEC | `peec` |

## Reference

- K. J. McCree, "The action spectrum, absorptance and quantum yield of
  photosynthesis in crop plants", *Agricultural Meteorology* **9**, 191 (1971/72)
- J. Goudriaan, *Crop Micrometeorology: A Simulation Study* (1977) — 散乱葉の
  消散係数 sqrt(1−ω) 補正
- CIE 1924 標準比視感度 V(λ)

## License

MIT
