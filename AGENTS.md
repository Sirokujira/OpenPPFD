# AGENTS.md — OpenPPFD

`CLAUDE.md` と同じ規約を、単独で読める形にまとめたもの
(`CLAUDE.md` / `.claude/` を読まないエージェント向け)。
**規約を変えたら両方直すこと。**

## このリポジトリは何か

植物工場照明 (PPFD / スペクトル / 輻射伝達) ソルバー (C 言語、
バイナリ `oppfd`)。[OpenFDTD](https://github.com/Sirokujira/OpenFDTD) の
姉妹プロジェクトで、ビルド規約・移植性規則を共有する。外部ライブラリ
(HDF5/LAPACK/BLAS) に依存しない。

処理内容: 直方体チャンバの直接光 (解析配光 / IES・EULUMDAT 実測配光) →
遮蔽物 (棚板等) の影 → 壁・遮蔽物のラジオシティ相互反射 →
群落 (葉群) の Beer-Lambert 減衰 → 測定面ごとの
PPFD / YPFD / ePAR / DLI / 照度 / R:FR / 均斉度。

### 前提として崩してはいけない 2 点

1. **波動光学ではない。** FDTD (OpenFDTD) や RCWA (OpenRCWA) が扱う
   回折・干渉ではなく、光線・輻射伝達の問題として解く。
2. **一般照明の測光量では評価できない。** 一般照明は lm / lx / lm/W
   (標準比視感度 V(λ)、緑 555nm ピーク) だが、光合成は 400-700nm の
   光子を数え (PPFD [µmol/m²/s])、量子収率は McCree 作用曲線
   (赤 600-680nm ピーク) に従う。赤主体の園芸用 LED は lm/W では低く
   µmol/J では高い。**両方を必ず並べて出力する。**

## ディレクトリ構成

| パス | 役割 |
|---|---|
| `include/ppfd.h` | 定数・構造体 (`ppfd_t`)・全プロトタイプ |
| `src/input_data.c` | `.ppfd` テキストのパース (2 パス: 分光グリッド確定 → 本体) |
| `src/spectrum.c` | 分光グリッド、McCree / V(λ)、PPFD・YPFD・lx への換算 |
| `src/photometry.c` | 実測配光ファイル (IES LM-63 / EULUMDAT) の読み込みと補間 |
| `src/canopy.c` | 群落の二流層カーネル、衝突源の預け入れ、拡散場の解 |
| `src/geometry.c` | 壁・遮蔽物のパッチ分割、測定面、群落スラブ、遮蔽判定 |
| `src/formfactor.c` | 微小面→多角形の閉形式形態係数、パッチ間形態係数、可視率 |
| `src/direct.c` | 光源→パッチ / 任意点の直接放射照度 |
| `src/radiosity.c` | 波長ビンごとのラジオシティ (Jacobi 反復) |
| `src/solve.c` | ドライバ、測定面の合成、エネルギー収支 |
| `src/output.c` | `ppfd.log` と CSV |
| `data/sample/` | サンプル `.ppfd` と検証スクリプト `ppfd_check.sh` |

## ビルド / テスト

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j"$(nproc)"

# 検証 (すべて解析解・SI 定義値との比較)
sh data/sample/ppfd_check.sh bin/oppfd /tmp/ppfd-check
```

依存は無し (OpenMP は任意)。Linux / macOS / Windows (MSVC + Ninja) の
3 OS を CI で検証する。

## 移植性の絶対規則 (OpenFDTD の Windows CI で実際に踏んだもの)

- **C99 VLA 禁止** (MSVC C2057/C2466)。`malloc` + 明示インデックスの
  フラット配列を使う。
- **OpenMP のループ変数は pragma の前に宣言する** (MSVC は OpenMP 2.0)。
  ループ変数は符号付き `int`。
- **float\*/double\* の取り違え禁止**: 配列の実型と読み出しポインタ型の
  不一致は Windows で 0xC0000005 クラッシュ (glibc は偶然耐える)。
  `plen` のみ float、それ以外の物理量は double。
- libm リンクは CMake の `MATH_LIB` 変数経由 (Windows には m.lib が無い)。
- MSVC フラグは CMakeLists の既存ブロックに従う
  (`/utf-8`, `_USE_MATH_DEFINES`, `_CRT_SECURE_NO_WARNINGS`)。
- 数学定数は `PI` / `C0` / `H_PLANCK` 等の `ppfd.h` の自前マクロを使う
  (`M_PI` に依存しない)。
- `strcasecmp` は MSVC に無い。`streq_ci()` (utils.c) を使う。
- `-Wall -Wextra -Wshadow -Wconversion -Wpedantic` で警告ゼロを維持する。
- **「厳密に 0 になるはず」の量を丸めに委ねない**。退化した配置 (受光点の
  接平面に乗った多角形など) は積和の打ち消しに頼ると FMA の有無で結果が
  変わる。Apple Silicon の macOS だけ落ちた実例があるので、退化は
  打ち消しではなく明示的な判定で落とす (`ff_point_poly` の冒頭)。
  ローカルで再現するには `-march=native -ffp-contract=fast` を付けてビルドする。

## 設計の規則

- グローバル変数は使わない。状態は `ppfd_t` コンテキスト構造体 1 個を
  main で確保して関数に渡す。
- **分光量はすべて「ビン積分値」で持つ** (`[W/m²]` であって `[W/m²/nm]`
  ではない)。光量子換算が重み付き総和になり、帯域の端の扱いが
  `band_weight()` 1 箇所に集約される。CSV へ出すときだけ `/dlam` する。
- 入力キー追加は `src/input_data.c` に、既定値は「キー省略時に従来動作と
  完全一致」になるよう初期化する (後方互換)。未知キーは無視 (前方互換)。
- 新機能には `data/sample/` の**解析解付き**検証ケースを追加し、
  `ppfd_check.sh` に判定を足す。解析解が取れない機能は、既存の厳密な
  性質 (エネルギー収支・単調性・対称性) で縛れないか先に考える。
- **エネルギー収支を壊さない**。群落が無い閉キャビティでは反射率にも
  形状にもよらず「壁の吸収 = 光源の放射束」が厳密に成り立つ。
  `closure error` が 1e-5 より悪化したら形態係数まわりの回帰を疑う
  (遮蔽物があるケースは別。影の境界でパッチ求積の被積分関数が不連続に
  なるので 1e-4 台が下限で、これは形態係数の問題ではない)。
- **形態係数は「受光側の法線」しか見ない**。放射側のパッチが裏を向いて
  いないかは呼び出し側で判定する (`ff_point_poly` の外)。凸キャビティ
  だけなら常に成り立つが、遮蔽物は同じ多角形を表裏 2 枚のパッチで
  共有するので、これを落とすと形態係数が二重計上になり行和が 1 を超え、
  ラジオシティが発散する。
- **精度を落とす最適化は必ずログに出す**。`quadrature` の自動選択のように
  暗黙にサンプル数を減らす箇所は、実際に使った値を `plog` する。
- OpenMP は任意依存。`#ifdef _OPENMP` でガードする。
- 外部ライブラリを追加しない。

## 物理・データの規則

- **組み込みの McCree 作用曲線は公表曲線の近似値** (誤差の目安 ±0.03、
  625nm で 1.0 に正規化)。`spectrum.c` のコメントにその旨を明記してある。
  ここを「公式データ」と書き換えない。精度が要る用途は入力キー
  `actionspectrum = file <path>` で上書きする。
- CIE 1924 V(λ) は標準表 (5nm, 380-780nm)。555nm でちょうど 1.0 という
  性質が `1 W @555nm → 683 lm` の検証に効いているので崩さない。
- 光量子換算係数 `PHOTON_K` は h, c, N_A の SI 定義値だけからなる厳密量
  (λ[nm] × 8.359346e-3 µmol/J)。丸めた定数に置き換えない。
- **群落の散乱は「素の遮断で運び、取り除いた分を配る」で組む**。現行の
  減衰係数の `√(1-ω)` は素の遮断率ではなく二流方程式の減衰固有値なので、
  これを据え置いたまま散乱を戻すと二重計上になる。`leafscatter = on` では
  ビームを `k = G a` (波長非依存) で運び、`√(1-ω)` は解の固有値として
  出力側に現れる。
- **実測配光ファイル (IES/LDT) からは配光の形だけを取る**。記載の cd/lm
  は測光量なので、SPD 抜きに放射束 [W] や PPF へは換算できない。
  放射束は入力ファイル側 (W 欄 / `ppf` / `lumens`) が決める。配光は
  `∫Î dΩ = 1` へ**厳密に**正規化する (台形則で近似するとその誤差が
  そのまま放射束の誤差になり `closure error` を汚す)。

## 検証ケース一覧 (`ppfd_check.sh`)

| ケース | 何を縛っているか |
|---|---|
| `unit.ppfd` | 光量子換算と光束の定義 (4.639437 µmol/s, 683 lm)、逆二乗+余弦則、DLI の算術 |
| `cavity.ppfd` | 閉キャビティのエネルギー収支 (ρ=0.9 と ρ=0)、平均壁面照度 |
| `canopy.ppfd` | Beer-Lambert 減衰 exp(−G·LAI) を 2 通りの G で |
| `panel.ppfd` | 面光源分割の収束 (矩形形態係数の閉形式) |
| `unit625.ppfd` | 作用曲線の正規化点で YPFD = PPFD (曲線の値に依存しない) |
| `photometry.ppfd` | 実測配光ファイル : 放射束の保存、等方/Lambert 配光の解析解との一致、C 面の向き、IES と EULUMDAT の一致、`rot` の作用、カタログ値 (`ppf`/`lumens`) からの換算 |
| `shelf.ppfd` | 遮蔽物込みのエネルギー収支と平均照度、仕切られた室が厳密に真っ暗になること、断面いっぱいの棚板なら可視率が厳密 (標本化に落ちる対が 0) |
| `shadow.ppfd` | 本影が厳密に 0、影の外が逆二乗則どおり、影があっても収支が閉じること |
| `rack2.ppfd` | 2 段ラック : 棚板で仕切った段どうしが厳密に独立 (下段の LED を消しても上段の分布が 1 ビットも変わらない) |
| 層カーネル (`--selftest`) | 幾何によらない代数的恒等式 (R+T+A=1、ω=0 で Beer-Lambert、ω=1 で無吸収、半無限アルベド) |
| `canopy.ppfd` 散乱 | ω=0 でのビット一致、群落込みの収支、ω=1 で吸収が厳密に 0 |
| `rack.ppfd` | 実用例の単調性 (群落を通ると PPFD と R:FR が下がる) |

## CI

`.github/workflows/ci.yml`: Linux / macOS / Windows (MSVC + Ninja) の 3 ジョブ。
検証スクリプトは 3 OS とも同一の POSIX sh を実行する (Windows は Git Bash)。
タグ `v*` push で Release にバイナリ添付。
