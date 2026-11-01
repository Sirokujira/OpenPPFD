#!/bin/sh
# ppfd_check.sh — OpenPPFD 検証 (CI 用)
#
# data/sample/ のケースを実行し、解析解と比較する。期待値の導出は
# 各 .ppfd ファイルのコメントを参照。すべて厳密解 (定義または閉形式) で、
# 数値誤差の出どころは形態係数の求積と面光源の分割のみ。
#
# 使い方 : ppfd_check.sh <oppfd 実行ファイル(絶対パス)> [作業ディレクトリ]

set -e

OPPFD="$1"
WORK="${2:-.}"
SRC="$(cd "$(dirname "$0")" && pwd)"

if [ -z "$OPPFD" ]; then
	echo "Usage: ppfd_check.sh <oppfd> [workdir]" >&2
	exit 2
fi

mkdir -p "$WORK"
cp "$SRC"/*.ppfd "$WORK/"
cp "$SRC"/*.ies "$SRC"/*.ldt "$WORK/"
status=0

# 555 nm の光子換算係数 [umol/J] = 555 * 1e-3 / (h c N_A)
KPH555=4.6394372

# run <input.ppfd>
run() {
	(cd "$WORK" && "$OPPFD" -n 2 "$1" > /dev/null)
	grep -q "normal end" "$WORK/ppfd.log"
}

# chk <label> <actual> <expected> <tol> : 相対誤差判定
chk() {
	awk -v a="$2" -v e="$3" -v tol="$4" -v lb="$1" 'BEGIN {
		d = (e == 0) ? (a - e) : ((a - e) / e); ad = (d < 0) ? -d : d;
		printf "%-28s actual=%.7g expected=%.7g -> %s %+.4f%%\n", lb, a, e, (ad <= tol) ? "OK" : "NG", d * 100;
		exit (ad <= tol) ? 0 : 1
	}' || status=1
}

# chkabs <label> <actual> <limit> : |actual| <= limit 判定
chkabs() {
	awk -v a="$2" -v l="$3" -v lb="$1" 'BEGIN {
		aa = (a < 0) ? -a : a;
		printf "%-28s actual=%.4g -> %s (|x| <= %g)\n", lb, a, (aa <= l) ? "OK" : "NG", l;
		exit (aa <= l) ? 0 : 1
	}' || status=1
}

# logval <ラベル> : ppfd.log の "<ラベル> = <値> <単位>" から値を取る
logval() {
	awk -F= -v k="$1" 'index($0, k) == 1 {gsub(/^[ \t]+/, "", $2); split($2, f, " "); print f[1]; exit}' "$WORK/ppfd.log"
}

# cell <map.csv> <ix> <iy> <列番号> : 指定セルの値
cell() {
	awk -F, -v ix="$2" -v iy="$3" -v c="$4" 'NR > 1 && $1 == ix && $2 == iy {print $c; exit}' "$WORK/$1"
}

# trunc : 鏡像の打ち切りで捨てた放射束 [W]。
# logval では使えない ("specular ff :" の行が先に当たり、その行に = が無い)。
trunc() {
	awk -F'truncated flux =' 'NF > 1 {split($2, f, " "); print f[1]; exit}' "$WORK/ppfd.log"
}

# sumval <列番号> : ppfd_summary.csv の 1 行目 (最初の target)
sumval() {
	awk -F, -v c="$1" 'NR == 2 {print $c; exit}' "$WORK/ppfd_summary.csv"
}

# sumrow <面の名前> <列番号> : ppfd_summary.csv の指定の面の値
sumrow() {
	awk -F, -v n="$1" -v c="$2" 'NR > 1 && $1 == n {print $c; exit}' "$WORK/ppfd_summary.csv"
}

# maxdiff <map1.csv> <map2.csv> [1] : PPFD 列の最大相対差。
# 第 3 引数 1 で 2 番目の表の (ix, iy) を入れ替えて比較する (90 度回転の判定用)。
maxdiff() {
	awk -F, -v tr="$3" '
		FNR == NR { if (FNR > 1) a[$1 "," $2] = $5; next }
		FNR > 1 {
			k = (tr == "1") ? ($2 "," $1) : ($1 "," $2)
			if (k in a) {
				d = a[k] - $5; if (d < 0) d = -d
				m = (a[k] > $5) ? a[k] : $5
				if (m > 0 && d / m > x) x = d / m
			}
		}
		END {printf "%.3e", x}
	' "$WORK/$1" "$WORK/$2"
}

echo "== (a) unit conversion / inverse square law =="
run unit.ppfd
chk "PPF (1W @555nm)"      "$(logval 'PPF (400-700nm)')" "$KPH555" 1e-4
chk "luminous flux (1W)"   "$(logval 'luminous flux')"   683.0     1e-4
chkabs "cavity closure (rho=0)" "$(logval 'closure error')" 1e-6

# 中心セル (ix=iy=10, x=y=0.5) : E = (Phi/pi) h^2/(h^2+r^2)^2, r=0, h=0.5
E0=$(awk -v k="$KPH555" 'BEGIN {h = 0.5; printf "%.9g", (1.0 / (4 * atan2(1, 1) * h * h)) * k}')
chk "PPFD at r=0"          "$(cell ppfd_map_floor.csv 10 10 5)" "$E0" 2e-3
chk "PPFD max = centre"    "$(sumval 7)"                        "$E0" 2e-3

# 斜め方向のセル (ix=15, iy=10) : cos^4 則
E1=$(awk -v k="$KPH555" 'BEGIN {
	pi = 4 * atan2(1, 1); h = 0.5; r = (15 + 0.5) / 21 - 0.5;
	printf "%.9g", (1.0 / pi) * h * h / ((h * h + r * r) ^ 2) * k
}')
chk "PPFD at r=0.238 (cos^4)" "$(cell ppfd_map_floor.csv 15 10 5)" "$E1" 2e-3

# DLI = PPFD_avg * photoperiod * 3600 * 1e-6
DLI=$(awk -v q="$(sumval 5)" 'BEGIN {printf "%.9g", q * 16 * 3600 * 1e-6}')
chk "DLI arithmetic"       "$(sumval 10)" "$DLI" 1e-5

echo
echo "== (b) closed cavity energy balance =="
run cavity.ppfd
chk "absorbed by walls (rho=0.9)" "$(logval 'absorbed by walls')" 1.0 5e-3
# 一様反射率なら壁面の平均入射放射照度 = Phi / (A_tot (1 - rho))
chk "mean wall irradiance"        "$(logval 'mean wall irrad.')" \
	"$(awk 'BEGIN {printf "%.9g", 1.0 / (3.24 * 0.1)}')" 5e-3
# 反射率 0 の変種 : 直接光ソルバー単体の収支
sed 's/gray 0.90*/gray 0.0/' "$SRC/cavity.ppfd" > "$WORK/cavity_black.ppfd"
run cavity_black.ppfd
chk "absorbed by walls (rho=0)"   "$(logval 'absorbed by walls')" 1.0 5e-3

echo
echo "== (c) Beer-Lambert canopy extinction =="
grep -v '^canopy' "$SRC/canopy.ppfd" > "$WORK/canopy_none.ppfd"
sed 's/^canopy = 0.8 0.4 2.0 0.5/canopy = 0.8 0.4 2.0 1.0/' "$SRC/canopy.ppfd" > "$WORK/canopy_g1.ppfd"

run canopy_none.ppfd
Q0=$(cell ppfd_map_below.csv 10 10 5)
# 群落なしの絶対値も解析解と比較 (h = 1.0 - 0.2 = 0.8)
E2=$(awk -v k="$KPH555" 'BEGIN {h = 0.8; printf "%.9g", (1.0 / (4 * atan2(1, 1) * h * h)) * k}')
chk "no canopy, PPFD at r=0" "$Q0" "$E2" 2e-3

run canopy.ppfd
Q1=$(cell ppfd_map_below.csv 10 10 5)
chk "transmittance G=0.5 LAI=2" "$(awk -v a="$Q1" -v b="$Q0" 'BEGIN {printf "%.9g", a / b}')" \
	"$(awk 'BEGIN {printf "%.9g", exp(-1.0)}')" 1e-2

run canopy_g1.ppfd
Q2=$(cell ppfd_map_below.csv 10 10 5)
chk "transmittance G=1.0 LAI=2" "$(awk -v a="$Q2" -v b="$Q0" 'BEGIN {printf "%.9g", a / b}')" \
	"$(awk 'BEGIN {printf "%.9g", exp(-2.0)}')" 1e-2

echo
echo "== (d) Lambertian rectangular panel, on-axis =="
run panel.ppfd
# E = Phi * F / A_s,  F = 4 * F_corner(X = a/h, Y = b/h)  (a,b = 半幅)
E3=$(awk -v k="$KPH555" 'BEGIN {
	pi = 4 * atan2(1, 1); a = 0.2; b = 0.15; h = 0.5; As = (2 * a) * (2 * b);
	X = a / h; Y = b / h; sx = sqrt(1 + X * X); sy = sqrt(1 + Y * Y);
	F = 4 * (X / sx * atan2(Y / sx, 1) + Y / sy * atan2(X / sy, 1)) / (2 * pi);
	printf "%.9g", (F / As) * k
}')
chk "panel PPFD on axis"   "$(cell ppfd_map_floor.csv 10 10 5)" "$E3" 1e-2

echo
echo "== (e) McCree action spectrum normalisation =="
# 625 nm は組み込み作用曲線の正規化点なので YPFD == PPFD (曲線の値によらない)
sed 's/mono 555/mono 625/' "$SRC/unit.ppfd" > "$WORK/unit625.ppfd"
run unit625.ppfd
chk "YPFD/PPFD at 625nm"   "$(awk -v y="$(cell ppfd_map_floor.csv 10 10 6)" -v q="$(cell ppfd_map_floor.csv 10 10 5)" \
	'BEGIN {printf "%.9g", y / q}')" 1.0 5e-3

echo
echo "== (f) plant factory rack (smoke) =="
run rack.ppfd
RTOP=$(sumval 5)
RBASE=$(awk -F, 'NR == 3 {print $5}' "$WORK/ppfd_summary.csv")
FRTOP=$(sumval 15)
FRBASE=$(awk -F, 'NR == 3 {print $15}' "$WORK/ppfd_summary.csv")
awk -v t="$RTOP" -v b="$RBASE" -v ft="$FRTOP" -v fb="$FRBASE" 'BEGIN {
	printf "%-28s canopy top=%.4g  base=%.4g umol/m2/s\n", "rack PPFD", t, b;
	printf "%-28s canopy top=%.4g  base=%.4g\n", "rack R:FR", ft, fb;
	ok = (b < t) && (fb < ft) && (t > 0);
	printf "%-28s -> %s (PPFD and R:FR must fall through the canopy)\n", "rack monotonicity", ok ? "OK" : "NG";
	exit ok ? 0 : 1
}' || status=1

echo
echo "== (g) photometric distribution files (IES LM-63 / EULUMDAT) =="
# 配光ファイルは形だけを与え、放射束は入力ファイルの W 欄から取る。
# したがって配光によらず PPF は保存される。
run photometry.ppfd
chk "IES flux conservation"   "$(logval 'PPF (400-700nm)')" "$KPH555" 1e-4
# 等方配光 : E = Phi/(4 pi h^2), h = 0.5  ("iso" オプションと同じ)
EI=$(awk -v k="$KPH555" 'BEGIN {h = 0.5; printf "%.9g", k / (16 * atan2(1, 1) * h * h)}')
chk "IES isotropic at r=0"    "$(cell ppfd_map_floor.csv 10 10 5)" "$EI" 2e-3
cp "$WORK/ppfd_map_floor.csv" "$WORK/map_ies_iso.csv"

# Lambert 配光 (1 deg 刻みの表) : "beam 1" の解析解と一致する
sed 's/ies iso\.ies/ies lambert.ies/' "$SRC/photometry.ppfd" > "$WORK/phot_lambert.ppfd"
run phot_lambert.ppfd
chk "IES Lambertian at r=0"   "$(cell ppfd_map_floor.csv 10 10 5)" "$E0" 2e-3

# 同じ配光を EULUMDAT (Isym=1) で書いたもの : IES と一致するはず
sed 's/ies iso\.ies/ldt iso.ldt/' "$SRC/photometry.ppfd" > "$WORK/phot_ldt.ppfd"
run phot_ldt.ppfd
chkabs "EULUMDAT == IES (iso)"  "$(maxdiff map_ies_iso.csv ppfd_map_floor.csv)" 1e-6

# 非対称配光 I ~ cos(gamma) (1 + 0.6 sin(gamma) cos 2C) : IES (C 面 0-90 の
# 4 分の 1 対称) と EULUMDAT (Isym=4) の一致、C=0 が +x であること、
# rot 90 で床面分布が転置されること、放射束が保存されること。
sed 's/ies iso\.ies/ies asym.ies/' "$SRC/photometry.ppfd" > "$WORK/phot_asym.ppfd"
run phot_asym.ppfd
cp "$WORK/ppfd_map_floor.csv" "$WORK/map_asym.csv"
# 中心から等距離の 2 セル (+x 方向と +y 方向) の比 = (1 + 0.6 s)/(1 - 0.6 s)
chk "asym C=0 (+x) vs C=90"   "$(awk -v a="$(cell map_asym.csv 20 10 5)" -v b="$(cell map_asym.csv 10 20 5)" \
	'BEGIN {printf "%.9g", a / b}')" \
	"$(awk 'BEGIN {r = 20.5 / 21 - 0.5; s = r / sqrt((r * r) + 0.25);
		printf "%.9g", (1 + 0.6 * s) / (1 - 0.6 * s)}')" 2e-3
# 配光ファイルの正規化 (∫ I dOmega = 1) が厳密でなければ収支がずれる
chkabs "IES closure (normalise)" "$(logval 'closure error')" 1e-6
sed 's/ies iso\.ies/ldt asym.ldt/' "$SRC/photometry.ppfd" > "$WORK/phot_asym_ldt.ppfd"
run phot_asym_ldt.ppfd
chkabs "EULUMDAT == IES (asym)"  "$(maxdiff map_asym.csv ppfd_map_floor.csv)" 1e-6
sed 's/ies iso\.ies/ies asym.ies rot 90/' "$SRC/photometry.ppfd" > "$WORK/phot_rot.ppfd"
run phot_rot.ppfd
chkabs "rot 90 transposes map"   "$(maxdiff map_asym.csv ppfd_map_floor.csv 1)" 1e-6

# カタログ値 (PPF [umol/s] / 光束 [lm]) からの放射束換算。
# 555 nm 単色なら lumens 683 がちょうど 1 W。
sed 's/input 1\.0/input 1.0 ppf 10/' "$SRC/unit.ppfd" > "$WORK/unit_ppf.ppfd"
run unit_ppf.ppfd
chk "flux from ppf"        "$(logval 'PPF (400-700nm)')" 10.0 1e-6
sed 's/input 1\.0/input 1.0 lumens 683/' "$SRC/unit.ppfd" > "$WORK/unit_lm.ppfd"
run unit_lm.ppfd
chk "flux from lumens"     "$(logval 'radiant flux')" 1.0 1e-6

echo
echo "== (h) occluders (shelves / baffles) =="
# (h-1) 断面いっぱいの棚板で 2 室に仕切る : 下室は真っ暗、収支は閉じたまま
run shelf.ppfd
chk "shelf: absorbed by walls"  "$(logval 'absorbed by walls')" 1.0 1e-4
# 全パッチ面積 8 m2 (チャンバ 6 + 棚板の表裏 2)、rho = 0.8
chk "shelf: mean irradiance"    "$(logval 'mean wall irrad.')" \
	"$(awk 'BEGIN {printf "%.9g", 1.0 / (0.2 * 8.0)}')" 1e-4
chkabs "shelf: sealed (lower=0)" "$(awk -F, 'NR == 3 {print $7}' "$WORK/ppfd_summary.csv")" 0
# 遮蔽が完全に 2 値なので形態係数は厳密なまま (標本化に落ちた対が 0)
chkabs "shelf: exact visibility" \
	"$(awk -F'=' '/partially shadowed/ {split($2, f, ","); print f[1]}' "$WORK/ppfd.log")" 0
chkabs "shelf: ff row sum"      "$(awk -F'=' '/row sum error/ {split($2, f, ","); print f[1]}' "$WORK/ppfd.log")" 1e-6

# (h-2) 点光源 + 小さな板の本影 : 影の中は厳密に 0、影の外は逆二乗則
run shadow.ppfd
chkabs "shadow: umbra is dark"  "$(cell ppfd_map_floor.csv 5 5 5)" 0
# 影の外のセル (ix=0, iy=5) : E = (Phi/pi) h^2/(h^2+r^2)^2, h=0.9, r=0.4545455
E4=$(awk -v k="$KPH555" 'BEGIN {
	pi = 4 * atan2(1, 1); h = 0.9; r = 0.5 - 0.5 / 11;
	printf "%.9g", (1.0 / pi) * h * h / ((h * h + r * r) ^ 2) * k
}')
chk "shadow: lit cell"          "$(cell ppfd_map_floor.csv 0 5 5)" "$E4" 2e-3
# 影に隠れた分は板が吸収するので収支は閉じる (影の境界の求積誤差だけ残る)
chkabs "shadow: closure"        "$(logval 'closure error')" 1e-3

# (h-3) 2 段ラック : 棚板で仕切った段どうしは完全に独立になる。
# 下段の LED を消しても上段の分布は 1 ビットも変わらず、下段は真っ暗になる。
run rack2.ppfd
cp "$WORK/ppfd_map_tier2.csv" "$WORK/map_tier2.csv"
awk -F, 'NR == 3 {printf "%-28s tier1 (lower, with LEDs) = %.4g umol/m2/s -> %s\n", \
	"rack2 lower tier lit", $5, ($5 > 0) ? "OK" : "NG"; exit ($5 > 0) ? 0 : 1}' \
	"$WORK/ppfd_summary.csv" || status=1
sed '/^array = 0.1 0.15 0.40/d' "$SRC/rack2.ppfd" > "$WORK/rack2_top.ppfd"
run rack2_top.ppfd
chkabs "rack2: tiers independent" "$(maxdiff map_tier2.csv ppfd_map_tier2.csv)" 1e-6
chkabs "rack2: lower tier dark"   "$(awk -F, 'NR == 3 {print $7}' "$WORK/ppfd_summary.csv")" 0

echo
echo "== (i) canopy two-stream layer kernel (algebraic identities) =="
# 幾何にも入力にも依存しない代数的性質。macOS (Apple Silicon) だけ FMA の
# 有無で結果が変わった実例があるので、この種の恒等式は 3 OS で常時縛る。
if "$OPPFD" --selftest > "$WORK/selftest.log" 2>&1; then
	awk 'END {printf "%-28s -> OK (%s)\n", "canopy kernel self test", $0}' "$WORK/selftest.log"
else
	echo "canopy kernel self test         -> NG" >&2
	grep -- '-> NG' "$WORK/selftest.log" >&2 || true
	status=1
fi

echo
echo "== (j) canopy scattering returned to the cavity =="
# (j-1) omega = 0 なら散乱は無いので、leafscatter = on でも出力がビット一致する
sed 's/^photoperiod/leafscatter = on\nphotoperiod/' "$SRC/canopy.ppfd" > "$WORK/scat0.ppfd"
run canopy.ppfd
cp "$WORK/ppfd_map_below.csv" "$WORK/map_noscat.csv"
run scat0.ppfd
chkabs "scatter: omega=0 identical" "$(maxdiff map_noscat.csv ppfd_map_below.csv)" 0

# (j-2) 散乱を戻すと群落の吸収は残差ではなく実計算になるので、収支が独立に閉じる
sed -e 's/^material = blackleaf gray 0.0 0.0/material = blackleaf gray 0.25 0.25/' \
    -e 's/^photoperiod/leafscatter = on\nphotoperiod/' "$SRC/canopy.ppfd" > "$WORK/scat_w.ppfd"
run scat_w.ppfd
chkabs "scatter: closure (omega=0.5)" "$(logval 'closure error')" 1e-5
# 戻した散乱光は正でなければならない
awk -v s="$(logval 'canopy scattered')" 'BEGIN {
	printf "%-28s actual=%.4g -> %s (must be > 0)\n", "scatter: flux returned", s, (s > 0) ? "OK" : "NG";
	exit (s > 0) ? 0 : 1}' || status=1

# (j-3) omega = 1 (無吸収葉) なら群落は 1 光子も吸収しない (構造的に厳密)
sed -e 's/^material = blackleaf gray 0.0 0.0/material = blackleaf gray 0.5 0.5/' \
    -e 's/^photoperiod/leafscatter = on\nphotoperiod/' "$SRC/canopy.ppfd" > "$WORK/scat_w1.ppfd"
run scat_w1.ppfd
chkabs "scatter: omega=1 no capture" "$(logval 'absorbed by canopy')" 0
chk "scatter: omega=1 walls take all" "$(logval 'absorbed by walls')" 1.0 2e-3

echo
echo "== (k) specular walls (mirror-image method) =="
run mirror.ppfd
# E = Phi/(4 pi) (1/d1^2 + rho_s/d2^2), d1 = 0.5, d2 = 1.0, rho_s = 0.8
EM=$(awk -v k="$KPH555" 'BEGIN {pi = 4 * atan2(1, 1); printf "%.9g", (1 / (4 * pi)) * ((1 / 0.25) + 0.8) * k}')
chk "mirror: direct + image"    "$(cell ppfd_map_mid.csv 10 10 5)" "$EM" 2e-3
# 鏡面が 1 面なら鏡像は自分自身の面で折り返さないので打ち切り誤差が無い
chkabs "mirror: no truncation"  "$(trunc)" 0
# 鏡面成分は鏡像が運ぶので、壁が吸収するのは (1 - rho_d - rho_s)。収支は閉じる
chkabs "mirror: closure"        "$(logval 'closure error')" 1e-6
chk "mirror: walls take all"    "$(logval 'absorbed by walls')" 1.0 1e-6
# 鏡面反射率を 0 にすると鏡像が消え、直接光だけの解析解に戻る
sed 's/^material = mirror gray 0.0 specular 0.8/material = mirror gray 0.0 specular 0.0/' \
	"$SRC/mirror.ppfd" > "$WORK/mirror_off.ppfd"
run mirror_off.ppfd
E5=$(awk -v k="$KPH555" 'BEGIN {pi = 4 * atan2(1, 1); printf "%.9g", k / (4 * pi * 0.25)}')
chk "mirror: rho_s=0 direct only" "$(cell ppfd_map_mid.csv 10 10 5)" "$E5" 2e-3

echo
echo "== (l) specular + diffuse interreflection (extended form factors) =="
# 拡散床の放射が鏡面天井で正反射して戻る経路 (拡張形態係数)。
# これが無いと rho_s x (鏡面壁への拡散入射) が消え、収支が数 % 破れる。
run mirror2.ppfd
cp "$WORK/ppfd_map_mid.csv" "$WORK/map_spec.csv"
chkabs "specdiff: closure (rho_s=1)" "$(logval 'closure error')" 1e-5
chkabs "specdiff: no truncation"     "$(trunc)" 0
# 二重箱等価 : 完全鏡面の折り返し = 鏡映した幾何 (鏡面コードを通らない参照解)
run mirror2ref.ppfd
chkabs "specdiff: folded equivalence" "$(maxdiff map_spec.csv ppfd_map_mid.csv)" 1e-6
# 部分鏡面 (rho_s = 0.8) でも収支が閉じる (修正前は -7.7% 破れていた)
sed 's/specular 1.0/specular 0.8/' "$SRC/mirror2.ppfd" > "$WORK/mirror2p.ppfd"
run mirror2p.ppfd
chkabs "specdiff: closure (rho_s=.8)" "$(logval 'closure error')" 1e-5
# 非物理な材料 (rho_d + rho_s > 1) は入力段で弾かれる
sed 's/^material = white  gray 0.6/material = white  gray 0.6 specular 0.5/' \
	"$SRC/mirror2.ppfd" > "$WORK/mirror2bad.ppfd"
if (cd "$WORK" && "$OPPFD" mirror2bad.ppfd > /dev/null 2>&1); then
	echo "specdiff: rho_d+rho_s>1 rejected -> NG (accepted invalid material)" >&2
	status=1
else
	echo "specdiff: rho_d+rho_s>1 rejected -> OK"
fi

echo
echo "== (m) specular side wall + canopy =="
# 側壁の鏡映は z を保つので群落スラブ (z だけで決まる) は鏡映で不変。
# 折り返した経路の z プロファイル = 鏡像への直線の z プロファイルとなり、
# 群落減衰は直線評価のまま厳密。完全鏡面 (rho_s=1) の側壁は
# 「x 方向に 2 倍にした幾何 + 光源の鏡像」と厳密に等価。
run mirror3.ppfd
cp "$WORK/ppfd_map_mid.csv" "$WORK/map_m3.csv"
M3AVG=$(sumval 5)
run mirror3ref.ppfd
chkabs "speccan: folded equivalence" "$(maxdiff map_m3.csv ppfd_map_mid.csv)" 1e-6
# 鏡を外すと床の平均 PPFD は下がる (鏡像の寄与の符号判定)
sed 's/specular 1.0/specular 0.0/' "$SRC/mirror3.ppfd" > "$WORK/mirror3off.ppfd"
run mirror3off.ppfd
awk -v m="$M3AVG" -v o="$(sumval 5)" 'BEGIN {
	printf "%-28s mirror=%.6g  off=%.6g -> %s (mirror must add light)\n", \
		"speccan: mirror adds PPFD", m, o, (m > o) ? "OK" : "NG";
	exit (m > o) ? 0 : 1
}' || status=1
# 側壁の鏡面は leafscatter と併用できる (検証 (v))。収支が閉じること
sed 's/^photoperiod/leafscatter = on\nphotoperiod/' "$SRC/mirror3.ppfd" > "$WORK/mirror3ls.ppfd"
run mirror3ls.ppfd
chkabs "speccan: +leafscatter closure" "$(logval 'closure error')" 1e-5

echo
echo "== (n) specular ceiling + canopy (folded path attenuation) =="
# z の鏡映は群落スラブを動かすので、鏡像への直線をそのまま測ると誤る。
# 鏡像への直線は箱型ビリヤードの展開なので、実経路の z は周期 2Lz の
# 三角波で厳密に得られる (canopy_path)。ztop = Lz にしてあるので効果は
# 最大 : この折り返しを落とすと下の等価判定が 8.2e-2 まで悪化する。
#
# (1) 直接光だけの解析解。壁をすべて黒にし測定面を 1 セル (中心 0.5,0.5) に
#     すると、残るのは実光源 (z=0.75) と鏡像 (z=1.25) の 2 本だけ :
#       E = Phi/(4 pi) * [ 1/0.70^2 * exp(-k*0.25) + 1/1.20^2 * exp(-k*0.75) ]
#     鏡像側の経路長 0.75 が折り返しの値 (素の直線評価なら 0.50)。
#     k = G * a * sqrt(1-omega) = 0.5 * (1.0/0.5) * sqrt(0.8) = 0.8944272
#     PPFD = E * 4.6394372 = 0.7335746 umol/m2/s
sed -e 's/^wall = zmin white/wall = zmin black/' \
    -e 's/^target      = base 0.05 12 6/target      = base 0.05 1 1/' \
	"$SRC/mirror4.ppfd" > "$WORK/mirror4d.ppfd"
run mirror4d.ppfd
chk "speccan-z: folded direct" "$(sumval 5)" 0.7335746 1e-6
# (2) 完全鏡面の天井は「z 方向に 2 倍にした幾何 + 光源の鏡像」と厳密に等価。
#     群落は自分の鏡像と合体して 1 枚のスラブになる (LAI は厚みに比例)。
run mirror4.ppfd
cp "$WORK/ppfd_map_base.csv" "$WORK/map_m4.csv"
M4AVG=$(sumval 5)
run mirror4ref.ppfd
chkabs "speccan-z: folded equivalence" "$(maxdiff map_m4.csv ppfd_map_base.csv)" 1e-6
# (3) 鏡を外すと床の平均 PPFD は下がる (鏡像の寄与の符号判定)
sed 's/specular 1.0/specular 0.0/' "$SRC/mirror4.ppfd" > "$WORK/mirror4off.ppfd"
run mirror4off.ppfd
awk -v m="$M4AVG" -v o="$(sumval 5)" 'BEGIN {
	printf "%-28s mirror=%.6g  off=%.6g -> %s (mirror must add light)\n", \
		"speccan-z: mirror adds PPFD", m, o, (m > o) ? "OK" : "NG";
	exit (m > o) ? 0 : 1
}' || status=1
# 散乱の還流だけは鏡像経路ぶんの層への預け入れが要るので、まだ弾かれる
# 鏡面の天井 + leafscatter (検証 (x))。群落の上面が鏡面の天井に乗る配置で収支が閉じる
sed 's/^photoperiod/leafscatter = on\nphotoperiod/' "$SRC/mirror4.ppfd" > "$WORK/mirror4ls.ppfd"
run mirror4ls.ppfd
chkabs "speccan-z: +leafscatter closure" "$(logval 'closure error')" 1e-6

echo
echo "== (o) stacked canopy slabs (multi-tier racks) =="
# canopy 行を複数書くと段ごとに群落を置ける。減衰は段ごとの寄与の和
#   tau = Σ_j G_j a_j s_j
# になり、鉛直光線 (s_j = 各段の厚み) では Σ_j G_j LAI_j に一致する。
grep -v '^canopy' "$SRC/canopy2.ppfd" > "$WORK/canopy2_none.ppfd"
run canopy2_none.ppfd
C0=$(cell ppfd_map_below.csv 10 10 5)
C1=$(cell ppfd_map_below.csv 15 10 5)
run canopy2.ppfd
# 中心セル : tau = 0.5*1.5 + 1.0*0.8 = 1.55
chk "stacked: vertical transmit" \
	"$(awk -v a="$(cell ppfd_map_below.csv 10 10 5)" -v b="$C0" 'BEGIN {printf "%.9g", a / b}')" \
	"$(awk 'BEGIN {printf "%.9g", exp(-1.55)}')" 1e-3
# 斜めセル : どちらの段も経路長が 1/cos(theta) 倍になる
chk "stacked: oblique transmit" \
	"$(awk -v a="$(cell ppfd_map_below.csv 15 10 5)" -v b="$C1" 'BEGIN {printf "%.9g", a / b}')" \
	"$(awk 'BEGIN {h = 0.8; r = 15.5 / 21 - 0.5; printf "%.9g", exp(-1.55 * sqrt(h * h + r * r) / h)}')" 1e-3
# 1 枚のスラブを 2 枚に割っても (葉面積密度が同じなら) 結果は変わらない。
# rack.ppfd は白壁 (rho = 0.9) + 実測相当の葉なので、直接光だけでなく
# パッチ間の群落透過率 (plen) と相互反射を通した経路も一緒に縛れる。
run rack.ppfd
cp "$WORK/ppfd_map_canopybase.csv" "$WORK/map_r1.csv"
sed 's/^canopy = 0.25 0.05 3.0 0.5 leaf/canopy = 0.25 0.15 1.5 0.5 leaf\ncanopy = 0.15 0.05 1.5 0.5 leaf/' \
	"$SRC/rack.ppfd" > "$WORK/rack_split.ppfd"
run rack_split.ppfd
chkabs "stacked: split == single" "$(maxdiff map_r1.csv ppfd_map_canopybase.csv)" 1e-12
# 重なったスラブは同じ葉を二重に数えるので入力段で弾かれる
sed 's/^canopy = 0.5 0.3 0.8 1.0 blackleaf/canopy = 0.8 0.3 0.8 1.0 blackleaf/' \
	"$SRC/canopy2.ppfd" > "$WORK/canopy2ov.ppfd"
if (cd "$WORK" && "$OPPFD" canopy2ov.ppfd > /dev/null 2>&1); then
	echo "stacked: overlap rejected -> NG (accepted overlapping slabs)" >&2
	status=1
else
	echo "stacked: overlap rejected -> OK"
fi
# 段ごとに違う葉材料は経路長を 1 個にまとめられないので弾かれる
sed -e 's/^material = black     gray 0.0/material = black     gray 0.0\nmaterial = greenleaf gray 0.1 0.1/' \
    -e 's/^canopy = 0.5 0.3 0.8 1.0 blackleaf/canopy = 0.5 0.3 0.8 1.0 greenleaf/' \
	"$SRC/canopy2.ppfd" > "$WORK/canopy2mat.ppfd"
if (cd "$WORK" && "$OPPFD" canopy2mat.ppfd > /dev/null 2>&1); then
	echo "stacked: mixed leaf rejected -> NG (accepted two leaf materials)" >&2
	status=1
else
	echo "stacked: mixed leaf rejected -> OK"
fi
# 棚板の無い複数スラブ + leafscatter : 向かい合う面どうしで拡散光をやりとりする (検証 (y))。
# 黒い壁でも反射の反復が回る (スラブどうしのやりとりを収束判定に入れてある)
sed 's/^photoperiod/leafscatter = on\nphotoperiod/' "$SRC/canopy2.ppfd" > "$WORK/canopy2ls.ppfd"
run canopy2ls.ppfd
chkabs "stacked: +leafscatter closure" "$(logval 'closure error')" 1e-6

echo
echo "== (p) two specular faces (commuting mirrors) =="
# 直交する面の鏡映は可換 (M_zmax M_xmax = M_xmax M_zmax) なので、折り返した
# 面の並びで鏡像を数えると角の鏡像を 2 回数えてしまう。軸ごとの展開指数で
# 数えるとこれが起きない。壁を黒くした変種で直接光の解析解と比べる :
#   鏡像は 3 個 (1.5,0.5,0.5) (0.5,0.5,1.5) (1.5,0.5,1.5)
#   E = Phi/(4 pi) Σ cos/d^2 -> PPFD = 2.7155297 umol/m2/s
# 面の並びで数えると角が 2 回入って 2.8170250 (+3.7%) になる。
sed -e 's/^wall = ymin white/wall = ymin black/' \
    -e 's/^target      = s 0.1 12 6/target      = s 0.1 1 1/' \
	"$SRC/mirror5.ppfd" > "$WORK/mirror5d.ppfd"
run mirror5d.ppfd
chk "twomirror: 3 images" "$(sumval 5)" 2.7155297 1e-6
# 反対側の面が黒なので 3 段目の経路は消える = 打ち切り誤差は厳密に 0
chkabs "twomirror: no truncation" "$(trunc)" 0
# 鏡面成分は鏡像が運ぶので、拡散壁を戻しても収支は閉じる
run mirror5.ppfd
chkabs "twomirror: closure" "$(logval 'closure error')" 1e-6
cp "$WORK/ppfd_map_s.csv" "$WORK/map_m5.csv"
# 完全鏡面の 2 面 = x と z の両方で 2 倍にした幾何 (角の変換も 1 対 1)
run mirror5ref.ppfd
chkabs "twomirror: folded equivalence" "$(maxdiff map_m5.csv ppfd_map_s.csv)" 1e-6

echo
echo "== (q) specular walls + occluders (folded shadow test) =="
# 鏡像への直線は箱型ビリヤードの展開なので、実経路は展開空間に並んだ
# チャンバの複製を順に通る。遮蔽判定は「通りうる複製それぞれにある棚板の
# 鏡像」と照合すれば厳密になる (occ_blocked_img)。mirror6.ppfd の板は
# 反射光の 1 本目の脚だけを遮り、鏡像への直線は実在の板に当たらないので、
# 素の直線判定では影が消える (+18.7%)。
# (1) 板と床を黒くした変種 : 残るのは実光源の 1 本だけ
sed -e 's/^wall = zmin white/wall = zmin black/' -e 's/ white div 4 4/ black div 4 4/' \
    -e 's/^target      = s 0.1 12 6/target      = s 0.1 1 1/' \
	"$SRC/mirror6.ppfd" > "$WORK/mirror6d.ppfd"
run mirror6d.ppfd
chk "specocc: reflected leg shadowed" "$(sumval 5)" 0.5267227 1e-6
# (2) 板にも拡散反射させて収支と二重箱等価 (板の鏡像を実在の 2 枚目として置く)。
# 遮蔽物があるときの収支は影の境界の求積で 1e-4 台が下限 (形態係数の問題ではない)
run mirror6.ppfd
chkabs "specocc: closure" "$(logval 'closure error')" 1e-3
cp "$WORK/ppfd_map_s.csv" "$WORK/map_m6.csv"
run mirror6ref.ppfd
chkabs "specocc: folded equivalence" "$(maxdiff map_m6.csv ppfd_map_s.csv)" 1e-6
# (3) 鏡面 2 面 (xmax, zmax) + 板 : 角の鏡像への直線は 4 つの複製を通るので、
# 参照解は 2 x 1 x 2 の箱に板 4 枚。mirror5 / mirror5ref に板を足して作る。
# 板の高さ 0.61 はかすめにくい配置として選んだもの。下の 0.57 の判定も参照。
sed 's/^wall = zmax mirror/wall = zmax mirror\noccluder = plate 0.7 0.9 0.3 0.7 0.61 0.61 white div 4 4/' \
	"$SRC/mirror5.ppfd" > "$WORK/mirror5o.ppfd"
sed 's/^wall = ymin white/wall = ymin white\noccluder = p00 0.7 0.9 0.3 0.7 0.61 0.61 white div 4 4\noccluder = p10 1.1 1.3 0.3 0.7 0.61 0.61 white div 4 4\noccluder = p01 0.7 0.9 0.3 0.7 1.39 1.39 white div 4 4\noccluder = p11 1.1 1.3 0.3 0.7 1.39 1.39 white div 4 4/' \
	"$SRC/mirror5ref.ppfd" > "$WORK/mirror5oref.ppfd"
run mirror5o.ppfd
chkabs "specocc2: closure" "$(logval 'closure error')" 1e-3
cp "$WORK/ppfd_map_s.csv" "$WORK/map_m5o.csv"
run mirror5oref.ppfd
chkabs "specocc2: folded equivalence" "$(maxdiff map_m5o.csv ppfd_map_s.csv)" 1e-6
# 板の高さ 0.57 : 直接光の影の境界判定を「5% 内側の 4 隅 + 中心」の 5 点
# (有理数の位置) でしていたころは、判定点からの光線が t = 0.4 でちょうど板の辺
# y = 0.3 に乗り、参照解の鏡映の複製どうしで細分の有無が割れて 3.5e-5 ずれた。
# 判定を粗い求積の Gauss 点 (√(3/5) を含む無理数の位置) で行うようにしてから 0。
sed 's/ 0\.61 0\.61 / 0.57 0.57 /' "$WORK/mirror5o.ppfd" > "$WORK/mirror5g.ppfd"
sed 's/ 0\.61 0\.61 / 0.57 0.57 /g; s/ 1\.39 1\.39 / 1.43 1.43 /g' "$WORK/mirror5oref.ppfd" > "$WORK/mirror5gref.ppfd"
run mirror5g.ppfd
cp "$WORK/ppfd_map_s.csv" "$WORK/map_m5g.csv"
run mirror5gref.ppfd
chkabs "specocc2: grazing height 0.57" "$(maxdiff map_m5g.csv ppfd_map_s.csv)" 1e-6

echo
echo "== (r) tilted panels (arbitrary parallelogram plates) =="
# tiltpanel.ppfd : 45 度傾けた黒い板。床の 3 点の直接光 (壁も黒なので間接光なし)
run tiltpanel.ppfd
# a : 光線が板の中心を貫く -> 厳密に 0
chkabs "panel: shadowed point" "$(sumval 5)" 0
# b : 光線は板の外接箱を通り抜けるが板の平面とは板の外で交わる -> 遮られない。
#     箱として扱う実装ならここが 0 になる
chk "panel: through bounding box" "$(awk -F, 'NR == 3 {print $5; exit}' "$WORK/ppfd_summary.csv")" 0.4127722 1e-6
chk "panel: clear point" "$(awk -F, 'NR == 4 {print $5; exit}' "$WORK/ppfd_summary.csv")" 0.5267227 1e-6
# 白い壁 (rho 0.5) と白い板 : 閉キャビティの収支。遮蔽物の影の境界は可視率の
# 標本化なので 1e-3 台が下限 (同じ位置の軸並行の板でも同じ桁 : 3.7e-4)
sed -e 's/^wall = all black/wall = all white/' -e 's/ black div 4 4/ white div 4 4/' \
    -e 's/^material = white gray 0.8/material = white gray 0.5/' \
    -e '/^target/d' -e 's/^photoperiod/target = floor 0.1 12 12\nphotoperiod/' \
	"$SRC/tiltpanel.ppfd" > "$WORK/panelw.ppfd"
run panelw.ppfd
chkabs "panel: closure (rho 0.5)" "$(logval 'closure error')" 2e-3
cp "$WORK/ppfd_map_floor.csv" "$WORK/map_pw.csv"
# 板を x = y で鏡映すると (光源も箱も x = y について対称) 床の分布は転置になる
sed 's/^panel = tilt 0.3 0.3 0.4  0.5 0.3 0.6  0.3 0.7 0.4/panel = tilt 0.3 0.3 0.4  0.3 0.5 0.6  0.7 0.3 0.4/' \
	"$WORK/panelw.ppfd" > "$WORK/panelwt.ppfd"
run panelwt.ppfd
chkabs "panel: x=y mirror -> transpose" "$(maxdiff map_pw.csv ppfd_map_floor.csv 1)" 1e-9
# 水平な panel は厚さ 0 の occluder と同じ (mirror6.ppfd の板を置き換え。鏡面の
# 折り返し判定 occ_blocked_img の板の側も通る)
run mirror6.ppfd
cp "$WORK/ppfd_map_s.csv" "$WORK/map_m6b.csv"
sed 's/^occluder = plate 0.7 0.9 0.3 0.7 0.57 0.57 white div 4 4/panel = plate 0.7 0.3 0.57  0.9 0.3 0.57  0.7 0.7 0.57 white div 4 4/' \
	"$SRC/mirror6.ppfd" > "$WORK/mirror6p.ppfd"
run mirror6p.ppfd
chkabs "panel: horizontal == occluder" "$(maxdiff map_m6b.csv ppfd_map_s.csv)" 1e-9
# 同じ板を 2 通りに書いても (panel / occluder) 結果がそろうこと。ちょうど中間の
# 高さ (z = 0.5) の板は、有理数の格子の標本だと床と天井の標本を結ぶ光線の族が
# 板の辺にまとめて乗り、交差判定の丸めの違いで可視率が跳んでいた (7.2e-4)。
# 標本位置を無理数 (受光側 : 黄金比、放射側 : √2) にずらしてからは 1e-5 台。
sed 's/^panel = tilt .*/panel = tilt 0.3 0.3 0.5  0.5 0.3 0.5  0.3 0.7 0.5  white div 4 4/' \
	"$WORK/panelw.ppfd" > "$WORK/panelmid.ppfd"
sed 's/^panel = tilt .*/occluder = tilt 0.3 0.5 0.3 0.7 0.5 0.5 white div 4 4/' \
	"$WORK/panelw.ppfd" > "$WORK/occmid.ppfd"
run panelmid.ppfd
cp "$WORK/ppfd_map_floor.csv" "$WORK/map_pmid.csv"
run occmid.ppfd
chkabs "visibility: panel == occluder" "$(maxdiff map_pmid.csv ppfd_map_floor.csv)" 1e-4
chkabs "visibility: mid-height closure" "$(logval 'closure error')" 1e-3
# 鏡面 + 斜めの板 : 板の鏡像 (x -> 2 - x) を実在の 2 枚目として置いた参照解と比べる
sed 's/^occluder = plate .*/panel = tilt 0.71 0.33 0.52  0.88 0.33 0.67  0.71 0.69 0.52  white div 4 4/' \
	"$SRC/mirror6.ppfd" > "$WORK/mirror6t.ppfd"
sed -e 's/^occluder = plate .*/panel = tilt 0.71 0.33 0.52  0.88 0.33 0.67  0.71 0.69 0.52  white div 4 4/' \
    -e 's/^occluder = plate2 .*/panel = tilt2 1.29 0.33 0.52  1.12 0.33 0.67  1.29 0.69 0.52  white div 4 4/' \
	"$SRC/mirror6ref.ppfd" > "$WORK/mirror6tref.ppfd"
run mirror6t.ppfd
chkabs "panel: specular closure" "$(logval 'closure error')" 1e-3
cp "$WORK/ppfd_map_s.csv" "$WORK/map_m6t.csv"
run mirror6tref.ppfd
chkabs "panel: specular folded equiv." "$(maxdiff map_m6t.csv ppfd_map_s.csv)" 1e-6
# 退化した板 (3 隅が一直線) とチャンバからはみ出す板は入力段で弾かれる
sed 's/^panel = tilt .*/panel = tilt 0.3 0.3 0.4  0.5 0.3 0.6  0.7 0.3 0.8  black/' "$SRC/tiltpanel.ppfd" > "$WORK/paneldeg.ppfd"
if (cd "$WORK" && "$OPPFD" paneldeg.ppfd > /dev/null 2>&1); then
	echo "panel: degenerate rejected -> NG (accepted a collinear panel)" >&2
	status=1
else
	echo "panel: degenerate rejected -> OK"
fi
sed 's/^panel = tilt .*/panel = tilt 0.3 0.3 0.4  0.5 0.3 0.6  0.3 1.2 0.4  black/' "$SRC/tiltpanel.ppfd" > "$WORK/panelout.ppfd"
if (cd "$WORK" && "$OPPFD" panelout.ppfd > /dev/null 2>&1); then
	echo "panel: outside chamber rejected -> NG (accepted a panel outside the box)" >&2
	status=1
else
	echo "panel: outside chamber rejected -> OK"
fi
# 鏡面を運ぶのは壁だけ。棚板・板に鏡面材料を付けると ρs が帳簿から消える
# (-1.4%) ので入力段で弾かれる (遮蔽物にもともとあった穴を同時に塞いだ)
sed 's/^occluder = plate 0.7 0.9 0.3 0.7 0.57 0.57 white div 4 4/occluder = plate 0.7 0.9 0.3 0.7 0.57 0.57 mirror div 4 4/' \
	"$SRC/mirror6.ppfd" > "$WORK/occspec.ppfd"
if (cd "$WORK" && "$OPPFD" occspec.ppfd > /dev/null 2>&1); then
	echo "plate: specular material rejected -> NG (accepted a specular occluder)" >&2
	status=1
else
	echo "plate: specular material rejected -> OK"
fi

echo "== (s) oriented sensors (vertical / downward measuring planes) =="
# sensor.ppfd : 黒いチャンバの中心の点光源。どの向きの面も直接光の解析解
#   E = Phi/(4 pi) (n . (S - x)) / d^3   (n . (S - x) > 0 のときだけ)
run sensor.ppfd
chk "sensor: vertical (-x)"      "$(sumrow wall 5)" 1.6510888 1e-6
chkabs "sensor: downward, below" "$(sumrow down 5)" 0
chk "sensor: downward, above"    "$(sumrow ceil 5)" 4.1021631 1e-6
# 同じ水平面を target と sensor で書いても分布がそろう。rack.ppfd (白壁・
# 棚板・群落) に leafscatter を入れ、群落の内部 (z = 0.06) で比べる
sed 's/^target = canopybase 0.06 25 13/target = canopybase 0.06 25 13\nsensor = cb2 0 0 0.06  1.2 0 0.06  0 0.6 0.06  25 13\nleafscatter = on/' \
	"$SRC/rack.ppfd" > "$WORK/rack_sensor.ppfd"
run rack_sensor.ppfd
chkabs "sensor: == target (in canopy)" "$(maxdiff ppfd_map_canopybase.csv ppfd_map_cb2.csv)" 1e-12
# 群落の内部の傾いた面 : 二流近似では各半球の放射輝度が一様なので
#   E = F_down (1 + nz)/2 + F_up (1 - nz)/2
# 光源を 3 つの面すべての平面上に置き (直接光は cos = 0 で 0)、壁を黒く
# すると、残るのは群落の拡散場だけ : 上向き = F_down、下向き = F_up、
# 鉛直 = その平均。要約 CSV は 6 桁なので 2e-6 で判定する
# (F_down だけを使う誤った実装なら 1e-3 ずれる)
sed -e 's/^material = black gray 0.0/material = black gray 0.0\nmaterial = leaf  gray 0.4 0.4/' \
    -e 's/^led = 0.5 0.5 0.5 1.0 na555 iso/led = 0.5 0.2 0.5 1.0 na555 iso\ncanopy = 0.7 0.3 2.0 0.5 leaf\nleafscatter = on/' \
    -e '/^sensor = /d' \
    -e 's/^target = floor 0.1 1 1/target = floor 0.1 1 1\nsensor = up   0.4 0.4 0.5  0.6 0.4 0.5  0.4 0.6 0.5  1 1\nsensor = dn   0.4 0.4 0.5  0.4 0.6 0.5  0.6 0.4 0.5  1 1\nsensor = vert 0.5 0.4 0.4  0.5 0.4 0.6  0.5 0.6 0.4  1 1/' \
	"$SRC/sensor.ppfd" > "$WORK/sensor_canopy.ppfd"
run sensor_canopy.ppfd
chk "sensor: tilted in canopy" "$(sumrow vert 5)" \
	"$(awk -v u="$(sumrow up 5)" -v d="$(sumrow dn 5)" 'BEGIN {printf "%.9g", (u + d) / 2}')" 2e-6
# はみ出す面・退化した面は入力段で弾かれる
sed 's/^sensor = wall .*/sensor = wall  0.9 0.4 0.2  1.3 0.4 0.4  0.9 0.6 0.2  1 1/' "$SRC/sensor.ppfd" > "$WORK/sensorout.ppfd"
if (cd "$WORK" && "$OPPFD" sensorout.ppfd > /dev/null 2>&1); then
	echo "sensor: outside chamber rejected -> NG (accepted a sensor outside the box)" >&2
	status=1
else
	echo "sensor: outside chamber rejected -> OK"
fi
sed 's/^sensor = wall .*/sensor = wall  0.9 0.4 0.2  0.9 0.4 0.4  0.9 0.4 0.6  1 1/' "$SRC/sensor.ppfd" > "$WORK/sensordeg.ppfd"
if (cd "$WORK" && "$OPPFD" sensordeg.ppfd > /dev/null 2>&1); then
	echo "sensor: degenerate rejected -> NG (accepted a collinear sensor)" >&2
	status=1
else
	echo "sensor: degenerate rejected -> OK"
fi

echo "== (t) diffuse transmittance through a canopy slab (2 E3) =="
# difftrans.ppfd : 側壁 4 面が完全鏡面の平行平板。天井 (白) の拡散光が
# 黒葉の群落を透過して黒い床に届く割合は 2 E3(tau0) で厳密。
# パッチ重心どうしの経路長で減衰させた以前の実装は -0.34%〜-1.47%
run difftrans.ppfd
chk "difftrans: tau0=1 floor"     "$(logval 'absorbed by walls')" 0.2193839 1e-4
chk "difftrans: tau0=1 target"    "$(sumval 14)" 0.2193839 1e-4
sed 's/^material = blackleaf gray 0.0 0.0/material = blackleaf gray 0.3 0.45/' \
	"$SRC/difftrans.ppfd" > "$WORK/difftrans05.ppfd"
run difftrans05.ppfd
chk "difftrans: tau0=0.5 floor"   "$(logval 'absorbed by walls')" 0.4432087 5e-4
sed 's/^canopy = 0.25 0.05 2.0/canopy = 0.25 0.05 6.0/' \
	"$SRC/difftrans.ppfd" > "$WORK/difftrans3.ppfd"
run difftrans3.ppfd
chk "difftrans: tau0=3 floor"     "$(logval 'absorbed by walls')" 0.0178613 5e-4
# leafscatter : 層への預け入れを透過率と同じ標本で作るので、白壁の
# キャビティでも収支が閉じる (canopy.ppfd を白壁・散乱葉にした変種)
sed -e 's/^material = blackleaf gray 0.0 0.0/material = blackleaf gray 0.25 0.25\nmaterial = white gray 0.7/' \
    -e 's/^wall     = all black/wall = all white/' \
    -e 's/^photoperiod/leafscatter = on\nphotoperiod/' "$SRC/canopy.ppfd" > "$WORK/scat_white.ppfd"
run scat_white.ppfd
chkabs "difftrans: scatter closure" "$(logval 'closure error')" 1e-5

echo "== (u) canopy over growing beds only (finite horizontal extent) =="
# canopybed.ppfd : x 0.3..0.7, y 0.3..0.7, z 0.4..0.8 の箱だけに黒葉の群落。
# 天井中心の光源から 3 点へ向かう光線の箱の中の長さを解析的に求め、
# canopy 行を消した変種との比 (透過率) で判定する
run canopybed.ppfd
bc_c=$(sumrow c 5); bc_e=$(sumrow e 5); bc_k=$(sumrow k 5)
sed '/^canopy/d' "$SRC/canopybed.ppfd" > "$WORK/canopybed0.ppfd"
run canopybed0.ppfd
chk "bed: vertical through"   "$(awk -v a="$bc_c" -v b="$(sumrow c 5)" 'BEGIN {printf "%.9g", a / b}')" 0.3678794 1e-5
chk "bed: exits side face"    "$(awk -v a="$bc_e" -v b="$(sumrow e 5)" 'BEGIN {printf "%.9g", a / b}')" 0.5717708 1e-5
chk "bed: exits edge"         "$(awk -v a="$bc_k" -v b="$(sumrow k 5)" 'BEGIN {printf "%.9g", a / b}')" 0.5420633 1e-5
# 範囲を断面全体にすると従来のスラブと同じ経路になる (ビット一致)。
# rack.ppfd (白壁・棚板) の群落を 2 枚のベッドに割っても結果が変わらない
run rack.ppfd
cp "$WORK/ppfd_map_canopybase.csv" "$WORK/map_bed0.csv"
sed 's/^canopy = 0.25 0.05 3.0 0.5 leaf/canopy = 0.25 0.05 3.0 0.5 leaf  0 1.2 0 0.6/' \
	"$SRC/rack.ppfd" > "$WORK/rack_bedfull.ppfd"
run rack_bedfull.ppfd
chkabs "bed: full extent == slab" "$(maxdiff map_bed0.csv ppfd_map_canopybase.csv)" 0
sed 's/^canopy = 0.25 0.05 3.0 0.5 leaf/canopy = 0.25 0.05 3.0 0.5 leaf  0 0.6 0 0.6\ncanopy = 0.25 0.05 3.0 0.5 leaf  0.6 1.2 0 0.6/' \
	"$SRC/rack.ppfd" > "$WORK/rack_beds.ppfd"
run rack_beds.ppfd
chkabs "bed: two beds == slab" "$(maxdiff map_bed0.csv ppfd_map_canopybase.csv)" 1e-9
# 鏡面 2 面 (xmax, zmax) + ベッド : 鏡像への直線を通りうる複製それぞれの
# 箱の鏡像と照合する。参照解は 2 x 1 x 2 の箱に光源 4 個と箱 4 個
run mirror7.ppfd
cp "$WORK/ppfd_map_s.csv" "$WORK/map_m7.csv"
run mirror7ref.ppfd
chkabs "bed: specular folded equiv." "$(maxdiff map_m7.csv ppfd_map_s.csv)" 1e-6
# 重なったベッド・はみ出すベッドは入力段で弾かれる (leafscatter との併用は検証 (aa))
for v in ov out; do
	case $v in
	ov)  sed 's/^canopy = .*/canopy = 0.8 0.4 2.0 0.5 blackleaf  0.3 0.7 0.3 0.7\ncanopy = 0.6 0.5 1.0 0.5 blackleaf  0.6 0.9 0.3 0.7/' \
		"$SRC/canopybed.ppfd" > "$WORK/canopybed_$v.ppfd" ;;
	out) sed 's/^canopy = .*/canopy = 0.8 0.4 2.0 0.5 blackleaf  0.3 1.2 0.3 0.7/' \
		"$SRC/canopybed.ppfd" > "$WORK/canopybed_$v.ppfd" ;;
	esac
	if (cd "$WORK" && "$OPPFD" "canopybed_$v.ppfd" > /dev/null 2>&1); then
		echo "bed: invalid ($v) rejected -> NG (accepted)" >&2
		status=1
	else
		echo "bed: invalid ($v) rejected -> OK"
	fi
done

echo "== (v) specular side walls + leafscatter =="
# mirror8.ppfd : 側壁 2 面 (x = 1, y = 1) が完全鏡面、散乱葉の群落。
# 鏡面のパッチで終わる光線は (1 - rho_s) の分だけ預け、rho_s の分は
# 鏡像経路が展開した全経路で預ける。群落の上下面も鏡面で折り返す
run mirror8.ppfd
chkabs "specscat: closure" "$(logval 'closure error')" 1e-6
m8c=$(logval 'absorbed by canopy')
for t in below above mid; do cp "$WORK/ppfd_map_$t.csv" "$WORK/map_m8$t.csv"; done
# 完全鏡面の 2 面 = x と y で 2 倍にした幾何 (光源 4 個、群落は 2 x 2 に広がる)
run mirror8ref.ppfd
chkabs "specscat: equiv. below canopy" "$(maxdiff map_m8below.csv ppfd_map_below.csv)" 1e-6
chkabs "specscat: equiv. above canopy" "$(maxdiff map_m8above.csv ppfd_map_above.csv)" 1e-6
chkabs "specscat: equiv. in canopy"    "$(maxdiff map_m8mid.csv ppfd_map_mid.csv)" 1e-6
chk "specscat: canopy absorption" "$m8c" "$(awk -v a="$(logval 'absorbed by canopy')" 'BEGIN {printf "%.9g", a / 4}')" 1e-5
# 部分的な鏡面 (rho_d 0.05 + rho_s 0.9) と白壁でも収支が閉じる
sed -e 's/^material = mirror gray 0.0 specular 1.0/material = mirror gray 0.05 specular 0.9/' \
    -e 's/^wall = all  black/wall = all  white/' "$SRC/mirror8.ppfd" > "$WORK/mirror8p.ppfd"
run mirror8p.ppfd
chkabs "specscat: partial mirror closure" "$(logval 'closure error')" 1e-5

echo "== (w) multi-tier racks with leafscatter (one two-stream column per tier) =="
# rack3.ppfd : 断面いっぱいの棚板で仕切った 2 段、各段に光源と散乱葉の群落。
# 段ごとに独立なカラムなので、各段を独立なチャンバとして解いた結果と一致する
run rack3.ppfd
r3c=$(logval 'absorbed by canopy'); r3w=$(logval 'absorbed by walls')
for t in a_top a_in b_top b_in; do cp "$WORK/ppfd_map_$t.csv" "$WORK/map_r3$t.csv"; done
run rack3a.ppfd
rac=$(logval 'absorbed by canopy'); raw=$(logval 'absorbed by walls')
chkabs "tiers: lower == alone (top)"  "$(maxdiff map_r3a_top.csv ppfd_map_a_top.csv)" 1e-6
chkabs "tiers: lower == alone (in)"   "$(maxdiff map_r3a_in.csv ppfd_map_a_in.csv)" 1e-6
sed -e 's/^wall = zmin tray/wall = zmin white/' -e 's/^led = .*/led = 0.6 0.4 0.45 0.8 white4000 beam 1/' \
    -e 's/^canopy = .*/canopy = 0.25 0.05 1.5 0.5 leaf/' \
    -e 's/^target = a_top  0.30/target = b_top  0.30/' -e 's/^target = a_in   0.15/target = b_in   0.15/' \
	"$SRC/rack3a.ppfd" > "$WORK/rack3b.ppfd"
run rack3b.ppfd
rbc=$(logval 'absorbed by canopy'); rbw=$(logval 'absorbed by walls')
chkabs "tiers: upper == alone (top)"  "$(maxdiff map_r3b_top.csv ppfd_map_b_top.csv)" 1e-6
chkabs "tiers: upper == alone (in)"   "$(maxdiff map_r3b_in.csv ppfd_map_b_in.csv)" 1e-6
chk "tiers: canopy abs = sum" "$r3c" "$(awk -v a="$rac" -v b="$rbc" 'BEGIN {printf "%.9g", a + b}')" 1e-5
chk "tiers: wall abs = sum"   "$r3w" "$(awk -v a="$raw" -v b="$rbw" 'BEGIN {printf "%.9g", a + b}')" 1e-5
# 段ごとに違う葉材料 (leafscatter = on のとき) : 上段だけ別の葉にしても、
# 上段は「その葉で単独に解いたもの」と一致する (ω は段ごとのカラムにしか現れない)
sed -e 's/^material = leaf  table .*/&\nmaterial = leaf2 table 400:0.10:0.05 550:0.25:0.20 660:0.08:0.04 730:0.40:0.40/' \
    -e 's/^canopy = 0.75 0.55 1.5 0.5 leaf/canopy = 0.75 0.55 1.5 0.5 leaf2/' "$SRC/rack3.ppfd" > "$WORK/rack3l2.ppfd"
run rack3l2.ppfd
for t in b_top b_in; do cp "$WORK/ppfd_map_$t.csv" "$WORK/map_r3l2$t.csv"; done
sed -e 's/^canopy = 0.25 0.05 1.5 0.5 leaf/canopy = 0.25 0.05 1.5 0.5 leaf2/' \
    -e 's/^material = leaf  table .*/&\nmaterial = leaf2 table 400:0.10:0.05 550:0.25:0.20 660:0.08:0.04 730:0.40:0.40/' \
	"$WORK/rack3b.ppfd" > "$WORK/rack3bl2.ppfd"
run rack3bl2.ppfd
chkabs "tiers: per-tier leaf (top)" "$(maxdiff map_r3l2b_top.csv ppfd_map_b_top.csv)" 1e-6
chkabs "tiers: per-tier leaf (in)"  "$(maxdiff map_r3l2b_in.csv ppfd_map_b_in.csv)" 1e-6
# 散乱なしでは段ごとの葉材料は経路長 1 本にまとめられないので弾く
sed '/^leafscatter/d' "$WORK/rack3l2.ppfd" > "$WORK/rack3l2ns.ppfd"
if (cd "$WORK" && "$OPPFD" rack3l2ns.ppfd > /dev/null 2>&1); then
	echo "tiers: mixed leaf w/o scatter rejected -> NG (accepted)" >&2
	status=1
else
	echo "tiers: mixed leaf w/o scatter rejected -> OK"
fi
# 棚板の上へ散乱光が漏れない : 上段の光源と群落を外すと上段の測定面は 0
sed -e '/^led = 0.6 0.4 0.95/d' -e '/^canopy = 0.75/d' "$SRC/rack3.ppfd" > "$WORK/rack3leak.ppfd"
run rack3leak.ppfd
chkabs "tiers: no leak above shelf" "$(sumrow b_top 5)" 0
# 棚板を外すと 2 段は向かい合う面どうしでやりとりする (検証 (y))。収支が閉じること
sed '/^occluder = shelf/d' "$SRC/rack3.ppfd" > "$WORK/rack3open.ppfd"
run rack3open.ppfd
chkabs "tiers: no shelf closure" "$(logval 'closure error')" 2e-4
# specbounce の上限は変換の数で決まる : 鏡面 4 面 (2 軸) なら 14 段まで入り、
# 打ち切りの収支誤差が 6 段の -1.4% から下がる。上限を超えた要求はログに出して下げる
sed -e '/^canopy/d' -e 's/^specbounce .*/specbounce  = 20/' "$SRC/difftrans.ppfd" > "$WORK/specdeep.ppfd"
run specdeep.ppfd
if grep -q "specbounce reduced 20 -> 14" "$WORK/ppfd.log"; then
	echo "specbounce: reduction logged -> OK"
else
	echo "specbounce: reduction logged -> NG" >&2
	status=1
fi
chkabs "specbounce: deep closure" "$(logval 'closure error')" 3e-3

echo "== (x) specular floor/ceiling + leafscatter (folded deposits, mirror boundary) =="
# mirror9.ppfd : 完全鏡面の床に乗った群落 = z で鏡映して 2 倍にした幾何
# (群落は鏡映面をはさんで 1 枚、光源の鏡像を実光源として置く)
run mirror9.ppfd
chkabs "zscat: closure" "$(logval 'closure error')" 1e-6
m9c=$(logval 'absorbed by canopy')
for t in above mid; do cp "$WORK/ppfd_map_$t.csv" "$WORK/map_m9$t.csv"; done
run mirror9ref.ppfd
chkabs "zscat: equiv. above canopy" "$(maxdiff map_m9above.csv ppfd_map_above.csv)" 1e-6
chkabs "zscat: equiv. in canopy"    "$(maxdiff map_m9mid.csv ppfd_map_mid.csv)" 1e-6
chk "zscat: canopy absorption" "$m9c" "$(awk -v a="$(logval 'absorbed by canopy')" 'BEGIN {printf "%.9g", a / 2}')" 1e-5
# 床と天井の両方が部分的な鏡面 (rho_d 0.05 + rho_s 0.9) で、群落と鏡面の間に
# すき間がある : 鏡像の連鎖は無限に続くので打ち切りが収支に出るが、段数を
# 上げると 0 へ収束する (4 段 -6.8e-3、16 段 -1.2e-4、32 段 -2.8e-6)
sed -e 's/^material = mirror gray 0.0 specular 1.0/material = mirror gray 0.05 specular 0.9/' \
    -e 's/^wall = all  black/wall = all  white/' -e 's/^wall = xmax mirror/wall = xmax white/' \
    -e 's/^wall = ymax mirror/wall = ymax white\nwall = zmax mirror\nwall = zmin mirror/' \
    -e 's/^specbounce .*/specbounce = 16/' "$SRC/mirror8.ppfd" > "$WORK/mirror8z.ppfd"
run mirror8z.ppfd
chkabs "zscat: two z-mirrors, 16 orders" "$(logval 'closure error')" 2e-4

echo "== (y) facing canopy slabs without a shelf (slab-to-slab diffuse exchange) =="
# 1 枚の群落を接する 2 枚に割っても結果が変わらない : 接した面どうしは
# 拡散光を全部やりとりし (cx = 1)、二流の adding と同じ解に収束する。
# rack.ppfd (白壁) の群落を 16 層 1 枚 / 8 層 x 2 枚 (層の厚さは同じ) で比べる
sed -e 's/^photoperiod/leafscatter = on\nleaflayers = 16\nphotoperiod/' \
    -e 's/^target = canopybase 0.06 25 13/target = canopybase 0.06 25 13\ntarget = upper 0.2 25 13/' \
	"$SRC/rack.ppfd" > "$WORK/rack_one.ppfd"
run rack_one.ppfd
for t in canopytop canopybase upper; do cp "$WORK/ppfd_map_$t.csv" "$WORK/map_one$t.csv"; done
r1c=$(logval 'absorbed by canopy')
sed -e 's/^canopy = 0.25 0.05 3.0 0.5 leaf/canopy = 0.25 0.15 1.5 0.5 leaf\ncanopy = 0.15 0.05 1.5 0.5 leaf/' \
    -e 's/^leaflayers = 16/leaflayers = 8/' "$WORK/rack_one.ppfd" > "$WORK/rack_two.ppfd"
run rack_two.ppfd
chkabs "facing: split == single (top)"   "$(maxdiff map_onecanopytop.csv ppfd_map_canopytop.csv)" 1e-6
chkabs "facing: split == single (upper)" "$(maxdiff map_oneupper.csv ppfd_map_upper.csv)" 1e-6
chkabs "facing: split == single (base)"  "$(maxdiff map_onecanopybase.csv ppfd_map_canopybase.csv)" 1e-6
chk "facing: split canopy abs" "$(logval 'absorbed by canopy')" "$r1c" 1e-6
# mirror10.ppfd : 鏡面の天井と群落の間にすき間 = z で 2 倍にした幾何
# (群落と鏡像の 2 枚が棚板なしで向かい合う)
run mirror10.ppfd
chkabs "facing: ceiling-gap closure" "$(logval 'closure error')" 1e-6
m10c=$(logval 'absorbed by canopy')
for t in above mid below; do cp "$WORK/ppfd_map_$t.csv" "$WORK/map_m10$t.csv"; done
run mirror10ref.ppfd
chkabs "facing: equiv. above canopy" "$(maxdiff map_m10above.csv ppfd_map_above.csv)" 1e-6
chkabs "facing: equiv. in canopy"    "$(maxdiff map_m10mid.csv ppfd_map_mid.csv)" 1e-6
chkabs "facing: equiv. below canopy" "$(maxdiff map_m10below.csv ppfd_map_below.csv)" 1e-6
chk "facing: canopy absorption" "$m10c" "$(awk -v a="$(logval 'absorbed by canopy')" 'BEGIN {printf "%.9g", a / 2}')" 1e-5

echo "== (z) cylinders and arcs (inscribed polygon of flat panels) =="
# cylinder.ppfd : 壁から壁まで渡した黒い 16 面の管の真上の点光源。
# 軸からの距離 0.50 R の光線は必ず遮られ、1.47 R の光線は必ず通る
run cylinder.ppfd
chkabs "cylinder: shadowed (d = 0.50 R)" "$(sumrow c 5)" 0
chk "cylinder: clear (d = 1.47 R)" "$(sumrow o 5)" 0.5477301 1e-6
if grep -q "cylinder     : rod .* 16 flat panels (inscribed polygon : sagitta" "$WORK/ppfd.log"; then
	echo "cylinder: sagitta logged -> OK"
else
	echo "cylinder: sagitta logged -> NG" >&2
	status=1
fi
# 4 面・45° 始まりの筒 = 同じ断面の occluder の箱 (白い壁・白い管。壁パッチの
# 境界が箱の辺に乗るよう patchdiv 10、面の分割もそろえる)。箱は slab 法、
# 板は平行四辺形の交差で遮蔽を判定するので、辺をかすめる光線の分だけ違う
sed -e 's/^patchdiv   = 6 6/patchdiv   = 10 10/' -e 's/^wall     = all black/wall     = all white/' \
    -e 's/^cylinder = .*/cylinder = rod  0.5 0.0 0.5  0.5 1.0 0.5  0.1414213562373095 4 white arc 45 405 div 10 10/' \
	"$SRC/cylinder.ppfd" > "$WORK/cylprism.ppfd"
run cylprism.ppfd
chkabs "cylinder: prism closure" "$(logval 'closure error')" 2e-3
cp "$WORK/ppfd_map_floor.csv" "$WORK/map_cylprism.csv"
sed -e 's/^patchdiv   = 6 6/patchdiv   = 10 10/' -e 's/^wall     = all black/wall     = all white/' \
    -e 's/^cylinder = .*/occluder = box 0.4 0.6 0.0 1.0 0.4 0.6 white div 10 10/' \
	"$SRC/cylinder.ppfd" > "$WORK/cylbox.ppfd"
run cylbox.ppfd
chkabs "cylinder: 4-gon == box" "$(maxdiff map_cylprism.csv ppfd_map_floor.csv)" 1e-3
# 半割りの反射フード (光源の上の円弧 180°、白) も板として収支に入る。
# 開口から中を見るパッチ対の可視率の標本化が残差を決める (-1.7e-3)
sed -e 's/^wall     = all black/wall     = all white/' \
    -e 's/^led = .*/led = 0.5 0.5 0.62 1.0 na555 iso/' \
    -e 's/^cylinder = .*/cylinder = hood 0.5 0.1 0.6  0.5 0.9 0.6  0.15 8 white arc -90 90/' \
	"$SRC/cylinder.ppfd" > "$WORK/cylhood.ppfd"
run cylhood.ppfd
chkabs "cylinder: half-pipe hood closure" "$(logval 'closure error')" 3e-3
# 退化した管 (半径 0 / 長さ 0) とチャンバからはみ出す管は入力段で弾く
for v in r0 out; do
	case $v in
	r0)  sed 's/^cylinder = .*/cylinder = rod  0.5 0.0 0.5  0.5 1.0 0.5  0.0 16 black/' "$SRC/cylinder.ppfd" > "$WORK/cyl_$v.ppfd" ;;
	out) sed 's/^cylinder = .*/cylinder = rod  0.5 0.0 0.5  0.5 1.0 0.5  0.6 16 black/' "$SRC/cylinder.ppfd" > "$WORK/cyl_$v.ppfd" ;;
	esac
	if (cd "$WORK" && "$OPPFD" "cyl_$v.ppfd" > /dev/null 2>&1); then
		echo "cylinder: invalid ($v) rejected -> NG (accepted)" >&2
		status=1
	else
		echo "cylinder: invalid ($v) rejected -> OK"
	fi
done

echo "== (aa) canopy beds with leafscatter (one two-stream column per bed) =="
# bedscat.ppfd : 断面全体の 1 枚 = x = 0.5 で割った対称な 2 つのベッド
run bedscat.ppfd
bsc=$(logval 'absorbed by canopy'); bss=$(logval 'canopy scattered')
for t in above mid below; do cp "$WORK/ppfd_map_$t.csv" "$WORK/map_bs$t.csv"; done
sed 's/^canopy = 0.4 0.2 2.0 0.5 leaf/canopy = 0.4 0.2 2.0 0.5 leaf  0 0.5 0 1\ncanopy = 0.4 0.2 2.0 0.5 leaf  0.5 1 0 1/' \
	"$SRC/bedscat.ppfd" > "$WORK/bedscat2.ppfd"
run bedscat2.ppfd
chkabs "bedscat: 2 beds == slab (above)" "$(maxdiff map_bsabove.csv ppfd_map_above.csv)" 1e-6
chkabs "bedscat: 2 beds == slab (in)"    "$(maxdiff map_bsmid.csv ppfd_map_mid.csv)" 1e-6
chkabs "bedscat: 2 beds == slab (below)" "$(maxdiff map_bsbelow.csv ppfd_map_below.csv)" 1e-6
chk "bedscat: canopy absorption" "$(logval 'absorbed by canopy')" "$bsc" 1e-6
chk "bedscat: scattered returned" "$(logval 'canopy scattered')" "$bss" 1e-6
# 通路のある 1 つのベッド (非対称) でも収支が閉じる
sed 's/^canopy = 0.4 0.2 2.0 0.5 leaf/canopy = 0.4 0.2 2.0 0.5 leaf  0.125 0.625 0.25 0.75/' \
	"$SRC/bedscat.ppfd" > "$WORK/bedscat1.ppfd"
run bedscat1.ppfd
chkabs "bedscat: one bed closure" "$(logval 'closure error')" 1e-6
# canopybed.ppfd (黒い壁のベッド) + leafscatter (以前は拒否)
sed 's/^photoperiod/leafscatter = on\nphotoperiod/' "$SRC/canopybed.ppfd" > "$WORK/canopybedls.ppfd"
run canopybedls.ppfd
chkabs "bedscat: canopybed closure" "$(logval 'closure error')" 1e-6
# mirror11.ppfd : 鏡面の側壁の隣のベッド = x 方向に 2 倍 (ベッドが 2 つ並ぶ)
run mirror11.ppfd
chkabs "bedscat: mirror closure" "$(logval 'closure error')" 1e-6
m11c=$(logval 'absorbed by canopy')
for t in above mid below; do cp "$WORK/ppfd_map_$t.csv" "$WORK/map_m11$t.csv"; done
run mirror11ref.ppfd
chkabs "bedscat: mirror equiv. (above)" "$(maxdiff map_m11above.csv ppfd_map_above.csv)" 1e-6
chkabs "bedscat: mirror equiv. (in)"    "$(maxdiff map_m11mid.csv ppfd_map_mid.csv)" 1e-6
chkabs "bedscat: mirror equiv. (below)" "$(maxdiff map_m11below.csv ppfd_map_below.csv)" 1e-6
chk "bedscat: mirror canopy abs" "$m11c" "$(awk -v a="$(logval 'absorbed by canopy')" 'BEGIN {printf "%.9g", a / 2}')" 1e-5
# 弾く組み合わせ : 床・天井の鏡面 / 高さの違うベッドの並び /
# 棚板なしで別の段にあるベッド (床に乗ったベッドは検証 (ad) で受け付ける)
for v in zmir height noshelf; do
	case $v in
	zmir)    sed 's/^wall = zmax white/wall = zmax mirror/' "$SRC/mirror11.ppfd" > "$WORK/bedls_$v.ppfd" ;;
	height)  sed 's/^canopy = 0.4 0.2 2.0 0.5 leaf  1.1 1.5 0.2 0.8/canopy = 0.45 0.25 2.0 0.5 leaf  1.1 1.5 0.2 0.8/' \
	             "$SRC/mirror11ref.ppfd" > "$WORK/bedls_$v.ppfd" ;;
	noshelf) sed 's/^canopy = 0.4 0.2 2.0 0.5 leaf  1.1 1.5 0.2 0.8/canopy = 0.9 0.7 2.0 0.5 leaf  1.1 1.5 0.2 0.8/' \
	             "$SRC/mirror11ref.ppfd" > "$WORK/bedls_$v.ppfd" ;;
	esac
	if (cd "$WORK" && "$OPPFD" "bedls_$v.ppfd" > /dev/null 2>&1); then
		echo "bedscat: invalid ($v) rejected -> NG (accepted)" >&2
		status=1
	else
		echo "bedscat: invalid ($v) rejected -> OK"
	fi
done

echo "== (ab) monte carlo photon tracing (general reflection, 3-D canopy scattering) =="
# mc.ppfd : 黒いチャンバの点光源。測定面全体の平均 = Φ Ω / (4π A) (統計誤差 0.097%)
run mc.ppfd
chk "mc: point source (solid angle)" "$(sumval 14)" 0.2087326 4e-3
chkabs "mc: black chamber closure" "$(logval 'closure error')" 1e-12
# スレッド数によらずビット一致 (光子の数は整数)
cp "$WORK/ppfd_map_floor.csv" "$WORK/map_mc2.csv"
(cd "$WORK" && "$OPPFD" -n 1 mc.ppfd > /dev/null)
chkabs "mc: 1 thread == 2 threads" "$(maxdiff map_mc2.csv ppfd_map_floor.csv)" 0
# 平行平板の拡散透過 2 E3(tau0) (difftrans.ppfd)。鏡像の打ち切りが無いので
# 決定論の解法より解析解に近い (統計誤差 : tau0 = 1 で 0.094%、3 で 0.37%)
sed 's/^photoperiod/montecarlo = 4000000 11\nphotoperiod/' "$SRC/difftrans.ppfd" > "$WORK/mcdt1.ppfd"
run mcdt1.ppfd
chk "mc: 2E3(1) floor"   "$(logval 'absorbed by walls')" 0.2193839 4e-3
sed -e 's/^photoperiod/montecarlo = 4000000 11\nphotoperiod/' -e 's/^canopy = 0.25 0.05 2.0/canopy = 0.25 0.05 1.0/' \
	"$SRC/difftrans.ppfd" > "$WORK/mcdt05.ppfd"
run mcdt05.ppfd
chk "mc: 2E3(0.5) floor" "$(logval 'absorbed by walls')" 0.4432087 3e-3
sed -e 's/^photoperiod/montecarlo = 4000000 11\nphotoperiod/' -e 's/^canopy = 0.25 0.05 2.0/canopy = 0.25 0.05 6.0/' \
	"$SRC/difftrans.ppfd" > "$WORK/mcdt3.ppfd"
run mcdt3.ppfd
chk "mc: 2E3(3) floor"   "$(logval 'absorbed by walls')" 0.0178613 1.5e-2
# 無吸収葉 (ω = 1) は 1 光子も吸収しない (構造的に厳密)。白壁なら収支も厳密
sed -e 's/^material = blackleaf gray 0.0 0.0/material = blackleaf gray 0.5 0.5\nmaterial = white gray 0.7/' \
    -e 's/^wall     = all black/wall = all white/' -e 's/^photoperiod/montecarlo = 1000000 5\nphotoperiod/' \
	"$SRC/canopy.ppfd" > "$WORK/mcw1.ppfd"
run mcw1.ppfd
chkabs "mc: omega=1 no capture" "$(logval 'absorbed by canopy')" 0
chkabs "mc: white cavity closure" "$(logval 'closure error')" 1e-12
# 白いキャビティ (cavity.ppfd) : 光子追跡 = ラジオシティ (ラジオシティの離散化の分)
run cavity.ppfd
cavd=$(sumval 14)
sed 's/^photoperiod/montecarlo = 1000000 5\nphotoperiod/' "$SRC/cavity.ppfd" > "$WORK/mccav.ppfd"
run mccav.ppfd
chk "mc: == radiosity (white cavity)" "$(sumval 14)" "$cavd" 1e-2
# 光沢反射 : ローブの指数を大きくすると鏡面に一致する (mirror5.ppfd の鏡面 2 面)
sed 's/^specbounce .*/montecarlo = 4000000 3/' "$SRC/mirror5.ppfd" > "$WORK/mcm5.ppfd"
run mcm5.ppfd
mcs=$(sumval 5)
sed 's/^material = mirror gray 0.0 specular 1.0/material = mirror gray 0.0 glossy 1.0 1e7/' "$WORK/mcm5.ppfd" > "$WORK/mcm5g.ppfd"
run mcm5g.ppfd
chk "mc: glossy n->inf == mirror" "$(sumval 5)" "$mcs" 3e-3
chkabs "mc: glossy closure" "$(logval 'closure error')" 1e-12
# 光沢の幅を広げると鏡面の像がぼける (同じ光子の数で、鏡像の向きの集中が下がる) :
# 収支は厳密のまま
sed 's/^material = mirror gray 0.0 specular 1.0/material = mirror gray 0.0 glossy 0.9 20/' "$WORK/mcm5.ppfd" > "$WORK/mcm5g20.ppfd"
run mcm5g20.ppfd
chkabs "mc: glossy n=20 closure" "$(logval 'closure error')" 1e-12
# 決定論の解法では光沢を受け付けない (ラジオシティは放射輝度が向きによらない前提)
sed 's/^material = mirror gray 0.0 specular 1.0/material = mirror gray 0.0 glossy 1.0 100/' "$SRC/mirror5.ppfd" > "$WORK/detglossy.ppfd"
if (cd "$WORK" && "$OPPFD" detglossy.ppfd > /dev/null 2>&1); then
	echo "mc: glossy w/o montecarlo rejected -> NG (accepted)" >&2
	status=1
else
	echo "mc: glossy w/o montecarlo rejected -> OK"
fi
# 棚板の鏡面 (決定論では弾く) も光子追跡なら扱える : 収支は厳密
sed -e 's/^occluder = plate .*/occluder = plate 0.7 0.9 0.3 0.7 0.57 0.57 mfilm div 4 4/' \
    -e 's/^material = mirror gray 0.0 specular 1.0/material = mirror gray 0.0 specular 1.0\nmaterial = mfilm gray 0.05 specular 0.9/' \
    -e 's/^specbounce .*/montecarlo = 1000000 3/' "$SRC/mirror6.ppfd" > "$WORK/mcocc.ppfd"
run mcocc.ppfd
chkabs "mc: specular shelf closure" "$(logval 'closure error')" 1e-12
# 群落の横方向の散乱 (mclateral.ppfd) : 通路の鉛直面はベッドの側面からの散乱光を
# 受ける (黒葉より明るい)。決定論の二流 (鉛直 1 次元) は横に出さない (黒葉と同じ)
run mclateral.ppfd
lat1=$(sumrow aisle 5)
sed 's/^material = leaf  gray 0.25 0.25/material = leaf  gray 0.0 0.0/' "$SRC/mclateral.ppfd" > "$WORK/mclat0.ppfd"
run mclat0.ppfd
lat0=$(sumrow aisle 5)
awk -v a="$lat1" -v b="$lat0" 'BEGIN {
	r = a / b;
	printf "%-28s scatter/black = %.4f -> %s (must be > 1.1)\n", "mc: lateral scattering to aisle", r, (r > 1.1) ? "OK" : "NG";
	exit (r > 1.1) ? 0 : 1}' || status=1
sed 's/^montecarlo.*//' "$WORK/mclat0.ppfd" > "$WORK/detlat0.ppfd"
run detlat0.ppfd
chk "mc: black leaf == deterministic" "$lat0" "$(sumrow aisle 5)" 1.5e-2

echo "== (ac) IES photometric types A / B =="
# phottype.ppfd : +x 向きの器具、「|V| <= 30° で一定」の帯。正規化の分母
# (帯の立体角) が型で違い、型 A の (H, V) = (60°, 25°) の向きは型 B では帯の外
run phottype.ppfd
chk "phottype: A on axis"      "$(sumrow axis 14)" 3.0870075 1e-5
chk "phottype: A (H,V)=(60,25)" "$(sumrow off 14)" 3.0870075 1e-5
sed 's/bandA\.ies/bandB.ies/' "$SRC/phottype.ppfd" > "$WORK/phottypeb.ppfd"
run phottypeb.ppfd
chk "phottype: B on axis"      "$(sumrow axis 14)" 4.5472841 1e-5
chkabs "phottype: B same dir = 0" "$(sumrow off 14)" 0
# 正規化は双線形の表の厳密な積分 : 白いキャビティでも壁の吸収 = 放射束
# (残差は帯の縁が壁パッチをまたぐ直接光の求積 : quadrature 4 で -2.7e-4、16 で -3.0e-5)
sed -e 's/^material = black gray 0.0/material = black gray 0.0\nmaterial = white gray 0.8/' \
    -e 's/^wall     = all black/wall     = all white/' -e 's/^photoperiod/quadrature = 16\nphotoperiod/' \
	"$WORK/phottypeb.ppfd" > "$WORK/phottypebw.ppfd"
run phottypebw.ppfd
chkabs "phottype: B white closure" "$(logval 'closure error')" 1e-4
# 光子追跡でも同じ配光 (棄却法で向きを引く) : 帯を受ける面 x = 0.95 の平均が決定論と一致
sed 's/^target = floor 0.0 4 4/sensor = front 0.95 1 0  0.95 0 0  0.95 1 1  40 40/' "$WORK/phottypeb.ppfd" > "$WORK/phottypebf.ppfd"
run phottypebf.ppfd
ptd=$(sumrow front 14)
sed 's/^photoperiod/montecarlo = 4000000 3\nphotoperiod/' "$WORK/phottypebf.ppfd" > "$WORK/phottypebmc.ppfd"
run phottypebmc.ppfd
chk "phottype: B monte carlo == det" "$(sumrow front 14)" "$ptd" 5e-3

echo "== (ad) canopies resting on a surface (floor / shelf) with leafscatter =="
# rest.ppfd : 床に乗ったベッド (下段) と棚板に乗った断面全体の群落 (上段)。
# 面に乗った面の散乱光はその下のパッチが重なりの割合で受け、面はほかの
# どこからも見えない。棚板で仕切った 2 段 = 段ごとの単独計算 (上段の棚板 =
# 単独計算の床)。以前は棚板に乗った群落の下面の光が棚板を素通りして下の段へ
# 漏れ (rack3 の上段を棚板に乗せると下段の群落の上 +2.4%)、床に乗ったベッドは弾いていた
run rest.ppfd
rsc=$(logval 'absorbed by canopy')
chkabs "rest: closure" "$(logval 'closure error')" 1e-4
for t in a_top a_in b_top b_in; do cp "$WORK/ppfd_map_$t.csv" "$WORK/map_rs$t.csv"; done
# 面に乗った面の分配は「重なりの割合」の総和 = 面積 (ベッドの縁は床のパッチの格子にそろえていない)
chk "rest: bed on floor plane sum" "$(awk '/^  slab 1 / {print $NF; exit}' "$WORK/ppfd.log")" 1 1e-4
sed -e 's/^canopy = .*/canopy = 0.20 0.0 2.0 0.5 leaf  0.2 0.7 0.15 0.85/' -e 's/^target = a_in   0.15/target = a_in   0.10/' \
	"$SRC/rack3a.ppfd" > "$WORK/restlo.ppfd"
run restlo.ppfd
rlc=$(logval 'absorbed by canopy')
chkabs "rest: lower == alone (top)" "$(maxdiff map_rsa_top.csv ppfd_map_a_top.csv)" 1e-6
chkabs "rest: lower == alone (in)"  "$(maxdiff map_rsa_in.csv ppfd_map_a_in.csv)" 1e-6
sed -e 's/^wall = zmin tray/wall = zmin white/' -e 's/^led = .*/led = 0.6 0.4 0.45 0.8 white4000 beam 1/' \
    -e 's/^canopy = .*/canopy = 0.20 0.0 1.5 0.5 leaf/' \
    -e 's/^target = a_top  0.30/target = b_top  0.30/' -e 's/^target = a_in   0.15/target = b_in   0.10/' \
	"$SRC/rack3a.ppfd" > "$WORK/restup.ppfd"
run restup.ppfd
ruc=$(logval 'absorbed by canopy')
chkabs "rest: upper == alone (top)" "$(maxdiff map_rsb_top.csv ppfd_map_b_top.csv)" 1e-6
chkabs "rest: upper == alone (in)"  "$(maxdiff map_rsb_in.csv ppfd_map_b_in.csv)" 1e-6
chk "rest: canopy abs = sum" "$rsc" "$(awk -v a="$rlc" -v b="$ruc" 'BEGIN {printf "%.9g", a + b}')" 1e-5
# 連続性 : 上段の群落を 1e-6 m 浮かせる (面はどこにも乗らず、形態係数で配る)
sed 's/^canopy = 0.70 0.50/canopy = 0.70 0.500001/' "$SRC/rest.ppfd" > "$WORK/restlift.ppfd"
run restlift.ppfd
chkabs "rest: lifted 1e-6 (top)" "$(maxdiff map_rsb_top.csv ppfd_map_b_top.csv)" 1e-5
chkabs "rest: lifted 1e-6 (in)"  "$(maxdiff map_rsb_in.csv ppfd_map_b_in.csv)" 1e-5
chkabs "rest: lifted 1e-6 (lower)" "$(maxdiff map_rsa_in.csv ppfd_map_a_in.csv)" 1e-5
# 床に乗ったベッドどうし : x = 0.5 で割った 2 つのベッド = 床に乗った断面全体の 1 枚
sed 's/^canopy = 0.4 0.2 2.0 0.5 leaf/canopy = 0.2 0.0 2.0 0.5 leaf/' "$SRC/bedscat.ppfd" > "$WORK/restf1.ppfd"
run restf1.ppfd
for t in above mid below; do cp "$WORK/ppfd_map_$t.csv" "$WORK/map_rf$t.csv"; done
sed 's/^canopy = 0.4 0.2 2.0 0.5 leaf/canopy = 0.2 0.0 2.0 0.5 leaf  0 0.5 0 1\ncanopy = 0.2 0.0 2.0 0.5 leaf  0.5 1 0 1/' \
	"$SRC/bedscat.ppfd" > "$WORK/restf2.ppfd"
run restf2.ppfd
chkabs "rest: 2 beds on floor == slab (above)" "$(maxdiff map_rfabove.csv ppfd_map_above.csv)" 1e-6
chkabs "rest: 2 beds on floor == slab (in)"    "$(maxdiff map_rfbelow.csv ppfd_map_below.csv)" 1e-6
# 鏡面の側壁の隣の、床に乗ったベッド = x 方向に 2 倍 (ベッドの縁は床の格子にそろえていない)
sed 's/^canopy = 0.4 0.2 2.0/canopy = 0.4 0.0 2.0/' "$SRC/mirror11.ppfd" > "$WORK/restm.ppfd"
run restm.ppfd
chkabs "rest: mirror bed on floor closure" "$(logval 'closure error')" 1e-6
rmc=$(logval 'absorbed by canopy')
for t in above mid below; do cp "$WORK/ppfd_map_$t.csv" "$WORK/map_rm$t.csv"; done
sed 's/^canopy = 0.4 0.2 2.0/canopy = 0.4 0.0 2.0/' "$SRC/mirror11ref.ppfd" > "$WORK/restmr.ppfd"
run restmr.ppfd
chkabs "rest: mirror equiv. (above)" "$(maxdiff map_rmabove.csv ppfd_map_above.csv)" 1e-6
chkabs "rest: mirror equiv. (in)"    "$(maxdiff map_rmmid.csv ppfd_map_mid.csv)" 1e-6
chk "rest: mirror canopy abs" "$rmc" "$(awk -v a="$(logval 'absorbed by canopy')" 'BEGIN {printf "%.9g", a / 2}')" 1e-5
# 弾く組み合わせ : 面の一部だけが乗る (覆われない部分だけを配る経路を持たない)
for v in shelf tray; do
	case $v in
	shelf) sed 's/^occluder = shelf 0 1 0 1/occluder = shelf 0 1 0 0.9/' "$SRC/rest.ppfd" > "$WORK/restbad_$v.ppfd" ;;
	tray)  sed -e 's/^occluder = shelf.*/&\noccluder = tray 0.3 0.6 0.3 0.6 0.0 0.1 tray/' \
	           -e 's/^canopy = 0.20 0.0 2.0/canopy = 0.30 0.1 2.0/' "$SRC/rest.ppfd" > "$WORK/restbad_$v.ppfd" ;;
	esac
	if (cd "$WORK" && "$OPPFD" "restbad_$v.ppfd" > /dev/null 2>&1); then
		echo "rest: partly resting ($v) rejected -> NG (accepted)" >&2
		status=1
	else
		echo "rest: partly resting ($v) rejected -> OK"
	fi
done

echo "== (ae) occluders in contact (grid aligned to contact lines, covered patches excluded) =="
# contact.ppfd : 側壁の格子線に乗らない z = 0.4 の棚板、光源は下だけ。
# 側壁の格子線を棚板の高さへ動かすので、棚板の上は厳密に 0 (以前は 0.077 漏れていた)
run contact.ppfd
chkabs "contact: no leak above shelf" "$(sumrow above 5)" 0
chkabs "contact: shelf closure" "$(logval 'closure error')" 1e-6
if grep -q 'aligned to contact lines : 4 lines moved' "$WORK/ppfd.log"; then
	echo "contact: grid alignment logged -> OK"
else
	echo "contact: grid alignment logged -> NG" >&2
	status=1
fi
# 斜めの板 (panel) の仕切りでも同じ : 板と壁の交線へ格子をそろえる。x = 0.4 の鉛直な
# 仕切り (光源は x < 0.4 だけ) の向こうは厳密に 0、同じ板を occluder で書いたものとビット一致
# (以前は panel だけ格子をそろえず 0.119 漏れていた)
for v in occ pan; do
	case $v in
	occ) o='occluder = shelf 0.4 0.4 0 1 0 1 white' ;;
	pan) o='panel = shelf 0.4 0 0  0.4 1 0  0.4 0 1 white' ;;
	esac
	sed -e "s/^occluder = .*/$o/" -e 's/^led = .*/led = 0.2 0.5 0.5 1.0 g iso/' \
	    -e 's/^target = above .*/target = left  0.3 4 8 0.0 0.4 0 1/' -e 's/^target = below .*/target = right 0.3 4 8 0.4 1.0 0 1/' \
		"$SRC/contact.ppfd" > "$WORK/contact_div$v.ppfd"
	run "contact_div$v.ppfd"
	cp "$WORK/ppfd_map_left.csv" "$WORK/map_cdiv$v.csv"
	chkabs "contact: $v divider no leak" "$(sumrow right 5)" 0
done
chkabs "contact: panel divider == occluder" "$(maxdiff map_cdivocc.csv map_cdivpan.csv)" 0
# 格子の軸に平行でない接触線 (斜めの仕切りと床・天井の交線) はそろえられないのでログで警告する
sed -e 's/^occluder = .*/panel = d 0.2 0 0  0.8 1 0  0.2 0 1 white/' "$SRC/contact.ppfd" > "$WORK/contact_diag.ppfd"
run contact_diag.ppfd
if grep -q 'warning : 2 panel contact lines are not parallel' "$WORK/ppfd.log"; then
	echo "contact: diagonal contact warned -> OK"
else
	echo "contact: diagonal contact warned -> NG" >&2
	status=1
fi
# 面どうしの接触 : 光源を上へ移し、棚板を接触する遮蔽物に置き換えた変種の収支
for v in block onshelf mat poster box; do
	case $v in
	block)   o='occluder = block 0 1 0 1 0.0 0.2 white' ; tol=1e-6 ;;
	onshelf) o='occluder = shelf 0 1 0 1 0.5 0.5 white\noccluder = block 0 1 0 1 0.5 0.7 white' ; tol=1e-6 ;;
	mat)     o='occluder = mat 0.2 0.7 0.15 0.85 0.0 0.0 white' ; tol=1e-6 ;;
	poster)  o='occluder = poster 0.0 0.0 0.2 0.7 0.3 0.8 black' ; tol=1e-6 ;;
	box)     o='occluder = tray 0.2 0.7 0.15 0.85 0.0 0.125 white' ; tol=3e-3 ;;
	esac
	sed -e "s/^occluder = .*/$o/" -e 's/^led = .*/led = 0.5 0.5 0.9 1.0 g iso/' "$SRC/contact.ppfd" > "$WORK/contact_$v.ppfd"
	run "contact_$v.ppfd"
	chkabs "contact: $v closure" "$(logval 'closure error')" "$tol"
done
# 床いっぱいの黒いマット = 黒い床 (覆われた床のパッチは外れ、マットの上面が床になる)
sed -e 's/^occluder = .*/occluder = mat 0 1 0 1 0.0 0.0 black/' -e 's/^led = .*/led = 0.5 0.5 0.9 1.0 g iso/' \
	"$SRC/contact.ppfd" > "$WORK/contact_fullmat.ppfd"
run contact_fullmat.ppfd
cp "$WORK/ppfd_map_below.csv" "$WORK/map_cmat.csv"
sed -e '/^occluder = /d' -e 's/^led = .*/led = 0.5 0.5 0.9 1.0 g iso/' -e 's/^wall = all  white/&\nwall = zmin black/' \
	"$SRC/contact.ppfd" > "$WORK/contact_blackfloor.ppfd"
run contact_blackfloor.ppfd
chkabs "contact: black mat == black floor" "$(maxdiff map_cmat.csv ppfd_map_below.csv)" 1e-6
# 光子追跡でも同じ : 同じ t に床とマットの上面があるとき、覆われた床ではなくマットに当たる
# (以前は床に当たり、白い床として反射していた : 床の近くの測定面 +26%)
sed 's/^solver/montecarlo = 200000 1\nsolver/' "$WORK/contact_fullmat.ppfd" > "$WORK/contact_fullmatmc.ppfd"
run contact_fullmatmc.ppfd
cp "$WORK/ppfd_map_below.csv" "$WORK/map_cmatmc.csv"
sed 's/^solver/montecarlo = 200000 1\nsolver/' "$WORK/contact_blackfloor.ppfd" > "$WORK/contact_blackfloormc.ppfd"
run contact_blackfloormc.ppfd
chkabs "contact: mc black mat == black floor" "$(maxdiff map_cmatmc.csv ppfd_map_below.csv)" 0
# 弾く置き方 : 同じ平面で重なる板 / 体積の重なる箱 / 床に寝かせた斜めの板
for v in plates boxes panel; do
	case $v in
	plates) o='occluder = a 0.2 0.7 0.15 0.85 0.3 0.3 white\noccluder = b 0.5 0.9 0.15 0.85 0.3 0.3 white' ;;
	boxes)  o='occluder = a 0.2 0.7 0.15 0.85 0.3 0.4 white\noccluder = b 0.5 0.9 0.15 0.85 0.3 0.4 white' ;;
	panel)  o='panel = p 0.2 0.15 0.0  0.7 0.15 0.0  0.2 0.85 0.0 white' ;;
	esac
	sed -e "s/^occluder = .*/$o/" "$SRC/contact.ppfd" > "$WORK/contactbad_$v.ppfd"
	if (cd "$WORK" && "$OPPFD" "contactbad_$v.ppfd" > /dev/null 2>&1); then
		echo "contact: overlap ($v) rejected -> NG (accepted)" >&2
		status=1
	else
		echo "contact: overlap ($v) rejected -> OK"
	fi
done

echo
if [ "$status" -ne 0 ]; then
	echo "*** PPFD validation FAILED" >&2
else
	echo "PPFD validation passed"
fi
exit $status
