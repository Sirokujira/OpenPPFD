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

# sumval <列番号> : ppfd_summary.csv の 1 行目 (最初の target)
sumval() {
	awk -F, -v c="$1" 'NR == 2 {print $c; exit}' "$WORK/ppfd_summary.csv"
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
chkabs "mirror: no truncation"  "$(awk '/^specular  / {print $9}' "$WORK/ppfd.log")" 0
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
chkabs "specdiff: no truncation"     "$(awk '/^specular  / {print $9}' "$WORK/ppfd.log")" 0
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
if [ "$status" -ne 0 ]; then
	echo "*** PPFD validation FAILED" >&2
else
	echo "PPFD validation passed"
fi
exit $status
