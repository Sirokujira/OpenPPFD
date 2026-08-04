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
if [ "$status" -ne 0 ]; then
	echo "*** PPFD validation FAILED" >&2
else
	echo "PPFD validation passed"
fi
exit $status
