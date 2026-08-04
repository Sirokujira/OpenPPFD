/*
spectrum.c

分光グリッド、組み込み作用曲線 (McCree) と標準比視感度 V(λ)、
分光放射照度 [W/m²] から各測光量 (PPFD / YPFD / ePAR / lx) への換算。

■ 分光量の持ち方
すべての分光量は等間隔グリッド lam[i] = lam0 + i*dlam 上の
「ビン積分値」として持つ。すなわち E[i] は幅 dlam のビン i に含まれる
放射照度 [W/m²] であって [W/m²/nm] ではない。こうすると光量子換算が
単なる重み付き総和になり、ビンの端の扱いも band_weight() に集約できる。

■ 光量子束密度
    PPFD = Σ_i E[i] * k(λ_i) * w_i     k(λ) = λ * 8.359346e-3 [µmol/J]
w_i は PAR 帯 (400-700 nm) とビン [λ_i-dλ/2, λ_i+dλ/2] の重なり割合。
k(λ) が λ の 1 次式なので、単色光をグリッドの隣接 2 点へ線形配分しても
PPFD は厳密に一致する (input_data.c の mono スペクトル)。
*/

#include "ppfd.h"

/*
CIE 1924 標準比視感度 V(λ) (明所視)。5 nm 間隔 380-780 nm。
555 nm でちょうど 1.0 (lm の定義 683 lm/W の基準)。
*/
static const double vlam_lam0 = 380.0, vlam_dlam = 5.0;
static const double vlam_tbl[] = {
	0.000039, 0.000064, 0.000120, 0.000217, 0.000396, 0.000640, 0.001210, 0.002180,
	0.004000, 0.007300, 0.011600, 0.016840, 0.023000, 0.029800, 0.038000, 0.048000,
	0.060000, 0.073900, 0.090980, 0.112600, 0.139020, 0.169300, 0.208020, 0.258600,
	0.323000, 0.407300, 0.503000, 0.608200, 0.710000, 0.793200, 0.862000, 0.914850,
	0.954000, 0.980300, 0.994950, 1.000000, 0.995000, 0.978600, 0.952000, 0.915400,
	0.870000, 0.816300, 0.757000, 0.694900, 0.631000, 0.566800, 0.503000, 0.441200,
	0.381000, 0.321000, 0.265000, 0.217000, 0.175000, 0.138200, 0.107000, 0.081600,
	0.061000, 0.044580, 0.032000, 0.023200, 0.017000, 0.011920, 0.008210, 0.005723,
	0.004102, 0.002929, 0.002091, 0.001484, 0.001047, 0.000740, 0.000520, 0.000361,
	0.000249, 0.000172, 0.000120, 0.000085, 0.000060, 0.000042, 0.000030, 0.000021,
	0.000015
};
#define NVLAM ((int)(sizeof(vlam_tbl) / sizeof(vlam_tbl[0])))

/*
McCree (1972) の相対光合成量子収率 (22 作物種の平均) の近似値。
25 nm 間隔 350-750 nm、625 nm で 1.0 に正規化。

【重要】公表曲線をディジタイズした近似値であって公式データではない
(誤差の目安 ±0.03)。厳密な作用曲線が要る場合は
    actionspectrum = file <path>
で「λ[nm] 相対量子収率」のテキストを与えて上書きすること。
YPFD の検証は 625 nm 単色で YPFD == PPFD になること (曲線の値に依存
しない構造的な性質) で行う。
*/
static const double mccree_lam0 = 350.0, mccree_dlam = 25.0;
static const double mccree_tbl[] = {
	0.30, 0.55, 0.66, 0.70, 0.73, 0.73, 0.71, 0.75,
	0.80, 0.88, 0.96, 1.00, 0.98, 0.90, 0.48, 0.13,
	0.03
};
#define NMCCREE ((int)(sizeof(mccree_tbl) / sizeof(mccree_tbl[0])))

/* 等間隔テーブルの線形補間 (範囲外は 0) */
static double interp_uniform(const double *tbl, int n, double lam0, double dlam, double lam)
{
	double f;
	int    i;

	if ((lam < lam0) || (lam > lam0 + (n - 1) * dlam)) return 0.0;

	f = (lam - lam0) / dlam;
	i = (int)f;
	if (i >= n - 1) return tbl[n - 1];
	f -= i;
	return (tbl[i] * (1.0 - f)) + (tbl[i + 1] * f);
}

/* 非等間隔テーブル (昇順) の線形補間。範囲外は端の値でクランプ */
double interp_table(const double *x, const double *y, int n, double xq)
{
	int i;

	if (n <= 0) return 0.0;
	if (n == 1) return y[0];
	if (xq <= x[0]) return y[0];
	if (xq >= x[n - 1]) return y[n - 1];

	for (i = 0; i < n - 1; i++) {
		if ((xq >= x[i]) && (xq <= x[i + 1])) {
			const double d = x[i + 1] - x[i];
			const double f = (d > EPS) ? (xq - x[i]) / d : 0.0;
			return (y[i] * (1.0 - f)) + (y[i + 1] * f);
		}
	}
	return y[n - 1];
}

/* 分光グリッドの確保 (input_data で lam0/dlam/nlam が決まったあと) */
void spectrum_grid(ppfd_t *p)
{
	int i;

	p->lam = (double *)xmalloc((size_t)p->nlam * sizeof(double));
	for (i = 0; i < p->nlam; i++) {
		p->lam[i] = p->lam0 + (i * p->dlam);
	}
}

/* 組み込み曲線をグリッドへ再標本化 (actionspectrum = file で上書き可) */
void spectrum_builtin(ppfd_t *p)
{
	int i;

	if (p->vlambda == NULL) {
		p->vlambda = (double *)xmalloc((size_t)p->nlam * sizeof(double));
		for (i = 0; i < p->nlam; i++) {
			p->vlambda[i] = interp_uniform(vlam_tbl, NVLAM, vlam_lam0, vlam_dlam, p->lam[i]);
		}
	}
	if (p->action == NULL) {
		p->action = (double *)xmalloc((size_t)p->nlam * sizeof(double));
		for (i = 0; i < p->nlam; i++) {
			p->action[i] = interp_uniform(mccree_tbl, NMCCREE, mccree_lam0, mccree_dlam, p->lam[i]);
		}
	}
}

/* 総和が 1 になるよう正規化 */
void spec_normalize(double *w, int n)
{
	int    i;
	double s = 0.0;

	for (i = 0; i < n; i++) s += w[i];
	if (s > EPS) {
		for (i = 0; i < n; i++) w[i] /= s;
	}
}

/* ビン i と帯域 [lam1, lam2] の重なり割合 (0..1) */
double band_weight(const ppfd_t *p, int i, double lam1, double lam2)
{
	const double b1 = p->lam[i] - (0.5 * p->dlam);
	const double b2 = p->lam[i] + (0.5 * p->dlam);
	double lo = (b1 > lam1) ? b1 : lam1;
	double hi = (b2 < lam2) ? b2 : lam2;

	if (hi <= lo) return 0.0;
	return (hi - lo) / p->dlam;
}

/* 帯域 [lam1,lam2] の光量子束密度 [µmol/m²/s] */
double calc_bandppf(const ppfd_t *p, const double *E, double lam1, double lam2)
{
	int    i;
	double s = 0.0;

	for (i = 0; i < p->nlam; i++) {
		const double w = band_weight(p, i, lam1, lam2);
		if (w > 0.0) s += E[i] * PHOTON_K(p->lam[i]) * w;
	}
	return s;
}

/* PPFD [µmol/m²/s] : PAR = 400-700 nm */
double calc_ppfd(const ppfd_t *p, const double *E)
{
	return calc_bandppf(p, E, 400.0, 700.0);
}

/* ePAR (拡張 PAR) [µmol/m²/s] : 400-750 nm (遠赤色を含む) */
double calc_epar(const ppfd_t *p, const double *E)
{
	return calc_bandppf(p, E, 400.0, 750.0);
}

/* YPFD [µmol/m²/s] : McCree 作用曲線で重み付けした PAR 光量子束密度 */
double calc_ypfd(const ppfd_t *p, const double *E)
{
	int    i;
	double s = 0.0;

	for (i = 0; i < p->nlam; i++) {
		const double w = band_weight(p, i, 400.0, 700.0);
		if (w > 0.0) s += E[i] * PHOTON_K(p->lam[i]) * p->action[i] * w;
	}
	return s;
}

/* 照度 [lx] : 一般照明タブとの比較用 */
double calc_lux(const ppfd_t *p, const double *E)
{
	int    i;
	double s = 0.0;

	for (i = 0; i < p->nlam; i++) {
		s += E[i] * p->vlambda[i];
	}
	return K_MAX * s;
}

/* 放射照度 [W/m²] (全グリッド積分) */
double calc_irradiance(const ppfd_t *p, const double *E)
{
	int    i;
	double s = 0.0;

	for (i = 0; i < p->nlam; i++) s += E[i];
	return s;
}
