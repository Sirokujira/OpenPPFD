/*
photometry.c

実測配光ファイル (IESNA LM-63 / EULUMDAT) の読み込みと補間。

■ 何を使い、何を捨てるか
配光ファイルから使うのは**配光の形だけ**で、記載されている光度の絶対値
(cd) と定格光束 (lm) は使わない。理由は本プロジェクトの前提そのもので、
lm は V(λ) 基準の測光量なので、SPD を知らずに放射束 [W] や PPF へは
換算できないため。放射束は従来どおり led / array キーの W 欄 (または
ppf / lumens オプション) で与え、ファイルは

    I(γ, C) = Φ * Î(γ, C)          ∫ Î dΩ = 1

の Î だけを供給する。こうすると配光ファイルを使っても放射束は厳密に
保存され、閉キャビティのエネルギー収支の検証がそのまま効く。

■ 正規化を厳密にやる
Î は (γ, C) の格子上の値を双線形補間したものと定義する。すると正規化
定数 ∫ I dΩ は補間関数の厳密な積分として閉形式で書ける :

    C 方向 : 区分線形なので台形則が厳密
    γ 方向 : J(γ) が区分線形のとき
        ∫_a^b (p + q(γ-a)) sinγ dγ
            = p (cos a - cos b) + q [ (sin b - sin a) - (b-a) cos b ]

台形則で近似すると O(Δγ²) の誤差がそのまま放射束の誤差になり
closure error を汚すので、ここは厳密にやる価値がある。

■ 角度の定義
γ = 配光軸 (dir) からの角度 [deg]。IES 光度分布型 C / EULUMDAT と同じく
γ = 0 が配光軸方向 (既定では真下)。C = 配光軸まわりの方位角で、
C = 0 が +x 方向、C = 90 が +y 方向 (既定の dir = 0 0 -1 のとき)。
`rot <deg>` で器具を配光軸まわりに回せる。
表の範囲外の γ は光度 0 (上向き成分の無い器具の表が 0-90 で終わるため)。
*/

#include "ppfd.h"

/* ---- 数値トークンの読み出し (空白 / 改行 / カンマ区切り) ----------- */
static int read_num(FILE *fp, double *v)
{
	char buf[128];
	int  n = 0;

	for (;;) {
		const int c = fgetc(fp);
		if (c == EOF) break;
		if ((c == ' ') || (c == '\t') || (c == '\r') || (c == '\n') || (c == ',')) {
			if (n > 0) break;
			continue;
		}
		if (n < (int)sizeof(buf) - 1) buf[n++] = (char)c;
	}
	if (n == 0) return 0;
	buf[n] = '\0';
	*v = atof(buf);
	return 1;
}

static int read_nums(FILE *fp, double *v, int n)
{
	int i;
	for (i = 0; i < n; i++) {
		if (!read_num(fp, &v[i])) return 0;
	}
	return 1;
}

/* 1 行読み捨て */
static int skip_line(FILE *fp, char *buf, size_t sz)
{
	return (fgets(buf, (int)sz, fp) != NULL);
}

/* ---- 正規化 : ∫ Î dΩ = 1 にする。戻り値 = 正規化前の ∫ I dΩ ------- */
/* ∫_a^b (Ka + q (x - a)) cos x dx (区分線形な関数と cos の積の閉形式) */
static double lin_cos(double a, double b, double ka, double kb)
{
	const double h = b - a;
	const double q = (h > 0.0) ? ((kb - ka) / h) : 0.0;
	return (ka * (sin(b) - sin(a))) + (q * ((h * sin(b)) + cos(b) - cos(a)));
}

/*
光度分布型 A / B の正規化。表は水平角 H (cpl) x 鉛直角 V (gam) の双線形。
  型 A : dΩ = cos V dV dH (H は鉛直軸まわりの方位角、V は仰角)
  型 B : dΩ = cos H dH dV (V は横軸まわりの面の回転角、H は横軸への緯度)
内側の積分は区分線形なので台形則が厳密、外側は lin_cos の閉形式。
*/
static double photdist_normalize_ab(photdist_t *d)
{
	const double deg = PI / 180.0;
	const int ng = d->ng, nc = d->nc;
	double total = 0.0;
	int    ig, ic;

	if (d->ptype == 3) {
		/* 型 A : K(V) = ∫ I dH、total = ∫ K(V) cos V dV */
		double *K = (double *)xcalloc((size_t)ng, sizeof(double));
		for (ig = 0; ig < ng; ig++) {
			if (nc <= 1) {
				K[ig] = PI * d->I[ig];                  /* H は -90..90 で一定 */
			}
			else {
				for (ic = 0; ic < nc - 1; ic++) {
					K[ig] += 0.5 * (d->cpl[ic + 1] - d->cpl[ic]) * deg
					       * (d->I[((size_t)ic * ng) + ig] + d->I[((size_t)(ic + 1) * ng) + ig]);
				}
			}
		}
		for (ig = 0; ig < ng - 1; ig++) {
			total += lin_cos(d->gam[ig] * deg, d->gam[ig + 1] * deg, K[ig], K[ig + 1]);
		}
		free(K);
	}
	else {
		/* 型 B : K(H) = ∫ I dV、total = ∫ K(H) cos H dH */
		const int nh = (nc <= 1) ? 1 : nc;
		double *K = (double *)xcalloc((size_t)nh, sizeof(double));
		for (ic = 0; ic < nh; ic++) {
			for (ig = 0; ig < ng - 1; ig++) {
				K[ic] += 0.5 * (d->gam[ig + 1] - d->gam[ig]) * deg
				       * (d->I[((size_t)ic * ng) + ig] + d->I[((size_t)ic * ng) + ig + 1]);
			}
		}
		if (nh == 1) {
			total = 2.0 * K[0];                          /* ∫ cos H dH (-90..90) = 2 */
		}
		else {
			for (ic = 0; ic < nh - 1; ic++) {
				total += lin_cos(d->cpl[ic] * deg, d->cpl[ic + 1] * deg, K[ic], K[ic + 1]);
			}
		}
		free(K);
	}

	if (total > EPS) {
		const size_t n = (size_t)ng * ((nc > 0) ? nc : 1);
		size_t k;
		for (k = 0; k < n; k++) d->I[k] /= total;
	}
	return total;
}

static double photdist_normalize(photdist_t *d)
{
	const double deg = PI / 180.0;
	double *J = (double *)xcalloc((size_t)d->ng, sizeof(double));
	double  total = 0.0;
	int     ig, ic;

	/* J(γ) = ∫ I dC : C 方向は区分線形なので台形則が厳密 */
	for (ig = 0; ig < d->ng; ig++) {
		if (d->nc <= 1) {
			J[ig] = 2.0 * PI * d->I[ig];
		}
		else {
			double s = 0.0;
			for (ic = 0; ic < d->nc - 1; ic++) {
				const double dc = (d->cpl[ic + 1] - d->cpl[ic]) * deg;
				s += 0.5 * dc * (d->I[((size_t)ic * d->ng) + ig] + d->I[((size_t)(ic + 1) * d->ng) + ig]);
			}
			J[ig] = s;
		}
	}

	/* ∫ J(γ) sinγ dγ : J が区分線形なので閉形式 */
	for (ig = 0; ig < d->ng - 1; ig++) {
		const double a = d->gam[ig] * deg;
		const double b = d->gam[ig + 1] * deg;
		const double h = b - a;
		double q;
		if (h <= 0.0) continue;
		q = (J[ig + 1] - J[ig]) / h;
		total += (J[ig] * (cos(a) - cos(b)))
		       + (q * ((sin(b) - sin(a)) - (h * cos(b))));
	}

	free(J);

	if (total > EPS) {
		const size_t n = (size_t)d->ng * ((d->nc > 0) ? d->nc : 1);
		size_t k;
		for (k = 0; k < n; k++) d->I[k] /= total;
	}
	return total;
}

/*
IES の水平角 (C 面) を対称性から 0..360 へ展開する。
入力 c[0..n-1] は昇順で c[0] = 0、末尾が 90 / 180 / 360 のいずれか。
戻り値 = 展開後の点数、*cout / *iout (各点が参照する元の面番号) を確保する。
*/
static int expand_cplanes(const double *c, int n, double **cout, int **iout)
{
	const double cmax = c[n - 1];
	double *ca;
	int    *ia;
	int     m = 0, i;

	if (n < 2) return 0;

	ca = (double *)xmalloc((size_t)(4 * n) * sizeof(double));
	ia = (int *)xmalloc((size_t)(4 * n) * sizeof(int));

	for (i = 0; i < n; i++) {
		ca[m] = c[i];
		ia[m] = i;
		m++;
	}
	if (fabs(cmax - 90.0) < 1e-6) {
		/* 4 分の 1 対称 : 90-180, 180-270, 270-360 を鏡像で埋める */
		for (i = n - 2; i >= 0; i--) { ca[m] = 180.0 - c[i]; ia[m] = i; m++; }
		for (i = 1; i < n; i++)      { ca[m] = 180.0 + c[i]; ia[m] = i; m++; }
		for (i = n - 2; i >= 0; i--) { ca[m] = 360.0 - c[i]; ia[m] = i; m++; }
	}
	else if (fabs(cmax - 180.0) < 1e-6) {
		/* 左右対称 : 180-360 を鏡像で埋める */
		for (i = n - 2; i >= 0; i--) { ca[m] = 360.0 - c[i]; ia[m] = i; m++; }
	}
	else if (cmax < (360.0 - 1e-6)) {
		/* 全周だが末尾が 360 で閉じていない (0..350 等) : 周期性で閉じる */
		ca[m] = 360.0;
		ia[m] = 0;
		m++;
	}
	/* cmax が 360 ならそのまま (非対称) */

	*cout = ca;
	*iout = ia;
	return m;
}

/* ---- IESNA LM-63 -------------------------------------------------- */
static int load_ies(photdist_t *d, const char *path)
{
	FILE   *fp = fopen(path, "r");
	char    buf[BUFSIZ];
	double  hdr[13], tilt[4];
	double *gam = NULL, *craw = NULL, *iraw = NULL;
	double *cexp = NULL;
	int    *cidx = NULL;
	int     ng, ncr, nce, i, j, ok = 0;
	double  scale;

	if (fp == NULL) return 1;

	/* TILT= の行まで読み飛ばす */
	for (;;) {
		if (!skip_line(fp, buf, sizeof(buf))) { fclose(fp); return 1; }
		if (!strncmp(buf, "TILT=", 5) || !strncmp(buf, "tilt=", 5)) break;
	}
	if ((strstr(buf, "INCLUDE") != NULL) || (strstr(buf, "include") != NULL)) {
		/* 傾斜補正データ : 本ソルバーは器具を傾けないので読み飛ばす */
		int nt;
		if (!read_nums(fp, tilt, 2)) { fclose(fp); return 1; }
		nt = (int)(tilt[1] + 0.5);
		for (i = 0; i < (2 * nt); i++) {
			double v;
			if (!read_num(fp, &v)) { fclose(fp); return 1; }
		}
	}

	if (!read_nums(fp, hdr, 10)) { fclose(fp); return 1; }
	if (!read_nums(fp, hdr + 10, 3)) { fclose(fp); return 1; }

	ng  = (int)(hdr[3] + 0.5);
	ncr = (int)(hdr[4] + 0.5);
	if ((ng < 2) || (ncr < 1) || (ng > 4096) || (ncr > 4096)) { fclose(fp); return 1; }

	/* 光度の倍率 : multiplier * ballast factor (正規化するので形にのみ効く) */
	scale = hdr[2] * hdr[10];
	if (scale == 0.0) scale = 1.0;

	gam  = (double *)xmalloc((size_t)ng * sizeof(double));
	craw = (double *)xmalloc((size_t)ncr * sizeof(double));
	iraw = (double *)xmalloc((size_t)ng * ncr * sizeof(double));

	if (read_nums(fp, gam, ng) && read_nums(fp, craw, ncr)) {
		ok = 1;
		for (j = 0; j < ncr; j++) {
			if (!read_nums(fp, iraw + ((size_t)j * ng), ng)) { ok = 0; break; }
		}
	}
	fclose(fp);
	if (!ok) { free(gam); free(craw); free(iraw); return 1; }

	for (i = 0; i < ng * ncr; i++) iraw[i] *= scale;

	/* 光度分布型 (hdr[5]) : 1 = type C、2 = type B、3 = type A */
	d->ptype = (int)(hdr[5] + 0.5);
	d->watt  = hdr[12];

	if ((d->ptype == 2) || (d->ptype == 3)) {
		/*
		型 A / B : 鉛直角 V (gam) と水平角 H (cpl) はどちらも -90..90。
		0..90 だけが書かれていれば、その角について対称に展開する
		(V : 水平面について、H : 鉛直面について)。
		*/
		const int vm = (gam[0] > -EPS) && (ng > 1);            /* V を鏡映して展開 */
		const int hm = (ncr > 1) && (craw[0] > -EPS);          /* H を鏡映して展開 */
		const int ngx = vm ? ((2 * ng) - ((fabs(gam[0]) <= EPS) ? 1 : 0)) : ng;
		const int ncx = (ncr == 1) ? 1 : (hm ? ((2 * ncr) - ((fabs(craw[0]) <= EPS) ? 1 : 0)) : ncr);
		double *gx = (double *)xmalloc((size_t)ngx * sizeof(double));
		double *cx = (ncx > 1) ? (double *)xmalloc((size_t)ncx * sizeof(double)) : NULL;
		double *ix = (double *)xmalloc((size_t)ngx * ncx * sizeof(double));
		int     k, m;
		/* 角度の並び (昇順) と、元の表の番号 */
		int    *gi = (int *)xmalloc((size_t)ngx * sizeof(int));
		int    *ci = (int *)xmalloc((size_t)ncx * sizeof(int));
		m = 0;
		if (vm) {
			for (k = ng - 1; k >= ((fabs(gam[0]) <= EPS) ? 1 : 0); k--) { gx[m] = -gam[k]; gi[m] = k; m++; }
		}
		for (k = 0; k < ng; k++) { gx[m] = gam[k]; gi[m] = k; m++; }
		m = 0;
		if (ncx > 1) {
			if (hm) {
				for (k = ncr - 1; k >= ((fabs(craw[0]) <= EPS) ? 1 : 0); k--) { cx[m] = -craw[k]; ci[m] = k; m++; }
			}
			for (k = 0; k < ncr; k++) { cx[m] = craw[k]; ci[m] = k; m++; }
		}
		else {
			ci[0] = 0;
		}
		for (j = 0; j < ncx; j++) {
			for (k = 0; k < ngx; k++) {
				ix[((size_t)j * ngx) + k] = iraw[((size_t)ci[j] * ng) + gi[k]];
			}
		}
		free(gi);
		free(ci);
		free(gam);
		free(craw);
		free(iraw);
		d->ng  = ngx;
		d->gam = gx;
		d->nc  = ncx;
		d->cpl = cx;
		d->I   = ix;
		d->lm  = photdist_normalize_ab(d);
		return 0;
	}

	if (ncr == 1) {
		d->nc  = 1;
		d->cpl = NULL;
		d->I   = iraw;
	}
	else {
		nce = expand_cplanes(craw, ncr, &cexp, &cidx);
		if (nce < 2) { free(gam); free(craw); free(iraw); return 1; }
		d->nc  = nce;
		d->cpl = cexp;
		d->I   = (double *)xmalloc((size_t)nce * ng * sizeof(double));
		for (j = 0; j < nce; j++) {
			memcpy(d->I + ((size_t)j * ng), iraw + ((size_t)cidx[j] * ng), (size_t)ng * sizeof(double));
		}
		free(cidx);
		free(iraw);
	}
	free(craw);

	d->ng  = ng;
	d->gam = gam;
	d->lm  = photdist_normalize(d);
	return 0;
}

/* ---- EULUMDAT (LDT) ------------------------------------------------ */
static int load_ldt(photdist_t *d, const char *path)
{
	FILE   *fp = fopen(path, "r");
	char    buf[BUFSIZ];
	double  v[16];
	double *gam = NULL, *craw = NULL, *iraw = NULL, *cexp = NULL;
	int     isym, mc, ng, mc1, nlamp, nce;
	int     i, j, ok;
	double  lampflux = 0.0;

	if (fp == NULL) return 1;

	/* 1-26 行目 : ヘッダ (1 行 1 項目) */
	if (!skip_line(fp, buf, sizeof(buf))) { fclose(fp); return 1; }   /* 1 : 会社名 */
	for (i = 1; i < 7; i++) {
		if (!skip_line(fp, buf, sizeof(buf))) { fclose(fp); return 1; }
		v[i] = atof(buf);
	}
	isym = (int)(v[2] + 0.5);
	mc   = (int)(v[3] + 0.5);
	ng   = (int)(v[5] + 0.5);
	for (i = 7; i < 26; i++) {
		if (!skip_line(fp, buf, sizeof(buf))) { fclose(fp); return 1; }
	}
	nlamp = (int)(atof(buf) + 0.5);   /* 26 行目 = 標準ランプ組数 */
	if (nlamp < 1) nlamp = 1;
	if ((mc < 1) || (ng < 2) || (mc > 4096) || (ng > 4096)) { fclose(fp); return 1; }

	/* ランプ諸元 : 1 組 6 行。3 行目が総光束 [lm] */
	for (i = 0; i < nlamp; i++) {
		for (j = 0; j < 6; j++) {
			if (!skip_line(fp, buf, sizeof(buf))) { fclose(fp); return 1; }
			if ((i == 0) && (j == 2)) lampflux = fabs(atof(buf));
		}
	}

	/* 直接比 DR (10 個) は使わない */
	for (i = 0; i < 10; i++) {
		double t;
		if (!read_num(fp, &t)) { fclose(fp); return 1; }
	}

	/*
	Isym : 0 = 非対称、1 = 回転対称、2 = C0-C180 面対称、
	       3 = C90-C270 面対称、4 = 両面対称
	格納されている C 面の数はこれで決まる。
	*/
	switch (isym) {
	case 1:  mc1 = 1;              break;
	case 2:  mc1 = (mc / 2) + 1;   break;
	case 3:  mc1 = (mc / 2) + 1;   break;
	case 4:  mc1 = (mc / 4) + 1;   break;
	default: mc1 = mc;             break;
	}
	if (mc1 < 1) mc1 = 1;

	craw = (double *)xmalloc((size_t)mc * sizeof(double));
	gam  = (double *)xmalloc((size_t)ng * sizeof(double));
	iraw = (double *)xmalloc((size_t)mc1 * ng * sizeof(double));

	ok = read_nums(fp, craw, mc) && read_nums(fp, gam, ng);
	if (ok) {
		for (j = 0; j < mc1; j++) {
			if (!read_nums(fp, iraw + ((size_t)j * ng), ng)) { ok = 0; break; }
		}
	}
	fclose(fp);
	if (!ok) { free(craw); free(gam); free(iraw); return 1; }

	/* 光度は cd/1000lm なので定格光束を掛けて cd にする (形は変わらない) */
	if (lampflux > 0.0) {
		for (i = 0; i < mc1 * ng; i++) iraw[i] *= lampflux * 1e-3;
	}

	if (isym == 1) {
		d->nc  = 1;
		d->cpl = NULL;
		d->I   = iraw;
	}
	else {
		/* 一様な C 面 0, Δ, ..., 360 (末尾は先頭の複製) へ対称性を展開する */
		const double dc = 360.0 / mc;
		nce = mc + 1;
		cexp = (double *)xmalloc((size_t)nce * sizeof(double));
		d->I = (double *)xmalloc((size_t)nce * ng * sizeof(double));
		for (j = 0; j < nce; j++) {
			int src;
			double a = (j % mc) * dc;
			cexp[j] = j * dc;
			switch (isym) {
			case 2: {                       /* C0-C180 面対称 : I(-x) = I(x) */
				double x = (a > 180.0) ? (360.0 - a) : a;
				src = (int)((x / dc) + 0.5);
				break;
			}
			case 3: {                       /* C90-C270 面対称 : 格納は C90..C270 */
				double x = a - 90.0;
				if (x < -180.0) x += 360.0;
				if (x > 180.0) x -= 360.0;
				src = (int)((fabs(x) / dc) + 0.5);
				break;
			}
			case 4: {                       /* 両面対称 : 格納は C0..C90 */
				double x = a;
				while (x >= 180.0) x -= 180.0;
				if (x > 90.0) x = 180.0 - x;
				src = (int)((x / dc) + 0.5);
				break;
			}
			default:
				src = j % mc;
				break;
			}
			if (src < 0) src = 0;
			if (src > mc1 - 1) src = mc1 - 1;
			memcpy(d->I + ((size_t)j * ng), iraw + ((size_t)src * ng), (size_t)ng * sizeof(double));
		}
		d->nc  = nce;
		d->cpl = cexp;
		free(iraw);
	}
	free(craw);

	d->ng    = ng;
	d->gam   = gam;
	d->ptype = 1;
	d->lm    = photdist_normalize(d);
	return 0;
}

/*
配光ファイルを読み込んで登録する。同じパスを 2 度読まない (array で
同じ器具を並べても 1 回で済む)。戻り値 = 配光番号、失敗は -1。
*/
int photdist_load(ppfd_t *p, const char *path, int isldt)
{
	photdist_t d;
	int        i, rc;

	for (i = 0; i < p->ndist; i++) {
		if (!strcmp(p->dist[i].path, path) && (p->dist[i].isldt == isldt)) return i;
	}

	memset(&d, 0, sizeof(d));
	strncpy(d.path, path, sizeof(d.path) - 1);
	d.isldt = isldt;

	rc = isldt ? load_ldt(&d, path) : load_ies(&d, path);
	if (rc) return -1;

	p->dist = (photdist_t *)realloc(p->dist, (size_t)(p->ndist + 1) * sizeof(photdist_t));
	if (p->dist == NULL) { fprintf(stderr, "*** out of memory\n"); exit(1); }
	p->dist[p->ndist] = d;
	return p->ndist++;
}

/* 型 A / B の Î(H, V) : 表の範囲外は 0 (双線形) */
static double photdist_value_hv(const photdist_t *d, double h, double v)
{
	int    ig, ic;
	double fg, fc, i0, i1;

	if ((v < d->gam[0]) || (v > d->gam[d->ng - 1])) return 0.0;
	ig = 0;
	while ((ig < d->ng - 2) && (v > d->gam[ig + 1])) ig++;
	{
		const double hh = d->gam[ig + 1] - d->gam[ig];
		fg = (hh > EPS) ? ((v - d->gam[ig]) / hh) : 0.0;
	}
	if (fg < 0.0) fg = 0.0;
	if (fg > 1.0) fg = 1.0;
	if (d->nc <= 1) {
		if ((h < -90.0) || (h > 90.0)) return 0.0;
		return (d->I[ig] * (1.0 - fg)) + (d->I[ig + 1] * fg);
	}
	if ((h < d->cpl[0]) || (h > d->cpl[d->nc - 1])) return 0.0;
	ic = 0;
	while ((ic < d->nc - 2) && (h > d->cpl[ic + 1])) ic++;
	{
		const double hh = d->cpl[ic + 1] - d->cpl[ic];
		fc = (hh > EPS) ? ((h - d->cpl[ic]) / hh) : 0.0;
	}
	if (fc < 0.0) fc = 0.0;
	if (fc > 1.0) fc = 1.0;
	i0 = (d->I[((size_t)ic * d->ng) + ig] * (1.0 - fg)) + (d->I[((size_t)ic * d->ng) + ig + 1] * fg);
	i1 = (d->I[((size_t)(ic + 1) * d->ng) + ig] * (1.0 - fg)) + (d->I[((size_t)(ic + 1) * d->ng) + ig + 1] * fg);
	return (i0 * (1.0 - fc)) + (i1 * fc);
}

/*
光源 e の向き u (単位ベクトル) への相対光度 Î [1/sr]。
  型 C (と EULUMDAT) : γ = 配光軸からの角、C = cx から cy へ向かう方位 - rot。
  型 A / B : 器具の座標 x = 配光軸 (dir)、z = 「上」(+z を配光軸に直交させた
      向き。配光軸が鉛直なら cx)、y = z x x を rot だけ配光軸まわりに回す。
      型 A : V = 仰角 asin(z)、H = 方位角 atan2(y, x)
      型 B : V = 横軸 y まわりの面の回転角 atan2(z, x)、H = 緯度 asin(y)
*/
double photdist_dir(const ppfd_t *p, const emitter_t *e, vec3_t u)
{
	const photdist_t *d = &p->dist[e->idist];

	if ((d->ptype != 2) && (d->ptype != 3)) {
		double cg = v_dot(u, e->dir);
		double ca;
		if (cg > 1.0) cg = 1.0;
		if (cg < -1.0) cg = -1.0;
		ca = atan2(v_dot(u, e->cy), v_dot(u, e->cx)) * (180.0 / PI) - e->crot;
		return photdist_value(d, acos(cg) * (180.0 / PI), ca);
	}
	{
		const vec3_t ex = e->dir;
		vec3_t ez = (fabs(ex.z) < 0.99) ? v_make(0.0, 0.0, 1.0) : e->cx;
		vec3_t ey;
		double x, y, z, h, v;
		const double cr = cos(e->crot * (PI / 180.0)), sr = sin(e->crot * (PI / 180.0));
		ez = v_sub(ez, v_scale(ex, v_dot(ez, ex)));
		ez = v_scale(ez, 1.0 / v_norm(ez));
		ey = v_cross(ez, ex);
		{
			/* rot : 器具を配光軸まわりに回す */
			const vec3_t ty = v_add(v_scale(ey, cr), v_scale(ez, sr));
			const vec3_t tz = v_sub(v_scale(ez, cr), v_scale(ey, sr));
			ey = ty;
			ez = tz;
		}
		x = v_dot(u, ex);
		y = v_dot(u, ey);
		z = v_dot(u, ez);
		if (y > 1.0) y = 1.0;
		if (y < -1.0) y = -1.0;
		if (z > 1.0) z = 1.0;
		if (z < -1.0) z = -1.0;
		if (d->ptype == 3) {
			v = asin(z) * (180.0 / PI);
			h = atan2(y, x) * (180.0 / PI);
		}
		else {
			v = atan2(z, x) * (180.0 / PI);
			h = asin(y) * (180.0 / PI);
		}
		return photdist_value_hv(d, h, v);
	}
}

/*
相対光度 Î(γ, C) [1/sr]。∫ Î dΩ = 1 に正規化済みなので、放射束 Φ の
器具の放射強度は I = Φ * Î。表の γ 範囲外は 0。
*/
double photdist_value(const photdist_t *d, double gam, double cang)
{
	int    ig, ic;
	double fg, fc, i0, i1;

	if ((gam < d->gam[0]) || (gam > d->gam[d->ng - 1])) return 0.0;

	/* γ の区間を探す (表は高々数百点なので線形探索で十分) */
	ig = 0;
	while ((ig < d->ng - 2) && (gam > d->gam[ig + 1])) ig++;
	{
		const double h = d->gam[ig + 1] - d->gam[ig];
		fg = (h > EPS) ? ((gam - d->gam[ig]) / h) : 0.0;
	}
	if (fg < 0.0) fg = 0.0;
	if (fg > 1.0) fg = 1.0;

	if (d->nc <= 1) {
		return (d->I[ig] * (1.0 - fg)) + (d->I[ig + 1] * fg);
	}

	cang = fmod(cang, 360.0);
	if (cang < 0.0) cang += 360.0;

	ic = 0;
	while ((ic < d->nc - 2) && (cang > d->cpl[ic + 1])) ic++;
	{
		const double h = d->cpl[ic + 1] - d->cpl[ic];
		fc = (h > EPS) ? ((cang - d->cpl[ic]) / h) : 0.0;
	}
	if (fc < 0.0) fc = 0.0;
	if (fc > 1.0) fc = 1.0;

	i0 = (d->I[((size_t)ic * d->ng) + ig] * (1.0 - fg))
	   + (d->I[((size_t)ic * d->ng) + ig + 1] * fg);
	i1 = (d->I[((size_t)(ic + 1) * d->ng) + ig] * (1.0 - fg))
	   + (d->I[((size_t)(ic + 1) * d->ng) + ig + 1] * fg);

	return (i0 * (1.0 - fc)) + (i1 * fc);
}
