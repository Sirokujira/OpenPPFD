/*
input_data.c

.ppfd 入力ファイルの読み込み (key = value 形式、OpenFDTD sol/input_data.c の方式)。

分光グリッド (wavelength キー) と作用曲線 (actionspectrum キー) は他の
すべてのキーの解釈に必要なので、ファイルを 2 回走査する
(1 回目 : グリッドの確定 / 2 回目 : 本体)。

キー省略時の既定値は「省略しても意味が変わらない」値に初期化し、
未知のキーは黙って読み飛ばす (前方互換)。
*/

#include "ppfd.h"

#define MAXTOKEN 512
#define K_BOLTZ (1.380649e-23)

/* 動的配列の伸長 (C99 VLA は MSVC 非対応のため malloc/realloc + フラット配列) */
#define APPEND(ptr, num, cap, type) \
	do { \
		if ((num) >= (cap)) { \
			(cap) = (cap) ? (2 * (cap)) : 16; \
			(ptr) = (type *)realloc((ptr), (size_t)(cap) * sizeof(type)); \
			if ((ptr) == NULL) { fprintf(stderr, "*** out of memory\n"); exit(1); } \
		} \
	} while (0)

/* 昇順に並べ替える (遮蔽物の範囲指定はどちら向きでもよい) */
#define SORT2(a, b) \
	do { \
		if ((a) > (b)) { const double t_ = (a); (a) = (b); (b) = t_; } \
	} while (0)

int find_mat(const ppfd_t *p, const char *name)
{
	int i;
	for (i = 0; i < p->nmat; i++) {
		if (!strcmp(p->mat[i].name, name)) return i;
	}
	return -1;
}

int find_spec(const ppfd_t *p, const char *name)
{
	int i;
	for (i = 0; i < p->nspec; i++) {
		if (!strcmp(p->spec[i].name, name)) return i;
	}
	return -1;
}

/* 面名 -> 面番号 (0:xmin 1:xmax 2:ymin 3:ymax 4:zmin 5:zmax) */
static int faceid(const char *s)
{
	if (streq_ci(s, "xmin")) return 0;
	if (streq_ci(s, "xmax")) return 1;
	if (streq_ci(s, "ymin")) return 2;
	if (streq_ci(s, "ymax")) return 3;
	if (streq_ci(s, "zmin") || streq_ci(s, "floor")) return 4;
	if (streq_ci(s, "zmax") || streq_ci(s, "ceiling")) return 5;
	return -1;
}

/* 行を 1 本読んで正規化する。EOF / "end" で 0 を返す */
static int readline_norm(FILE *fp, char *buf, size_t sz)
{
	while (fgets(buf, (int)sz, fp) != NULL) {
		size_t n;
		if (strlen(buf) <= 1) continue;
		if (buf[0] == '#') continue;
		n = strlen(buf);
		while ((n > 0) && ((buf[n - 1] == '\n') || (buf[n - 1] == '\r'))) {
			buf[--n] = '\0';
		}
		if (!strncmp(buf, "end", 3) && ((buf[3] == '\0') || (buf[3] == ' ') || (buf[3] == '\t'))) return 0;
		if (n == 0) continue;
		return 1;
	}
	return 0;
}

/* "λ:v" または "λ:v1:v2" を分解する。戻り値 = 取れた値の数 (1 or 2)、失敗は 0 */
static int parse_pair(const char *tok, double *lam, double *v1, double *v2)
{
	char   tmp[256];
	char  *c1, *c2;

	strncpy(tmp, tok, sizeof(tmp) - 1);
	tmp[sizeof(tmp) - 1] = '\0';
	c1 = strchr(tmp, ':');
	if (c1 == NULL) return 0;
	*c1++ = '\0';
	c2 = strchr(c1, ':');
	if (c2 != NULL) *c2++ = '\0';

	*lam = atof(tmp);
	*v1 = atof(c1);
	if (c2 != NULL) {
		*v2 = atof(c2);
		return 2;
	}
	return 1;
}

/*
分光テーブルをグリッドへ再標本化する。
tok[0..ntok-1] は "λ:v[:v2]" 形式、または file <path> のときは読み込み済みの配列。
*/
static void resample(const ppfd_t *p, const double *xl, const double *y, int n, double *out)
{
	int i;
	for (i = 0; i < p->nlam; i++) {
		out[i] = interp_table(xl, y, n, p->lam[i]);
	}
}

/* "λ v [v2]" 形式のテキストを読む。戻り値 = 行数 (0 = 失敗) */
static int read_tablefile(const char *path, double **xl, double **y1, double **y2)
{
	FILE   *fp = fopen(path, "r");
	char    buf[BUFSIZ];
	int     n = 0, cap = 0;
	double *a = NULL, *b = NULL, *c = NULL;

	if (fp == NULL) return 0;

	while (fgets(buf, sizeof(buf), fp) != NULL) {
		double lam = 0, v1 = 0, v2 = 0;
		int    nv;
		if ((buf[0] == '#') || (strlen(buf) <= 1)) continue;
		nv = sscanf(buf, "%lf %lf %lf", &lam, &v1, &v2);
		if (nv < 2) continue;
		if (n >= cap) {
			cap = cap ? (2 * cap) : 64;
			a = (double *)realloc(a, (size_t)cap * sizeof(double));
			b = (double *)realloc(b, (size_t)cap * sizeof(double));
			c = (double *)realloc(c, (size_t)cap * sizeof(double));
			if ((a == NULL) || (b == NULL) || (c == NULL)) { fprintf(stderr, "*** out of memory\n"); exit(1); }
		}
		a[n] = lam;
		b[n] = v1;
		c[n] = (nv >= 3) ? v2 : 0.0;
		n++;
	}
	fclose(fp);

	*xl = a;
	*y1 = b;
	*y2 = c;
	return n;
}

/*
材料の分光特性を組み立てる。tok[k] 以降が mode + 引数。
戻り値 0 = 成功。
*/
static int build_material(ppfd_t *p, pmat_t *m, char **tok, int ntok, int k)
{
	int i;

	m->rho = (double *)xcalloc((size_t)p->nlam, sizeof(double));
	m->tau = (double *)xcalloc((size_t)p->nlam, sizeof(double));
	m->rhos = 0.0;

	/*
	鏡面反射率は末尾の "specular <v>" で与える (モードによらず共通)。
	波長非依存のスカラ : 反射フィルムや研磨アルミはほぼ中性なので、
	分光にせず 1 個の値で持つ (README に明記)。
	*/
	for (i = k; i + 1 < ntok; i++) {
		if (streq_ci(tok[i], "specular") || streq_ci(tok[i], "mirror")) {
			m->rhos = atof(tok[i + 1]);
			if (m->rhos < 0.0) m->rhos = 0.0;
			if (m->rhos > 1.0) m->rhos = 1.0;
			ntok = i;               /* 以降はモード引数から外す */
			break;
		}
	}

	if (k >= ntok) return 1;

	if (streq_ci(tok[k], "gray") || streq_ci(tok[k], "grey")) {
		const double r = (k + 1 < ntok) ? atof(tok[k + 1]) : 0.0;
		const double t = (k + 2 < ntok) ? atof(tok[k + 2]) : 0.0;
		for (i = 0; i < p->nlam; i++) {
			m->rho[i] = r;
			m->tau[i] = t;
		}
		return 0;
	}
	else if (streq_ci(tok[k], "table")) {
		double *xl, *y1, *y2;
		int     n = 0, j;
		const int nmax = ntok - k - 1;
		if (nmax <= 0) return 1;
		xl = (double *)xmalloc((size_t)nmax * sizeof(double));
		y1 = (double *)xmalloc((size_t)nmax * sizeof(double));
		y2 = (double *)xmalloc((size_t)nmax * sizeof(double));
		for (j = k + 1; j < ntok; j++) {
			double lam = 0, v1 = 0, v2 = 0;
			if (parse_pair(tok[j], &lam, &v1, &v2) == 0) continue;
			xl[n] = lam;
			y1[n] = v1;
			y2[n] = v2;
			n++;
		}
		if (n == 0) { free(xl); free(y1); free(y2); return 1; }
		resample(p, xl, y1, n, m->rho);
		resample(p, xl, y2, n, m->tau);
		free(xl); free(y1); free(y2);
		return 0;
	}
	else if (streq_ci(tok[k], "file")) {
		double *xl = NULL, *y1 = NULL, *y2 = NULL;
		int     n;
		if (k + 1 >= ntok) return 1;
		n = read_tablefile(tok[k + 1], &xl, &y1, &y2);
		if (n == 0) return 1;
		resample(p, xl, y1, n, m->rho);
		resample(p, xl, y2, n, m->tau);
		free(xl); free(y1); free(y2);
		return 0;
	}

	return 1;
}

/* スペクトル (ビン重み、総和 1) を組み立てる。戻り値 0 = 成功 */
static int build_spectrum(ppfd_t *p, spec_t *s, char **tok, int ntok, int k)
{
	int i;

	s->w = (double *)xcalloc((size_t)p->nlam, sizeof(double));
	if (k >= ntok) return 1;

	if (streq_ci(tok[k], "mono")) {
		/*
		単色光。グリッド点に乗らない波長は隣接 2 点へ線形配分する。
		光量子換算係数 k(λ) は λ の 1 次式なので、この配分で PPFD は厳密。
		*/
		double f;
		int    i0;
		if (k + 1 >= ntok) return 1;
		f = (atof(tok[k + 1]) - p->lam0) / p->dlam;
		if ((f < -EPS) || (f > (p->nlam - 1) + EPS)) return 2;   /* グリッド外 */
		i0 = (int)f;
		if (i0 >= p->nlam - 1) {
			s->w[p->nlam - 1] = 1.0;
		}
		else {
			const double fr = f - i0;
			s->w[i0]     = 1.0 - fr;
			s->w[i0 + 1] = fr;
		}
		return 0;
	}
	else if (streq_ci(tok[k], "gauss")) {
		const double c = (k + 1 < ntok) ? atof(tok[k + 1]) : 0.0;
		const double fw = (k + 2 < ntok) ? atof(tok[k + 2]) : 0.0;
		if (fw <= 0.0) return 1;
		for (i = 0; i < p->nlam; i++) {
			const double t = (p->lam[i] - c) / fw;
			s->w[i] = exp(-4.0 * log(2.0) * t * t);
		}
		spec_normalize(s->w, p->nlam);
		return 0;
	}
	else if (streq_ci(tok[k], "sum")) {
		/* 3 つ組 (中心, 半値全幅, 重み) の繰り返し : 白色 LED = 青ピーク + 蛍光体 */
		int j;
		int npeak = 0;
		for (j = k + 1; j + 2 < ntok; j += 3) {
			const double c = atof(tok[j]);
			const double fw = atof(tok[j + 1]);
			const double a = atof(tok[j + 2]);
			if (fw <= 0.0) continue;
			for (i = 0; i < p->nlam; i++) {
				const double t = (p->lam[i] - c) / fw;
				s->w[i] += a * exp(-4.0 * log(2.0) * t * t);
			}
			npeak++;
		}
		if (npeak == 0) return 1;
		spec_normalize(s->w, p->nlam);
		return 0;
	}
	else if (streq_ci(tok[k], "blackbody")) {
		const double T = (k + 1 < ntok) ? atof(tok[k + 1]) : 0.0;
		if (T <= 0.0) return 1;
		for (i = 0; i < p->nlam; i++) {
			const double lm = p->lam[i] * 1e-9;
			const double x = (H_PLANCK * C0) / (lm * K_BOLTZ * T);
			s->w[i] = (x < 700.0) ? (1.0 / (pow(lm, 5.0) * (exp(x) - 1.0))) : 0.0;
		}
		spec_normalize(s->w, p->nlam);
		return 0;
	}
	else if (streq_ci(tok[k], "table") || streq_ci(tok[k], "file")) {
		double *xl = NULL, *y1 = NULL, *y2 = NULL;
		int     n = 0;
		if (streq_ci(tok[k], "file")) {
			if (k + 1 >= ntok) return 1;
			n = read_tablefile(tok[k + 1], &xl, &y1, &y2);
		}
		else {
			int j;
			const int nmax = ntok - k - 1;
			if (nmax <= 0) return 1;
			xl = (double *)xmalloc((size_t)nmax * sizeof(double));
			y1 = (double *)xmalloc((size_t)nmax * sizeof(double));
			y2 = (double *)xmalloc((size_t)nmax * sizeof(double));
			for (j = k + 1; j < ntok; j++) {
				double lam = 0, v1 = 0, v2 = 0;
				if (parse_pair(tok[j], &lam, &v1, &v2) == 0) continue;
				xl[n] = lam;
				y1[n] = v1;
				y2[n] = v2;
				n++;
			}
		}
		if (n == 0) { free(xl); free(y1); free(y2); return 1; }
		/*
		テーブルは分光放射束密度 [W/nm 相対]。グリッド中心で標本化して
		ビン重みにする (どのみち総和 1 へ正規化するので dλ 倍は不要)。
		範囲外は interp_table が端値でクランプするため、テーブルの外側で
		0 にしたい場合は両端に 0 の点を明示すること。
		*/
		resample(p, xl, y1, n, s->w);
		for (i = 0; i < p->nlam; i++) {
			if (s->w[i] < 0.0) s->w[i] = 0.0;
		}
		spec_normalize(s->w, p->nlam);
		free(xl); free(y1); free(y2);
		return 0;
	}

	return 1;
}

/*
スペクトルの裾を切って有効ビン範囲 [i0, i1] を決める。

gauss / blackbody は数学的にはどこまでも裾を引くが、ピークの 1e-10 以下は
物理的にも数値的にも無意味で、内側ループ長 (= 実行時間) を 3 倍近く
悪化させる。切ったあと総和 1 へ再正規化するので放射束は厳密に保存する。
mono はこれで 1〜2 ビンになる。
*/
static void spec_trim(spec_t *s, int nlam)
{
	const double thresh = 1e-10;
	double wmax = 0.0;
	int    i;

	for (i = 0; i < nlam; i++) {
		if (s->w[i] > wmax) wmax = s->w[i];
	}
	for (i = 0; i < nlam; i++) {
		if (s->w[i] < (thresh * wmax)) s->w[i] = 0.0;
	}
	spec_normalize(s->w, nlam);

	s->i0 = 0;
	s->i1 = nlam - 1;
	while ((s->i0 < nlam - 1) && (s->w[s->i0] == 0.0)) s->i0++;
	while ((s->i1 > s->i0) && (s->w[s->i1] == 0.0)) s->i1--;
}

/*
dir から正規直交系を 2 組作る。

  ax, ay : 面光源のローカル軸 (向きは任意でよい)
  cx, cy : 実測配光の C 面基準。dir = (0,0,-1) のとき cx = +x, cy = +y に
           なるよう +x を dir に直交化して作る (C=0 が +x、C=90 が +y)。
*/
static void make_frame(emitter_t *e)
{
	vec3_t up = v_make(0, 0, 1);
	vec3_t t;
	double s;

	if (fabs(v_dot(e->dir, up)) > 0.9) up = v_make(1, 0, 0);
	e->ax = v_unit(v_cross(up, e->dir));
	e->ay = v_unit(v_cross(e->dir, e->ax));

	t = v_make(1, 0, 0);
	if (fabs(v_dot(e->dir, t)) > 0.9) t = v_make(0, 1, 0);
	s = v_dot(t, e->dir);
	e->cx = v_unit(v_sub(t, v_scale(e->dir, s)));
	e->cy = v_unit(v_cross(e->cx, e->dir));
}

/* スペクトル既知の光源で PPF [µmol/s] から放射束 [W] へ */
static double flux_from_ppf(const ppfd_t *p, int is, double ppf)
{
	double k = 0.0;
	int    i;

	for (i = 0; i < p->nlam; i++) {
		const double w = band_weight(p, i, 400.0, 700.0);
		if (w > 0.0) k += p->spec[is].w[i] * PHOTON_K(p->lam[i]) * w;
	}
	return (k > EPS) ? (ppf / k) : 0.0;
}

/* スペクトル既知の光源で 光束 [lm] から放射束 [W] へ */
static double flux_from_lm(const ppfd_t *p, int is, double lm)
{
	double k = 0.0;
	int    i;

	for (i = 0; i < p->nlam; i++) {
		k += p->spec[is].w[i] * p->vlambda[i];
	}
	k *= K_MAX;
	return (k > EPS) ? (lm / k) : 0.0;
}

/* led / array 共通のオプション引数を読む。戻り値 0 = 成功 */
static int emitter_options(ppfd_t *p, emitter_t *e, char **tok, int ntok, int k)
{
	while (k < ntok) {
		if (streq_ci(tok[k], "beam") && (k + 1 < ntok)) {
			e->mexp = atof(tok[k + 1]);
			k += 2;
		}
		else if (streq_ci(tok[k], "iso")) {
			e->mexp = -1.0;
			k += 1;
		}
		else if (streq_ci(tok[k], "dir") && (k + 3 < ntok)) {
			e->dir = v_unit(v_make(atof(tok[k + 1]), atof(tok[k + 2]), atof(tok[k + 3])));
			k += 4;
		}
		else if (streq_ci(tok[k], "size") && (k + 2 < ntok)) {
			e->wsize = atof(tok[k + 1]);
			e->hsize = atof(tok[k + 2]);
			k += 3;
			if ((k + 1 < ntok) && (tok[k][0] >= '0') && (tok[k][0] <= '9')) {
				e->nu = atoi(tok[k]);
				e->nv = atoi(tok[k + 1]);
				k += 2;
			}
		}
		else if (streq_ci(tok[k], "input") && (k + 1 < ntok)) {
			e->watt = atof(tok[k + 1]);
			k += 2;
		}
		else if ((streq_ci(tok[k], "ies") || streq_ci(tok[k], "ldt")) && (k + 1 < ntok)) {
			e->idist = photdist_load(p, tok[k + 1], streq_ci(tok[k], "ldt"));
			if (e->idist < 0) {
				fprintf(stderr, "*** cannot read photometric file : %s\n", tok[k + 1]);
				return 1;
			}
			k += 2;
		}
		else if (streq_ci(tok[k], "rot") && (k + 1 < ntok)) {
			e->crot = atof(tok[k + 1]);
			k += 2;
		}
		else if (streq_ci(tok[k], "ppf") && (k + 1 < ntok)) {
			/* 放射束をカタログの PPF [µmol/s] で与える (W 欄は無視される) */
			e->flux = flux_from_ppf(p, e->ispec, atof(tok[k + 1]));
			k += 2;
		}
		else if (streq_ci(tok[k], "lumens") && (k + 1 < ntok)) {
			e->flux = flux_from_lm(p, e->ispec, atof(tok[k + 1]));
			k += 2;
		}
		else {
			return 1;   /* 未知のオプション */
		}
	}
	make_frame(e);
	if (e->nu < 1) e->nu = 1;
	if (e->nv < 1) e->nv = 1;
	return 0;
}

int input_data(FILE *fp, ppfd_t *p)
{
	int    nline = 0, ierr = 0;
	char   strline[BUFSIZ], strsave[BUFSIZ];
	char  *token[MAXTOKEN];
	const char sep[] = " \t";
	int    cmat = 0, cspec = 0, cemit = 0, ctarget = 0, cband = 0, cocc = 0;
	int    i;
	double *ax_lam = NULL, *ax_v = NULL, *ax_v2 = NULL;
	int    nax = 0;

	/* 既定値 */
	memset(p, 0, sizeof(ppfd_t));
	p->lam0 = 380.0;
	p->dlam = 5.0;
	p->nlam = 81;                   /* 380..780 nm */
	p->Lx = p->Ly = p->Lz = 1.0;
	p->ndivu = p->ndivv = 6;
	p->photoperiod = 16.0;
	p->maxiter = 200;
	p->converg = 1e-6;
	p->canopy.k0 = 0.5;
	p->canopy.nlayer = 8;
	p->specbounce = 2;

	/* ---- 1 回目 : 分光グリッドと作用曲線を確定する ---------------- */
	while (readline_norm(fp, strline, sizeof(strline))) {
		int ntok;
		strcpy(strsave, strline);
		ntok = tokenize(strline, sep, token, MAXTOKEN);
		if ((nline == 0) && (ntok >= 1) && !strcmp(token[0], PROGRAM)) {
			nline++;
			continue;
		}
		nline++;
		if ((ntok < 3) || strcmp(token[1], "=")) continue;

		if (!strcmp(token[0], "wavelength") && (ntok >= 5)) {
			const double l1 = atof(token[2]);
			const double l2 = atof(token[3]);
			const double dl = atof(token[4]);
			if ((dl <= 0.0) || (l2 <= l1)) {
				fprintf(stderr, "*** invalid wavelength data\n");
				ierr = 1;
			}
			else {
				p->lam0 = l1;
				p->dlam = dl;
				p->nlam = (int)floor(((l2 - l1) / dl) + 0.5) + 1;
				if (p->nlam < 2) p->nlam = 2;
			}
		}
		else if (!strcmp(token[0], "actionspectrum") && (ntok >= 4) && streq_ci(token[2], "file")) {
			nax = read_tablefile(token[3], &ax_lam, &ax_v, &ax_v2);
			if (nax == 0) {
				fprintf(stderr, "*** cannot read actionspectrum file : %s\n", token[3]);
				ierr = 1;
			}
		}
	}
	if (ierr) return ierr;

	spectrum_grid(p);
	if (nax > 0) {
		p->action = (double *)xmalloc((size_t)p->nlam * sizeof(double));
		resample(p, ax_lam, ax_v, nax, p->action);
		free(ax_lam); free(ax_v); free(ax_v2);
	}
	spectrum_builtin(p);

	/* 材料 0 = 完全吸収 (黒)。wall 未指定の面はこれになる */
	APPEND(p->mat, p->nmat, cmat, pmat_t);
	memset(&p->mat[0], 0, sizeof(pmat_t));
	strcpy(p->mat[0].name, "black");
	p->mat[0].rho = (double *)xcalloc((size_t)p->nlam, sizeof(double));
	p->mat[0].tau = (double *)xcalloc((size_t)p->nlam, sizeof(double));
	p->nmat = 1;

	/* ---- 2 回目 : 本体 -------------------------------------------- */
	rewind(fp);
	nline = 0;
	while (readline_norm(fp, strline, sizeof(strline))) {
		int ntok;
		strcpy(strsave, strline);
		ntok = tokenize(strline, sep, token, MAXTOKEN);
		if ((nline == 0) && (ntok >= 1) && !strcmp(token[0], PROGRAM)) {
			nline++;
			continue;
		}
		nline++;
		if (ntok < 3) continue;
		if (strcmp(token[1], "=")) continue;

		if (!strcmp(token[0], "title")) {
			const char *q = strstr(strsave, "=");
			strncpy(p->title, (q != NULL) ? (q + 1) : "", sizeof(p->title) - 1);
			p->title[sizeof(p->title) - 1] = '\0';
		}
		else if (!strcmp(token[0], "wavelength") || !strcmp(token[0], "actionspectrum")) {
			/* 1 回目で処理済み */
		}
		else if (!strcmp(token[0], "chamber")) {
			if (ntok < 5) { fprintf(stderr, "*** invalid chamber data\n"); ierr = 1; continue; }
			p->Lx = atof(token[2]);
			p->Ly = atof(token[3]);
			p->Lz = atof(token[4]);
			if ((p->Lx <= 0) || (p->Ly <= 0) || (p->Lz <= 0)) {
				fprintf(stderr, "*** invalid chamber data\n");
				ierr = 1;
			}
		}
		else if (!strcmp(token[0], "patchdiv")) {
			if (ntok < 4) { fprintf(stderr, "*** invalid patchdiv data\n"); ierr = 1; continue; }
			p->ndivu = atoi(token[2]);
			p->ndivv = atoi(token[3]);
			if ((p->ndivu < 1) || (p->ndivv < 1)) {
				fprintf(stderr, "*** invalid patchdiv data\n");
				ierr = 1;
			}
		}
		else if (!strcmp(token[0], "material")) {
			pmat_t m;
			memset(&m, 0, sizeof(m));
			strncpy(m.name, token[2], NAMELEN - 1);
			if (build_material(p, &m, token, ntok, 3)) {
				fprintf(stderr, "*** invalid material data : %s\n", strsave);
				ierr = 1;
				continue;
			}
			APPEND(p->mat, p->nmat, cmat, pmat_t);
			p->mat[p->nmat++] = m;
		}
		else if (!strcmp(token[0], "spectrum")) {
			spec_t s;
			int    rc;
			memset(&s, 0, sizeof(s));
			strncpy(s.name, token[2], NAMELEN - 1);
			rc = build_spectrum(p, &s, token, ntok, 3);
			if (rc == 2) {
				fprintf(stderr, "*** spectrum \"%s\" : wavelength outside the [%g, %g] nm grid\n",
					s.name, p->lam0, p->lam0 + ((p->nlam - 1) * p->dlam));
				ierr = 1;
				continue;
			}
			else if (rc) {
				fprintf(stderr, "*** invalid spectrum data : %s\n", strsave);
				ierr = 1;
				continue;
			}
			spec_trim(&s, p->nlam);
			APPEND(p->spec, p->nspec, cspec, spec_t);
			p->spec[p->nspec++] = s;
		}
		else if (!strcmp(token[0], "wall")) {
			int im, f;
			if (ntok < 4) { fprintf(stderr, "*** invalid wall data\n"); ierr = 1; continue; }
			im = find_mat(p, token[3]);
			if (im < 0) {
				fprintf(stderr, "*** unknown material : %s\n", token[3]);
				ierr = 1;
				continue;
			}
			if (streq_ci(token[2], "all")) {
				for (f = 0; f < 6; f++) p->wallmat[f] = im;
			}
			else if (streq_ci(token[2], "side")) {
				for (f = 0; f < 4; f++) p->wallmat[f] = im;
			}
			else {
				f = faceid(token[2]);
				if (f < 0) {
					fprintf(stderr, "*** unknown face : %s\n", token[2]);
					ierr = 1;
					continue;
				}
				p->wallmat[f] = im;
			}
		}
		else if (!strcmp(token[0], "led")) {
			emitter_t e;
			/* led = x y z W spectrum [options] */
			if (ntok < 7) { fprintf(stderr, "*** invalid led data\n"); ierr = 1; continue; }
			memset(&e, 0, sizeof(e));
			e.pos = v_make(atof(token[2]), atof(token[3]), atof(token[4]));
			e.flux = atof(token[5]);
			e.ispec = find_spec(p, token[6]);
			if (e.ispec < 0) {
				fprintf(stderr, "*** unknown spectrum : %s\n", token[6]);
				ierr = 1;
				continue;
			}
			e.dir = v_make(0, 0, -1);
			e.mexp = 1.0;
			e.idist = -1;
			e.nu = e.nv = 4;
			if (emitter_options(p, &e, token, ntok, 7)) {
				fprintf(stderr, "*** invalid led option : %s\n", strsave);
				ierr = 1;
				continue;
			}
			APPEND(p->emit, p->nemit, cemit, emitter_t);
			p->emit[p->nemit++] = e;
		}
		else if (!strcmp(token[0], "array")) {
			/* array = x0 y0 z nx ny px py W_each spectrum [options] */
			emitter_t e0;
			int    ax, ay, nxa, nya;
			double x0, y0, z0, px, py;
			if (ntok < 11) { fprintf(stderr, "*** invalid array data\n"); ierr = 1; continue; }
			x0 = atof(token[2]);
			y0 = atof(token[3]);
			z0 = atof(token[4]);
			nxa = atoi(token[5]);
			nya = atoi(token[6]);
			px = atof(token[7]);
			py = atof(token[8]);
			memset(&e0, 0, sizeof(e0));
			e0.flux = atof(token[9]);
			e0.ispec = find_spec(p, token[10]);
			if (e0.ispec < 0) {
				fprintf(stderr, "*** unknown spectrum : %s\n", token[10]);
				ierr = 1;
				continue;
			}
			if ((nxa < 1) || (nya < 1)) { fprintf(stderr, "*** invalid array data\n"); ierr = 1; continue; }
			e0.dir = v_make(0, 0, -1);
			e0.mexp = 1.0;
			e0.idist = -1;
			e0.nu = e0.nv = 4;
			if (emitter_options(p, &e0, token, ntok, 11)) {
				fprintf(stderr, "*** invalid array option : %s\n", strsave);
				ierr = 1;
				continue;
			}
			for (ay = 0; ay < nya; ay++) {
				for (ax = 0; ax < nxa; ax++) {
					emitter_t e = e0;
					e.pos = v_make(x0 + (ax * px), y0 + (ay * py), z0);
					APPEND(p->emit, p->nemit, cemit, emitter_t);
					p->emit[p->nemit++] = e;
				}
			}
		}
		else if (!strcmp(token[0], "target")) {
			target_t t;
			if (ntok < 6) { fprintf(stderr, "*** invalid target data\n"); ierr = 1; continue; }
			memset(&t, 0, sizeof(t));
			strncpy(t.name, token[2], NAMELEN - 1);
			t.z = atof(token[3]);
			t.nx = atoi(token[4]);
			t.ny = atoi(token[5]);
			if ((t.nx < 1) || (t.ny < 1)) { fprintf(stderr, "*** invalid target data\n"); ierr = 1; continue; }
			if (ntok >= 10) {
				t.x0 = atof(token[6]);
				t.x1 = atof(token[7]);
				t.y0 = atof(token[8]);
				t.y1 = atof(token[9]);
			}
			else {
				t.x0 = 0.0;
				t.x1 = -1.0;   /* setup_target でチャンバ全面に展開 */
				t.y0 = 0.0;
				t.y1 = -1.0;
			}
			APPEND(p->target, p->ntarget, ctarget, target_t);
			p->target[p->ntarget++] = t;
		}
		else if (!strcmp(token[0], "occluder")) {
			/* occluder = name x0 x1 y0 y1 z0 z1 material [div nu nv] */
			occluder_t o;
			int im;
			if (ntok < 10) { fprintf(stderr, "*** invalid occluder data\n"); ierr = 1; continue; }
			memset(&o, 0, sizeof(o));
			strncpy(o.name, token[2], NAMELEN - 1);
			o.x0 = atof(token[3]);  o.x1 = atof(token[4]);
			o.y0 = atof(token[5]);  o.y1 = atof(token[6]);
			o.z0 = atof(token[7]);  o.z1 = atof(token[8]);
			SORT2(o.x0, o.x1);
			SORT2(o.y0, o.y1);
			SORT2(o.z0, o.z1);
			im = find_mat(p, token[9]);
			if (im < 0) {
				fprintf(stderr, "*** unknown material : %s\n", token[9]);
				ierr = 1;
				continue;
			}
			o.imat = im;
			if ((ntok >= 13) && streq_ci(token[10], "div")) {
				o.ndivu = atoi(token[11]);
				o.ndivv = atoi(token[12]);
			}
			if (((o.x1 - o.x0) <= 0.0) && ((o.y1 - o.y0) <= 0.0)) {
				fprintf(stderr, "*** invalid occluder data (degenerate in 2 axes)\n");
				ierr = 1;
				continue;
			}
			if (((o.x1 - o.x0) <= 0.0) && ((o.z1 - o.z0) <= 0.0)) {
				fprintf(stderr, "*** invalid occluder data (degenerate in 2 axes)\n");
				ierr = 1;
				continue;
			}
			if (((o.y1 - o.y0) <= 0.0) && ((o.z1 - o.z0) <= 0.0)) {
				fprintf(stderr, "*** invalid occluder data (degenerate in 2 axes)\n");
				ierr = 1;
				continue;
			}
			APPEND(p->occ, p->nocc, cocc, occluder_t);
			p->occ[p->nocc++] = o;
		}
		else if (!strcmp(token[0], "canopy")) {
			int im;
			/* canopy = ztop zbot LAI G leafmaterial */
			if (ntok < 7) { fprintf(stderr, "*** invalid canopy data\n"); ierr = 1; continue; }
			p->canopy.ztop = atof(token[2]);
			p->canopy.zbot = atof(token[3]);
			p->canopy.lai  = atof(token[4]);
			p->canopy.k0   = atof(token[5]);
			im = find_mat(p, token[6]);
			if (im < 0) {
				fprintf(stderr, "*** unknown material : %s\n", token[6]);
				ierr = 1;
				continue;
			}
			p->canopy.imat = im;
			if (p->canopy.ztop <= p->canopy.zbot) {
				fprintf(stderr, "*** invalid canopy data (ztop <= zbot)\n");
				ierr = 1;
				continue;
			}
			p->canopy.a = p->canopy.lai / (p->canopy.ztop - p->canopy.zbot);
			p->canopy.on = 1;
		}
		else if (!strcmp(token[0], "band")) {
			band_t b;
			if (ntok < 5) { fprintf(stderr, "*** invalid band data\n"); ierr = 1; continue; }
			memset(&b, 0, sizeof(b));
			strncpy(b.name, token[2], NAMELEN - 1);
			b.lam1 = atof(token[3]);
			b.lam2 = atof(token[4]);
			APPEND(p->band, p->nband, cband, band_t);
			p->band[p->nband++] = b;
		}
		else if (!strcmp(token[0], "photoperiod")) {
			p->photoperiod = atof(token[2]);
		}
		else if (!strcmp(token[0], "solver")) {
			if (ntok >= 3) p->maxiter = atoi(token[2]);
			if (ntok >= 4) p->converg = atof(token[3]);
			if (p->maxiter < 1) p->maxiter = 1;
		}
		else if (!strcmp(token[0], "leafscatter")) {
			/*
			群落の散乱光をキャビティへ戻すか。既定 off は v1 と完全に同じ挙動。
			キー名を canopy* にしないのは ppfd_check.sh が群落なしの変種を
			grep -v '^canopy' で作っているため (巻き添えで消える)。
			*/
			p->canopy.scatter = (streq_ci(token[2], "on") || streq_ci(token[2], "twostream")) ? 1 : 0;
		}
		else if (!strcmp(token[0], "leaflayers")) {
			p->canopy.nlayer = atoi(token[2]);
			if (p->canopy.nlayer < 1) p->canopy.nlayer = 1;
			if (p->canopy.nlayer > 256) p->canopy.nlayer = 256;
		}
		else if (!strcmp(token[0], "specbounce")) {
			p->specbounce = atoi(token[2]);
			if (p->specbounce < 0) p->specbounce = 0;
			if (p->specbounce > 6) p->specbounce = 6;
		}
		else if (!strcmp(token[0], "quadrature")) {
			p->msub = atoi(token[2]);
			if (p->msub < 0) p->msub = 0;
		}
		/* 未知のキーは黙って無視 (前方互換) */
	}
	if (ierr) return ierr;

	/* 帯域の既定値 (band キーが 1 つも無いとき) */
	if (p->nband == 0) {
		static const char *nm[] = {"UV-A", "blue", "green", "red", "far-red"};
		static const double l1[] = {315.0, 400.0, 500.0, 600.0, 700.0};
		static const double l2[] = {400.0, 500.0, 600.0, 700.0, 750.0};
		for (i = 0; i < 5; i++) {
			band_t b;
			memset(&b, 0, sizeof(b));
			strncpy(b.name, nm[i], NAMELEN - 1);
			b.lam1 = l1[i];
			b.lam2 = l2[i];
			APPEND(p->band, p->nband, cband, band_t);
			p->band[p->nband++] = b;
		}
	}

	if (p->nemit == 0) {
		fprintf(stderr, "*** no light source (led / array)\n");
		return 1;
	}
	if (p->ntarget == 0) {
		fprintf(stderr, "*** no target plane (target)\n");
		return 1;
	}

	/*
	鏡像法は「折り返した経路を直線で表す」ので、群落や遮蔽物があると
	経路長も遮蔽判定も合わない。黙って誤差を出すより弾く。
	*/
	{
		int f, spec = 0;
		for (f = 0; f < 6; f++) {
			if (p->mat[p->wallmat[f]].rhos > 0.0) spec = 1;
		}
		/* 反射率の合計が 1 を超える材料は非物理 (収支の検証も壊れる) */
		if (spec) {
			int im2, il2;
			for (im2 = 0; im2 < p->nmat; im2++) {
				if (p->mat[im2].rhos <= 0.0) continue;
				for (il2 = 0; il2 < p->nlam; il2++) {
					if (p->mat[im2].rho[il2] + p->mat[im2].rhos > 1.0 + 1e-12) {
						fprintf(stderr, "*** material \"%s\" : rho_d + rho_s > 1\n",
							p->mat[im2].name);
						return 1;
					}
				}
			}
		}
		if (spec && (p->nocc > 0)) {
			fprintf(stderr, "*** specular walls cannot be combined with occluders\n");
			fprintf(stderr, "    (the mirror-image method folds the path; the shadow test\n");
			fprintf(stderr, "     on the straight image ray would be wrong)\n");
			return 1;
		}
		if (spec && p->canopy.on) {
			/*
			側壁 (x/y 面) の鏡映は z を保つので群落スラブは鏡映で不変。
			床・天井の鏡映はスラブを動かすが、鏡像への直線の z を周期 2Lz の
			三角波で折り返せば実経路の z プロファイルが厳密に得られるので、
			canopy_path (geometry.c) がその折り返しで経路長を測る。どちらも
			弧長は鏡映で保たれるため、群落減衰は鏡像への直線のまま厳密。
			散乱の還流だけは鏡像経路ぶんの層への預け入れが要るため未対応。
			*/
			if (p->canopy.scatter) {
				fprintf(stderr, "*** specular walls cannot be combined with leafscatter yet\n");
				return 1;
			}
		}
	}

	return 0;
}
