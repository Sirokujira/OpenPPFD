/*
canopy.c

群落の散乱光をキャビティへ戻す (衝突源 + 層別二流)。

■ なぜ書き直したか
v1 の群落は「減衰のみ」で、葉が反射・透過した光を捨てていた。捨てた分は
エネルギー収支の残りとして `absorbed by canopy` に計上されるので収支は
帳尻が合うが、物理的には過大な吸収であり、群落下部の PPFD と R:FR、
白い反射壁の効きを系統的に過小評価する。

■ 二重計上をどう避けるか (ここが要)
v1 の減衰係数 μ = G a √(1-ω) の √(1-ω) は「素の遮断率」ではない。
後方散乱比 β = 1/2 の二流方程式

    dF↓/dτ = -(1-ω+ωβ) F↓ + ωβ F↑
   -dF↑/dτ = -(1-ω+ωβ) F↑ + ωβ F↓        (τ = 素の遮断の光学深さ)

の減衰固有値が u = √(a² - s²) = √((1-ω/2)² - (ω/2)²) = √(1-ω) であって、
すでに前方散乱を含んだ「総流束の減衰」を表している。したがって
「μ を据え置いたまま遮断分の ω を散乱として戻す」と実効消散がもう一度
√(1-ω) 倍だけ弱くなり、二重計上になる。

本実装の回避は構造的である :

  (1) 非衝突ビームの減衰は **素の遮断 k = G a** だけを使う (√(1-ω) は
      ビーム側に一切現れない。k は波長によらない)。
  (2) ビームから取り除いた流束はその場の層に預け、(1-ω) を吸収へ、
      ω を拡散場の湧き出しへ配る (取り除いた量 = 吸収 + 湧き出しは恒等式)。
  (3) 拡散場を二流で解く。その解の固有値として √(1-ω) が **出力側に**
      再登場する。1 光子あたり √(1-ω) が現れるのはちょうど 1 回。

つまり v1 が「仮定」していた Goudriaan の √(1-ω) を、新モデルは
「導出結果」として再現する。data/sample/canopy_deep.ppfd の判定はこの
性質を直接縛っている。

■ 層カーネル (閉形式、ω→1 でも特異にならない形)
厚さ t (素の遮断の光学厚さ) の層について u = √(1-ω), x = u t,
sinhc(x) = sinh(x)/x, f(x) = 2(cosh x - 1)/x² とおくと

    D  = 2 cosh(x) + (2-ω) t sinhc(x)
    T  = 2 / D                                  (拡散透過率)
    R  = ω t sinhc(x) / D                       (拡散反射率)
    A  = 2(1-ω) [ t² f(x)/2 + t sinhc(x) ] / D  (拡散吸収率)
    E  = [ t f(x)/2 + sinhc(x) ] / D            (層内一様等方源の片側脱出率)

R + T + A = 1 は代数的恒等式 ((1-ω)t²f = x²f = 2(cosh x - 1) を使う)。
3 つを独立に評価しているので、この恒等式はそのまま実装の検証に使える。
R は因子 ω を、A は因子 (1-ω) を構造的に持つので、ω = 0 で R が厳密に
0.0、ω = 1 で A が厳密に 0.0 になる (打ち消しに頼らない)。

極限 :
    ω → 0  : R = 0, T = e^{-t}, A = 1 - e^{-t}     (Beer-Lambert に厳密一致)
    ω → 1  : A = 0, T = 2/(2+t), R = t/(2+t)
    t → ∞  : R → (1-u)/(1+u)                        (半無限群落のアルベド)
    一様源 : ω = 0 で E = (1-e^{-t})/(2t)            (解析積分に厳密一致)

■ 1 次元近似の限界 (README にも明記)
拡散場は水平一様な 1 本のカラムとして解く。群落自体が水平一様
(canopy キーは x-y の広がりを持たない) で、散乱光は空間的に滑らかな
補正項なのでこの近似は妥当だが、群落の高さ範囲に埋没した側壁パッチとの
横方向のやりとりは表現できない。埋没パッチの枚数と面積比はログに出す。
*/

#include "ppfd.h"

/*
層カーネル。omega = 葉の散乱係数、t = 素の遮断の光学厚さ (= G a Δz)。
R/T/A/E をそれぞれ独立に (1 - 他の 2 つ、という作り方をせずに) 返す。

厚い層で cosh/sinh がオーバーフローしないよう、全体を e^{-x} でスケールした
量で組む (x = u t)。恒等式は同じ係数で割るだけなので保存される :

    ch = cosh(x) e^{-x}  = (1 + e^{-2x})/2
    sc = sinhc(x) e^{-x} = -expm1(-2x)/(2x)
    fc = f(x) e^{-x}     = (expm1(-x))² / x²
    D  = 2 ch + (2-ω) t sc
    T  = 2 e^{-x} / D        R = ω t sc / D
    A  = 2(1-ω)[t² fc/2 + t sc] / D        E = [t fc/2 + sc] / D

expm1 を使うので x → 0 でも桁落ちしない。x = 0 (ω = 1 または t = 0) だけ
極限値 ch = sc = fc = 1 を直接入れる。
*/
void canopy_kernel(double omega, double t, double *R, double *T, double *A, double *E)
{
	const double u = sqrt((omega < 1.0) ? (1.0 - omega) : 0.0);
	const double x = u * t;
	double ch, sc, fc, ex, den;

	if (t <= 0.0) {                         /* 厚さ 0 の縮退層 */
		*R = 0.0;
		*T = 1.0;
		*A = 0.0;
		*E = 0.5;
		return;
	}
	if (x > 0.0) {
		const double e1 = expm1(-x);        /* e^{-x} - 1 (小さい x で桁落ちしない) */
		ex = exp(-x);                       /* 1 + e1 で作ると x が大きいとき桁落ちする */
		ch = 0.5 * (1.0 + (ex * ex));
		sc = -expm1(-2.0 * x) / (2.0 * x);
		fc = (e1 * e1) / (x * x);
	}
	else {                                  /* x = 0 : ω = 1 の極限 */
		ex = 1.0;
		ch = 1.0;
		sc = 1.0;
		fc = 1.0;
	}

	den = (2.0 * ch) + ((2.0 - omega) * t * sc);
	*T = (2.0 * ex) / den;
	*R = (omega * t * sc) / den;
	*A = (2.0 * (1.0 - omega) * ((0.5 * t * t * fc) + (t * sc))) / den;
	*E = ((0.5 * t * fc) + sc) / den;
}

/*
層カーネルの自己検証 (oppfd --selftest)。

すべて幾何にも入力にも依存しない代数的性質なので、3 OS の CI で同じ値が
出なければならない。macOS (Apple Silicon) だけ FMA の有無で結果が変わった
実例があるため、この種の恒等式は CI で常時縛る。戻り値 = 失敗数。
*/
static int st_chk(const char *lbl, double a, double e, double tol)
{
	const double d = (e != 0.0) ? fabs((a - e) / e) : fabs(a - e);
	const int ok = (d <= tol);

	printf("%-40s %.17g vs %.17g -> %s (%.2e)\n", lbl, a, e, ok ? "OK" : "NG", d);
	return ok ? 0 : 1;
}

int canopy_selftest(void)
{
	double R, T, A, E;
	double worst;
	int    bad = 0, i, j;

	printf("=== canopy two-stream layer kernel self test ===\n");

	/* (1) R + T + A = 1 : 代数的恒等式 (3 量は独立に評価している) */
	worst = 0.0;
	for (i = 0; i <= 100; i++) {
		for (j = 0; j <= 60; j++) {
			const double om = i / 100.0;
			const double t = pow(10.0, -6.0 + (j / 6.0));
			double s;
			canopy_kernel(om, t, &R, &T, &A, &E);
			s = fabs((R + T + A) - 1.0);
			if (s > worst) worst = s;
			if ((R < 0.0) || (T < 0.0) || (A < 0.0) || (E < 0.0)) bad++;
		}
	}
	bad += st_chk("R+T+A = 1 (101 x 61 grid)", worst, 0.0, 1e-13);

	/* (2) ω = 0 : Beer-Lambert に一致し、反射は厳密に 0 */
	for (i = 0; i < 5; i++) {
		const double t = 0.25 * pow(2.0, i);
		canopy_kernel(0.0, t, &R, &T, &A, &E);
		bad += st_chk("omega=0 : T = exp(-t)", T, exp(-t), 1e-13);
		bad += st_chk("omega=0 : R = 0 exactly", R, 0.0, 0.0);
		bad += st_chk("omega=0 : escape = (1-e^-t)/2t", E, (1.0 - exp(-t)) / (2.0 * t), 1e-13);
	}

	/* (3) ω = 1 : 吸収は厳密に 0 */
	for (i = 0; i < 3; i++) {
		const double t = 0.25 * pow(4.0, i);
		canopy_kernel(1.0, t, &R, &T, &A, &E);
		bad += st_chk("omega=1 : A = 0 exactly", A, 0.0, 0.0);
		bad += st_chk("omega=1 : T = 2/(2+t)", T, 2.0 / (2.0 + t), 1e-14);
	}

	/* (4) 厚い群落のアルベド -> 半無限二流アルベド (1-u)/(1+u) */
	for (i = 1; i <= 4; i++) {
		const double om = 0.2 * i;
		const double u = sqrt(1.0 - om);
		canopy_kernel(om, 400.0, &R, &T, &A, &E);
		bad += st_chk("deep canopy : R = (1-u)/(1+u)", R, (1.0 - u) / (1.0 + u), 1e-12);
	}

	/* (5) 恒等式 A = 2(1-ω) t E */
	worst = 0.0;
	for (i = 0; i <= 40; i++) {
		for (j = 1; j <= 40; j++) {
			const double om = i / 40.0;
			const double t = j * 0.25;
			double d;
			canopy_kernel(om, t, &R, &T, &A, &E);
			d = fabs(A - (2.0 * (1.0 - om) * t * E));
			if (d > worst) worst = d;
		}
	}
	bad += st_chk("A = 2(1-omega) t E", worst, 0.0, 1e-14);

	/* (6) 層の積み上げ : 厚さ t の 2 枚 = 厚さ 2t の 1 枚 */
	{
		double R1, T1, A1, E1, R2, T2, A2, E2, g;
		canopy_kernel(0.6, 0.7, &R1, &T1, &A1, &E1);
		canopy_kernel(0.6, 1.4, &R2, &T2, &A2, &E2);
		g = 1.0 / (1.0 - (R1 * R1));
		bad += st_chk("adding : 2 x t == 1 x 2t (R)", R1 + (T1 * T1 * R1 * g), R2, 1e-12);
		bad += st_chk("adding : 2 x t == 1 x 2t (T)", T1 * T1 * g, T2, 1e-12);
	}

	printf("%s (%d failures)\n", bad ? "*** canopy self test FAILED" : "canopy self test passed", bad);
	return bad;
}

/* 葉の散乱係数 ω = ρ + τ (クランプはここ 1 箇所に集約する) */
double canopy_omega(const ppfd_t *p, int il)
{
	const pmat_t *m = &p->mat[p->canopy.imat];
	double omega = m->rho[il] + m->tau[il];

	if (omega < 0.0) omega = 0.0;
	if (omega > 1.0) omega = 1.0;
	return omega;
}

/*
群落の消散係数 [1/m]。
散乱モデルが有効なら素の遮断 k = G a (波長非依存)、
無効なら従来どおり Goudriaan 補正込みの G a √(1-ω)。
*/
double canopy_kext(const ppfd_t *p, int il)
{
	const canopy_t *c = &p->canopy;

	if (!c->on) return 0.0;
	if (c->scatter) return c->k0 * c->a;
	return c->k0 * sqrt(1.0 - canopy_omega(p, il)) * c->a;
}

/*
線分 a→b が層 m を通る長さ [m] を seg[0..nlayer-1] に加算する。
群落スラブの外にある部分は無視する。戻り値 = スラブ内の全長 [m]
(= canopy_path と一致する)。
*/
double canopy_segments(const ppfd_t *p, vec3_t a, vec3_t b, double *seg)
{
	const canopy_t *c = &p->canopy;
	const vec3_t d = v_sub(b, a);
	const double len = v_norm(d);
	const double dz = c->ztop - c->zbot;
	double total = 0.0;
	int    m;

	for (m = 0; m < c->nlayer; m++) seg[m] = 0.0;
	if (!c->on || (len < EPS)) return 0.0;

	if (fabs(d.z) < (EPS * (len + 1.0))) {
		/* ほぼ水平な線分 : 1 つの層に丸ごと入る */
		if ((a.z < c->zbot) || (a.z > c->ztop)) return 0.0;
		m = (int)(((c->ztop - a.z) / dz) * c->nlayer);
		if (m < 0) m = 0;
		if (m > c->nlayer - 1) m = c->nlayer - 1;
		seg[m] = len;
		return len;
	}

	/* 層 m は z ∈ [ztop - (m+1)h, ztop - m h] (m = 0 が最上層) */
	for (m = 0; m < c->nlayer; m++) {
		const double h = dz / c->nlayer;
		const double zhi = c->ztop - (m * h);
		const double zlo = c->ztop - ((m + 1) * h);
		double t0 = (zlo - a.z) / d.z;
		double t1 = (zhi - a.z) / d.z;
		double lo = (t0 < t1) ? t0 : t1;
		double hi = (t0 < t1) ? t1 : t0;
		if (lo < 0.0) lo = 0.0;
		if (hi > 1.0) hi = 1.0;
		if (hi <= lo) continue;
		seg[m] = (hi - lo) * len;
		total += seg[m];
	}
	return total;
}

/*
線分 a→b に沿って単位放射束が層ごとに落とす「遮断割合」を dep[] に加算する。
w は光線が運ぶ重み (放射束 [W] など)。

    dep[m] += w * exp(-k s_before) * (1 - exp(-k Δs_m))

k は波長によらないので、この係数は分光ループの外で 1 回だけ作れる。
Σ_m dep[m] + w exp(-k s_total) = w が望遠鏡和で厳密に成り立つ。
戻り値 = 全経路の透過率 exp(-k s_total)。
*/
double canopy_deposit(const ppfd_t *p, vec3_t a, vec3_t b, double w, double *dep, double *seg)
{
	const canopy_t *c = &p->canopy;
	const double k = c->k0 * c->a;
	double tr = 1.0;
	int    m;

	if (canopy_segments(p, a, b, seg) <= 0.0) return 1.0;

	/* a に近い層から順に処理する (光線の向きで層の順序が変わる) */
	if (b.z <= a.z) {
		for (m = 0; m < c->nlayer; m++) {
			const double e = exp(-k * seg[m]);
			dep[m] += w * tr * (1.0 - e);
			tr *= e;
		}
	}
	else {
		for (m = c->nlayer - 1; m >= 0; m--) {
			const double e = exp(-k * seg[m]);
			dep[m] += w * tr * (1.0 - e);
			tr *= e;
		}
	}
	return tr;
}

/* 群落の上面 (iz = 1) / 下面 (iz = 0) の矩形 */
static void canopy_plane(const ppfd_t *p, int istop, vec3_t *poly)
{
	const double z = istop ? p->canopy.ztop : p->canopy.zbot;

	poly[0] = v_make(0.0,    0.0,    z);
	poly[1] = v_make(p->Lx,  0.0,    z);
	poly[2] = v_make(p->Lx,  p->Ly,  z);
	poly[3] = v_make(0.0,    p->Ly,  z);
}

/*
散乱モデルの前処理 (setup_ff のあと)。

  cffup[i] / cffdn[i] : パッチ i から群落上面 / 下面への形態係数。
      上面より上のパッチだけが上面を、下面より下のパッチだけが下面を見る。
      群落の高さ範囲に重心があるパッチ (埋没パッチ) はどちらとも交換しない
      — 1 次元カラム近似では横方向の輸送が定義されないため。数はログに出す。

      面から出た放射束が過不足なくパッチへ渡るよう Σ_i A_i F_i = A_c へ
      正規化する。正規化前の比 (1 が理想) を診断として出す。求積誤差と
      埋没パッチのぶんだけ 1 から外れるので、近似の強さがそのまま見える。

  cdepw[i][m] : パッチ i の放射発散度 1 単位あたり、層 m へ預けられる
      放射束 [m²]。A_i Σ_j F_ij (i→j の経路が層 m で落とす割合)。
      素の遮断 k = G a は波長によらないので、この係数は 1 回作れば足りる。
*/
void canopy_setup(ppfd_t *p)
{
	const canopy_t *c = &p->canopy;
	const int n = p->npatch;
	const int nl = c->nlayer;
	vec3_t ptop[4], pbot[4];
	double sup = 0.0, sdn = 0.0;
	int    i, m;

	p->carea = p->Lx * p->Ly;
	p->cffup = (double *)xcalloc((size_t)n, sizeof(double));
	p->cffdn = (double *)xcalloc((size_t)n, sizeof(double));
	p->cdepw = (double *)xcalloc((size_t)n * nl, sizeof(double));
	p->cdep  = (double *)xcalloc((size_t)nl * p->nlam, sizeof(double));
	p->cfdn  = (double *)xcalloc((size_t)(nl + 1) * p->nlam, sizeof(double));
	p->cfup  = (double *)xcalloc((size_t)(nl + 1) * p->nlam, sizeof(double));
	p->cbtop = (double *)xcalloc((size_t)p->nlam, sizeof(double));
	p->cbbot = (double *)xcalloc((size_t)p->nlam, sizeof(double));
	p->cabs  = (double *)xcalloc((size_t)p->nlam, sizeof(double));

	canopy_plane(p, 1, ptop);
	canopy_plane(p, 0, pbot);

	p->cnburied = 0;
	for (i = 0; i < n; i++) {
		const patch_t *q = &p->patch[i];
		/* 上面は「群落より上」の求積点だけ、下面は「群落より下」だけが見る */
		p->cffup[i] = ff_patch_poly(q, ptop, 4, 4, ptop[0], v_make(0, 0, 1));
		p->cffdn[i] = ff_patch_poly(q, pbot, 4, 4, pbot[0], v_make(0, 0, -1));
		sup += q->area * p->cffup[i];
		sdn += q->area * p->cffdn[i];
		if ((p->cffup[i] <= 0.0) && (p->cffdn[i] <= 0.0)
		 && (q->c.z < c->ztop) && (q->c.z > c->zbot)) {
			p->cnburied++;
		}
	}

	/* 面から出た放射束を過不足なく配るための正規化 */
	p->cnorm_up = (p->carea > 0.0) ? (sup / p->carea) : 0.0;
	p->cnorm_dn = (p->carea > 0.0) ? (sdn / p->carea) : 0.0;
	if (sup > EPS) {
		for (i = 0; i < n; i++) p->cffup[i] *= p->carea / sup;
	}
	if (sdn > EPS) {
		for (i = 0; i < n; i++) p->cffdn[i] *= p->carea / sdn;
	}

	/* 壁パッチ間のビームが層へ預ける係数 */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
	for (i = 0; i < n; i++) {
		double *dep = &p->cdepw[(size_t)i * nl];
		double *seg = (double *)xmalloc((size_t)nl * sizeof(double));
		int     j;
		for (j = 0; j < n; j++) {
			const double f = p->ff[((size_t)i * n) + j];
			if (f <= 0.0) continue;
			canopy_deposit(p, p->patch[i].c, p->patch[j].c, f, dep, seg);
		}
		for (j = 0; j < nl; j++) dep[j] *= p->patch[i].area;
		free(seg);
	}
	(void)m;
}

/*
波長ビン il について、壁の放射発散度 B[] から群落の拡散場を解き、
上面 / 下面から出る流束と層界面の流束を更新する。
*aout に群落が吸収した放射束 [W] を返す。

湧き出し Q_m = ω (光源ビームの預け入れ + Σ_i cdepw[i][m] B_i) / A_c
*/
void canopy_field(ppfd_t *p, int il, const double *B, double *aout, double *sout)
{
	const canopy_t *c = &p->canopy;
	const int nl = c->nlayer;
	const double omega = canopy_omega(p, il);
	double *Q = (double *)xcalloc((size_t)nl, sizeof(double));
	double *fdn = &p->cfdn[(size_t)il * (nl + 1)];
	double *fup = &p->cfup[(size_t)il * (nl + 1)];
	double  etop = 0.0, ebot = 0.0, adiff = 0.0, dsum = 0.0;
	int     i, m;

	for (m = 0; m < nl; m++) {
		/* 総和の順序を固定する (スレッド数によらず同じ値にするため) */
		double d = p->cdep[((size_t)m * p->nlam) + il];
		for (i = 0; i < p->npatch; i++) {
			d += p->cdepw[((size_t)i * nl) + m] * B[(size_t)i * p->nlam + il];
		}
		dsum += d;
		Q[m] = (omega * d) / p->carea;
	}

	canopy_column(p, omega, Q, fdn, fup, &etop, &ebot, &adiff);

	p->cbtop[il] = etop;
	p->cbbot[il] = ebot;
	/* 遮断分の (1-ω) が直接吸収、拡散場の吸収が adiff (どちらも W) */
	*aout = ((1.0 - omega) * dsum) + (adiff * p->carea);
	*sout = (etop + ebot) * p->carea;

	free(Q);
}

/*
群落内部の高さ z における上向き受光面の拡散放射照度 [W/m²]。
層界面の下向き流束を線形内挿する (測定面は上向きの量子センサ)。
*/
double canopy_interior(const ppfd_t *p, int il, double z)
{
	const canopy_t *c = &p->canopy;
	const int nl = c->nlayer;
	const double *fdn = &p->cfdn[(size_t)il * (nl + 1)];
	double f;
	int    m;

	if ((z >= c->ztop) || (z <= c->zbot)) return 0.0;
	f = ((c->ztop - z) / (c->ztop - c->zbot)) * nl;
	m = (int)f;
	if (m < 0) m = 0;
	if (m > nl - 1) m = nl - 1;
	f -= m;
	return (fdn[m] * (1.0 - f)) + (fdn[m + 1] * f);
}

/*
1 本のカラムの層別二流を解く (外部からの拡散入射は無い : 外から来る光は
すべてビームとして dep に預けられている)。

  Q[m]      層 m の湧き出し [W/m²] (= ω × 預け入れ / 水平面積)
  *etop     上面 (ztop) から出る上向き拡散流束 [W/m²]
  *ebot     下面 (zbot) から出る下向き拡散流束 [W/m²]
  fdn/fup   層界面 (nlayer+1 点) の下向き / 上向き流束 [W/m²]
  *aabs     拡散場が層内で吸収した総量 [W/m²]

下から上への adding sweep で (R_b, S_b) を積み上げ、上から下へ流束を復元する。
R < 1 なので 1 - R R_b > 0 が保証される。
*/
void canopy_column(const ppfd_t *p, double omega, const double *Q,
                   double *fdn, double *fup, double *etop, double *ebot, double *aabs)
{
	const canopy_t *c = &p->canopy;
	const int n = c->nlayer;
	const double h = (c->ztop - c->zbot) / n;
	const double t = c->k0 * c->a * h;
	double R, T, A, E;
	double *Rb = (double *)xmalloc((size_t)(n + 1) * sizeof(double));
	double *Sb = (double *)xmalloc((size_t)(n + 1) * sizeof(double));
	int    m;

	canopy_kernel(omega, t, &R, &T, &A, &E);

	/* 下面より下は空 : 反射も湧き出しも無い */
	Rb[n] = 0.0;
	Sb[n] = 0.0;
	for (m = n - 1; m >= 0; m--) {
		const double g = 1.0 / (1.0 - (R * Rb[m + 1]));
		const double src = E * Q[m];
		Rb[m] = R + (T * T * Rb[m + 1] * g);
		Sb[m] = src + (T * g * ((Rb[m + 1] * ((R * Sb[m + 1]) + src)) + Sb[m + 1]));
	}

	/* 上面には外から拡散光は入らない (d_0 = 0) */
	fdn[0] = 0.0;
	fup[0] = Sb[0];
	for (m = 0; m < n; m++) {
		const double g = 1.0 / (1.0 - (R * Rb[m + 1]));
		const double src = E * Q[m];
		/* 層 m を下へ抜ける流束 : 下側からの戻りを含めて解いた形 */
		fdn[m + 1] = g * ((T * fdn[m]) + (R * Sb[m + 1]) + src);
		fup[m + 1] = (Rb[m + 1] * fdn[m + 1]) + Sb[m + 1];
	}
	/* 上向きの復元 (層 m の上界面) */
	for (m = n - 1; m >= 0; m--) {
		const double src = E * Q[m];
		fup[m] = (R * fdn[m]) + (T * fup[m + 1]) + src;
	}

	*etop = fup[0];
	*ebot = fdn[n];

	/*
	吸収 = 入れたもの - 出たもの。カラムの収支そのものなので、層カーネルが
	R+T+A=1 を満たす限りこれは厳密 (打ち消しではなく帳簿の残り)。
	ω = 1 は「吸収が無い」ことが構造的に決まっているので明示的に 0 にする。
	*/
	{
		double qsum = 0.0;
		for (m = 0; m < n; m++) qsum += Q[m];
		*aabs = (omega >= 1.0) ? 0.0 : (qsum - *etop - *ebot);
		if (*aabs < 0.0) *aabs = 0.0;
	}

	free(Rb);
	free(Sb);
}
