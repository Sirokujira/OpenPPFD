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

	/*
	(7) カラムの解き方 (adding sweep) : 対称な湧き出しなら上下の流出は等しい、
	完全反射の境界 (crbot = 1) は「鏡映して 2 倍にしたカラム」と等しい。
	吸収を残差で数えるので収支の判定では見えない (以前の Sb の余分な
	(1 + R Rb) はこれで見つかった)。
	*/
	{
		ppfd_t   q;
		cslab_t  sl[1];
		double   rt[1] = { 0.0 }, rb[1] = { 0.0 };
		double   Q[16], fdn[17], fup[17], et, eb, ab, et2, eb2, ab2;
		int      m;

		memset(&q, 0, sizeof(q));
		memset(sl, 0, sizeof(sl));
		sl[0].ztop = 1.0;
		sl[0].zbot = 0.0;
		sl[0].k0 = 0.5;
		sl[0].a = 3.0;
		q.slab = sl;
		q.nslab = 1;
		q.crtop = rt;
		q.crbot = rb;

		/* 対称な湧き出し (8 層) */
		q.canopy.nlayer = 8;
		for (m = 0; m < 8; m++) Q[m] = 1.0 + (0.3 * ((m < 4) ? m : (7 - m)));
		canopy_column(&q, 0, 0.6, Q, fdn, fup, &et, &eb, &ab);
		bad += st_chk("column : symmetric source -> top = bottom", et, eb, 1e-13);

		/* 完全反射の下端 (4 層) = 鏡映した 8 層のカラム (厚さも 2 倍) */
		q.canopy.nlayer = 4;
		sl[0].ztop = 0.5;
		rb[0] = 1.0;
		for (m = 0; m < 4; m++) Q[m] = 1.0 + (0.5 * m);
		canopy_column(&q, 0, 0.6, Q, fdn, fup, &et, &eb, &ab);
		q.canopy.nlayer = 8;
		sl[0].ztop = 1.0;
		rb[0] = 0.0;
		for (m = 0; m < 4; m++) Q[7 - m] = Q[m];
		canopy_column(&q, 0, 0.6, Q, fdn, fup, &et2, &eb2, &ab2);
		bad += st_chk("column : mirror boundary == doubled (top)", et, et2, 1e-13);
		bad += st_chk("column : mirror boundary == doubled (abs)", 2.0 * ab, ab2, 1e-12);
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
線分 a→b が層 g を通る長さ [m] を seg[0..cnl-1] に入れる。層の通し番号は
スラブ j (上から) の層 m (上から) で g = j * nlayer + m なので、g が増える
向きが z の下向き。散乱モデルでないときはスラブ 0 の層だけ。群落スラブの
外にある部分は無視する。戻り値 = スラブ内の全長 [m]。
*/
double canopy_segments(const ppfd_t *p, vec3_t a, vec3_t b, double *seg)
{
	const canopy_t *c = &p->canopy;
	const vec3_t d = v_sub(b, a);
	const double len = v_norm(d);
	const int ns = (c->scatter && (p->nslab > 0)) ? p->nslab : 1;
	const int nl = c->nlayer;
	double total = 0.0;
	int    m, j;

	for (m = 0; m < ns * nl; m++) seg[m] = 0.0;
	if (!c->on || (len < EPS)) return 0.0;

	for (j = 0; j < ns; j++) {
		const double ztop = (p->nslab > 0) ? p->slab[j].ztop : c->ztop;
		const double zbot = (p->nslab > 0) ? p->slab[j].zbot : c->zbot;
		const double dz = ztop - zbot;
		double *sj = &seg[(size_t)j * nl];

		if (fabs(d.z) < (EPS * (len + 1.0))) {
			/* ほぼ水平な線分 : 1 つの層に丸ごと入る */
			if ((a.z < zbot) || (a.z > ztop)) continue;
			m = (int)(((ztop - a.z) / dz) * nl);
			if (m < 0) m = 0;
			if (m > nl - 1) m = nl - 1;
			sj[m] = len;
			return len;
		}

		/* 層 m は z ∈ [ztop - (m+1)h, ztop - m h] (m = 0 が最上層) */
		for (m = 0; m < nl; m++) {
			const double h = dz / nl;
			const double zhi = ztop - (m * h);
			const double zlo = ztop - ((m + 1) * h);
			double t0 = (zlo - a.z) / d.z;
			double t1 = (zhi - a.z) / d.z;
			double lo = (t0 < t1) ? t0 : t1;
			double hi = (t0 < t1) ? t1 : t0;
			if (lo < 0.0) lo = 0.0;
			if (hi > 1.0) hi = 1.0;
			if (hi <= lo) continue;
			sj[m] = (hi - lo) * len;
			total += sj[m];
		}
	}
	return total;
}

/*
線分 a→b に沿って単位放射束が層ごとに落とす「遮断割合」を dep[] に加算する。
w は光線が運ぶ重み (放射束 [W] など)。

    dep[g] += w * exp(-Σ_{前の層} k s) * (1 - exp(-k_g Δs_g))

k_g = G_j a_j (層 g が属するスラブ j の素の遮断) は波長によらないので、
この係数は分光ループの外で 1 回だけ作れる。
Σ_g dep[g] + w exp(-Σ k s) = w が望遠鏡和で厳密に成り立つ。
光線は z について単調なので、通る順は g の昇順 (下向き) か降順 (上向き)。
戻り値 = 全経路の透過率。
*/
static double deposit_straight(const ppfd_t *p, vec3_t a, vec3_t b, double w, double *dep, double *seg)
{
	const canopy_t *c = &p->canopy;
	const int ns = (c->scatter && (p->nslab > 0)) ? p->nslab : 1;
	const int nl = c->nlayer;
	double tr = 1.0;
	int    g;

	if (canopy_segments(p, a, b, seg) <= 0.0) return 1.0;

	/* a に近い層から順に処理する (光線の向きで層の順序が変わる) */
	if (b.z <= a.z) {
		for (g = 0; g < ns * nl; g++) {
			const cslab_t *sl = &p->slab[g / nl];
			const double k = (p->nslab > 0) ? (sl->k0 * sl->a) : (c->k0 * c->a);
			const double e = exp(-k * seg[g]);
			dep[g] += w * tr * (1.0 - e);
			tr *= e;
		}
	}
	else {
		for (g = (ns * nl) - 1; g >= 0; g--) {
			const cslab_t *sl = &p->slab[g / nl];
			const double k = (p->nslab > 0) ? (sl->k0 * sl->a) : (c->k0 * c->a);
			const double e = exp(-k * seg[g]);
			dep[g] += w * tr * (1.0 - e);
			tr *= e;
		}
	}
	return tr;
}

/*
床・天井が鏡面のとき、箱の外へ出る線分は鏡像への直線 (展開した経路)。
実経路は z = k Lz の平面で折り返すたびに層を通る順序が入れ替わる (下りの
あとに上り) ので、線分をその平面で区切り、区間ごとに z を三角波 fold_z で
実チャンバへ戻してから、経路に沿った順に預ける (透過率は区間をまたいで
掛け継ぐ)。鏡映は弧長を保ち、各区間の中で実経路の z は単調なので、これで
展開しない実経路を 1 本ずつたどったのと同じ預け入れになる。
床・天井が鏡面でない、または線分が箱の中に収まるなら従来どおり 1 区間。
*/
double canopy_deposit(const ppfd_t *p, vec3_t a, vec3_t b, double w, double *dep, double *seg)
{
	const int zmir = (p->mat[p->wallmat[4]].rhos > 0.0) || (p->mat[p->wallmat[5]].rhos > 0.0);
	const vec3_t d = v_sub(b, a);
	double tcut[66];
	double tr = 1.0;
	int    nc = 0, k, k0, k1, i;

	if (!zmir || ((a.z >= 0.0) && (a.z <= p->Lz) && (b.z >= 0.0) && (b.z <= p->Lz))) {
		return deposit_straight(p, a, b, w, dep, seg);
	}

	/* 折り返しの平面 z = k Lz で区切る (区間の端の t を昇順に並べる) */
	tcut[nc++] = 0.0;
	if (fabs(d.z) > 0.0) {
		const double zlo = (a.z < b.z) ? a.z : b.z;
		const double zhi = (a.z < b.z) ? b.z : a.z;
		k0 = (int)floor(zlo / p->Lz) + 1;
		k1 = (int)ceil(zhi / p->Lz) - 1;
		for (k = k0; (k <= k1) && (nc < 65); k++) {
			const double t = ((k * p->Lz) - a.z) / d.z;
			if ((t > 0.0) && (t < 1.0)) tcut[nc++] = t;
		}
		/* d.z < 0 なら t は降順に積まれているので並べ直す (挿入ソート) */
		for (i = 1; i < nc; i++) {
			const double t = tcut[i];
			int j = i - 1;
			while ((j >= 0) && (tcut[j] > t)) { tcut[j + 1] = tcut[j]; j--; }
			tcut[j + 1] = t;
		}
	}
	tcut[nc++] = 1.0;

	for (i = 0; i + 1 < nc; i++) {
		const vec3_t u = v_add(a, v_scale(d, tcut[i]));
		const vec3_t v = v_add(a, v_scale(d, tcut[i + 1]));
		const vec3_t ur = v_make(u.x, u.y, fold_z(p, u.z));
		const vec3_t vr = v_make(v.x, v.y, fold_z(p, v.z));
		tr *= deposit_straight(p, ur, vr, w * tr, dep, seg);
	}
	return tr;
}

/*
2 点則で持った経路長 pl = (s1, s2, w1) の透過率 (canopy_plen を参照)。
s2 = 0 なら群落を通らない (減衰なし)。
*/
double plen_trans(const float *pl, double k)
{
	const double w1 = pl[2];

	if (pl[1] <= 0.0f) return 1.0;
	if (w1 >= 1.0) return exp(-k * pl[0]);
	return (w1 * exp(-k * pl[0])) + ((1.0 - w1) * exp(-k * pl[1]));
}

/* スラブ j の上面 (istop = 1) / 下面 (istop = 0) の矩形 */
static void canopy_plane(const ppfd_t *p, int j, int istop, vec3_t *poly)
{
	const double z = istop ? p->slab[j].ztop : p->slab[j].zbot;

	poly[0] = v_make(0.0,    0.0,    z);
	poly[1] = v_make(p->Lx,  0.0,    z);
	poly[2] = v_make(p->Lx,  p->Ly,  z);
	poly[3] = v_make(0.0,    p->Ly,  z);
}

/*
群落の上面 / 下面 (またはその鏡像) の多角形 poly4 への、パッチ q からの
形態係数。半空間 gn·(x - gp) > 0 にある求積点だけが見る (gn が零ベクトル
なら判定しない)。
tn は面が鏡像のときの展開指数 (実在なら NULL)。

遮蔽物が無ければ従来どおり多角形 1 枚の閉形式。遮蔽物があるときは、
面を patchdiv の格子の小矩形に割って、小矩形ごとの閉形式の形態係数に
可視率 (vis_patch_patch) を掛けて足す。面は断面全体の大きな多角形なので、
ここで可視率を見ないと棚板の向こうのパッチにも群落の散乱光が届いてしまう
(断面いっぱいの棚板の上の測定面に光が漏れていた)。
*/
static double plane_ff(const ppfd_t *p, const patch_t *q, const vec3_t *poly4, vec3_t gp,
                       vec3_t gn, const int *tn)
{
	const int nu = p->ndivu, nv = p->ndivv;
	const vec3_t eu = v_sub(poly4[1], poly4[0]);
	const vec3_t ev = v_sub(poly4[3], poly4[0]);
	double f = 0.0;
	int    iu, iv, v;

	if (!occ_any(p)) return ff_patch_poly(q, poly4, 4, 4, gp, gn);

	for (iv = 0; iv < nv; iv++) {
		for (iu = 0; iu < nu; iu++) {
			patch_t sub;
			double  fs;
			sub.p[0] = v_add(poly4[0], v_add(v_scale(eu, (double)iu / nu), v_scale(ev, (double)iv / nv)));
			sub.p[1] = v_add(sub.p[0], v_scale(eu, 1.0 / nu));
			sub.p[2] = v_add(sub.p[1], v_scale(ev, 1.0 / nv));
			sub.p[3] = v_add(sub.p[0], v_scale(ev, 1.0 / nv));
			fs = ff_patch_poly(q, sub.p, 4, 4, gp, gn);
			if (fs <= 0.0) continue;
			sub.n = gn;
			sub.c = v_make(0.0, 0.0, 0.0);
			for (v = 0; v < 4; v++) sub.c = v_add(sub.c, v_scale(sub.p[v], 0.25));
			sub.area = 0.0;
			sub.imat = 0;
			sub.iface = -1;
			f += fs * vis_patch_patch(p, q, &sub, tn, NULL);
		}
	}
	return f;
}

/*
測定面の点 x (法線 nz) から群落の上面 / 下面 (またはその鏡像) への形態係数。
遮蔽物があれば plane_ff と同じく小矩形ごとに可視率を掛ける。
*/
double canopy_plane_ff_point(const ppfd_t *p, vec3_t x, vec3_t nz, const vec3_t *poly4, const int *tn)
{
	const int nu = p->ndivu, nv = p->ndivv;
	const vec3_t eu = v_sub(poly4[1], poly4[0]);
	const vec3_t ev = v_sub(poly4[3], poly4[0]);
	double f = 0.0;
	int    iu, iv, v;

	if (!occ_any(p)) return ff_point_poly(x, nz, poly4, 4);

	for (iv = 0; iv < nv; iv++) {
		for (iu = 0; iu < nu; iu++) {
			patch_t sub;
			double  fs;
			sub.p[0] = v_add(poly4[0], v_add(v_scale(eu, (double)iu / nu), v_scale(ev, (double)iv / nv)));
			sub.p[1] = v_add(sub.p[0], v_scale(eu, 1.0 / nu));
			sub.p[2] = v_add(sub.p[1], v_scale(ev, 1.0 / nv));
			sub.p[3] = v_add(sub.p[0], v_scale(ev, 1.0 / nv));
			fs = ff_point_poly(x, nz, sub.p, 4);
			if (fs <= 0.0) continue;
			sub.n = v_make(0.0, 0.0, 0.0);
			sub.c = v_make(0.0, 0.0, 0.0);
			for (v = 0; v < 4; v++) sub.c = v_add(sub.c, v_scale(sub.p[v], 0.25));
			sub.area = 0.0;
			sub.imat = 0;
			sub.iface = -1;
			f += fs * vis_point_patch(p, x, &sub, tn, NULL);
		}
	}
	return f;
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
      i→j の経路は setup_ff が透過率を作ったのと同じ標本 (canopy_pair_samples)
      で運ぶので、対ごとに「透過 + 各層への預け入れ = 1」が成り立つ。
      素の遮断 k = G a は波長によらないので、この係数は 1 回作れば足りる。

■ 側壁の鏡面との併用
側壁 (x/y 面) の鏡映は z を保つので、鏡像への直線の z プロファイルは
実経路のものと同じで、層への預け入れはそのまま展開した直線で測れる。
帳簿は「光線は止まった所で預ける」に揃える :

  - 鏡面のパッチ j に届いた光のうち、止まる (吸収または拡散反射) のは
    (1 - ρs_j) の分。ρs_j の分は j の鏡像を通る拡張形態係数の経路が
    展開した全経路で運び、預け入れもそちらで数える。そこで i → j の
    預け入れは (1 - ρs_j) 倍する (直接光の光線も同じ。direct.c)。
    鏡面が無ければ係数は厳密に 1。
  - 群落の上下面から出た放射も鏡面で折り返す : 面の鏡像への形態係数を
    Π ρs の重みで足す。面から出た放射束が過不足なく止まるよう
    Σ_i A_i (1 - ρs_i) Σ_t w_t F(i → T_t(面)) = A_c へ正規化する。
    鏡面が無ければ従来の Σ_i A_i F_i = A_c そのもの。

完全鏡面 (ρs = 1) ならこれは「鏡映した幾何を実在させた二重箱」の帳簿と
項ごとに一致する (鏡面のパッチは何も預けず、鏡像経路が全経路を運ぶ)。
床・天井の鏡面は層の順序を入れ替える (預け入れの順序が折り返す) ので
併用させない (入力段で弾く)。
*/
void canopy_setup(ppfd_t *p)
{
	const canopy_t *c = &p->canopy;
	const int n = p->npatch;
	const int nsl = p->nslab;
	const int nl = nsl * c->nlayer;          /* 全スラブの層の通し数 */
	strans_t *tr;
	int    nt;
	int    i, m, js;

	p->cnl   = nl;
	p->carea = p->Lx * p->Ly;
	p->cffup = (double *)xcalloc((size_t)nsl * n, sizeof(double));
	p->cffdn = (double *)xcalloc((size_t)nsl * n, sizeof(double));
	p->cdepw = (double *)xcalloc((size_t)n * nl, sizeof(double));
	p->cdep  = (double *)xcalloc((size_t)nl * p->nlam, sizeof(double));
	p->cfdn  = (double *)xcalloc((size_t)nsl * (c->nlayer + 1) * p->nlam, sizeof(double));
	p->cfup  = (double *)xcalloc((size_t)nsl * (c->nlayer + 1) * p->nlam, sizeof(double));
	p->cbtop = (double *)xcalloc((size_t)nsl * p->nlam, sizeof(double));
	p->cbbot = (double *)xcalloc((size_t)nsl * p->nlam, sizeof(double));
	p->cabs  = (double *)xcalloc((size_t)p->nlam, sizeof(double));
	p->cnorm_up = (double *)xcalloc((size_t)nsl, sizeof(double));
	p->cnorm_dn = (double *)xcalloc((size_t)nsl, sizeof(double));
	p->crtop = (double *)xcalloc((size_t)nsl, sizeof(double));
	p->crbot = (double *)xcalloc((size_t)nsl, sizeof(double));

	tr = (strans_t *)xmalloc((size_t)MAXSPECT * sizeof(strans_t));
	nt = spec_transforms(p, tr, MAXSPECT);

	p->cnburied = 0;
	for (js = 0; js < nsl; js++) {
		const cslab_t *sl = &p->slab[js];
		double *fu = &p->cffup[(size_t)js * n];
		double *fd = &p->cffdn[(size_t)js * n];
		const vec3_t gup = v_make(0.0, 0.0, 1.0), gdn = v_make(0.0, 0.0, -1.0);
		const vec3_t zero = v_make(0.0, 0.0, 0.0);
		const double tol = 1e-9 * p->Lz;
		/* 面が床・天井に接している (面がそのまま壁に乗る) */
		const int ontop = (sl->ztop >= p->Lz - tol), onbot = (sl->zbot <= tol);
		vec3_t ptop[4], pbot[4];
		double sup = 0.0, sdn = 0.0, rtop = 0.0, rbot = 0.0;

		canopy_plane(p, js, 1, ptop);
		canopy_plane(p, js, 0, pbot);
		for (i = 0; i < n; i++) {
			const patch_t *q = &p->patch[i];
			/* 止まる割合 (鏡面のパッチでは正反射の分は先へ運ばれる) */
			const double stop = 1.0 - p->mat[q->imat].rhos;
			int t, v;
			/* 上面は「群落より上」の求積点だけ、下面は「群落より下」だけが見る */
			fu[i] = plane_ff(p, q, ptop, ptop[0], gup, NULL);
			fd[i] = plane_ff(p, q, pbot, pbot[0], gdn, NULL);
			/*
			鏡面で折り返して見える上下面 (面の鏡像)。側壁の鏡映 (n2 = 0) は
			z を保つ。上面は天井で 1 回だけ (n2 = +1)、下面は床で 1 回だけ
			(n2 = -1) 折り返した経路まで : それより先は必ず自分の面を通って
			群落へ戻る (それは境界の反射率 crtop / crbot が運ぶ)。
			*/
			for (t = 0; t < nt; t++) {
				vec3_t it[4], ib[4];
				for (v = 0; v < 4; v++) {
					it[v] = spec_apply(p, &tr[t], ptop[v]);
					ib[v] = spec_apply(p, &tr[t], pbot[v]);
				}
				if ((tr[t].n[2] == 0) || (tr[t].n[2] == 1)) {
					fu[i] += tr[t].w * plane_ff(p, q, it, ptop[0], gup, tr[t].n);
				}
				if ((tr[t].n[2] == 0) || (tr[t].n[2] == -1)) {
					fd[i] += tr[t].w * plane_ff(p, q, ib, pbot[0], gdn, tr[t].n);
				}
			}
			/* 面が壁に乗っているとき、その壁のパッチは面の放射を全部受ける */
			if (ontop && (q->iface == 5)) fu[i] = 1.0;
			if (onbot && (q->iface == 4)) fd[i] = 1.0;
			sup += (q->area * stop) * fu[i];
			sdn += (q->area * stop) * fd[i];
			if ((fu[i] <= 0.0) && (fd[i] <= 0.0)
			 && (q->c.z < sl->ztop) && (q->c.z > sl->zbot)) {
				p->cnburied++;
			}
		}

		/*
		境界の反射率 : 面から出た放射のうち、床・天井の鏡面で折り返して
		そのまま同じ面へ戻る割合。戻った光は二流の境界からの拡散入射として
		カラムへ入る (canopy_column)。面を小矩形の受光面に割って、面の鏡像
		(n2 = +1 / -1、側壁の鏡映との合成も) への形態係数を足す。面が鏡面の
		壁に乗っているなら、出た光はただちにその壁で正反射して戻る (ρs)。
		*/
		{
			const int nu = p->ndivu, nv = p->ndivv;
			int iu, iv, t, v, side;
			for (side = 0; side < 2; side++) {
				const vec3_t *pl = (side == 0) ? ptop : pbot;
				const int want = (side == 0) ? 1 : -1;
				const int onw = (side == 0) ? ontop : onbot;
				const double rsw = p->mat[p->wallmat[(side == 0) ? 5 : 4]].rhos;
				double r = 0.0;
				if (onw) {
					r = rsw;
				}
				else if (rsw > 0.0) {
					const vec3_t eu = v_sub(pl[1], pl[0]);
					const vec3_t ev = v_sub(pl[3], pl[0]);
					for (iv = 0; iv < nv; iv++) {
						for (iu = 0; iu < nu; iu++) {
							patch_t sub;
							sub.p[0] = v_add(pl[0], v_add(v_scale(eu, (double)iu / nu), v_scale(ev, (double)iv / nv)));
							sub.p[1] = v_add(sub.p[0], v_scale(eu, 1.0 / nu));
							sub.p[2] = v_add(sub.p[1], v_scale(ev, 1.0 / nv));
							sub.p[3] = v_add(sub.p[0], v_scale(ev, 1.0 / nv));
							if (side == 1) {
								/* 下面の受光面は下向き (頂点の順を逆にして法線を -z に) */
								const vec3_t tmp = sub.p[1];
								sub.p[1] = sub.p[3];
								sub.p[3] = tmp;
							}
							sub.n = (side == 0) ? gup : gdn;
							sub.c = v_make(0.0, 0.0, 0.0);
							for (v = 0; v < 4; v++) sub.c = v_add(sub.c, v_scale(sub.p[v], 0.25));
							sub.area = p->carea / (nu * nv);
							sub.imat = 0;
							sub.iface = -1;
							for (t = 0; t < nt; t++) {
								vec3_t im[4];
								if (tr[t].n[2] != want) continue;
								for (v = 0; v < 4; v++) im[v] = spec_apply(p, &tr[t], pl[v]);
								r += (sub.area / p->carea) * tr[t].w * plane_ff(p, &sub, im, zero, zero, tr[t].n);
							}
						}
					}
				}
				if (side == 0) rtop = r;
				else           rbot = r;
			}
		}

		/*
		面から出た放射束を過不足なく配るための正規化 : 止まる分 (パッチ) と
		戻る分 (境界の反射) の和が A_c になるようにそろえて縮める。鏡面の
		床・天井が無ければ戻る分は 0 で、従来の Σ_i A_i F_i = A_c そのもの。
		*/
		sup += p->carea * rtop;
		sdn += p->carea * rbot;
		p->cnorm_up[js] = (p->carea > 0.0) ? (sup / p->carea) : 0.0;
		p->cnorm_dn[js] = (p->carea > 0.0) ? (sdn / p->carea) : 0.0;
		if (sup > EPS) {
			for (i = 0; i < n; i++) fu[i] *= p->carea / sup;
			rtop *= p->carea / sup;
		}
		if (sdn > EPS) {
			for (i = 0; i < n; i++) fd[i] *= p->carea / sdn;
			rbot *= p->carea / sdn;
		}
		p->crtop[js] = rtop;
		p->crbot[js] = rbot;
	}

	/* 壁パッチ間のビームが層へ預ける係数 */
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
	for (i = 0; i < n; i++) {
		double *dep = &p->cdepw[(size_t)i * nl];
		double *seg = (double *)xmalloc((size_t)nl * sizeof(double));
		vec3_t *xa = (vec3_t *)xmalloc((size_t)CANOPY_MAXSAMP * sizeof(vec3_t));
		vec3_t *yb = (vec3_t *)xmalloc((size_t)CANOPY_MAXSAMP * sizeof(vec3_t));
		double *W  = (double *)xmalloc((size_t)CANOPY_MAXSAMP * sizeof(double));
		int     j, t;
		for (j = 0; j < n; j++) {
			/* 止まる分だけ預ける (鏡面のパッチ j では ρs の分を鏡像経路が運ぶ) */
			const double f = p->ff[((size_t)i * n) + j] * (1.0 - p->mat[p->patch[j].imat].rhos);
			int ns, a;
			if (f <= 0.0) continue;
			/* setup_ff が透過率を作ったのと同じ標本で預ける (対ごとに透過 + 預け入れ = 1) */
			ns = canopy_pair_samples(p, &p->patch[i], p->patch[j].p, p->patch[j].n, NULL, xa, yb, W);
			for (a = 0; a < ns; a++) canopy_deposit(p, xa[a], yb[a], f * W[a], dep, seg);
		}
		/* 鏡像経路 : 展開した直線のまま預ける (側壁の鏡映は z を保つ) */
		for (t = 0; t < p->nst; t++) {
			for (j = 0; j < n; j++) {
				const patch_t *qj = &p->patch[j];
				const double f = p->ffs[((size_t)t * n * n) + ((size_t)i * n) + j]
				               * (1.0 - p->mat[qj->imat].rhos);
				vec3_t q4[4];
				int ns, a, v;
				if (f <= 0.0) continue;
				for (v = 0; v < 4; v++) q4[v] = spec_apply(p, &tr[t], qj->p[v]);
				ns = canopy_pair_samples(p, &p->patch[i], q4, spec_apply_dir(&tr[t], qj->n),
				                         tr[t].n, xa, yb, W);
				for (a = 0; a < ns; a++) canopy_deposit(p, xa[a], yb[a], f * W[a], dep, seg);
			}
		}
		for (j = 0; j < nl; j++) dep[j] *= p->patch[i].area;
		free(seg);
		free(xa);
		free(yb);
		free(W);
	}
	free(tr);
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
	double  atot = 0.0, stot = 0.0;
	int     i, m, js;

	/* スラブごとに独立なカラム (棚板で仕切られているので互いに見えない) */
	for (js = 0; js < p->nslab; js++) {
		double *fdn = &p->cfdn[(((size_t)js * p->nlam) + il) * (nl + 1)];
		double *fup = &p->cfup[(((size_t)js * p->nlam) + il) * (nl + 1)];
		double  etop = 0.0, ebot = 0.0, adiff = 0.0, dsum = 0.0;

		for (m = 0; m < nl; m++) {
			/* 総和の順序を固定する (スレッド数によらず同じ値にするため) */
			const int g = (js * nl) + m;
			double d = p->cdep[((size_t)g * p->nlam) + il];
			for (i = 0; i < p->npatch; i++) {
				d += p->cdepw[((size_t)i * p->cnl) + g] * B[(size_t)i * p->nlam + il];
			}
			dsum += d;
			Q[m] = (omega * d) / p->carea;
		}

		canopy_column(p, js, omega, Q, fdn, fup, &etop, &ebot, &adiff);

		p->cbtop[((size_t)js * p->nlam) + il] = etop;
		p->cbbot[((size_t)js * p->nlam) + il] = ebot;
		/* 遮断分の (1-ω) が直接吸収、拡散場の吸収が adiff (どちらも W) */
		atot += ((1.0 - omega) * dsum) + (adiff * p->carea);
		/* 返した放射束は正味 (境界で同じ面へ戻った分は返していない) */
		stot += ((etop - fdn[0]) + (ebot - fup[nl])) * p->carea;
	}
	*aout = atot;
	*sout = stot;

	free(Q);
}

/*
群落内部の高さ z で、単位法線の z 成分が nz の受光面が受ける拡散放射照度
[W/m²]。層界面の流束を線形内挿する。

二流近似では各半球の放射輝度が一様 (下向き L↓ = F↓/π、上向き L↑ = F↑/π)
なので、天頂から角 β だけ傾いた面 (cos β = nz) が見る上半球の割合は
(1 + cos β)/2、下半球は (1 - cos β)/2 になる :

    E = F↓ (1 + nz)/2 + F↑ (1 - nz)/2

上向き (nz = 1) なら F↓、下向きなら F↑、鉛直なら両者の平均。nz = 1 では
係数が厳密に 1 と 0 なので、従来の上向きの値と 1 ビットも変わらない。
*/
double canopy_interior(const ppfd_t *p, int js, int il, double z, double nz)
{
	const canopy_t *c = &p->canopy;
	const cslab_t *sl = &p->slab[js];
	const int nl = c->nlayer;
	const double *fdn = &p->cfdn[(((size_t)js * p->nlam) + il) * (nl + 1)];
	const double *fup = &p->cfup[(((size_t)js * p->nlam) + il) * (nl + 1)];
	double f, dn, up;
	int    m;

	if ((z >= sl->ztop) || (z <= sl->zbot)) return 0.0;
	f = ((sl->ztop - z) / (sl->ztop - sl->zbot)) * nl;
	m = (int)f;
	if (m < 0) m = 0;
	if (m > nl - 1) m = nl - 1;
	f -= m;
	dn = (fdn[m] * (1.0 - f)) + (fdn[m + 1] * f);
	up = (fup[m] * (1.0 - f)) + (fup[m + 1] * f);
	return (dn * (0.5 * (1.0 + nz))) + (up * (0.5 * (1.0 - nz)));
}

/*
スラブ js の 1 本のカラムの層別二流を解く (外部からの拡散入射は無い : 外から来る光は
すべてビームとして dep に預けられている)。

  Q[m]      層 m の湧き出し [W/m²] (= ω × 預け入れ / 水平面積)
  *etop     上面 (ztop) から出る上向き拡散流束 [W/m²]
  *ebot     下面 (zbot) から出る下向き拡散流束 [W/m²]
  fdn/fup   層界面 (nlayer+1 点) の下向き / 上向き流束 [W/m²]
  *aabs     拡散場が層内で吸収した総量 [W/m²]

下から上への adding sweep で (R_b, S_b) を積み上げ、上から下へ流束を復元する。
R < 1 なので 1 - R R_b > 0 が保証される。
*/
void canopy_column(const ppfd_t *p, int js, double omega, const double *Q,
                   double *fdn, double *fup, double *etop, double *ebot, double *aabs)
{
	const canopy_t *c = &p->canopy;
	const cslab_t *sl = &p->slab[js];
	const int n = c->nlayer;
	const double h = (sl->ztop - sl->zbot) / n;
	const double t = sl->k0 * sl->a * h;
	double R, T, A, E;
	double *Rb = (double *)xmalloc((size_t)(n + 1) * sizeof(double));
	double *Sb = (double *)xmalloc((size_t)(n + 1) * sizeof(double));
	int    m;

	canopy_kernel(omega, t, &R, &T, &A, &E);

	/* 下面より下 : 鏡面の床で戻る分だけ反射する (無ければ 0) */
	Rb[n] = p->crbot[js];
	Sb[n] = 0.0;
	for (m = n - 1; m >= 0; m--) {
		const double g = 1.0 / (1.0 - (R * Rb[m + 1]));
		const double src = E * Q[m];
		Rb[m] = R + (T * T * Rb[m + 1] * g);
		/*
		上界面の上向き流束 (層 m の上から入射が無いとき) :
		  u_b = Rb d + Sb,  d = R u_b + src,  u_t = T u_b + src
		  -> u_t = src + T g (Rb src + Sb)
		(以前は Sb に余分な (1 + R Rb) が掛かっていた。吸収を残差で数えるので
		収支には出ず、対称な湧き出しで上下の流出が割れることで見つかった)
		*/
		Sb[m] = src + (T * g * ((Rb[m + 1] * src) + Sb[m + 1]));
	}

	/*
	上面 : 鏡面の天井で同じ面へ戻る分 crtop だけ拡散入射がある
	(d_0 = crtop u_0, u_0 = Rb_0 d_0 + Sb_0)。戻りが無ければ d_0 = 0。
	*/
	{
		const double rt = p->crtop[js];
		const double u0 = Sb[0] / (1.0 - (Rb[0] * rt));
		fdn[0] = rt * u0;
		fup[0] = Sb[0] + (Rb[0] * fdn[0]);
	}
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
		/* 出たもの = 正味 (境界で戻った分は差し引く) */
		*aabs = (omega >= 1.0) ? 0.0 : (qsum - (*etop - fdn[0]) - (*ebot - fup[n]));
		if (*aabs < 0.0) *aabs = 0.0;
	}

	free(Rb);
	free(Sb);
}
