/*
formfactor.c

微小面 -> 多角形の形態係数 (閉形式) と、それを受光面上で数値積分した
パッチ間形態係数。

■ 微小面 dA (位置 p、単位法線 n) から平面多角形 P への形態係数
    F = (1/2π) | Σ_i (n · û_i) θ_i |
        û_i = (R_i × R_{i+1}) / |R_i × R_{i+1}|
        θ_i = ∠(R_i, R_{i+1})           R_i = P_i - p
これは多角形が dA の上半空間に完全に含まれるときの厳密解なので、
先に多角形を平面 {x : n·(x-p) ≥ 0} で Sutherland-Hodgman クリップする。

■ パッチ間
    F_ij = (1/A_i) ∫_{A_i} F_{dA→P_j} dA
を受光パッチ上の Gauss-Legendre 求積で評価する (放射側は閉形式なので
誤差は受光側の求積のみ)。近接パッチは被積分関数の変化が急なので
複合求積の分割数を上げる。

同一平面のパッチどうしは互いに見えないので明示的に 0 とする。

■ 遮蔽
遮蔽物 (occluder キー) が無ければキャビティは凸なので遮蔽判定は要らず、
F_ij は上の求積の精度で厳密になる。遮蔽物があるときは可視率

    V_ij = (両パッチ上の標本点対のうち遮られない割合)

を掛ける。まず 2×2 × 2×2 = 16 本の線分で判定し、全部同じ (完全に見える
/ 完全に隠れる) ならそこで確定する — 棚板がキャビティを仕切るような
配置ではすべての対がこれで確定し、遮蔽があっても F は厳密なままになる。
判定が割れた対だけ 4×4 × 4×4 = 256 本へ細分する。この細分された対の
数はログに出す (標本化に頼った = 厳密でない箇所がどれだけあったか)。
*/

#include "ppfd.h"

/* Gauss-Legendre 3 点 (区間 [0,1] へ写像済み、重みの総和 = 1) */
static const double gx3[3] = {
	0.5 * (1.0 - 0.7745966692414834), 0.5, 0.5 * (1.0 + 0.7745966692414834)
};
static const double gw3[3] = {
	0.5 * 0.5555555555555556, 0.5 * 0.8888888888888888, 0.5 * 0.5555555555555556
};

/* 平面 {x : n·(x-p) >= 0} でのクリップ。戻り値 = 出力頂点数 */
static int clip_halfspace(const vec3_t *in, int nin, vec3_t p, vec3_t n, vec3_t *out, int maxout)
{
	int    i, nout = 0;
	double dprev;
	vec3_t vprev;

	if (nin < 3) return 0;

	vprev = in[nin - 1];
	dprev = v_dot(n, v_sub(vprev, p));

	for (i = 0; i < nin; i++) {
		const vec3_t vcur = in[i];
		const double dcur = v_dot(n, v_sub(vcur, p));

		if (dcur >= 0.0) {
			if (dprev < 0.0) {
				const double t = dprev / (dprev - dcur);
				if (nout < maxout) out[nout++] = v_add(vprev, v_scale(v_sub(vcur, vprev), t));
			}
			if (nout < maxout) out[nout++] = vcur;
		}
		else if (dprev >= 0.0) {
			const double t = dprev / (dprev - dcur);
			if (nout < maxout) out[nout++] = v_add(vprev, v_scale(v_sub(vcur, vprev), t));
		}
		vprev = vcur;
		dprev = dcur;
	}
	return nout;
}

double ff_point_poly(vec3_t p, vec3_t n, const vec3_t *poly, int nv)
{
	vec3_t buf[12];
	int    m, i;
	double sum = 0.0;
	double dmax = 0.0, rmax = 0.0;

	m = clip_halfspace(poly, nv, p, n, buf, 12);
	if (m < 3) return 0.0;

	/*
	クリップ後の多角形が受光点の接平面上に退化していないか。
	多角形が平面 n·(x-p) = 0 に完全に乗っているとき立体角は 0 だが、
	クリップは平面上の頂点 (d = 0) を残すので、線分に潰れた多角形が
	そのまま下の総和に流れ込む。頂点がほぼ同一直線上に並ぶと
	外積の丸め誤差がそのまま角度になるので、値は 0 ではなく数値ごみに
	なる (FMA を使うと顕在化し、棚板に接する壁パッチとの間に 1e-9
	程度の形態係数が生えて、仕切ったはずの部屋に光が漏れた)。
	*/
	for (i = 0; i < m; i++) {
		const vec3_t r = v_sub(buf[i], p);
		const double d = v_dot(n, r);
		const double a = v_norm(r);
		if (d > dmax) dmax = d;
		if (a > rmax) rmax = a;
	}
	if (dmax <= (1e-12 * rmax)) return 0.0;

	for (i = 0; i < m; i++) {
		const vec3_t r0 = v_sub(buf[i], p);
		const vec3_t r1 = v_sub(buf[(i + 1) % m], p);
		const vec3_t cr = v_cross(r0, r1);
		const double lc = v_norm(cr);
		double a0, a1, ct;

		if (lc < 1e-300) continue;
		a0 = v_norm(r0);
		a1 = v_norm(r1);
		if ((a0 < 1e-300) || (a1 < 1e-300)) continue;
		ct = v_dot(r0, r1) / (a0 * a1);
		if (ct > 1.0) ct = 1.0;
		if (ct < -1.0) ct = -1.0;
		sum += v_dot(n, v_scale(cr, 1.0 / lc)) * acos(ct);
	}

	sum = fabs(sum) / (2.0 * PI);
	return (sum > 0.0) ? sum : 0.0;
}

/*
矩形パッチ上の複合 Gauss 求積点を作る (msub × msub の小矩形 × 3 点 Gauss)。
xs / ws は最大 (msub*3)^2 要素。重みの総和は 1 (面積平均)。
*/
static int rect_quad(const patch_t *q, int msub, vec3_t *xs, double *ws)
{
	const vec3_t eu = v_sub(q->p[1], q->p[0]);
	const vec3_t ev = v_sub(q->p[3], q->p[0]);
	const double h = 1.0 / msub;
	int    su, sv, gu, gv, k = 0;

	for (sv = 0; sv < msub; sv++) {
		for (su = 0; su < msub; su++) {
			for (gv = 0; gv < 3; gv++) {
				for (gu = 0; gu < 3; gu++) {
					const double u = (su + gx3[gu]) * h;
					const double v = (sv + gx3[gv]) * h;
					xs[k] = v_add(q->p[0], v_add(v_scale(eu, u), v_scale(ev, v)));
					ws[k] = gw3[gu] * gw3[gv] * h * h;
					k++;
				}
			}
		}
	}
	return k;
}

/*
パッチ q から平面多角形 poly への形態係数 (受光側を複合 Gauss 求積)。
群落の上下面のような「パッチでない多角形」を相手にするとき用。

gn が零ベクトルでなければ、半空間 gn·(x - gp) > 0 にある求積点だけを
数える。群落面をまたぐ壁パッチが面の裏側からも寄与を拾わないようにする
(ff_point_poly は受光側の法線しか見ないので、これは呼び出し側の責任)。
*/
double ff_patch_poly(const patch_t *q, const vec3_t *poly, int nv, int msub,
                     vec3_t gp, vec3_t gn)
{
	vec3_t xs[144];
	double ws[144];
	const int gate = (v_dot(gn, gn) > 0.0);
	double f = 0.0;
	int    nq, g;

	if (msub < 1) msub = 1;
	if (msub > 4) msub = 4;
	nq = rect_quad(q, msub, xs, ws);
	for (g = 0; g < nq; g++) {
		if (gate && (v_dot(gn, v_sub(xs[g], gp)) <= 0.0)) continue;
		f += ws[g] * ff_point_poly(xs[g], q->n, poly, nv);
	}
	return f;
}

/* パッチの代表長さ (対角長) */
static double patch_size(const patch_t *q)
{
	return v_norm(v_sub(q->p[2], q->p[0]));
}

/*
可視率の標本点の 1 次元座標 (4 点の層化標本)。

小矩形の中心 (k + 1/2)/4 のような有理数の格子だと、2 桁の小数で書いた板の
辺に「光線の族」がまとめて乗る (例 : 板の高さ 0.53、辺 0.31 で
0.47 * 17/48 + 0.53 * 13/48 = 0.31 ちょうど)。乗った光線の当たり外れは
丸めで決まり、16 本単位で可視率が跳ぶ。そこで各層の中の位置を無理数に
ずらす。

  - どちらの集合も u -> 1 - u で閉じている : パッチを鏡映しても (角の順序が
    変わっても) 同じ点の集合になるので、鏡面の二重箱等価が保たれる。
  - 受光側 (vis_a : 黄金比) と放射側 (vis_b : √2) で別の集合を使う : 同じ
    集合だと s と 1 - s の組の中点が有理数に戻り、ちょうど中間の高さに
    置いた板の辺に乗る。√5 と √2 は有理数上で独立なので、2 点を結ぶ線分が
    有理数の座標を通ることはない。
*/
static const double vis_a[4] = {
	0.12732200375003502,        /* (2 - φ) / 3 */
	0.3819660112501051,        /* 2 - φ */
	0.6180339887498949,        /* φ - 1 = 1 - (2 - φ) */
	0.872677996249965         /* 1 - (2 - φ) / 3 */
};
static const double vis_b[4] = {
	0.13807118745769839,        /* (√2 - 1) / 3 */
	0.41421356237309515,        /* √2 - 1 */
	0.58578643762690485,        /* 2 - √2 = 1 - (√2 - 1) */
	0.86192881254230158         /* 1 - (√2 - 1) / 3 */
};

/* パッチ上の層化標本点 (4 x 4)。set = 0 : 受光側 (vis_a)、1 : 放射側 (vis_b) */
static vec3_t patch_sample(const patch_t *q, int iu, int iv, int set)
{
	const vec3_t eu = v_sub(q->p[1], q->p[0]);
	const vec3_t ev = v_sub(q->p[3], q->p[0]);
	const double *t = (set == 0) ? vis_a : vis_b;

	return v_add(q->p[0], v_add(v_scale(eu, t[iu]), v_scale(ev, t[iv])));
}

/* 点の並びの軸並行の外接箱 */
static void bbox_add(vec3_t *lo, vec3_t *hi, vec3_t a)
{
	if (a.x < lo->x) lo->x = a.x;
	if (a.y < lo->y) lo->y = a.y;
	if (a.z < lo->z) lo->z = a.z;
	if (a.x > hi->x) hi->x = a.x;
	if (a.y > hi->y) hi->y = a.y;
	if (a.z > hi->z) hi->z = a.z;
}

/*
可視率の推定 (vis_point_patch / vis_patch_patch) の組み立て :

  1. 両端 (点とパッチ、またはパッチとパッチ) の外接箱がどの遮蔽物とも
     重ならなければ、どの光線も当たりえないので厳密に 1。
  2. 重なれば、常に細かい標本 (点 -> パッチは 4 x 4 = 16 本、
     パッチ -> パッチは 4 x 4 x 4 x 4 = 256 本) で数える。

以前は 2 x 2 の粗い標本が全一致ならそこで打ち切っていた。それだと影の
境界が粗い標本のすき間を通ると見落とし、1 本でも判定が割れた瞬間に細かい
推定へ跳ぶので、規則格子の光線が板の辺にちょうど乗る配置では丸めの違い
1 本で F が ~15% 動いた (slab 法と平面交差で結果が割れた実例あり)。
いまの組み立てでは、1 の判定の境目 (外接箱が触れ始める所) で細かい推定も
「全光線が通る = 1」になるので値が跳ばない。残る段差は光線 1 本ぶん
(1/16 または 1/256) だけ。

*refined は可視率が 0 と 1 の間 (部分的な影) になった対で 1 にする (診断)。
n は q が鏡像パッチのときの展開指数 (実在のパッチなら NULL)。
*/
double vis_point_patch(const ppfd_t *p, vec3_t x, const patch_t *q, const int *n, int *refined)
{
	const int nsub = 4;
	vec3_t lo = x, hi = x;
	int    nvis = 0, ju, jv, v;

	if (!occ_any(p)) return 1.0;
	for (v = 0; v < 4; v++) bbox_add(&lo, &hi, q->p[v]);
	if (!occ_aabb_overlap(p, lo, hi, n)) return 1.0;

	for (jv = 0; jv < nsub; jv++) {
		for (ju = 0; ju < nsub; ju++) {
			if (!occ_blocked_img(p, x, patch_sample(q, ju, jv, 1), n)) nvis++;
		}
	}
	if ((nvis > 0) && (nvis < nsub * nsub) && (refined != NULL)) *refined = 1;
	return (double)nvis / (nsub * nsub);
}

/* パッチ i からパッチ j の可視率 (0..1)。組み立ては vis_point_patch と同じ */
double vis_patch_patch(const ppfd_t *p, const patch_t *qi, const patch_t *qj, const int *n, int *refined)
{
	const int nsub = 4;
	const int ntot = nsub * nsub * nsub * nsub;
	vec3_t lo = qi->p[0], hi = qi->p[0];
	int    nvis = 0, iu, iv, ju, jv, v;

	if (!occ_any(p)) return 1.0;
	for (v = 1; v < 4; v++) bbox_add(&lo, &hi, qi->p[v]);
	for (v = 0; v < 4; v++) bbox_add(&lo, &hi, qj->p[v]);
	if (!occ_aabb_overlap(p, lo, hi, n)) return 1.0;

	for (iv = 0; iv < nsub; iv++) {
		for (iu = 0; iu < nsub; iu++) {
			const vec3_t si = patch_sample(qi, iu, iv, 0);
			for (jv = 0; jv < nsub; jv++) {
				for (ju = 0; ju < nsub; ju++) {
					if (!occ_blocked_img(p, si, patch_sample(qj, ju, jv, 1), n)) nvis++;
				}
			}
		}
	}
	if ((nvis > 0) && (nvis < ntot) && (refined != NULL)) *refined = 1;
	return (double)nvis / ntot;
}

/* 2 つのパッチが同一平面にあるか (互いに見えない) */
static int coplanar(const patch_t *a, const patch_t *b)
{
	if (fabs(v_dot(a->n, b->n)) < (1.0 - 1e-12)) return 0;
	return (fabs(v_dot(a->n, v_sub(b->c, a->c))) < 1e-9);
}

void setup_ff(ppfd_t *p)
{
	const int n = p->npatch;
	int    i;
	int    npartial = 0;

	p->ff = (double *)xcalloc((size_t)n * n, sizeof(double));
	if (p->canopy.on) {
		p->plen = (float *)xcalloc((size_t)n * n, sizeof(float));
	}

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) reduction(+:npartial)
#endif
	for (i = 0; i < n; i++) {
		int    j;
		vec3_t xs[144];
		double ws[144];
		const patch_t *qi = &p->patch[i];
		const double si = patch_size(qi);

		for (j = 0; j < n; j++) {
			const patch_t *qj = &p->patch[j];
			double dist, sj, f, vis;
			int    msub, nq, g, refined = 0;

			if (i == j) continue;
			if (qi->iface == qj->iface) continue;   /* 同一面 : 互いに見えない */
			if (coplanar(qi, qj)) continue;         /* 同一平面 (遮蔽物と壁の接触等) */

			if (p->plen != NULL) {
				p->plen[((size_t)i * n) + j] = (float)canopy_path(p, qi->c, qj->c);
			}

			vis = vis_patch_patch(p, qi, qj, NULL, &refined);
			npartial += refined;
			if (vis <= 0.0) continue;

			sj = patch_size(qj);
			dist = v_norm(v_sub(qj->c, qi->c));
			/* 近接パッチほど被積分関数が急峻なので分割を増やす */
			msub = (dist < (0.75 * (si + sj))) ? 4 : ((dist < (2.0 * (si + sj))) ? 2 : 1);
			nq = rect_quad(qi, msub, xs, ws);

			f = 0.0;
			for (g = 0; g < nq; g++) {
				/*
				放射側の向きの判定 : 求積点が qj の裏側にあるなら qj の
				表面は見えない。ff_point_poly は受光側の法線しか見ないので
				ここで落とす必要がある (凸キャビティだけなら常に真だが、
				遮蔽物の表裏の面は同じ多角形を共有するため必須)。
				*/
				if (v_dot(qj->n, v_sub(xs[g], qj->p[0])) <= 0.0) continue;
				f += ws[g] * ff_point_poly(xs[g], qi->n, qj->p, 4);
			}
			p->ff[((size_t)i * n) + j] = f * vis;
		}
	}
	p->ff_npartial = npartial;

	/*
	診断 : 閉じた面の集合なら行和 = 1、相反則 A_i F_ij = A_j F_ji。
	遮蔽物が壁と接して完全に囲まれたパッチ (何も見えない) は行和 0 に
	なるが、そこには光も入らないので収支には効かない。数だけ数えて
	行和の判定からは外す。
	*/
	{
		double emax = 0.0, rmax = 0.0;
		int    j, nblind = 0;
		for (i = 0; i < n; i++) {
			double s = 0.0;
			for (j = 0; j < n; j++) s += p->ff[((size_t)i * n) + j];
			if (s <= 0.0) nblind++;
			else if (fabs(s - 1.0) > emax) emax = fabs(s - 1.0);
			for (j = i + 1; j < n; j++) {
				const double a = p->patch[i].area * p->ff[((size_t)i * n) + j];
				const double b = p->patch[j].area * p->ff[((size_t)j * n) + i];
				const double m = (fabs(a) > fabs(b)) ? fabs(a) : fabs(b);
				if (m > 1e-12) {
					const double e = fabs(a - b) / m;
					if (e > rmax) rmax = e;
				}
			}
		}
		p->ff_rowsum_err = emax;
		p->ff_recip_err = rmax;
		p->ff_nblind = nblind;
	}

	/*
	■ 鏡面壁の拡張形態係数
	拡散光が鏡面壁で正反射して届く経路を「放射側パッチの鏡像」として
	加算する。軸並行の直方体では、キャビティ内の点から鏡像パッチへの
	直線は必ず鏡面の矩形内で平面を横切る (鏡映が各座標を箱の範囲の
	鏡像へ写すため) ので、空のチャンバではこの拡張は厳密になる。

	これを入れないと、拡散壁からの放射が鏡面壁に当たったとき ρs 分が
	どこにも運ばれずに消え、収支が ρs × (鏡面壁への拡散入射) だけ
	破れる (実測 : 拡散床 ρd=0.6 + 鏡面天井 ρs=0.8 で closure -7.7%)。

	放射側パッチが最初の折り返し面の上にある場合は除外する (鏡映で
	自分自身に写り、直接項と二重計上になるため)。行和・相反則の診断は
	拡張前の基本形態係数に対して計算してある (拡張後の行和は 1 を
	超えるのが正しい : 鏡面壁への到達分は全量で数え、その続きの
	正反射分も別項として数える帳簿だから)。

	群落・遮蔽物との併用は入力段で弾かれているので、ここでは減衰も
	可視率も考えなくてよい。
	*/
	{
		strans_t *tr = (strans_t *)xmalloc((size_t)MAXSPECT * sizeof(strans_t));
		const int nt = spec_transforms(p, tr, MAXSPECT);

		/*
		群落があるときは ff へ畳み込まず分離保持する。鏡像経路の群落内
		経路長は直接経路と異なるので、減衰を掛けるには項を分けておく
		必要がある。経路長は canopy_path が厳密に測る : 側壁の鏡映は
		z を保つのでそのまま、床・天井の鏡映は鏡像への直線の z を周期
		2Lz の三角波で折り返すと実経路の z プロファイルになる。
		群落が無ければ従来どおり畳み込む。
		*/
		if ((nt > 0) && p->canopy.on) {
			if ((size_t)nt * (size_t)n * (size_t)n > (size_t)32 * 1024 * 1024) {
				fprintf(stderr, "*** specular + canopy : too many transform entries "
					"(reduce patchdiv or specbounce)\n");
				exit(1);
			}
			p->nst = nt;
			p->ffs = (double *)xcalloc((size_t)nt * n * n, sizeof(double));
			p->plens = (float *)xcalloc((size_t)nt * n * n, sizeof(float));
		}

		if (nt > 0) {
			int is;
			int npart_s = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) reduction(+:npart_s)
#endif
			for (is = 0; is < n; is++) {
				vec3_t xs[144];
				double ws[144];
				const patch_t *qi = &p->patch[is];
				const double si = patch_size(qi);
				int it, j;

				for (it = 0; it < nt; it++) {
					for (j = 0; j < n; j++) {
						const patch_t *qj = &p->patch[j];
						vec3_t q4[4], nj, cj;
						double dist, sj, f;
						int    msub, nq, g, v;

						if (spec_on_mirror(p, &tr[it], qj)) continue;
						for (v = 0; v < 4; v++) q4[v] = spec_apply(p, &tr[it], qj->p[v]);
						nj = spec_apply_dir(&tr[it], qj->n);
						cj = spec_apply(p, &tr[it], qj->c);

						sj = patch_size(qj);
						dist = v_norm(v_sub(cj, qi->c));
						msub = (dist < (0.75 * (si + sj))) ? 4 : ((dist < (2.0 * (si + sj))) ? 2 : 1);
						nq = rect_quad(qi, msub, xs, ws);

						f = 0.0;
						for (g = 0; g < nq; g++) {
							/* 鏡像パッチの表側にある求積点だけが受け取る */
							if (v_dot(nj, v_sub(xs[g], q4[0])) <= 0.0) continue;
							f += ws[g] * ff_point_poly(xs[g], qi->n, q4, 4);
						}
						if ((f > 0.0) && occ_any(p)) {
							/*
							棚板の影 : 鏡像パッチへの直線を折り返して、通りうる
							チャンバの複製それぞれにある棚板の鏡像と照合する
							(occ_blocked_img)。鏡像パッチは角を写した仮のパッチで表す。
							*/
							patch_t qm = *qj;
							int     refined = 0;
							for (v = 0; v < 4; v++) qm.p[v] = q4[v];
							qm.n = nj;
							qm.c = cj;
							f *= vis_patch_patch(p, qi, &qm, tr[it].n, &refined);
							npart_s += refined;
						}
						if (p->nst > 0) {
							const size_t k = ((size_t)it * n * n) + ((size_t)is * n) + j;
							p->ffs[k] = tr[it].w * f;
							p->plens[k] = (float)canopy_path(p, qi->c, cj);
						}
						else {
							p->ff[((size_t)is * n) + j] += tr[it].w * f;
						}
					}
				}
			}
			if (occ_any(p)) {
				plog(p, "specular ff : partially shadowed mirrored patch pairs = %d\n", npart_s);
			}
			plog(p, "specular ff : %d mirror transforms %s\n",
				nt, (p->nst > 0) ? "kept separate (canopy attenuation per folded path)"
				                 : "folded into the form factors");
		}
		free(tr);
	}
}
