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

を掛ける。両パッチの外接箱がどの遮蔽物とも重ならなければ厳密に 1、
重なれば 4×4 × 4×4 = 256 本の線分で数える (下の vis_patch_patch)。
棚板がキャビティを仕切るような配置ではすべての対が 0 か 1 になり、遮蔽が
あっても F は厳密なままになる。0 と 1 の間になった対の数はログに出す
(標本化に頼った = 厳密でない箇所がどれだけあったか)。

■ 群落
群落があるときは対ごとの経路長の分布を 2 点 Gauss 則に縮めて持つ
(canopy_plen の説明を参照)。
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

/*
■ パッチ間の群落透過率

形態係数の重みで平均した透過率は

    τ_ij(k) = ∫∫ K e^{-k s(x,y)} dA_x dA_y / ∫∫ K dA_x dA_y,
    K = cos θ_x cos θ_y / r²

で、経路長 s が対の中でばらつくので重心どうしの 1 本 (exp(-k s_c)) では
決まらない (指数の凸性で必ず過小評価になる。平行平板の間の拡散透過
2E₃(τ₀) で -0.35%)。k は葉の散乱係数を通じて波長ごとに変わるので、
対ごとに 1 個の透過率を持つことはできない。そこで経路長の分布 (K で
重み付けした標本) を 2 点の Gauss 則

    τ(k) ≈ w₁ e^{-k s₁} + (1 - w₁) e^{-k s₂}

に縮めて持つ。2 点は分布の 0〜3 次のモーメントを再現するので、
τ(k) を k = 0 のまわりで 3 次まで正しく与え、分布が 2 点以下の原子なら
(全部が同じ長さ、または「群落を通る / 通らない」が半々など) 厳密になる。

標本は両パッチ上の Gauss-Legendre 点の直積 (近接対は 4×4 × 4×4、遠い対は
2×2 × 2×2)。s(x,y) は端点について連続なので、影の判定のような有理数の
格子の問題は起きない (遮蔽物の影になる標本対は、可視率と同じく通らない
光線として重みから外す)。

leafscatter = on では k = G a が波長によらないので、標本の平均 τ を
そのまま 1 個の等価経路長 -ln τ / k で持つ (2 点に縮めない)。層への
預け入れ (canopy_setup) も同じ標本で作るので、対ごとに「透過 + 預け入れ
= 1」がそのまま成り立つ。
*/
static const double gx2[2] = {
	0.5 * (1.0 - 0.5773502691896258), 0.5 * (1.0 + 0.5773502691896258)
};
static const double gw2[2] = { 0.5, 0.5 };
static const double gx4[4] = {
	0.5 * (1.0 - 0.8611363115940526), 0.5 * (1.0 - 0.3399810435848563),
	0.5 * (1.0 + 0.3399810435848563), 0.5 * (1.0 + 0.8611363115940526)
};
static const double gw4[4] = {
	0.5 * 0.3478548451374538, 0.5 * 0.6521451548625461,
	0.5 * 0.6521451548625461, 0.5 * 0.3478548451374538
};

/*
受光側 (点 xr[0..nr-1]、重み wr、法線 nrm) と放射側 (角 q4、法線 nj) の
経路の標本を作る。放射側は m × m の Gauss 点。tn は放射側が鏡像パッチの
ときの展開指数 (実在なら NULL)。重みは総和 1 に規格化して返す。
戻り値 = 標本数 (0 : 群落を通りえない、または有効な標本が無い)。
*/
static int canopy_samples(const ppfd_t *p, const vec3_t *xr, const double *wr, int nr, vec3_t nrm,
                          const vec3_t *q4, vec3_t nj, int m, const int *tn,
                          vec3_t *xa, vec3_t *yb, double *W)
{
	const vec3_t eu = v_sub(q4[1], q4[0]);
	const vec3_t ev = v_sub(q4[3], q4[0]);
	const double *gx = (m == 4) ? gx4 : gx2;
	const double *gw = (m == 4) ? gw4 : gw2;
	vec3_t lo = q4[0], hi = q4[0];
	double wsum = 0.0;
	int    ns = 0, occl, a, iu, iv, v;

	for (v = 1; v < 4; v++) bbox_add(&lo, &hi, q4[v]);
	for (a = 0; a < nr; a++) bbox_add(&lo, &hi, xr[a]);
	if (!canopy_may_cross(p, lo.z, hi.z)) return 0;
	occl = occ_aabb_overlap(p, lo, hi, tn);

	for (a = 0; a < nr; a++) {
		for (iv = 0; iv < m; iv++) {
			for (iu = 0; iu < m; iu++) {
				const vec3_t y = v_add(q4[0], v_add(v_scale(eu, gx[iu]), v_scale(ev, gx[iv])));
				const vec3_t d = v_sub(y, xr[a]);
				const double r2 = v_dot(d, d);
				const double ci = v_dot(nrm, d);
				const double cj = -v_dot(nj, d);
				double w;
				if ((ci <= 0.0) || (cj <= 0.0) || (r2 <= 0.0)) continue;
				if (occl && occ_blocked_img(p, xr[a], y, tn)) continue;
				w = wr[a] * gw[iu] * gw[iv] * ci * cj / (r2 * r2);
				xa[ns] = xr[a];
				yb[ns] = y;
				W[ns] = w;
				wsum += w;
				ns++;
			}
		}
	}
	if ((ns == 0) || (wsum <= 0.0)) return 0;
	for (a = 0; a < ns; a++) W[a] /= wsum;
	return ns;
}

/* パッチ q の m × m Gauss 点 (重みの総和 1) */
static int patch_gauss(const patch_t *q, int m, vec3_t *xs, double *ws)
{
	const vec3_t eu = v_sub(q->p[1], q->p[0]);
	const vec3_t ev = v_sub(q->p[3], q->p[0]);
	const double *gx = (m == 4) ? gx4 : gx2;
	const double *gw = (m == 4) ? gw4 : gw2;
	int iu, iv, k = 0;

	for (iv = 0; iv < m; iv++) {
		for (iu = 0; iu < m; iu++) {
			xs[k] = v_add(q->p[0], v_add(v_scale(eu, gx[iu]), v_scale(ev, gx[iv])));
			ws[k] = gw[iu] * gw[iv];
			k++;
		}
	}
	return k;
}

/*
パッチ qi (受光側) とパッチ qj (放射側。鏡像なら角 q4・法線 nj・展開指数
tn) の経路の標本。近接対ほど経路長のばらつきが大きいので細かく取る。
setup_ff と canopy_setup (層への預け入れ) が同じ標本を使うための窓口。
*/
int canopy_pair_samples(const ppfd_t *p, const patch_t *qi, const vec3_t *q4, vec3_t nj,
                        const int *tn, vec3_t *xa, vec3_t *yb, double *W)
{
	const vec3_t cj = v_scale(v_add(q4[0], q4[2]), 0.5);
	const double si = patch_size(qi);
	const double sj = v_norm(v_sub(q4[2], q4[0]));
	const double dist = v_norm(v_sub(cj, qi->c));
	const int    m = (dist < (2.0 * (si + sj))) ? 4 : 2;
	vec3_t xs[16];
	double ws[16];
	const int nr = patch_gauss(qi, m, xs, ws);

	return canopy_samples(p, xs, ws, nr, qi->n, q4, nj, m, tn, xa, yb, W);
}

/* 点 x (法線 nrm) とパッチ (角 q4・法線 nj・展開指数 tn) の経路の標本 */
int canopy_point_samples(const ppfd_t *p, vec3_t x, vec3_t nrm, const vec3_t *q4, vec3_t nj,
                         const int *tn, vec3_t *xa, vec3_t *yb, double *W)
{
	const vec3_t cj = v_scale(v_add(q4[0], q4[2]), 0.5);
	const double sj = v_norm(v_sub(q4[2], q4[0]));
	const double dist = v_norm(v_sub(cj, x));
	const int    m = (dist < (2.0 * sj)) ? 4 : 2;
	const double w1 = 1.0;

	return canopy_samples(p, &x, &w1, 1, nrm, q4, nj, m, tn, xa, yb, W);
}

/*
標本 (xa, yb, W) の経路長の分布を 2 点 Gauss 則 (s₁, s₂, w₁) に縮めて
out[0..2] へ書く。ns = 0 なら経路長 0 (減衰なし)。leafscatter のときは
平均の透過率を 1 個の等価経路長で持つ (上の説明)。
*/
void canopy_plen(const ppfd_t *p, int ns, const vec3_t *xa, const vec3_t *yb, const double *W,
                 float *out)
{
	double s[CANOPY_MAXSAMP];
	double mu = 0.0, var = 0.0, m3 = 0.0;
	int    a;

	out[0] = 0.0f;
	out[1] = 0.0f;
	out[2] = 1.0f;
	if (ns <= 0) return;

	for (a = 0; a < ns; a++) {
		s[a] = canopy_path(p, xa[a], yb[a]);
		mu += W[a] * s[a];
	}
	if (p->canopy.scatter) {
		const double k = p->canopy.k0 * p->canopy.a;
		double tr = 0.0;
		for (a = 0; a < ns; a++) tr += W[a] * exp(-k * s[a]);
		if ((k > 0.0) && (tr > 0.0)) out[0] = (float)(-log(tr) / k);
		else                          out[0] = (float)mu;
		out[1] = out[0];
		return;
	}
	for (a = 0; a < ns; a++) {
		const double y = s[a] - mu;
		var += W[a] * y * y;
		m3  += W[a] * y * y * y;
	}
	if (var <= (1e-14 * ((mu * mu) + 1e-300))) {
		out[0] = (float)mu;
		out[1] = (float)mu;
		return;
	}
	{
		/* 平均まわりの 2 次の直交多項式 y² - c y - σ² の根 */
		const double c = m3 / var;
		const double r = sqrt((c * c) + (4.0 * var));
		const double y1 = 0.5 * (c - r);
		const double y2 = 0.5 * (c + r);
		double s1 = mu + y1;
		const double s2 = mu + y2;
		if (s1 < 0.0) s1 = 0.0;
		out[0] = (float)s1;
		out[1] = (float)s2;
		out[2] = (float)(y2 / (y2 - y1));
	}
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
		p->plen = (float *)xcalloc((size_t)n * n * PLEN_N, sizeof(float));
	}

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) reduction(+:npartial)
#endif
	for (i = 0; i < n; i++) {
		int    j;
		vec3_t xs[144];
		double ws[144];
		vec3_t xa[CANOPY_MAXSAMP], yb[CANOPY_MAXSAMP];
		double W[CANOPY_MAXSAMP];
		const patch_t *qi = &p->patch[i];
		const double si = patch_size(qi);

		for (j = 0; j < n; j++) {
			const patch_t *qj = &p->patch[j];
			double dist, sj, f, vis;
			int    msub, nq, g, refined = 0;

			if (i == j) continue;
			if (p->pcov[i] || p->pcov[j]) continue; /* 接して覆われたパッチ : 光が届かない */
			if (qi->iface == qj->iface) continue;   /* 同一面 : 互いに見えない */
			if (coplanar(qi, qj)) continue;         /* 同一平面 (遮蔽物と壁の接触等) */

			vis = vis_patch_patch(p, qi, qj, NULL, &refined);
			npartial += refined;
			if (vis <= 0.0) continue;

			if (p->plen != NULL) {
				const int ns = canopy_pair_samples(p, qi, qj->p, qj->n, NULL, xa, yb, W);
				canopy_plen(p, ns, xa, yb, W, &p->plen[(((size_t)i * n) + j) * PLEN_N]);
			}

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

	群落と併用するときは鏡像経路ごとの減衰を後で掛けるので分離保持し
	(下記)、遮蔽物があれば鏡像パッチへの直線で可視率を掛ける。
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
			p->plens = (float *)xcalloc((size_t)nt * n * n * PLEN_N, sizeof(float));
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
				vec3_t xa[CANOPY_MAXSAMP], yb[CANOPY_MAXSAMP];
				double W[CANOPY_MAXSAMP];
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
						if (p->pcov[is] || p->pcov[j]) continue;
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
							if (f > 0.0) {
								/* 鏡像パッチへの直線の z を折り返して測る (canopy_path) */
								const int ns = canopy_pair_samples(p, qi, q4, nj, tr[it].n, xa, yb, W);
								canopy_plen(p, ns, xa, yb, W, &p->plens[k * PLEN_N]);
							}
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
