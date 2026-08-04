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

キャビティが凸なので遮蔽判定は不要。同一面のパッチどうしは互いに
見えないので明示的に 0 とする。
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

	m = clip_halfspace(poly, nv, p, n, buf, 12);
	if (m < 3) return 0.0;

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

/* パッチの代表長さ (対角長) */
static double patch_size(const patch_t *q)
{
	return v_norm(v_sub(q->p[2], q->p[0]));
}

void setup_ff(ppfd_t *p)
{
	const int n = p->npatch;
	int    i;

	p->ff = (double *)xcalloc((size_t)n * n, sizeof(double));
	if (p->canopy.on) {
		p->plen = (float *)xcalloc((size_t)n * n, sizeof(float));
	}

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
	for (i = 0; i < n; i++) {
		int    j;
		vec3_t xs[144];
		double ws[144];
		const patch_t *qi = &p->patch[i];
		const double si = patch_size(qi);

		for (j = 0; j < n; j++) {
			const patch_t *qj = &p->patch[j];
			double dist, sj, f;
			int    msub, nq, g;

			if (i == j) continue;
			if (qi->iface == qj->iface) continue;   /* 同一平面 : 互いに見えない */

			if (p->plen != NULL) {
				p->plen[((size_t)i * n) + j] = (float)canopy_path(p, qi->c, qj->c);
			}

			sj = patch_size(qj);
			dist = v_norm(v_sub(qj->c, qi->c));
			/* 近接パッチほど被積分関数が急峻なので分割を増やす */
			msub = (dist < (0.75 * (si + sj))) ? 4 : ((dist < (2.0 * (si + sj))) ? 2 : 1);
			nq = rect_quad(qi, msub, xs, ws);

			f = 0.0;
			for (g = 0; g < nq; g++) {
				f += ws[g] * ff_point_poly(xs[g], qi->n, qj->p, 4);
			}
			p->ff[((size_t)i * n) + j] = f;
		}
	}

	/* 診断 : 閉キャビティなら行和 = 1、相反則 A_i F_ij = A_j F_ji */
	{
		double emax = 0.0, rmax = 0.0;
		int    j;
		for (i = 0; i < n; i++) {
			double s = 0.0;
			for (j = 0; j < n; j++) s += p->ff[((size_t)i * n) + j];
			if (fabs(s - 1.0) > emax) emax = fabs(s - 1.0);
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
	}
}
