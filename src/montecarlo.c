/*
montecarlo.c

モンテカルロ光子追跡 (montecarlo キー)。決定論の解法 (ラジオシティ + 鏡像法 +
層別二流) の外にある 2 つを扱う :

  1. 一般の反射 (BRDF) : Lambert 拡散 + 鏡面 + 光沢 (Phong ローブ) を、壁・
     棚板・斜めの板のどの面でも。ラジオシティは放射輝度が向きによらない前提
     なので光沢を表せないが、光子は反射のたびに向きを引き直すだけでよい。
  2. 群落の 3 次元の散乱 : 光子は群落の中で素の遮断 k = G a で止まり、葉の
     散乱係数 ω(λ) = ρ + τ の確率で等方に散乱する (β = 1/2 の二流が近似して
     いる過程そのもの)。横方向 (側壁・ベッドの側面・通路) のやりとりも
     そのまま入る。

■ 光子の重みと再現性
全光子を同じ重み w = Φ_total / N にする (光源は放射束の比で、波長ビンは
スペクトルの重みで確率的に選ぶ)。すると集計は「光子の数」(整数) になり、
足す順序によらず厳密に同じ値になるので、スレッド数によらず出力がビット
一致する。光子 i の乱数は (seed, i) から作る (splitmix64)。

■ 集計
  - 測定面 : 光子の自由行程が受光面を表側から横切ったら、そのセルに 1 個
    (面は光を遮らない。決定論の測定面と同じ)。セル平均の放射照度
    E = (個数) w / (セル面積)。
  - 吸収 : 壁 (棚板・板を含む) / 群落 / 打ち切り (MC_MAXEV 回の相互作用を
    超えた光子)。光子は必ずどれかに数えられるので、収支は打ち切り分を除いて
    厳密に閉じる。

決定論の解法との違い (README に明記) :
  - 測定面の値はセル平均 (決定論はセル中心の値)。
  - 群落は常に散乱を含む (leafscatter の有無によらない)。散乱なしの
    決定論の √(1-ω) は、この過程の近似。
  - 鏡面は実際に反射させる (鏡像も specbounce の打ち切りも無い)。
  - 統計誤差 : セルの相対誤差 ~ 1/√(個数)。ログに測定面ごとの個数を出す。
*/

#include "ppfd.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define MC_MAXEV 100000             /* 1 光子あたりの相互作用の上限 (超えたら打ち切りに数える) */

typedef unsigned long long u64;

static u64 sm64(u64 *s)
{
	u64 z = (*s += 0x9E3779B97F4A7C15ULL);
	z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
	z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
	return z ^ (z >> 31);
}

/* [0, 1) の一様乱数 (53 bit) */
static double urand(u64 *s)
{
	return (double)(sm64(s) >> 11) * (1.0 / 9007199254740992.0);
}

static vec3_t unitv(vec3_t a)
{
	return v_scale(a, 1.0 / v_norm(a));
}

/* 単位ベクトル n に直交する正規直交基底 */
static void frame(vec3_t n, vec3_t *t1, vec3_t *t2)
{
	const vec3_t a = (fabs(n.x) < 0.9) ? v_make(1.0, 0.0, 0.0) : v_make(0.0, 1.0, 0.0);
	*t1 = unitv(v_cross(n, a));
	*t2 = v_cross(n, *t1);
}

/* 軸 n のまわりに cos θ = ct の向き (方位角は一様) */
static vec3_t around(vec3_t n, double ct, u64 *s)
{
	vec3_t t1, t2;
	const double st = sqrt((ct < 1.0) ? (1.0 - (ct * ct)) : 0.0);
	const double ph = 2.0 * PI * urand(s);

	frame(n, &t1, &t2);
	return v_add(v_scale(n, ct), v_add(v_scale(t1, st * cos(ph)), v_scale(t2, st * sin(ph))));
}

static vec3_t sphere(u64 *s)
{
	return around(v_make(0.0, 0.0, 1.0), (2.0 * urand(s)) - 1.0, s);
}

/* 光源 e の配光にしたがう放射の向き。dmax = 実測配光の最大値 (棄却法) */
static vec3_t emit_dir(const ppfd_t *p, const emitter_t *e, double dmax, u64 *s)
{
	if (e->idist >= 0) {
		/* 実測配光 : 一様な球面の向きを Î / Î_max で受理する (∫Î dΩ = 1 に正規化済み) */
		int k;
		for (k = 0; k < 100000; k++) {
			const vec3_t u = sphere(s);
			double cg = v_dot(u, e->dir), ca, iv;
			if (cg > 1.0) cg = 1.0;
			if (cg < -1.0) cg = -1.0;
			ca = atan2(v_dot(u, e->cy), v_dot(u, e->cx)) * (180.0 / PI) - e->crot;
			iv = photdist_value(&p->dist[e->idist], acos(cg) * (180.0 / PI), ca);
			if (urand(s) * dmax < iv) return u;
		}
		return e->dir;
	}
	if (e->mexp < 0.0) return sphere(s);
	/* I ∝ cos^m θ (半球) : cos θ = u^{1/(m+1)} */
	return around(e->dir, pow(urand(s), 1.0 / (e->mexp + 1.0)), s);
}

/* 面との交点 */
typedef struct {
	double t;
	vec3_t n;                       /* 入ってきた側を向く単位法線 */
	int    imat;
} hit_t;

static void hit_walls(const ppfd_t *p, vec3_t o, vec3_t d, hit_t *h)
{
	const double L[3] = { p->Lx, p->Ly, p->Lz };
	const double oc[3] = { o.x, o.y, o.z };
	const double dc[3] = { d.x, d.y, d.z };
	int i;

	for (i = 0; i < 3; i++) {
		double t;
		int    f;
		if (dc[i] > 0.0)      { t = (L[i] - oc[i]) / dc[i]; f = (2 * i) + 1; }
		else if (dc[i] < 0.0) { t = -oc[i] / dc[i];         f = 2 * i; }
		else continue;
		if (t < 0.0) t = 0.0;
		if (t < h->t) {
			h->t = t;
			h->imat = p->wallmat[f];
			h->n = v_make((i == 0) ? ((f & 1) ? -1.0 : 1.0) : 0.0,
			              (i == 1) ? ((f & 1) ? -1.0 : 1.0) : 0.0,
			              (i == 2) ? ((f & 1) ? -1.0 : 1.0) : 0.0);
		}
	}
}

static void hit_occluders(const ppfd_t *p, vec3_t o, vec3_t d, double teps, hit_t *h)
{
	const double oc[3] = { o.x, o.y, o.z };
	const double dc[3] = { d.x, d.y, d.z };
	int io, i;

	for (io = 0; io < p->nocc; io++) {
		const occluder_t *b = &p->occ[io];
		const double lo[3] = { b->x0, b->y0, b->z0 };
		const double hi[3] = { b->x1, b->y1, b->z1 };
		double tn = -1e300, tf = 1e300;
		int    an = -1, ok = 1;
		for (i = 0; (i < 3) && ok; i++) {
			if (fabs(dc[i]) < 1e-300) {
				if ((oc[i] < lo[i]) || (oc[i] > hi[i])) ok = 0;
			}
			else {
				double t1 = (lo[i] - oc[i]) / dc[i];
				double t2 = (hi[i] - oc[i]) / dc[i];
				if (t1 > t2) { const double t = t1; t1 = t2; t2 = t; }
				if (t1 > tn) { tn = t1; an = i; }
				if (t2 < tf) tf = t2;
				if (tf < tn) ok = 0;
			}
		}
		if (!ok || (an < 0) || (tn <= teps) || (tn >= h->t)) continue;
		h->t = tn;
		h->imat = b->imat;
		h->n = v_make((an == 0) ? ((dc[0] > 0.0) ? -1.0 : 1.0) : 0.0,
		              (an == 1) ? ((dc[1] > 0.0) ? -1.0 : 1.0) : 0.0,
		              (an == 2) ? ((dc[2] > 0.0) ? -1.0 : 1.0) : 0.0);
	}
}

static void hit_panels(const ppfd_t *p, vec3_t o, vec3_t d, double teps, hit_t *h)
{
	int ip;

	for (ip = 0; ip < p->npanel; ip++) {
		const panel_t *pn = &p->panel[ip];
		const vec3_t N = v_cross(pn->eu, pn->ev);
		const double nn = v_dot(N, N);
		const double dn = v_dot(d, N);
		double t, u, v;
		vec3_t w;
		if (fabs(dn) <= (1e-15 * sqrt(nn))) continue;
		t = v_dot(v_sub(pn->org, o), N) / dn;
		if ((t <= teps) || (t >= h->t)) continue;
		w = v_sub(v_add(o, v_scale(d, t)), pn->org);
		u = v_dot(v_cross(w, pn->ev), N) / nn;
		v = v_dot(v_cross(pn->eu, w), N) / nn;
		if ((u < 0.0) || (u > 1.0) || (v < 0.0) || (v > 1.0)) continue;
		h->t = t;
		h->imat = pn->imat;
		h->n = (dn < 0.0) ? pn->n : v_scale(pn->n, -1.0);
	}
}

/*
自由行程 [0, tmax] の中で群落に止められるか。止められたら距離を返す
(*js = スラブ)。止められなければ -1。光学深さ τ = -ln u を、光線が
スラブ (ベッド) の箱の中にある区間に沿って k_j = G_j a_j で消費する。
*/
static double canopy_collide(const ppfd_t *p, vec3_t o, vec3_t d, double tmax, u64 *s, int *js)
{
	double tau = -log(1.0 - urand(s));
	double t0[64], t1[64];
	int    jj[64], n = 0, j, k, i;
	const vec3_t b = v_add(o, v_scale(d, tmax));

	if (!p->canopy.on || (tmax <= 0.0)) return -1.0;
	for (j = 0; j < p->nslab; j++) {
		const cslab_t *sl = &p->slab[j];
		double u0[8], u1[8];
		const int ni = canopy_box_intervals(p, sl->x0, sl->x1, sl->y0, sl->y1, sl->zbot, sl->ztop,
		                                    o, b, u0, u1, 8);
		for (k = 0; (k < ni) && (n < 64); k++) {
			t0[n] = u0[k];
			t1[n] = u1[k];
			jj[n] = j;
			n++;
		}
	}
	/* 光線に沿った順 */
	for (i = 1; i < n; i++) {
		const double a0 = t0[i], a1 = t1[i];
		const int    aj = jj[i];
		k = i - 1;
		while ((k >= 0) && (t0[k] > a0)) {
			t0[k + 1] = t0[k];
			t1[k + 1] = t1[k];
			jj[k + 1] = jj[k];
			k--;
		}
		t0[k + 1] = a0;
		t1[k + 1] = a1;
		jj[k + 1] = aj;
	}
	for (i = 0; i < n; i++) {
		const cslab_t *sl = &p->slab[jj[i]];
		const double kk = sl->k0 * sl->a;
		const double len = (t1[i] - t0[i]) * tmax;
		if (kk <= 0.0) continue;
		if (tau <= kk * len) {
			*js = jj[i];
			return (t0[i] * tmax) + (tau / kk);
		}
		tau -= kk * len;
	}
	return -1.0;
}

/* 測定面のセルを横切った光子を数える (表側から入るもの) */
static void tally_targets(const ppfd_t *p, vec3_t o, vec3_t e, int il, long long *cnt, const int *off)
{
	int it;

	for (it = 0; it < p->ntarget; it++) {
		const target_t *t = &p->target[it];
		const double s0 = v_dot(t->n, v_sub(o, t->org));
		const double s1 = v_dot(t->n, v_sub(e, t->org));
		vec3_t x, w, a, b;
		double u, v;
		int    ix, iy;
		if (!((s0 > 0.0) && (s1 <= 0.0))) continue;
		x = v_add(o, v_scale(v_sub(e, o), s0 / (s0 - s1)));
		w = v_sub(x, t->org);
		/* w = u eu + v ev (eu と ev は直交とは限らない) */
		a = v_cross(t->ev, t->n);
		b = v_cross(t->n, t->eu);
		u = v_dot(w, a) / v_dot(t->eu, a);
		v = v_dot(w, b) / v_dot(t->ev, b);
		if ((u < 0.0) || (u >= 1.0) || (v < 0.0) || (v >= 1.0)) continue;
		ix = (int)(u * t->nx);
		iy = (int)(v * t->ny);
		if (ix > t->nx - 1) ix = t->nx - 1;
		if (iy > t->ny - 1) iy = t->ny - 1;
		cnt[(size_t)off[it] + ((((size_t)iy * t->nx) + ix) * p->nlam) + il]++;
	}
}

void mc_solve(ppfd_t *p)
{
	const int nl = p->nlam;
	const int nreal = p->nemit;
	const int N = p->mc_n;
	const double teps = 1e-10 * (p->Lx + p->Ly + p->Lz);
	double *ecdf = (double *)xcalloc((size_t)nreal, sizeof(double));
	double *dmax = (double *)xcalloc((size_t)nreal, sizeof(double));
	double *scdf = (double *)xcalloc((size_t)p->nspec * nl, sizeof(double));
	int    *off = (int *)xcalloc((size_t)p->ntarget + 1, sizeof(int));
	double  ftot = 0.0, w;
	int     nth = 1, ie, il, it, i, th;
	size_t  ncnt;
	long long *cnt, *wab, *cab, *lst;

	/* 光源の選択 (放射束の比)、実測配光の最大値、スペクトルの累積 */
	for (ie = 0; ie < nreal; ie++) {
		ftot += p->emit[ie].flux;
		ecdf[ie] = ftot;
		if (p->emit[ie].idist >= 0) {
			const photdist_t *d = &p->dist[p->emit[ie].idist];
			int k;
			for (k = 0; k < d->ng * d->nc; k++) {
				if (d->I[k] > dmax[ie]) dmax[ie] = d->I[k];
			}
		}
	}
	for (i = 0; i < p->nspec; i++) {
		double c = 0.0;
		for (il = 0; il < nl; il++) {
			c += p->spec[i].w[il];
			scdf[((size_t)i * nl) + il] = c;
		}
	}
	for (it = 0; it < p->ntarget; it++) {
		off[it + 1] = off[it] + (p->target[it].nx * p->target[it].ny * nl);
	}
	ncnt = (size_t)off[p->ntarget];
	w = (N > 0) ? (ftot / N) : 0.0;

#ifdef _OPENMP
	nth = omp_get_max_threads();
#endif
	/* スレッドごとの個数 (整数なので足す順序によらず同じ値) */
	cnt = (long long *)xcalloc(((size_t)nth * ncnt) + 1, sizeof(long long));
	wab = (long long *)xcalloc((size_t)nth * nl, sizeof(long long));
	cab = (long long *)xcalloc((size_t)nth * nl, sizeof(long long));
	lst = (long long *)xcalloc((size_t)nth * nl, sizeof(long long));

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 4096) private(th)
#endif
	for (i = 0; i < N; i++) {
		u64 s = (p->mc_seed * 0xD1B54A32D192ED03ULL) ^ ((u64)i * 0x9E3779B97F4A7C15ULL);
		const emitter_t *e;
		vec3_t o, d;
		double r;
		int    lo = 0, hi = nreal - 1, ev, lb;
#ifdef _OPENMP
		th = omp_get_thread_num();
#else
		th = 0;
#endif
		sm64(&s);
		/* 光源 */
		r = urand(&s) * ftot;
		while (lo < hi) {
			const int mid = (lo + hi) / 2;
			if (ecdf[mid] > r) hi = mid;
			else lo = mid + 1;
		}
		e = &p->emit[lo];
		o = e->pos;
		if (e->nu * e->nv > 1) {
			/* 面光源 : nu x nv の小光源 (決定論と同じ分け方) から 1 つ */
			const int iu = (int)(urand(&s) * e->nu), iv = (int)(urand(&s) * e->nv);
			const double du = ((((iu < e->nu) ? iu : (e->nu - 1)) + 0.5) / e->nu - 0.5) * e->wsize;
			const double dv = ((((iv < e->nv) ? iv : (e->nv - 1)) + 0.5) / e->nv - 0.5) * e->hsize;
			o = v_add(o, v_add(v_scale(e->ax, du), v_scale(e->ay, dv)));
		}
		/* 波長ビン */
		{
			const double *c = &scdf[(size_t)e->ispec * nl];
			const double rr = urand(&s) * c[nl - 1];
			int a = 0, b = nl - 1;
			while (a < b) {
				const int mid = (a + b) / 2;
				if (c[mid] > rr) b = mid;
				else a = mid + 1;
			}
			lb = a;
		}
		d = emit_dir(p, e, dmax[lo], &s);

		for (ev = 0; ev < MC_MAXEV; ev++) {
			hit_t h;
			double tc;
			int    js = -1;
			h.t = 1e300;
			h.imat = 0;
			h.n = v_make(0.0, 0.0, 1.0);
			hit_walls(p, o, d, &h);
			hit_occluders(p, o, d, teps, &h);
			hit_panels(p, o, d, teps, &h);
			tc = canopy_collide(p, o, d, h.t, &s, &js);
			if (tc >= 0.0) {
				/* 葉に止められた : ω の確率で等方に散乱、残りは吸収 */
				const pmat_t *lm = &p->mat[p->slab[js].imat];
				const vec3_t x = v_add(o, v_scale(d, tc));
				tally_targets(p, o, x, lb, &cnt[(size_t)th * ncnt], off);
				o = x;
				if (urand(&s) < (lm->rho[lb] + lm->tau[lb])) {
					d = sphere(&s);
					continue;
				}
				cab[((size_t)th * nl) + lb]++;
				break;
			}
			{
				const vec3_t x = v_add(o, v_scale(d, h.t));
				const pmat_t *m = &p->mat[h.imat];
				const double u = urand(&s);
				tally_targets(p, o, x, lb, &cnt[(size_t)th * ncnt], off);
				o = x;
				if (u < m->rho[lb]) {
					d = around(h.n, sqrt(urand(&s)), &s);            /* Lambert */
					continue;
				}
				if (u < m->rho[lb] + m->rhos) {
					d = v_sub(d, v_scale(h.n, 2.0 * v_dot(d, h.n)));  /* 鏡面 */
					continue;
				}
				if (u < m->rho[lb] + m->rhos + m->rhog) {
					/* 光沢 : 鏡面の向きのまわりの cos^n ローブ。面の下へ出たら吸収 */
					const vec3_t mr = v_sub(d, v_scale(h.n, 2.0 * v_dot(d, h.n)));
					const vec3_t g = around(mr, pow(urand(&s), 1.0 / (m->gexp + 1.0)), &s);
					if (v_dot(g, h.n) > 0.0) {
						d = g;
						continue;
					}
				}
				wab[((size_t)th * nl) + lb]++;
				break;
			}
		}
		if (ev >= MC_MAXEV) lst[((size_t)th * nl) + lb]++;
	}

	/* 集計 : 整数の和 -> 放射照度・放射束 */
	for (it = 0; it < p->ntarget; it++) {
		target_t *t = &p->target[it];
		const double acell = v_norm(v_cross(t->eu, t->ev)) / (t->nx * t->ny);
		const int    nc = t->nx * t->ny;
		long long    tot = 0;
		int          ic;
		for (ic = 0; ic < nc; ic++) {
			for (il = 0; il < nl; il++) {
				long long c = 0;
				for (th = 0; th < nth; th++) {
					c += cnt[((size_t)th * ncnt) + (size_t)off[it] + ((size_t)ic * nl) + il];
				}
				tot += c;
				t->E[((size_t)ic * nl) + il] = ((double)c * w) / acell;
			}
		}
		plog(p, "monte carlo : target \"%s\" %lld photons (%.1f per cell, cell 1-sigma ~ %.2g)\n",
			t->name, tot, (double)tot / nc, (tot > 0) ? (1.0 / sqrt((double)tot / nc)) : 0.0);
	}
	{
		long long a = 0, c = 0, l = 0;
		double *Ew = (double *)xcalloc((size_t)nl, sizeof(double));
		double *Ec = (double *)xcalloc((size_t)nl, sizeof(double));
		for (il = 0; il < nl; il++) {
			long long sa = 0, sc = 0;
			for (th = 0; th < nth; th++) {
				sa += wab[((size_t)th * nl) + il];
				sc += cab[((size_t)th * nl) + il];
				l  += lst[((size_t)th * nl) + il];
			}
			a += sa;
			c += sc;
			Ew[il] = (double)sa * w;
			Ec[il] = (double)sc * w;
		}
		p->absorb_wall = (double)a * w;
		p->absorb_canopy = (double)c * w;
		p->ppf_wall = calc_ppfd(p, Ew);
		p->canopy_abs = calc_ppfd(p, Ec);
		p->mc_lost = (double)l * w;
		free(Ew);
		free(Ec);
	}

	free(ecdf);
	free(dmax);
	free(scdf);
	free(off);
	free(cnt);
	free(wab);
	free(cab);
	free(lst);
}
