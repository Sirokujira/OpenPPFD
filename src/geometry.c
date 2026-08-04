/*
geometry.c

チャンバ壁面のパッチ分割、測定面グリッドの確保、群落 (Beer-Lambert
スラブ) の経路長と分光透過率。

■ 群落モデル
群落を z ∈ [zbot, ztop] の一様な葉群スラブとみなす。葉面積密度
a = LAI/(ztop-zbot) [1/m]、葉群投影係数 G = k0 (球状分布 0.5 / 水平葉 1.0)
のとき、経路長 s の線分に沿う光学的厚さは

    τ_opt = G * a * s * sqrt(1 - ω_λ)

ω_λ = ρ_leaf + τ_leaf は葉の散乱係数。sqrt(1-ω) は散乱葉に対する
Goudriaan の近似で、ω = 0 (黒体葉) では厳密な Beer-Lambert 則
exp(-G·LAI/cosθ) に一致する (鉛直入射 s = LAI/a のとき exp(-G·LAI))。

赤・青は ω が小さく (~0.15) 強く減衰し、緑 (~0.5) と遠赤色 (~0.9) は
深くまで透過する — 植物工場で群落下部のスペクトルが変わる主因。
*/

#include "ppfd.h"

void setup_patch(ppfd_t *p)
{
	/* 面 : 原点 org と 2 辺 ue, ve。ue × ve が内向き法線を向くよう並べる */
	vec3_t org[6], ue[6], ve[6], nrm[6];
	int    f, iu, iv, k;

	org[0] = v_make(0,     0,     0);      ue[0] = v_make(0, p->Ly, 0); ve[0] = v_make(0, 0, p->Lz); nrm[0] = v_make( 1, 0, 0);
	org[1] = v_make(p->Lx, 0,     0);      ue[1] = v_make(0, 0, p->Lz); ve[1] = v_make(0, p->Ly, 0); nrm[1] = v_make(-1, 0, 0);
	org[2] = v_make(0,     0,     0);      ue[2] = v_make(0, 0, p->Lz); ve[2] = v_make(p->Lx, 0, 0); nrm[2] = v_make(0,  1, 0);
	org[3] = v_make(0,     p->Ly, 0);      ue[3] = v_make(p->Lx, 0, 0); ve[3] = v_make(0, 0, p->Lz); nrm[3] = v_make(0, -1, 0);
	org[4] = v_make(0,     0,     0);      ue[4] = v_make(p->Lx, 0, 0); ve[4] = v_make(0, p->Ly, 0); nrm[4] = v_make(0, 0,  1);
	org[5] = v_make(0,     0,     p->Lz);  ue[5] = v_make(0, p->Ly, 0); ve[5] = v_make(p->Lx, 0, 0); nrm[5] = v_make(0, 0, -1);

	p->npatch = 6 * p->ndivu * p->ndivv;
	if (p->npatch > MAXPATCH) {
		fprintf(stderr, "*** too many patches : %d > %d (reduce patchdiv)\n", p->npatch, MAXPATCH);
		exit(1);
	}
	p->patch = (patch_t *)xcalloc((size_t)p->npatch, sizeof(patch_t));

	k = 0;
	for (f = 0; f < 6; f++) {
		for (iv = 0; iv < p->ndivv; iv++) {
			for (iu = 0; iu < p->ndivu; iu++) {
				patch_t *q = &p->patch[k++];
				const double a0 = (double)iu / p->ndivu;
				const double a1 = (double)(iu + 1) / p->ndivu;
				const double b0 = (double)iv / p->ndivv;
				const double b1 = (double)(iv + 1) / p->ndivv;
				q->p[0] = v_add(org[f], v_add(v_scale(ue[f], a0), v_scale(ve[f], b0)));
				q->p[1] = v_add(org[f], v_add(v_scale(ue[f], a1), v_scale(ve[f], b0)));
				q->p[2] = v_add(org[f], v_add(v_scale(ue[f], a1), v_scale(ve[f], b1)));
				q->p[3] = v_add(org[f], v_add(v_scale(ue[f], a0), v_scale(ve[f], b1)));
				q->n = nrm[f];
				q->c = v_scale(v_add(v_add(q->p[0], q->p[1]), v_add(q->p[2], q->p[3])), 0.25);
				q->area = v_norm(v_cross(v_sub(q->p[1], q->p[0]), v_sub(q->p[3], q->p[0])));
				q->imat = p->wallmat[f];
				q->iface = f;
			}
		}
	}
}

void setup_target(ppfd_t *p)
{
	int it;

	for (it = 0; it < p->ntarget; it++) {
		target_t *t = &p->target[it];
		if (t->x1 < t->x0) {
			t->x0 = 0.0;
			t->x1 = p->Lx;
		}
		if (t->y1 < t->y0) {
			t->y0 = 0.0;
			t->y1 = p->Ly;
		}
		t->E = (double *)xcalloc((size_t)t->nx * t->ny * p->nlam, sizeof(double));
	}
}

/* 線分 a-b のうち群落スラブ内にある部分の長さ [m] */
double canopy_path(const ppfd_t *p, vec3_t a, vec3_t b)
{
	const canopy_t *cp = &p->canopy;
	const vec3_t d = v_sub(b, a);
	const double len = v_norm(d);
	double t0, t1, lo, hi;

	if (!cp->on || (len < EPS)) return 0.0;

	if (fabs(d.z) < (EPS * (len + 1.0))) {
		/* ほぼ水平な線分 : 全体がスラブ内かどうか */
		return ((a.z >= cp->zbot) && (a.z <= cp->ztop)) ? len : 0.0;
	}

	t0 = (cp->zbot - a.z) / d.z;
	t1 = (cp->ztop - a.z) / d.z;
	lo = (t0 < t1) ? t0 : t1;
	hi = (t0 < t1) ? t1 : t0;
	if (lo < 0.0) lo = 0.0;
	if (hi > 1.0) hi = 1.0;
	if (hi <= lo) return 0.0;

	return (hi - lo) * len;
}

/* 経路長 s [m] に対する群落の分光透過率 tr[0..nlam-1] */
void canopy_trans(const ppfd_t *p, double s, double *tr)
{
	const canopy_t *cp = &p->canopy;
	int i;

	if (!cp->on || (s <= 0.0)) {
		for (i = 0; i < p->nlam; i++) tr[i] = 1.0;
		return;
	}

	for (i = 0; i < p->nlam; i++) {
		double omega = p->mat[cp->imat].rho[i] + p->mat[cp->imat].tau[i];
		double kext;
		if (omega < 0.0) omega = 0.0;
		if (omega > 0.999) omega = 0.999;
		kext = cp->k0 * sqrt(1.0 - omega) * cp->a;
		tr[i] = exp(-kext * s);
	}
}
