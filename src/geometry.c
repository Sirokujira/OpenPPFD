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

■ 遮蔽物 (棚板・トレイ・バッフル)
軸並行直方体として与える。6 面 (厚さ 0 なら 2 面) をパッチに分割して
ラジオシティに参加させ、線分との交差判定 (スラブ法) で影を落とす。
チャンバ + 遮蔽物は閉じた面の集合なので、形態係数の行和 = 1 と
「全放射束が表面に吸収される」はそのまま成り立つ。

面が壁と同一平面で接する置き方 (床の上に厚さのある箱を直に置く等) は
避けること。同一平面のパッチどうしは互いに見えないものとして落とすので、
その隙間に光が消える。棚板は厚さ 0 の板として与えるのが確実。
*/

#include "ppfd.h"

/* 矩形面を nu × nv のパッチに分割して追加する。戻り値 = 次のパッチ番号 */
static int add_face(patch_t *pt, int k, vec3_t org, vec3_t ue, vec3_t ve, vec3_t nrm,
                    int nu, int nv, int imat, int iface)
{
	int iu, iv;

	for (iv = 0; iv < nv; iv++) {
		for (iu = 0; iu < nu; iu++) {
			patch_t *q = &pt[k++];
			const double a0 = (double)iu / nu;
			const double a1 = (double)(iu + 1) / nu;
			const double b0 = (double)iv / nv;
			const double b1 = (double)(iv + 1) / nv;
			q->p[0] = v_add(org, v_add(v_scale(ue, a0), v_scale(ve, b0)));
			q->p[1] = v_add(org, v_add(v_scale(ue, a1), v_scale(ve, b0)));
			q->p[2] = v_add(org, v_add(v_scale(ue, a1), v_scale(ve, b1)));
			q->p[3] = v_add(org, v_add(v_scale(ue, a0), v_scale(ve, b1)));
			q->n = nrm;
			q->c = v_scale(v_add(v_add(q->p[0], q->p[1]), v_add(q->p[2], q->p[3])), 0.25);
			q->area = v_norm(v_cross(v_sub(q->p[1], q->p[0]), v_sub(q->p[3], q->p[0])));
			q->imat = imat;
			q->iface = iface;
		}
	}
	return k;
}

/*
遮蔽物の面 f (0:xmin 1:xmax 2:ymin 3:ymax 4:zmin 5:zmax) の原点と 2 辺。
法線は外向き (キャビティ側)。ue × ve が法線を向くよう並べる。
厚さ 0 の板では退化する 4 面 (面積 0) を呼び出し側で落とす。
*/
static void occ_face(const occluder_t *o, int f, vec3_t *org, vec3_t *ue, vec3_t *ve, vec3_t *nrm)
{
	const double dx = o->x1 - o->x0;
	const double dy = o->y1 - o->y0;
	const double dz = o->z1 - o->z0;

	switch (f) {
	case 0:
		*org = v_make(o->x0, o->y0, o->z0); *ue = v_make(0, 0, dz); *ve = v_make(0, dy, 0);
		*nrm = v_make(-1, 0, 0); break;
	case 1:
		*org = v_make(o->x1, o->y0, o->z0); *ue = v_make(0, dy, 0); *ve = v_make(0, 0, dz);
		*nrm = v_make(1, 0, 0); break;
	case 2:
		*org = v_make(o->x0, o->y0, o->z0); *ue = v_make(dx, 0, 0); *ve = v_make(0, 0, dz);
		*nrm = v_make(0, -1, 0); break;
	case 3:
		*org = v_make(o->x0, o->y1, o->z0); *ue = v_make(0, 0, dz); *ve = v_make(dx, 0, 0);
		*nrm = v_make(0, 1, 0); break;
	case 4:
		*org = v_make(o->x0, o->y0, o->z0); *ue = v_make(0, dy, 0); *ve = v_make(dx, 0, 0);
		*nrm = v_make(0, 0, -1); break;
	default:
		*org = v_make(o->x0, o->y0, o->z1); *ue = v_make(dx, 0, 0); *ve = v_make(0, dy, 0);
		*nrm = v_make(0, 0, 1); break;
	}
}

/* 遮蔽物 io の面 f を生成するか (厚さ 0 の板では 2 面だけになる) */
static int occ_face_used(const occluder_t *o, int f)
{
	const double dx = o->x1 - o->x0;
	const double dy = o->y1 - o->y0;
	const double dz = o->z1 - o->z0;
	const double e = 1e-12;

	if (dx <= e) return (f < 2);
	if (dy <= e) return ((f >= 2) && (f < 4));
	if (dz <= e) return (f >= 4);
	return 1;
}

void setup_patch(ppfd_t *p)
{
	/* 面 : 原点 org と 2 辺 ue, ve。ue × ve が内向き法線を向くよう並べる */
	vec3_t org[6], ue[6], ve[6], nrm[6];
	int    f, k, io, n;

	org[0] = v_make(0,     0,     0);      ue[0] = v_make(0, p->Ly, 0); ve[0] = v_make(0, 0, p->Lz); nrm[0] = v_make( 1, 0, 0);
	org[1] = v_make(p->Lx, 0,     0);      ue[1] = v_make(0, 0, p->Lz); ve[1] = v_make(0, p->Ly, 0); nrm[1] = v_make(-1, 0, 0);
	org[2] = v_make(0,     0,     0);      ue[2] = v_make(0, 0, p->Lz); ve[2] = v_make(p->Lx, 0, 0); nrm[2] = v_make(0,  1, 0);
	org[3] = v_make(0,     p->Ly, 0);      ue[3] = v_make(p->Lx, 0, 0); ve[3] = v_make(0, 0, p->Lz); nrm[3] = v_make(0, -1, 0);
	org[4] = v_make(0,     0,     0);      ue[4] = v_make(p->Lx, 0, 0); ve[4] = v_make(0, p->Ly, 0); nrm[4] = v_make(0, 0,  1);
	org[5] = v_make(0,     0,     p->Lz);  ue[5] = v_make(0, p->Ly, 0); ve[5] = v_make(p->Lx, 0, 0); nrm[5] = v_make(0, 0, -1);

	n = 6 * p->ndivu * p->ndivv;
	for (io = 0; io < p->nocc; io++) {
		occluder_t *o = &p->occ[io];
		if (o->ndivu < 1) o->ndivu = p->ndivu;
		if (o->ndivv < 1) o->ndivv = p->ndivv;
		for (f = 0; f < 6; f++) {
			if (occ_face_used(o, f)) n += o->ndivu * o->ndivv;
		}
	}

	p->npatch = n;
	if (p->npatch > MAXPATCH) {
		fprintf(stderr, "*** too many patches : %d > %d (reduce patchdiv / occluder div)\n",
			p->npatch, MAXPATCH);
		exit(1);
	}
	p->patch = (patch_t *)xcalloc((size_t)p->npatch, sizeof(patch_t));

	k = 0;
	for (f = 0; f < 6; f++) {
		k = add_face(p->patch, k, org[f], ue[f], ve[f], nrm[f],
			p->ndivu, p->ndivv, p->wallmat[f], f);
	}
	for (io = 0; io < p->nocc; io++) {
		const occluder_t *o = &p->occ[io];
		for (f = 0; f < 6; f++) {
			vec3_t fo, fu, fv, fn;
			if (!occ_face_used(o, f)) continue;
			occ_face(o, f, &fo, &fu, &fv, &fn);
			/* iface は「同一平面どうしは見えない」判定に使うので面ごとに一意 */
			k = add_face(p->patch, k, fo, fu, fv, fn,
				o->ndivu, o->ndivv, o->imat, 6 + (6 * io) + f);
		}
	}
}

/*
線分 a-b が遮蔽物に遮られるか (スラブ法)。端点そのものは当たり判定から
外す (パッチ自身の面から出る光線が自分を遮ると判定されないように)。
厚さ 0 の板 (lo == hi) もそのまま扱える。
*/
static int slab(double a, double d, double lo, double hi, double *t0, double *t1)
{
	const double geps = 1e-12;
	double ta, tb;

	if (fabs(d) < 1e-15) {
		/* 軸に平行 : スラブ内を通っていれば他軸の判定に委ねる */
		return (a >= (lo - geps)) && (a <= (hi + geps));
	}
	ta = (lo - a) / d;
	tb = (hi - a) / d;
	if (ta > tb) {
		const double t = ta;
		ta = tb;
		tb = t;
	}
	if (ta > *t0) *t0 = ta;
	if (tb < *t1) *t1 = tb;
	return (*t0 <= *t1);
}

int occ_blocked(const ppfd_t *p, vec3_t a, vec3_t b)
{
	const double teps = 1e-9;
	const vec3_t d = v_sub(b, a);
	int    io;

	for (io = 0; io < p->nocc; io++) {
		const occluder_t *o = &p->occ[io];
		double t0 = teps, t1 = 1.0 - teps;
		if (!slab(a.x, d.x, o->x0, o->x1, &t0, &t1)) continue;
		if (!slab(a.y, d.y, o->y0, o->y1, &t0, &t1)) continue;
		if (!slab(a.z, d.z, o->z0, o->z1, &t0, &t1)) continue;
		return 1;
	}
	return 0;
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

	(void)cp;
	for (i = 0; i < p->nlam; i++) {
		tr[i] = exp(-canopy_kext(p, i) * s);
	}
}
