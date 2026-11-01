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

■ 斜めの板 (panel)
任意の向きの平行四辺形 (厚さ 0 の両面板)。表裏 2 面をパッチにし、遮蔽は
「線分と平行四辺形の交差」で判定する (外接箱では代用しない)。

■ 接触線への格子のそろえ (snap_cuts)
遮蔽物が壁や別の遮蔽物の面に接する置き方 (床に直に置いた箱、側壁に接する
棚板、棚板に載せた箱) では、面のパッチが「接触していて光が届かない部分」と
「キャビティに面した部分」をまたぐことがある。パッチの放射発散度は一様なので、
またいだパッチは接触部分のぶんまで放射したことになり、その光はどこへも
届かず消える (床に置いた箱で収支 -1.7%)。また側壁のパッチが断面いっぱいの
棚板をまたぐと、そのパッチ 1 枚を通して上下の段の間に光が漏れる。
そこで、面に接する遮蔽物の外形の線 (接触線) に面の格子線をそろえる。
最寄りの格子線を接触線へ動かし (半セル以内)、動かせなければ挿入する。
接触線が無い面の格子は従来と 1 ビットも変わらない。動かした / 挿入した
本数はログに出す。
*/

#include "ppfd.h"

/*
矩形面を nu × nv のパッチに分割して追加する。cu[0..nu] / cv[0..nv] は辺に
沿った分割点の割合 (等分なら k/n)。戻り値 = 次のパッチ番号
*/
static int add_face(patch_t *pt, int k, vec3_t org, vec3_t ue, vec3_t ve, vec3_t nrm,
                    int nu, int nv, const double *cu, const double *cv, int imat, int iface)
{
	int iu, iv;

	for (iv = 0; iv < nv; iv++) {
		for (iu = 0; iu < nu; iu++) {
			patch_t *q = &pt[k++];
			const double a0 = cu[iu];
			const double a1 = cu[iu + 1];
			const double b0 = cv[iv];
			const double b1 = cv[iv + 1];
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

static double vcomp(vec3_t v, int a)
{
	return (a == 0) ? v.x : ((a == 1) ? v.y : v.z);
}

/* 軸並行なベクトルの軸 (0..2) */
static int vaxis(vec3_t v)
{
	if (fabs(v.x) > 0.0) return 0;
	if (fabs(v.y) > 0.0) return 1;
	return 2;
}

/*
軸並行な矩形面 (原点 org、2 辺 ue, ve) に接する遮蔽物の外形の線を、ue / ve
方向の割合 (0, 1) として cu / cv に足す。接する = 遮蔽物の法線方向の範囲が
面の平面を含み (接するのも含む)、面内の範囲が面と重なるか接する。owner は
面の持ち主の遮蔽物 (自分自身は数えない。壁の面なら -1)。
*/
static void face_contacts(const ppfd_t *p, int owner, vec3_t org, vec3_t ue, vec3_t ve,
                          double *cu, int *ncu, double *cv, int *ncv, int *ndiag)
{
	const double tol = 1e-9 * (p->Lx + p->Ly + p->Lz);
	const int au = vaxis(ue), av = vaxis(ve), an = 3 - au - av;
	const double ou = vcomp(org, au), ov = vcomp(org, av), w = vcomp(org, an);
	const double lu = vcomp(ue, au), lv = vcomp(ve, av);
	int io, s;

	/*
	斜めの板 : 板と面の平面の交わり (板が面に接する辺、または面を横切る線) を
	求める。面内で軸に平行な線なら、その座標と両端を分割点に足す (板で仕切った
	部屋の壁のパッチが仕切りをまたがない)。斜めの線は格子をそろえられないので
	数えてログに出す (*ndiag)。
	*/
	for (io = 0; io < p->npanel; io++) {
		const panel_t *pn = &p->panel[io];
		vec3_t c[4], x[8];
		double sd[4], umin = 1e300, umax = -1e300, vmin = 1e300, vmax = -1e300;
		int    k, nx = 0;
		c[0] = pn->org;
		c[1] = v_add(pn->org, pn->eu);
		c[2] = v_add(c[1], pn->ev);
		c[3] = v_add(pn->org, pn->ev);
		for (k = 0; k < 4; k++) sd[k] = vcomp(c[k], an) - w;
		for (k = 0; k < 4; k++) {
			const int k2 = (k + 1) % 4;
			if (fabs(sd[k]) <= tol) {
				x[nx++] = c[k];
			}
			else if ((fabs(sd[k2]) > tol) && ((sd[k] > 0.0) != (sd[k2] > 0.0))) {
				x[nx++] = v_add(c[k], v_scale(v_sub(c[k2], c[k]), sd[k] / (sd[k] - sd[k2])));
			}
		}
		if (nx < 2) continue;
		for (k = 0; k < nx; k++) {
			const double u = vcomp(x[k], au), v = vcomp(x[k], av);
			if (u < umin) umin = u;
			if (u > umax) umax = u;
			if (v < vmin) vmin = v;
			if (v > vmax) vmax = v;
		}
		if ((umax - umin <= tol) && (vmax - vmin <= tol)) continue;      /* 1 点で触れるだけ */
		/* 面の矩形と重ならない交わりは関係ない */
		if ((umax < ou - tol) || (umin > ou + lu + tol) || (vmax < ov - tol) || (vmin > ov + lv + tol)) continue;
		if ((umax - umin > tol) && (vmax - vmin > tol)) {
			(*ndiag)++;
			continue;
		}
		for (s = 0; s < 2; s++) {
			const double fu = (((s == 0) ? umin : umax) - ou) / lu;
			const double fv = (((s == 0) ? vmin : vmax) - ov) / lv;
			if ((fu * lu > tol) && ((1.0 - fu) * lu > tol)) cu[(*ncu)++] = fu;
			if ((fv * lv > tol) && ((1.0 - fv) * lv > tol)) cv[(*ncv)++] = fv;
		}
	}

	for (io = 0; io < p->nocc; io++) {
		const occluder_t *o = &p->occ[io];
		const vec3_t lo = v_make(o->x0, o->y0, o->z0), hi = v_make(o->x1, o->y1, o->z1);
		const double lou = vcomp(lo, au), hiu = vcomp(hi, au);
		const double lov = vcomp(lo, av), hiv = vcomp(hi, av);
		if (io == owner) continue;
		if ((vcomp(lo, an) > w + tol) || (vcomp(hi, an) < w - tol)) continue;
		if ((hiu < ou - tol) || (lou > ou + lu + tol) || (hiv < ov - tol) || (lov > ov + lv + tol)) continue;
		for (s = 0; s < 2; s++) {
			const double fu = (((s == 0) ? lou : hiu) - ou) / lu;
			const double fv = (((s == 0) ? lov : hiv) - ov) / lv;
			if ((fu * lu > tol) && ((1.0 - fu) * lu > tol)) cu[(*ncu)++] = fu;
			if ((fv * lv > tol) && ((1.0 - fv) * lv > tol)) cv[(*ncv)++] = fv;
		}
	}
}

/*
n 等分の分割点 cut[0..n] (= k/n) を接触線 c[0..nc-1] にそろえる。最寄りの
内側の分割点が半セル以内でまだ動かしていなければそこへ動かし、そうでなければ
挿入する。すでに分割点に乗っている接触線は何もしない (等分の値のまま)。
cut の大きさは n + nc + 1 以上。戻り値 = 分割数。
*/
static int snap_cuts(double *cut, int n, const double *c, int nc, double tol, int *nmoved, int *nins)
{
	int *moved = (int *)xcalloc((size_t)(n + 1), sizeof(int));
	int  m = n + 1, k, i, j;

	for (k = 0; k <= n; k++) cut[k] = (double)k / n;
	for (i = 0; i < nc; i++) {
		int on = 0;
		for (j = 0; j < m; j++) {
			if (fabs(cut[j] - c[i]) <= tol) on = 1;
		}
		if (on) continue;
		k = (int)floor((c[i] * n) + 0.5);
		if ((k >= 1) && (k <= n - 1) && !moved[k] && (fabs(((double)k / n) - c[i]) <= (0.5 / n) + tol)) {
			cut[k] = c[i];
			moved[k] = 1;
			(*nmoved)++;
		}
		else {
			cut[m++] = c[i];
			(*nins)++;
		}
	}
	/* 挿入した点を並べ直す (件数は少ないので挿入ソート) */
	for (i = 1; i < m; i++) {
		const double t = cut[i];
		for (j = i - 1; (j >= 0) && (cut[j] > t); j--) cut[j + 1] = cut[j];
		cut[j + 1] = t;
	}
	free(moved);
	return m - 1;
}

void setup_patch(ppfd_t *p)
{
	/* 面 : 原点 org と 2 辺 ue, ve。ue × ve が内向き法線を向くよう並べる */
	vec3_t org[6], ue[6], ve[6], nrm[6];
	int    f, k, io, n, g, nface, cap;
	double *cuts, *cl;
	int    *nfu, *nfv;

	org[0] = v_make(0,     0,     0);      ue[0] = v_make(0, p->Ly, 0); ve[0] = v_make(0, 0, p->Lz); nrm[0] = v_make( 1, 0, 0);
	org[1] = v_make(p->Lx, 0,     0);      ue[1] = v_make(0, 0, p->Lz); ve[1] = v_make(0, p->Ly, 0); nrm[1] = v_make(-1, 0, 0);
	org[2] = v_make(0,     0,     0);      ue[2] = v_make(0, 0, p->Lz); ve[2] = v_make(p->Lx, 0, 0); nrm[2] = v_make(0,  1, 0);
	org[3] = v_make(0,     p->Ly, 0);      ue[3] = v_make(p->Lx, 0, 0); ve[3] = v_make(0, 0, p->Lz); nrm[3] = v_make(0, -1, 0);
	org[4] = v_make(0,     0,     0);      ue[4] = v_make(p->Lx, 0, 0); ve[4] = v_make(0, p->Ly, 0); nrm[4] = v_make(0, 0,  1);
	org[5] = v_make(0,     0,     p->Lz);  ue[5] = v_make(0, p->Ly, 0); ve[5] = v_make(p->Lx, 0, 0); nrm[5] = v_make(0, 0, -1);

	/* 面ごとの分割数 (wall の div。省略時は patchdiv) */
	n = 0;
	for (f = 0; f < 6; f++) {
		if (p->wdivu[f] < 1) p->wdivu[f] = p->ndivu;
		if (p->wdivv[f] < 1) p->wdivv[f] = p->ndivv;
	}
	for (io = 0; io < p->nocc; io++) {
		occluder_t *o = &p->occ[io];
		if (o->ndivu < 1) o->ndivu = p->ndivu;
		if (o->ndivv < 1) o->ndivv = p->ndivv;
	}

	/*
	面ごとの分割点 (壁 6 面 + 遮蔽物 6 面ずつ)。接触線にそろえる (snap_cuts)。
	面 g の分割点は cuts[g] の先頭 (nfu[g] + 1 個) と cuts[g] + cap の先頭。
	*/
	nface = 6 + (6 * p->nocc);
	cap = (2 * (p->nocc + p->npanel)) + 1;
	for (f = 0; f < 6; f++) {
		if (p->wdivu[f] + 1 > cap) cap = p->wdivu[f] + 1;
		if (p->wdivv[f] + 1 > cap) cap = p->wdivv[f] + 1;
	}
	for (io = 0; io < p->nocc; io++) {
		if (p->occ[io].ndivu + 1 > cap) cap = p->occ[io].ndivu + 1;
		if (p->occ[io].ndivv + 1 > cap) cap = p->occ[io].ndivv + 1;
	}
	cap += (2 * (p->nocc + p->npanel)) + 1;
	cuts = (double *)xmalloc((size_t)nface * 2 * (size_t)cap * sizeof(double));
	nfu = (int *)xcalloc((size_t)nface, sizeof(int));
	nfv = (int *)xcalloc((size_t)nface, sizeof(int));
	cl  = (double *)xmalloc((size_t)(2 * (2 * (p->nocc + p->npanel) + 1)) * sizeof(double));
	p->nsnap_moved = 0;
	p->nsnap_ins = 0;
	p->nsnap_diag = 0;
	for (g = 0; g < nface; g++) {
		const double tol = 1e-12;
		double *cu = &cuts[(size_t)g * 2 * cap];
		double *cv = cu + cap;
		double *clu = cl, *clv = cl + (2 * (p->nocc + p->npanel) + 1);
		int     ncu = 0, ncv = 0, du, dv, owner;
		vec3_t  fo, fu, fv, fn;
		if (g < 6) {
			fo = org[g]; fu = ue[g]; fv = ve[g];
			du = p->wdivu[g]; dv = p->wdivv[g];
			owner = -1;
		}
		else {
			const occluder_t *o;
			owner = (g - 6) / 6;
			o = &p->occ[owner];
			if (!occ_face_used(o, (g - 6) % 6)) continue;
			occ_face(o, (g - 6) % 6, &fo, &fu, &fv, &fn);
			du = o->ndivu; dv = o->ndivv;
		}
		face_contacts(p, owner, fo, fu, fv, clu, &ncu, clv, &ncv, &p->nsnap_diag);
		nfu[g] = snap_cuts(cu, du, clu, ncu, tol, &p->nsnap_moved, &p->nsnap_ins);
		nfv[g] = snap_cuts(cv, dv, clv, ncv, tol, &p->nsnap_moved, &p->nsnap_ins);
		n += nfu[g] * nfv[g];
	}
	free(cl);

	for (io = 0; io < p->npanel; io++) {
		panel_t *pn = &p->panel[io];
		if (pn->ndivu < 1) pn->ndivu = p->ndivu;
		if (pn->ndivv < 1) pn->ndivv = p->ndivv;
		n += 2 * pn->ndivu * pn->ndivv;        /* 表と裏 */
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
		const double *cu = &cuts[(size_t)f * 2 * cap];
		k = add_face(p->patch, k, org[f], ue[f], ve[f], nrm[f],
			nfu[f], nfv[f], cu, cu + cap, p->wallmat[f], f);
	}
	for (io = 0; io < p->nocc; io++) {
		const occluder_t *o = &p->occ[io];
		for (f = 0; f < 6; f++) {
			const int gf = 6 + (6 * io) + f;
			const double *cu = &cuts[(size_t)gf * 2 * cap];
			vec3_t fo, fu, fv, fn;
			if (!occ_face_used(o, f)) continue;
			occ_face(o, f, &fo, &fu, &fv, &fn);
			/* iface は「同一平面どうしは見えない」判定に使うので面ごとに一意 */
			k = add_face(p->patch, k, fo, fu, fv, fn,
				nfu[gf], nfv[gf], cu, cu + cap, o->imat, gf);
		}
	}
	for (io = 0; io < p->npanel; io++) {
		const panel_t *pn = &p->panel[io];
		const int base = 6 + (6 * p->nocc) + (2 * io);
		double *pu = (double *)xmalloc((size_t)(pn->ndivu + 1) * sizeof(double));
		double *pv = (double *)xmalloc((size_t)(pn->ndivv + 1) * sizeof(double));
		int     j;
		for (j = 0; j <= pn->ndivu; j++) pu[j] = (double)j / pn->ndivu;
		for (j = 0; j <= pn->ndivv; j++) pv[j] = (double)j / pn->ndivv;
		/* 表 : eu x ev = n。裏は 2 辺を入れ替えて ev x eu = -n (格子は表と同じ) */
		k = add_face(p->patch, k, pn->org, pn->eu, pn->ev, pn->n,
			pn->ndivu, pn->ndivv, pu, pv, pn->imat, base);
		k = add_face(p->patch, k, pn->org, pn->ev, pn->eu, v_scale(pn->n, -1.0),
			pn->ndivv, pn->ndivu, pv, pu, pn->imat, base + 1);
		free(pu);
		free(pv);
	}
	free(cuts);
	free(nfu);
	free(nfv);

	/*
	覆われたパッチ : 同一平面で向かい合う面 (壁・遮蔽物の面) の矩形の中に
	重心があるパッチ。格子を接触線にそろえてあるので、パッチは全部覆われるか
	全部覆われないかのどちらか。覆われたパッチには光が届かず、どこへも放射
	しない (直接光・形態係数・群落の面の交換から外す)。厚い箱なら光線が箱の
	内部を通るので遮られるが、厚さ 0 の板に載せた箱の底面や、床に寝かせた板の
	下の床は、光線の端が板の面に乗るので遮られない : 明示的に外す。
	*/
	p->pcov = (unsigned char *)xcalloc((size_t)p->npatch, 1);
	p->npcov = 0;
	for (k = 0; k < p->npatch; k++) {
		const patch_t *q = &p->patch[k];
		const double tol = 1e-9 * (p->Lx + p->Ly + p->Lz);
		const int an = (fabs(q->n.x) > 0.5) ? 0 : ((fabs(q->n.y) > 0.5) ? 1 : 2);
		const double sg = vcomp(q->n, an);
		const double w = vcomp(q->c, an);
		const int a1 = (an + 1) % 3, a2 = (an + 2) % 3;
		const double c1 = vcomp(q->c, a1), c2 = vcomp(q->c, a2);
		const int owner = ((q->iface >= 6) && (q->iface < 6 + (6 * p->nocc))) ? ((q->iface - 6) / 6) : -1;
		int cov = 0;
		if (q->iface >= 6 + (6 * p->nocc)) continue;     /* 斜めの板は覆わない・覆われない */
		/* 壁 : 内向き法線 -sg の面が平面 w にあるか (壁の矩形はチャンバの断面全体) */
		if ((sg < 0.0) && (fabs(w) <= tol)) cov = 1;
		if ((sg > 0.0) && (fabs(w - vcomp(v_make(p->Lx, p->Ly, p->Lz), an)) <= tol)) cov = 1;
		if (owner < 0) cov = 0;                          /* 壁のパッチは壁に覆われない */
		for (io = 0; (io < p->nocc) && !cov; io++) {
			const occluder_t *o = &p->occ[io];
			const vec3_t lo = v_make(o->x0, o->y0, o->z0), hi = v_make(o->x1, o->y1, o->z1);
			/* 外向き法線が -sg の面 : sg > 0 なら min 側の面 (lo)、sg < 0 なら max 側の面 */
			const double wf = (sg > 0.0) ? vcomp(lo, an) : vcomp(hi, an);
			if (io == owner) continue;
			if (fabs(wf - w) > tol) continue;
			if ((c1 > vcomp(lo, a1) + tol) && (c1 < vcomp(hi, a1) - tol)
			 && (c2 > vcomp(lo, a2) + tol) && (c2 < vcomp(hi, a2) - tol)) cov = 1;
		}
		if (cov) {
			p->pcov[k] = 1;
			p->npcov++;
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

/*
線分 a + t d (t ∈ (teps, 1-teps)) が平行四辺形 org + u eu + v ev
(u, v ∈ [0, 1]) を横切るか。N = eu x ev (正規化しない)。
平面に平行な線分 (面内を走るものを含む) は遮らない扱い : 厚さ 0 の板を
真横から見る光線なので、遮る面積は 0。
*/
static int panel_hit(vec3_t a, vec3_t d, vec3_t org, vec3_t eu, vec3_t ev, double teps)
{
	const double geps = 1e-12;
	const vec3_t N = v_cross(eu, ev);
	const double nn = v_dot(N, N);
	const double dn = v_dot(d, N);
	double t, u, v;
	vec3_t w;

	if (fabs(dn) <= (1e-15 * sqrt(nn) * (v_norm(d) + 1e-300))) return 0;
	t = v_dot(v_sub(org, a), N) / dn;
	if ((t <= teps) || (t >= (1.0 - teps))) return 0;
	w = v_sub(v_add(a, v_scale(d, t)), org);
	/* w = u eu + v ev を外積で解く (eu と ev は直交とは限らない) */
	u = v_dot(v_cross(w, ev), N) / nn;
	v = v_dot(v_cross(eu, w), N) / nn;
	return (u >= -geps) && (u <= (1.0 + geps)) && (v >= -geps) && (v <= (1.0 + geps));
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
	for (io = 0; io < p->npanel; io++) {
		const panel_t *pn = &p->panel[io];
		if (panel_hit(a, d, pn->org, pn->eu, pn->ev, teps)) return 1;
	}
	return 0;
}

/* 遮るものが 1 つでもあるか (棚板 / 斜めの板) */
int occ_any(const ppfd_t *p)
{
	return (p->nocc + p->npanel) > 0;
}

/* 複製 (k0, k1, k2) への点の写像 (鏡像の展開。spec_axis_map を軸ごとに) */
static vec3_t map_copy(const ppfd_t *p, vec3_t c, int k0, int k1, int k2)
{
	return v_make(spec_axis_map(c.x, p->Lx, k0), spec_axis_map(c.y, p->Ly, k1),
	              spec_axis_map(c.z, p->Lz, k2));
}


/*
鏡像への直線 a-b が遮蔽物に当たるか。n = 鏡像の展開指数
(NULL または全 0 なら実在の線分なので occ_blocked と同じ)。

鏡像への直線は箱型ビリヤードの展開なので、実経路は展開空間に並んだ
チャンバの複製を順に通る。各軸の座標は直線上で単調なので、複製番号も
軸ごとに 0 から n_a まで単調に進み、通りうる複製は添字の箱
[0..n0] x [0..n1] x [0..n2] に収まる。複製 k の中で実経路が遮蔽物 O に
当たる ⟺ 直線が O の複製 T_k(O) に当たる。T_k(O) は複製 k の内部に
あるので、直線が通らない複製の T_k(O) に誤って当たることはない。
軸並行の箱は軸ごとの鏡映で軸並行の箱に移るので、判定は slab 法のまま厳密。
*/
int occ_blocked_img(const ppfd_t *p, vec3_t a, vec3_t b, const int *n)
{
	const double teps = 1e-9;
	const vec3_t d = v_sub(b, a);
	int    m0, m1, m2, j0, j1, j2, io;

	if ((n == NULL) || ((n[0] == 0) && (n[1] == 0) && (n[2] == 0))) {
		return occ_blocked(p, a, b);
	}
	if (!occ_any(p)) return 0;

	m0 = (n[0] < 0) ? -n[0] : n[0];
	m1 = (n[1] < 0) ? -n[1] : n[1];
	m2 = (n[2] < 0) ? -n[2] : n[2];
	for (j0 = 0; j0 <= m0; j0++) {
		const int k0 = (n[0] < 0) ? -j0 : j0;
		for (j1 = 0; j1 <= m1; j1++) {
			const int k1 = (n[1] < 0) ? -j1 : j1;
			for (j2 = 0; j2 <= m2; j2++) {
				const int k2 = (n[2] < 0) ? -j2 : j2;
				for (io = 0; io < p->nocc; io++) {
					const occluder_t *o = &p->occ[io];
					const double xa = spec_axis_map(o->x0, p->Lx, k0), xb = spec_axis_map(o->x1, p->Lx, k0);
					const double ya = spec_axis_map(o->y0, p->Ly, k1), yb = spec_axis_map(o->y1, p->Ly, k1);
					const double za = spec_axis_map(o->z0, p->Lz, k2), zb = spec_axis_map(o->z1, p->Lz, k2);
					double t0 = teps, t1 = 1.0 - teps;
					if (!slab(a.x, d.x, (xa < xb) ? xa : xb, (xa < xb) ? xb : xa, &t0, &t1)) continue;
					if (!slab(a.y, d.y, (ya < yb) ? ya : yb, (ya < yb) ? yb : ya, &t0, &t1)) continue;
					if (!slab(a.z, d.z, (za < zb) ? za : zb, (za < zb) ? zb : za, &t0, &t1)) continue;
					return 1;
				}
				for (io = 0; io < p->npanel; io++) {
					/* 板の鏡像も平行四辺形 : 3 隅を写して辺を作り直す */
					const panel_t *pn = &p->panel[io];
					const vec3_t o0 = map_copy(p, pn->org, k0, k1, k2);
					const vec3_t o1 = map_copy(p, v_add(pn->org, pn->eu), k0, k1, k2);
					const vec3_t o3 = map_copy(p, v_add(pn->org, pn->ev), k0, k1, k2);
					if (panel_hit(a, d, o0, v_sub(o1, o0), v_sub(o3, o0), teps)) return 1;
				}
			}
		}
	}
	return 0;
}

/* 区間 [a0,a1] と [b0,b1] が重なるか (端で接するのも重なりに数える) */
static int span_overlap(double a0, double a1, double b0, double b1)
{
	const double e = 1e-12;
	return (a0 <= (b1 + e)) && (b0 <= (a1 + e));
}

/*
軸並行の箱 [lo, hi] がどれかの遮蔽物 (棚板 / 斜めの板。n が非 NULL なら
通りうる複製それぞれにあるその鏡像) と重なるか。

可視率の推定で「1 本も当たりえない」ことを幾何で先に確定させるための判定。
保守的 (重なっても実際には当たらないことがある) でよい : 重なると判定
された対は全光線で数えるので、当たらなければ厳密に 1 が返る。
*/
int occ_aabb_overlap(const ppfd_t *p, vec3_t lo, vec3_t hi, const int *n)
{
	int m0, m1, m2, j0, j1, j2, io;

	if (!occ_any(p)) return 0;

	m0 = (n == NULL) ? 0 : ((n[0] < 0) ? -n[0] : n[0]);
	m1 = (n == NULL) ? 0 : ((n[1] < 0) ? -n[1] : n[1]);
	m2 = (n == NULL) ? 0 : ((n[2] < 0) ? -n[2] : n[2]);
	for (j0 = 0; j0 <= m0; j0++) {
		const int k0 = ((n != NULL) && (n[0] < 0)) ? -j0 : j0;
		for (j1 = 0; j1 <= m1; j1++) {
			const int k1 = ((n != NULL) && (n[1] < 0)) ? -j1 : j1;
			for (j2 = 0; j2 <= m2; j2++) {
				const int k2 = ((n != NULL) && (n[2] < 0)) ? -j2 : j2;
				for (io = 0; io < p->nocc; io++) {
					const occluder_t *o = &p->occ[io];
					const vec3_t a = map_copy(p, v_make(o->x0, o->y0, o->z0), k0, k1, k2);
					const vec3_t b = map_copy(p, v_make(o->x1, o->y1, o->z1), k0, k1, k2);
					if (span_overlap(lo.x, hi.x, (a.x < b.x) ? a.x : b.x, (a.x < b.x) ? b.x : a.x)
					 && span_overlap(lo.y, hi.y, (a.y < b.y) ? a.y : b.y, (a.y < b.y) ? b.y : a.y)
					 && span_overlap(lo.z, hi.z, (a.z < b.z) ? a.z : b.z, (a.z < b.z) ? b.z : a.z)) {
						return 1;
					}
				}
				for (io = 0; io < p->npanel; io++) {
					const panel_t *pn = &p->panel[io];
					vec3_t c[4], plo, phi;
					int    v;
					c[0] = pn->org;
					c[1] = v_add(pn->org, pn->eu);
					c[2] = v_add(c[1], pn->ev);
					c[3] = v_add(pn->org, pn->ev);
					for (v = 0; v < 4; v++) c[v] = map_copy(p, c[v], k0, k1, k2);
					plo = c[0];
					phi = c[0];
					for (v = 1; v < 4; v++) {
						if (c[v].x < plo.x) plo.x = c[v].x;
						if (c[v].y < plo.y) plo.y = c[v].y;
						if (c[v].z < plo.z) plo.z = c[v].z;
						if (c[v].x > phi.x) phi.x = c[v].x;
						if (c[v].y > phi.y) phi.y = c[v].y;
						if (c[v].z > phi.z) phi.z = c[v].z;
					}
					if (span_overlap(lo.x, hi.x, plo.x, phi.x) && span_overlap(lo.y, hi.y, plo.y, phi.y)
					 && span_overlap(lo.z, hi.z, plo.z, phi.z)) {
						return 1;
					}
				}
			}
		}
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
		if (!t->sensor) {
			/* 水平で上向きの面も向きの情報を持たせておく (出力・群落内の扱い用) */
			t->org = v_make(t->x0, t->y0, t->z);
			t->eu  = v_make(t->x1 - t->x0, 0.0, 0.0);
			t->ev  = v_make(0.0, t->y1 - t->y0, 0.0);
			t->n   = v_make(0.0, 0.0, 1.0);
		}
		t->E = (double *)xcalloc((size_t)t->nx * t->ny * p->nlam, sizeof(double));
	}
}

/*
水平な平面多角形 (頂点 nv <= 8) を xy の矩形 [x0,x1] x [y0,y1] で切り取った
面積 (xy への射影)。Sutherland-Hodgman で 4 本の半平面に順に切る (頂点は
1 本ごとに高々 1 個増えるので 8 + 4 = 12 個で足りる)。群落の面がそのまま
棚板・床などに乗っているとき、面と重なる部分の割合を求めるのに使う。
*/
double poly_rect_area(const vec3_t *poly, int nv, double x0, double x1, double y0, double y1)
{
	double ax[16], ay[16], bx[16], by[16];
	int    na, nb, k, e, i;
	double a = 0.0;

	if ((nv < 3) || (nv > 8)) return 0.0;
	na = nv;
	for (i = 0; i < nv; i++) {
		ax[i] = poly[i].x;
		ay[i] = poly[i].y;
	}
	for (e = 0; e < 4; e++) {
		nb = 0;
		for (k = 0; k < na; k++) {
			const int k2 = (k + 1) % na;
			/* 半平面の内側なら d >= 0 (e = 0..3 : x >= x0, x <= x1, y >= y0, y <= y1) */
			const double d1 = (e == 0) ? (ax[k] - x0) : (e == 1) ? (x1 - ax[k])
			                : (e == 2) ? (ay[k] - y0) : (y1 - ay[k]);
			const double d2 = (e == 0) ? (ax[k2] - x0) : (e == 1) ? (x1 - ax[k2])
			                : (e == 2) ? (ay[k2] - y0) : (y1 - ay[k2]);
			if (d1 >= 0.0) {
				bx[nb] = ax[k];
				by[nb] = ay[k];
				nb++;
			}
			if (((d1 >= 0.0) && (d2 < 0.0)) || ((d1 < 0.0) && (d2 >= 0.0))) {
				const double t = d1 / (d1 - d2);
				bx[nb] = ax[k] + (t * (ax[k2] - ax[k]));
				by[nb] = ay[k] + (t * (ay[k2] - ay[k]));
				nb++;
			}
		}
		na = nb;
		for (k = 0; k < na; k++) {
			ax[k] = bx[k];
			ay[k] = by[k];
		}
		if (na < 3) return 0.0;
	}
	for (k = 0; k < na; k++) {
		const int k2 = (k + 1) % na;
		a += (ax[k] * ay[k2]) - (ax[k2] * ay[k]);
	}
	return 0.5 * fabs(a);
}

/*
鏡像へ引いた直線の z 座標を実経路の z 座標へ写す「折り返し」。

床・天井が鏡面のとき、鏡像への直線は箱型ビリヤードの展開そのものなので、
実経路の z は周期 P = 2Lz の三角波

    fold(z) = w        (w = z mod P,  w <= Lz)
            = P - w    (w > Lz)

で与えられる (床・天井で何回反射しても成り立つ)。両端が箱の内側にある
実在の線分では恒等写像になるため、鏡面が無い入力の結果は変わらない。
*/
double fold_z(const ppfd_t *p, double z)
{
	const double period = 2.0 * p->Lz;
	double w;

	if ((z >= 0.0) && (z <= p->Lz)) return z;
	if (period <= 0.0) return z;

	w = z - (period * floor(z / period));
	if (w < 0.0) w = 0.0;               /* 丸めで負に振れたときの保険 */
	if (w > period) w = period;
	return (w <= p->Lz) ? w : (period - w);
}

/* 区間 [a0,a1] と [b0,b1] の重なりの長さ */
static double overlap1d(double a0, double a1, double b0, double b1)
{
	const double lo = (a0 > b0) ? a0 : b0;
	const double hi = (a1 < b1) ? a1 : b1;

	return (hi > lo) ? (hi - lo) : 0.0;
}

/*
直線の z 範囲 [zlo,zhi] のうち fold(z) が群落スラブに入る z の測度。
fold の逆像は周期 P ごとに現れる 2 本の区間の和 :

    A_k = [zbot + kP, ztop + kP]      (上りの枝)
    B_k = [kP - ztop, kP - zbot]      (下りの枝)

zbot = 0 や ztop = Lz でも 2 本が接するだけなので二重計上にならない。
*/
static double fold_measure(const ppfd_t *p, const cslab_t *c, double zlo, double zhi)
{
	const double period = 2.0 * p->Lz;
	double sum = 0.0;
	int    k, k0, k1;

	if (period <= 0.0) return overlap1d(zlo, zhi, c->zbot, c->ztop);

	k0 = (int)floor((zlo - c->ztop) / period);
	k1 = (int)ceil((zhi + c->ztop) / period);
	for (k = k0; k <= k1; k++) {
		const double o = (double)k * period;
		sum += overlap1d(zlo, zhi, c->zbot + o, c->ztop + o);
		sum += overlap1d(zlo, zhi, o - c->ztop, o - c->zbot);
	}
	return sum;
}

/*
線分 a-b の群落内の経路長 [m]。

スラブが複数あるときは「先頭スラブ基準に規格化した等価経路長」
Σ_j w_j s_j (w_j = G_j a_j / (G_0 a_0)) を返す。消散係数 canopy_kext は
先頭スラブの k = G_0 a_0 √(1-ω) のままなので、その積は段ごとの諸元を
そのまま反映した Σ_j G_j a_j √(1-ω) s_j になる。w_0 は厳密に 1.0 なので
スラブが 1 枚だけの入力は従来と 1 ビットも変わらない。

床・天井が鏡面のとき、端点が箱の外にある線分は鏡像への直線なので
折り返した実経路で測る。鏡映は弧長を保つので「直線の z がスラブの逆像に
入っている割合 × 全長」が実経路の群落内長さと厳密に一致する。
水平範囲を限った群落 (栽培ベッド) は box_path で箱との交差長を測る。
*/
/* 線分 a + t d (t ∈ [0,1]、長さ len) と軸並行の箱 [lo, hi] の交差長 (slab 法) */
static double box_chord(vec3_t a, vec3_t d, double len, vec3_t lo, vec3_t hi)
{
	const double ac[3] = { a.x, a.y, a.z };
	const double dc[3] = { d.x, d.y, d.z };
	const double lc[3] = { lo.x, lo.y, lo.z };
	const double hc[3] = { hi.x, hi.y, hi.z };
	double t0 = 0.0, t1 = 1.0;
	int    i;

	for (i = 0; i < 3; i++) {
		if (fabs(dc[i]) < (EPS * (len + 1.0))) {
			/*
			軸に平行 : 範囲内なら全区間、外なら交わらない。範囲は半開区間
			[lo, hi) にする : 接して並べたベッドの共有面に乗る光線を両方の
			箱で数えると二重計上になるため (並べたベッド = 1 枚のスラブ)。
			*/
			if ((ac[i] < lc[i]) || (ac[i] >= hc[i])) return 0.0;
		}
		else {
			double ta = (lc[i] - ac[i]) / dc[i];
			double tb = (hc[i] - ac[i]) / dc[i];
			if (ta > tb) { const double t = ta; ta = tb; tb = t; }
			if (ta > t0) t0 = ta;
			if (tb < t1) t1 = tb;
			if (t1 <= t0) return 0.0;
		}
	}
	return (t1 - t0) * len;
}

/*
水平範囲を限った群落 (栽培ベッドの上だけ) の線分 a-b の経路長。

鏡像への直線は箱型ビリヤードの展開なので、実経路は展開空間に並んだ
チャンバの複製を順に通る (occ_blocked_img と同じ考え方)。複製 k の中の
実経路が群落の箱 C を通る長さ = 直線と C の複製 T_k(C) の交差長で、
鏡映は弧長を保つので総和が実経路の群落内長さに厳密に等しい。T_k(C) は
複製 k の内部にあるので、直線が通らない複製の T_k(C) に誤って当たる
ことはない。複製の番号は端点が入っている複製から決める (各軸の座標は
直線上で単調なので、その間の番号をすべて調べれば足りる)。鏡面の無い
軸は折り返さないので複製 0 だけ。
*/
static double box_path(const ppfd_t *p, const cslab_t *c, vec3_t a, vec3_t b, vec3_t d, double len)
{
	const double L[3] = { p->Lx, p->Ly, p->Lz };
	const double ac[3] = { a.x, a.y, a.z };
	const double bc[3] = { b.x, b.y, b.z };
	int    klo[3], khi[3], i, k0, k1, k2;
	double sum = 0.0;

	for (i = 0; i < 3; i++) {
		const int mir = (p->mat[p->wallmat[2 * i]].rhos > 0.0) ||
		                (p->mat[p->wallmat[(2 * i) + 1]].rhos > 0.0);
		if (mir && (L[i] > 0.0)) {
			const int ka = (int)floor(ac[i] / L[i]);
			const int kb = (int)floor(bc[i] / L[i]);
			klo[i] = (ka < kb) ? ka : kb;
			khi[i] = (ka < kb) ? kb : ka;
		}
		else {
			klo[i] = 0;
			khi[i] = 0;
		}
	}
	for (k0 = klo[0]; k0 <= khi[0]; k0++) {
		for (k1 = klo[1]; k1 <= khi[1]; k1++) {
			for (k2 = klo[2]; k2 <= khi[2]; k2++) {
				const vec3_t u = map_copy(p, v_make(c->x0, c->y0, c->zbot), k0, k1, k2);
				const vec3_t v = map_copy(p, v_make(c->x1, c->y1, c->ztop), k0, k1, k2);
				const vec3_t lo = v_make((u.x < v.x) ? u.x : v.x, (u.y < v.y) ? u.y : v.y,
				                         (u.z < v.z) ? u.z : v.z);
				const vec3_t hi = v_make((u.x < v.x) ? v.x : u.x, (u.y < v.y) ? v.y : u.y,
				                         (u.z < v.z) ? v.z : u.z);
				sum += box_chord(a, d, len, lo, hi);
			}
		}
	}
	return sum;
}

/*
線分 a + t (b - a) (t ∈ [0,1]) が軸並行の箱 [x0,x1]x[y0,y1]x[z0,z1] (と、
鏡面の軸ではその複製) の中にある区間 [t0, t1] を列挙する (box_path と同じ
複製の数え方)。範囲を限った群落 (栽培ベッド) の層ごとの預け入れを、光線に
沿った順に並べるために使う。戻り値 = 区間の数 (maxn で打ち切り)。
*/
int canopy_box_intervals(const ppfd_t *p, double x0, double x1, double y0, double y1, double z0, double z1,
                         vec3_t a, vec3_t b, double *t0, double *t1, int maxn)
{
	const double L[3] = { p->Lx, p->Ly, p->Lz };
	const double ac[3] = { a.x, a.y, a.z };
	const double bc[3] = { b.x, b.y, b.z };
	const vec3_t d = v_sub(b, a);
	const double len = v_norm(d);
	int    klo[3], khi[3], i, k0, k1, k2, n = 0;

	if (len <= 0.0) return 0;
	for (i = 0; i < 3; i++) {
		const int mir = (p->mat[p->wallmat[2 * i]].rhos > 0.0) ||
		                (p->mat[p->wallmat[(2 * i) + 1]].rhos > 0.0);
		if (mir && (L[i] > 0.0)) {
			const int ka = (int)floor(ac[i] / L[i]);
			const int kb = (int)floor(bc[i] / L[i]);
			klo[i] = (ka < kb) ? ka : kb;
			khi[i] = (ka < kb) ? kb : ka;
		}
		else {
			klo[i] = 0;
			khi[i] = 0;
		}
	}
	for (k0 = klo[0]; k0 <= khi[0]; k0++) {
		for (k1 = klo[1]; k1 <= khi[1]; k1++) {
			for (k2 = klo[2]; k2 <= khi[2]; k2++) {
				const vec3_t u = map_copy(p, v_make(x0, y0, z0), k0, k1, k2);
				const vec3_t v = map_copy(p, v_make(x1, y1, z1), k0, k1, k2);
				const double lo[3] = { (u.x < v.x) ? u.x : v.x, (u.y < v.y) ? u.y : v.y, (u.z < v.z) ? u.z : v.z };
				const double hi[3] = { (u.x < v.x) ? v.x : u.x, (u.y < v.y) ? v.y : u.y, (u.z < v.z) ? v.z : u.z };
				const double dc[3] = { d.x, d.y, d.z };
				double ta = 0.0, tb = 1.0;
				int    ok = 1;
				for (i = 0; (i < 3) && ok; i++) {
					if (fabs(dc[i]) < (EPS * (len + 1.0))) {
						/* 軸に平行 : 半開区間 [lo, hi) (box_chord と同じ) */
						if ((ac[i] < lo[i]) || (ac[i] >= hi[i])) ok = 0;
					}
					else {
						double s0 = (lo[i] - ac[i]) / dc[i];
						double s1 = (hi[i] - ac[i]) / dc[i];
						if (s0 > s1) { const double t = s0; s0 = s1; s1 = t; }
						if (s0 > ta) ta = s0;
						if (s1 < tb) tb = s1;
						if (tb <= ta) ok = 0;
					}
				}
				if (ok && (n < maxn)) {
					t0[n] = ta;
					t1[n] = tb;
					n++;
				}
			}
		}
	}
	return n;
}

double canopy_path(const ppfd_t *p, vec3_t a, vec3_t b)
{
	const vec3_t d = v_sub(b, a);
	const double len = v_norm(d);
	/*
	折り返すのは「床・天井が鏡面のときに箱の外へ出る線分」= 鏡像への直線
	だけ。z の鏡面が無ければ従来どおり直線で測る (鏡面を使わない入力の
	結果を変えないため。箱の外に置かれた光源等もそのまま扱う)。
	*/
	const int zmir = (p->mat[p->wallmat[4]].rhos > 0.0) ||
	                 (p->mat[p->wallmat[5]].rhos > 0.0);
	const int fold = zmir && !((a.z >= 0.0) && (a.z <= p->Lz) &&
	                           (b.z >= 0.0) && (b.z <= p->Lz));
	const int horiz = (fabs(d.z) < (EPS * (len + 1.0)));
	const double zfold = (fold && horiz) ? fold_z(p, a.z) : a.z;
	const double zlo = (a.z < b.z) ? a.z : b.z;
	const double zhi = (a.z < b.z) ? b.z : a.z;
	double sum = 0.0;
	int    j;

	if (!p->canopy.on || (len < EPS)) return 0.0;

	for (j = 0; j < p->nslab; j++) {
		const cslab_t *c = &p->slab[j];
		double t0, t1, lo, hi;

		if (!c->full) {
			sum += c->w * box_path(p, c, a, b, d, len);
			continue;
		}
		if (horiz) {
			/* ほぼ水平な線分 : 全体がスラブ内かどうか */
			if ((zfold >= c->zbot) && (zfold <= c->ztop)) sum += c->w * len;
			continue;
		}
		if (fold) {
			sum += c->w * ((fold_measure(p, c, zlo, zhi) / fabs(d.z)) * len);
			continue;
		}

		t0 = (c->zbot - a.z) / d.z;
		t1 = (c->ztop - a.z) / d.z;
		lo = (t0 < t1) ? t0 : t1;
		hi = (t0 < t1) ? t1 : t0;
		if (lo < 0.0) lo = 0.0;
		if (hi > 1.0) hi = 1.0;
		if (hi <= lo) continue;
		sum += c->w * ((hi - lo) * len);
	}
	return sum;
}

/*
z の範囲 [zlo, zhi] に収まる線分が群落を通りうるか。範囲が箱の外へ出る
(鏡像への直線で折り返しうる) ときは保守的に 1 を返す。0 なら、その範囲の
どの線分でも canopy_path は厳密に 0 になる (経路長の標本化を省くための
判定なので、閉区間で判定して境界に乗る線分も「通りうる」側に入れる)。
*/
int canopy_may_cross(const ppfd_t *p, double zlo, double zhi)
{
	int j;

	if (!p->canopy.on) return 0;
	if ((zlo < 0.0) || (zhi > p->Lz)) return 1;
	for (j = 0; j < p->nslab; j++) {
		if ((zhi >= p->slab[j].zbot) && (zlo <= p->slab[j].ztop)) return 1;
	}
	return 0;
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
