/*
direct.c

光源から壁パッチ / 測定点への直接放射照度。

■ 配光
放射束 Φ の光源の放射強度は
    I(θ) = Φ (m+1)/(2π) cos^m θ      (θ = 配光軸からの角度、上半空間のみ)
m = 1 が Lambert 配光 I = (Φ/π) cosθ、m = 0 が半球一様、m < 0 は等方
    I = Φ/(4π)
を表す。∫I dΩ = Φ は m によらず満たされる。

実測配光ファイル (IES LM-63 / EULUMDAT) を与えた場合は
    I(γ, C) = Φ Î(γ, C)          ∫ Î dΩ = 1
の Î を photometry.c から引く。こちらも放射束は厳密に保存される。

■ 面光源
寸法 w × h の平面 Lambert 光源は nu × nv 個の点光源に分割し、各々に
Φ/(nu·nv) を与える。分割数 → ∞ で軸上放射照度は解析解
    E = Φ F_{dA→A_s} / A_s
に収束する (data/sample/panel.ppfd の検証ケース)。

■ 群落
光源から受光点までの線分が群落スラブを横切る長さ s を求め、分光透過率
exp(-G a s sqrt(1-ω_λ)) を掛ける (geometry.c)。

■ 遮蔽物
光源と受光点を結ぶ線分が遮蔽物 (棚板等) に当たればその寄与を落とす。
直接光の遮蔽は 2 値なので近似は入らない。
*/

#include "ppfd.h"

/*
点 x (単位法線 n) における全光源からの直接分光放射照度 [W/m²] を
E[0..nlam-1] に加算する (E は呼び出し側でゼロ初期化しておく)。
*/
void direct_point_ex(const ppfd_t *p, vec3_t x, vec3_t n, double *E, double *dep, double *seg)
{
	int    ie, i;
	double *acc = NULL;
	double  tr[512];
	const int nl = p->canopy.nlayer;
	const int use_acc = (!p->canopy.on) && (p->nspec <= 4096);

	if (use_acc) {
		acc = (double *)xcalloc((size_t)p->nspec, sizeof(double));
	}

	for (ie = 0; ie < p->nemit; ie++) {
		const emitter_t *e = &p->emit[ie];
		const int nsub = e->nu * e->nv;
		const double fsub = e->flux / nsub;
		int iu, iv;

		for (iv = 0; iv < e->nv; iv++) {
			for (iu = 0; iu < e->nu; iu++) {
				vec3_t sp = e->pos;
				vec3_t r;
				double d, d2, cr, cs, inten, es;

				if (nsub > 1) {
					const double du = ((iu + 0.5) / e->nu - 0.5) * e->wsize;
					const double dv = ((iv + 0.5) / e->nv - 0.5) * e->hsize;
					sp = v_add(sp, v_add(v_scale(e->ax, du), v_scale(e->ay, dv)));
				}

				r = v_sub(x, sp);
				d2 = v_dot(r, r);
				if (d2 < 1e-18) continue;
				d = sqrt(d2);

				cr = -v_dot(n, r) / d;          /* 受光面が光源を向く向き */
				if (cr <= 0.0) continue;
				if (occ_blocked(p, sp, x)) continue;   /* 棚板等の影 */

				if (e->idist >= 0) {
					/* 実測配光 : Î は ∫Î dΩ = 1 に正規化済みなので I = Φ Î */
					const vec3_t u = v_scale(r, 1.0 / d);
					double cg = v_dot(u, e->dir);
					double ca;
					if (cg > 1.0) cg = 1.0;
					if (cg < -1.0) cg = -1.0;
					ca = atan2(v_dot(u, e->cy), v_dot(u, e->cx)) * (180.0 / PI) - e->crot;
					inten = fsub * photdist_value(&p->dist[e->idist], acos(cg) * (180.0 / PI), ca);
					if (inten <= 0.0) continue;
				}
				else if (e->mexp < 0.0) {
					inten = fsub / (4.0 * PI);
				}
				else {
					cs = v_dot(e->dir, r) / d;
					if (cs <= 0.0) continue;
					inten = fsub * (e->mexp + 1.0) / (2.0 * PI) * pow(cs, e->mexp);
				}

				es = inten * cr / d2;

				if (dep != NULL) {
					/* この光線が層へ預ける分 (幾何係数は波長によらない) */
					canopy_deposit(p, sp, x, es, &dep[(size_t)e->ispec * nl], seg);
				}

				if (use_acc) {
					acc[e->ispec] += es;
				}
				else {
					const spec_t *sc = &p->spec[e->ispec];
					const double s = canopy_path(p, sp, x);
					canopy_trans(p, s, tr);
					for (i = sc->i0; i <= sc->i1; i++) {
						E[i] += es * sc->w[i] * tr[i];
					}
				}
			}
		}
	}

	if (use_acc) {
		int is;
		for (is = 0; is < p->nspec; is++) {
			if (acc[is] == 0.0) continue;
			for (i = p->spec[is].i0; i <= p->spec[is].i1; i++) {
				E[i] += acc[is] * p->spec[is].w[i];
			}
		}
		free(acc);
	}
}

void direct_point(const ppfd_t *p, vec3_t x, vec3_t n, double *E)
{
	direct_point_ex(p, x, n, E, NULL, NULL);
}

/*
各パッチの面積平均直接放射照度 Ed[ip*nlam + i]。

複合 3 点 Gauss (msub x msub の小矩形 x 3x3 点) で面積平均する。計算量は
nemit * npatch * (3*msub)^2 * (有効ビン数) に比例するので、光源が多い
ケースでは msub を落とさないと現実的な時間に収まらない。既定は光源数
から自動で決め、実際に使った値をログに出す (quadrature キーで上書き可)。

光源が少ないほど被積分関数が急峻なので、この既定は精度を落とさない
方向に働く : 検証ケース (光源 1 個) は常に msub = 3 になる。

遮蔽物があると影の境界で被積分関数が不連続になり、複合 Gauss の収束が
O(1/msub) まで落ちる (エネルギー収支の closure error に効く)。境界が
またぐパッチだけ msub を 8 倍 (上限 32) に上げ、細分したパッチ数を
ログに出す。それでも収支は 1e-4 程度が限界で、これは形態係数ではなく
直接光のパッチ求積の限界。
*/

/*
影の境界がまたぐパッチか (光源からパッチの 5 点への遮蔽判定が割れるか)。
面光源は中心で代表させる (細分の要否の判定に使うだけなので粗くてよい)。

判定点は角そのものではなく少し内側に取る。棚板の面にちょうど接する壁
パッチは、角が遮蔽物の平面上に乗るせいで「角だけ遮られない」と出るが、
内部は一様に影なので細分しても無駄になる (2 段ラックでは側壁の 64 枚が
これに当たり、直接光の計算時間が数倍になっていた)。
*/
static int patch_on_shadow_edge(const ppfd_t *p, const patch_t *q)
{
	const double t = 0.05;          /* 内側への寄せ量 (パッチ辺の割合) */
	const vec3_t eu = v_sub(q->p[1], q->p[0]);
	const vec3_t ev = v_sub(q->p[3], q->p[0]);
	vec3_t xs[5];
	int    ie, k;

	if (p->nocc == 0) return 0;

	for (k = 0; k < 4; k++) {
		const double u = ((k == 0) || (k == 3)) ? t : (1.0 - t);
		const double v = (k < 2) ? t : (1.0 - t);
		xs[k] = v_add(q->p[0], v_add(v_scale(eu, u), v_scale(ev, v)));
	}
	xs[4] = q->c;

	for (ie = 0; ie < p->nemit; ie++) {
		int nb = 0;
		for (k = 0; k < 5; k++) {
			nb += occ_blocked(p, p->emit[ie].pos, xs[k]);
		}
		if ((nb > 0) && (nb < 5)) return 1;
	}
	return 0;
}

void direct_patch(ppfd_t *p)
{
	const int n = p->npatch;
	const int nl = p->canopy.nlayer;
	const int ns = p->nspec;
	double *cd = NULL;
	int ip;
	int msub = p->msub;
	int nedge = 0;

	if (msub <= 0) {
		msub = (p->nemit <= 8) ? 3 : ((p->nemit <= 64) ? 2 : 1);
	}
	plog(p, "direct quadrature : msub = %d (%d points/patch)\n", msub, 9 * msub * msub);

	p->Ed = (double *)xcalloc((size_t)n * p->nlam, sizeof(double));
	if (p->canopy.scatter) {
		/* パッチごとに分けて溜め、あとで固定順に足す (スレッド数不変のため) */
		cd = (double *)xcalloc((size_t)n * ns * nl, sizeof(double));
	}

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) reduction(+:nedge)
#endif
	for (ip = 0; ip < n; ip++) {
		const patch_t *q = &p->patch[ip];
		const vec3_t eu = v_sub(q->p[1], q->p[0]);
		const vec3_t ev = v_sub(q->p[3], q->p[0]);
		double *E = &p->Ed[(size_t)ip * p->nlam];
		double *tmp = (double *)xcalloc((size_t)p->nlam, sizeof(double));
		const int edge = patch_on_shadow_edge(p, q);
		const int msub_i = edge ? ((msub * 8 < 32) ? (msub * 8) : 32) : msub;
		const double h = 1.0 / msub_i;
		double *dp = (cd != NULL) ? &cd[(size_t)ip * ns * nl] : NULL;
		double *dtmp = (cd != NULL) ? (double *)xcalloc((size_t)ns * nl, sizeof(double)) : NULL;
		double *seg = (cd != NULL) ? (double *)xmalloc((size_t)nl * sizeof(double)) : NULL;
		int su, sv, gu, gv, i;

		nedge += edge;
		/*
		群落散乱が有効なら、このパッチへ向かう光束が層へ預ける分を集める。
		パッチ i に届くはずの放射束は A_i * (減衰前の面積平均放射照度) なので、
		面積平均の求積重みをそのまま使い、最後に A_i を掛ければよい。
		Σ_m 預け入れ + 透過 = 1 が望遠鏡和で厳密に成り立つ。
		*/
		for (sv = 0; sv < msub_i; sv++) {
			for (su = 0; su < msub_i; su++) {
				for (gv = 0; gv < 3; gv++) {
					for (gu = 0; gu < 3; gu++) {
						static const double gx[3] = {
							0.5 * (1.0 - 0.7745966692414834), 0.5, 0.5 * (1.0 + 0.7745966692414834)
						};
						static const double gw[3] = {
							0.5 * 0.5555555555555556, 0.5 * 0.8888888888888888, 0.5 * 0.5555555555555556
						};
						const double u = (su + gx[gu]) * h;
						const double v = (sv + gx[gv]) * h;
						const double w = gw[gu] * gw[gv] * h * h;
						const vec3_t x = v_add(q->p[0], v_add(v_scale(eu, u), v_scale(ev, v)));

						for (i = 0; i < p->nlam; i++) tmp[i] = 0.0;
						if (dtmp != NULL) {
							for (i = 0; i < ns * nl; i++) dtmp[i] = 0.0;
						}
						direct_point_ex(p, x, q->n, tmp, dtmp, seg);
						for (i = 0; i < p->nlam; i++) E[i] += w * tmp[i];
						if (dtmp != NULL) {
							for (i = 0; i < ns * nl; i++) dp[i] += w * dtmp[i];
						}
					}
				}
			}
		}
		if (dp != NULL) {
			/* 面積平均 -> 放射束 [W] */
			for (i = 0; i < ns * nl; i++) dp[i] *= q->area;
			free(dtmp);
			free(seg);
		}
		free(tmp);
	}

	if (cd != NULL) {
		/* 波長へ展開する (加算順序はパッチ -> スペクトル -> 層で固定) */
		int is, m, il;
		for (ip = 0; ip < n; ip++) {
			for (is = 0; is < ns; is++) {
				const spec_t *sc = &p->spec[is];
				for (m = 0; m < nl; m++) {
					const double d = cd[(((size_t)ip * ns) + is) * nl + m];
					if (d == 0.0) continue;
					for (il = sc->i0; il <= sc->i1; il++) {
						p->cdep[((size_t)m * p->nlam) + il] += d * sc->w[il];
					}
				}
			}
		}
		free(cd);
	}

	if (nedge > 0) {
		plog(p, "  shadow edges : %d patches refined to msub = %d\n",
			nedge, (msub * 8 < 32) ? (msub * 8) : 32);
	}
}
