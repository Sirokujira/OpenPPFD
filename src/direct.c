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
■ 鏡面反射 (鏡像法)
チャンバは軸並行の直方体なので、鏡面反射する壁は「光源の鏡像」と等価に
なる。壁 f の鏡面反射率を rho_s とすると、f に対する鏡像光源は放射束
rho_s * Phi を持ち、その直線が折り返し経路をそのまま表す。直接光に
ついてはこれで厳密 (近似は段数の打ち切りだけ)。

鏡像は幅優先で作る。直前に折り返した面へは戻らない (直方体では同じ面で
2 回続けて折り返すと元に戻るため)。段数は specbounce で決め、打ち切りで
捨てた放射束をログに出す。

壁の吸収は (1 - rho_diffuse - rho_specular) * Einc になる。鏡面成分は
鏡像が運ぶので、収支はこれで閉じる。

【制限】鏡像からの直線は折り返した経路を表すので、群落や遮蔽物があると
経路長も遮蔽判定も合わない。両立は入力段で明示的に弾く。
*/

/* 面 f (0:xmin 1:xmax 2:ymin 3:ymax 4:zmin 5:zmax) の平面座標と軸 */
static void face_plane(const ppfd_t *p, int f, int *axis, double *val)
{
	*axis = f / 2;
	switch (f) {
	case 0: *val = 0.0;    break;
	case 1: *val = p->Lx;  break;
	case 2: *val = 0.0;    break;
	case 3: *val = p->Ly;  break;
	case 4: *val = 0.0;    break;
	default: *val = p->Lz; break;
	}
}

/* 軸 axis の座標 val の平面で点を鏡映する */
static vec3_t mirror_point(vec3_t a, int axis, double val)
{
	if (axis == 0) a.x = (2.0 * val) - a.x;
	else if (axis == 1) a.y = (2.0 * val) - a.y;
	else a.z = (2.0 * val) - a.z;
	return a;
}

/* 方向ベクトルの鏡映 (平面の位置によらず軸成分の符号反転) */
static vec3_t mirror_dir(vec3_t a, int axis)
{
	if (axis == 0) a.x = -a.x;
	else if (axis == 1) a.y = -a.y;
	else a.z = -a.z;
	return a;
}

/*
鏡像光源を p->emit の末尾に追加する。p->nemit は実光源のみを数えたまま
にはできないので、実光源数を nreal に控えて全体を nemit とする。
*/
void setup_images(ppfd_t *p)
{
	const int nreal = p->nemit;
	int    lev, f, i;
	int    cap = p->nemit;
	int    beg = 0, end = nreal;
	int   *lastface = NULL;
	double rs[6];
	int    nspecular = 0;

	p->nimage = 0;
	p->spec_lost = 0.0;
	for (f = 0; f < 6; f++) {
		rs[f] = p->mat[p->wallmat[f]].rhos;
		if (rs[f] > 0.0) nspecular++;
	}
	if ((nspecular == 0) || (p->specbounce <= 0)) return;

	lastface = (int *)xmalloc((size_t)(nreal + 1) * sizeof(int));
	for (i = 0; i < nreal; i++) lastface[i] = -1;

	for (lev = 0; lev < p->specbounce; lev++) {
		const int b = beg, e = end;
		for (i = b; i < e; i++) {
			for (f = 0; f < 6; f++) {
				emitter_t im;
				int    axis;
				double val;
				if (rs[f] <= 0.0) continue;
				if (lastface[i] == f) continue;   /* 同じ面で 2 回続けては折り返さない */
				face_plane(p, f, &axis, &val);
				im = p->emit[i];
				im.flux = p->emit[i].flux * rs[f];
				im.watt = 0.0;                    /* 消費電力は実光源のみ */
				if (im.flux <= 0.0) continue;
				im.pos = mirror_point(im.pos, axis, val);
				im.dir = mirror_dir(im.dir, axis);
				im.ax  = mirror_dir(im.ax, axis);
				im.ay  = mirror_dir(im.ay, axis);
				im.cx  = mirror_dir(im.cx, axis);
				im.cy  = mirror_dir(im.cy, axis);
				if (p->nemit >= cap) {
					cap = cap ? (2 * cap) : 16;
					p->emit = (emitter_t *)realloc(p->emit, (size_t)cap * sizeof(emitter_t));
					lastface = (int *)realloc(lastface, (size_t)cap * sizeof(int));
					if ((p->emit == NULL) || (lastface == NULL)) {
						fprintf(stderr, "*** out of memory\n");
						exit(1);
					}
				}
				if (lev == p->specbounce - 1) {
					/*
					次の段は作らないので、この鏡像が生むはずだった分を数える。
					自分が折り返した面 f へは戻らないので、それ以外の鏡面だけ。
					鏡面が 1 面しかなければ 0 になる (打ち切り誤差なし)。
					*/
					int g;
					for (g = 0; g < 6; g++) {
						if ((g != f) && (rs[g] > 0.0)) p->spec_lost += im.flux * rs[g];
					}
				}
				lastface[p->nemit] = f;
				p->emit[p->nemit++] = im;
				p->nimage++;
			}
		}
		beg = e;
		end = p->nemit;
		if (beg >= end) break;
	}
	free(lastface);
}

/*
鏡面反射の変換列を列挙する (拡散光の輸送用)。

face[0] が放射側に最も近い折り返し、face[nmir-1] が受光側に最も近い
折り返し。同じ面で 2 回続けては折り返さない (恒等になる)。列挙の規則は
setup_images の鏡像光源と同一なので、直接光と拡散光で経路の数え方が
食い違わない。戻り値 = 変換の数 (maxn で打ち切り)。
*/
int spec_transforms(const ppfd_t *p, strans_t *tr, int maxn)
{
	double rs[6];
	int    f, lev, k, n = 0, beg, end;

	for (f = 0; f < 6; f++) rs[f] = p->mat[p->wallmat[f]].rhos;

	for (f = 0; f < 6; f++) {
		if ((rs[f] <= 0.0) || (n >= maxn)) continue;
		tr[n].nmir = 1;
		tr[n].face[0] = f;
		tr[n].w = rs[f];
		n++;
	}
	beg = 0;
	end = n;
	for (lev = 2; lev <= p->specbounce; lev++) {
		for (k = beg; k < end; k++) {
			const int last = tr[k].face[tr[k].nmir - 1];
			for (f = 0; f < 6; f++) {
				if ((rs[f] <= 0.0) || (f == last) || (n >= maxn)) continue;
				tr[n] = tr[k];
				tr[n].face[tr[n].nmir] = f;
				tr[n].nmir++;
				tr[n].w *= rs[f];
				n++;
			}
		}
		beg = end;
		end = n;
	}
	return n;
}

/* 変換列を点に適用する (face[0] から順に鏡映) */
vec3_t spec_apply(const ppfd_t *p, const strans_t *t, vec3_t a)
{
	int    k, axis;
	double val;

	for (k = 0; k < t->nmir; k++) {
		face_plane(p, t->face[k], &axis, &val);
		a = mirror_point(a, axis, val);
	}
	return a;
}

/* 変換列を方向ベクトルに適用する */
vec3_t spec_apply_dir(const strans_t *t, vec3_t a)
{
	int k;

	for (k = 0; k < t->nmir; k++) {
		a = mirror_dir(a, t->face[k] / 2);
	}
	return a;
}

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
