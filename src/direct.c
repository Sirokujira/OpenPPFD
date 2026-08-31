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

鏡像は「折り返した面の並び」ではなく軸ごとの展開指数で数える。直交する
面の鏡映は可換なので、並びで数えると同じ像を何度も数えてしまう
(xmax→zmax と zmax→xmax は同じ点)。

【制限】鏡像からの直線は折り返した経路を表すので、遮蔽物があると遮蔽
判定が合わない。両立は入力段で明示的に弾く (群落の経路長は canopy_path
が折り返して測るので併用できる)。
*/

/* 軸 a の展開指数 n に対応する重み (交差する面の Π ρs)。0 なら経路なし */
static double axis_weight(const double *rs, int a, int n)
{
	const int fmin = 2 * a, fmax = (2 * a) + 1;
	const int m = (n < 0) ? -n : n;
	double w = 1.0;
	int    j;

	/* n > 0 は max 面から、n < 0 は min 面から交互に折り返す */
	for (j = 1; j <= m; j++) {
		const int odd = (j % 2) != 0;
		w *= rs[(n > 0) ? (odd ? fmax : fmin) : (odd ? fmin : fmax)];
		if (w <= 0.0) return 0.0;
	}
	return w;
}

/* 展開指数 n における軸 a の座標写像 (箱の長さ L) */
static double axis_map(double c, double L, int n)
{
	/* C の % は負数で符号が残るので、奇数判定はそのまま使える */
	return ((n % 2) == 0) ? (c + ((double)n * L)) : (((double)(n + 1) * L) - c);
}

/*
鏡面反射の変換を列挙する。

軸ごとの展開指数 (n0, n1, n2) を |n0|+|n1|+|n2| = 1..maxord の範囲で
すべて回す。直交する鏡映が可換なせいで起きる二重計上は、この数え方
では原理的に起こらない (1 つの鏡像 = 1 つの (n0,n1,n2))。
戻り値 = 変換の数 (maxn で打ち切り)。
*/
static int enum_transforms(const ppfd_t *p, strans_t *tr, int maxn, int maxord)
{
	double rs[6];
	int    f, n0, n1, n2, n = 0;

	for (f = 0; f < 6; f++) rs[f] = p->mat[p->wallmat[f]].rhos;
	if (maxord < 1) return 0;

	for (n0 = -maxord; n0 <= maxord; n0++) {
		const double w0 = axis_weight(rs, 0, n0);
		const int    m0 = (n0 < 0) ? -n0 : n0;
		if (w0 <= 0.0) continue;
		for (n1 = -(maxord - m0); n1 <= (maxord - m0); n1++) {
			const double w1 = w0 * axis_weight(rs, 1, n1);
			const int    m1 = m0 + ((n1 < 0) ? -n1 : n1);
			if (w1 <= 0.0) continue;
			for (n2 = -(maxord - m1); n2 <= (maxord - m1); n2++) {
				const double w2 = w1 * axis_weight(rs, 2, n2);
				const int    m2 = m1 + ((n2 < 0) ? -n2 : n2);
				if ((m2 < 1) || (w2 <= 0.0)) continue;
				if (n >= maxn) return n;
				tr[n].n[0] = n0;
				tr[n].n[1] = n1;
				tr[n].n[2] = n2;
				tr[n].nmir = m2;
				tr[n].w = w2;
				n++;
			}
		}
	}
	return n;
}

/* 拡散光の輸送用 : specbounce 段までの変換 (setup_images と同一の数え方) */
int spec_transforms(const ppfd_t *p, strans_t *tr, int maxn)
{
	return enum_transforms(p, tr, maxn, p->specbounce);
}

/* 変換を点に適用する */
vec3_t spec_apply(const ppfd_t *p, const strans_t *t, vec3_t a)
{
	a.x = axis_map(a.x, p->Lx, t->n[0]);
	a.y = axis_map(a.y, p->Ly, t->n[1]);
	a.z = axis_map(a.z, p->Lz, t->n[2]);
	return a;
}

/* 変換を方向ベクトルに適用する (折り返し回数が奇数の軸だけ符号反転) */
vec3_t spec_apply_dir(const strans_t *t, vec3_t a)
{
	if ((t->n[0] % 2) != 0) a.x = -a.x;
	if ((t->n[1] % 2) != 0) a.y = -a.y;
	if ((t->n[2] % 2) != 0) a.z = -a.z;
	return a;
}

/*
パッチ q が変換 t の折り返し面の上に載っているか。

載っていると鏡像がパッチ自身と同じ平面に来て経路が退化する (自分が
置かれている鏡で自分を折り返すことになる) ので、拡張形態係数の計算から
落とす。「その軸で折り返していて (n != 0)、かつ座標が動かない」が条件。
*/
int spec_on_mirror(const ppfd_t *p, const strans_t *t, const patch_t *q)
{
	int    a;
	double L, c;

	if (q->iface >= 6) return 0;            /* 遮蔽物の面 (鏡面とは併用しない) */
	a = q->iface / 2;
	if (t->n[a] == 0) return 0;
	L = (a == 0) ? p->Lx : ((a == 1) ? p->Ly : p->Lz);
	c = (a == 0) ? q->p[0].x : ((a == 1) ? q->p[0].y : q->p[0].z);
	return (axis_map(c, L, t->n[a]) == c);
}

/*
鏡像光源を p->emit の末尾に追加する。p->nemit は実光源のみを数えたまま
にはできないので、実光源数を nreal に控えて全体を nemit とする。

打ち切り誤差 (spec_lost) は「1 段深い変換が運ぶはずだった放射束」=
実光源の総放射束 × Σ(order = specbounce+1 の変換の重み)。鏡面が 1 面
だけなら 2 段目以降は必ず反対側の面 (ρs = 0) を通るので厳密に 0 になる。
*/
void setup_images(ppfd_t *p)
{
	const int nreal = p->nemit;
	strans_t *tr;
	double flux_real = 0.0;
	int    nt, it, i, f, nspecular = 0, nimg = 0;

	p->nimage = 0;
	p->spec_lost = 0.0;
	for (f = 0; f < 6; f++) {
		if (p->mat[p->wallmat[f]].rhos > 0.0) nspecular++;
	}
	if ((nspecular == 0) || (p->specbounce <= 0)) return;

	for (i = 0; i < nreal; i++) flux_real += p->emit[i].flux;

	tr = (strans_t *)xmalloc((size_t)MAXSPECT * sizeof(strans_t));
	nt = enum_transforms(p, tr, MAXSPECT, p->specbounce + 1);

	/* 鏡像の総数はあらかじめ分かるので 1 回で確保する */
	for (it = 0; it < nt; it++) {
		if (tr[it].nmir <= p->specbounce) nimg += nreal;
	}
	p->emit = (emitter_t *)realloc(p->emit, (size_t)(nreal + nimg) * sizeof(emitter_t));
	if (p->emit == NULL) { fprintf(stderr, "*** out of memory\n"); exit(1); }

	for (it = 0; it < nt; it++) {
		if (tr[it].nmir > p->specbounce) {
			p->spec_lost += flux_real * tr[it].w;
			continue;
		}
		for (i = 0; i < nreal; i++) {
			emitter_t im = p->emit[i];
			im.flux = p->emit[i].flux * tr[it].w;
			im.watt = 0.0;                    /* 消費電力は実光源のみ */
			im.pos = spec_apply(p, &tr[it], im.pos);
			im.dir = spec_apply_dir(&tr[it], im.dir);
			im.ax  = spec_apply_dir(&tr[it], im.ax);
			im.ay  = spec_apply_dir(&tr[it], im.ay);
			im.cx  = spec_apply_dir(&tr[it], im.cx);
			im.cy  = spec_apply_dir(&tr[it], im.cy);
			p->emit[p->nemit++] = im;
			p->nimage++;
		}
	}
	free(tr);
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
