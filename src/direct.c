/*
direct.c

光源から壁パッチ / 測定点への直接放射照度。

■ 配光
放射束 Φ の光源の放射強度は
    I(θ) = Φ (m+1)/(2π) cos^m θ      (θ = 配光軸からの角度、上半空間のみ)
m = 1 が Lambert 配光 I = (Φ/π) cosθ、m = 0 が半球一様、m < 0 は等方
    I = Φ/(4π)
を表す。∫I dΩ = Φ は m によらず満たされる。

■ 面光源
寸法 w × h の平面 Lambert 光源は nu × nv 個の点光源に分割し、各々に
Φ/(nu·nv) を与える。分割数 → ∞ で軸上放射照度は解析解
    E = Φ F_{dA→A_s} / A_s
に収束する (data/sample/panel.ppfd の検証ケース)。

■ 群落
光源から受光点までの線分が群落スラブを横切る長さ s を求め、分光透過率
exp(-G a s sqrt(1-ω_λ)) を掛ける (geometry.c)。
*/

#include "ppfd.h"

/*
点 x (単位法線 n) における全光源からの直接分光放射照度 [W/m²] を
E[0..nlam-1] に加算する (E は呼び出し側でゼロ初期化しておく)。
*/
void direct_point(const ppfd_t *p, vec3_t x, vec3_t n, double *E)
{
	int    ie, i;
	double *acc = NULL;
	double  tr[512];
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

				if (e->mexp < 0.0) {
					inten = fsub / (4.0 * PI);
				}
				else {
					cs = v_dot(e->dir, r) / d;
					if (cs <= 0.0) continue;
					inten = fsub * (e->mexp + 1.0) / (2.0 * PI) * pow(cs, e->mexp);
				}

				es = inten * cr / d2;

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

/*
各パッチの面積平均直接放射照度 Ed[ip*nlam + i]。

複合 3 点 Gauss (msub x msub の小矩形 x 3x3 点) で面積平均する。計算量は
nemit * npatch * (3*msub)^2 * (有効ビン数) に比例するので、光源が多い
ケースでは msub を落とさないと現実的な時間に収まらない。既定は光源数
から自動で決め、実際に使った値をログに出す (quadrature キーで上書き可)。

光源が少ないほど被積分関数が急峻なので、この既定は精度を落とさない
方向に働く : 検証ケース (光源 1 個) は常に msub = 3 になる。
*/
void direct_patch(ppfd_t *p)
{
	const int n = p->npatch;
	int ip;
	int msub = p->msub;

	if (msub <= 0) {
		msub = (p->nemit <= 8) ? 3 : ((p->nemit <= 64) ? 2 : 1);
	}
	plog(p, "direct quadrature : msub = %d (%d points/patch)\n", msub, 9 * msub * msub);

	p->Ed = (double *)xcalloc((size_t)n * p->nlam, sizeof(double));

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
	for (ip = 0; ip < n; ip++) {
		const patch_t *q = &p->patch[ip];
		const vec3_t eu = v_sub(q->p[1], q->p[0]);
		const vec3_t ev = v_sub(q->p[3], q->p[0]);
		double *E = &p->Ed[(size_t)ip * p->nlam];
		double *tmp = (double *)xcalloc((size_t)p->nlam, sizeof(double));
		const double h = 1.0 / msub;
		int su, sv, gu, gv, i;

		for (sv = 0; sv < msub; sv++) {
			for (su = 0; su < msub; su++) {
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
						direct_point(p, x, q->n, tmp);
						for (i = 0; i < p->nlam; i++) E[i] += w * tmp[i];
					}
				}
			}
		}
		free(tmp);
	}
}
