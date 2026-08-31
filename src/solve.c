/*
solve.c

解析の流れ (OpenFDTD sol/solve.c に相当するメインドライバ) :

  1. setup_patch      チャンバ壁面をパッチに分割
  2. setup_target     測定面グリッドの確保
  3. setup_ff         パッチ間形態係数 (+ 群落経路長)
  4. direct_patch     光源 -> パッチの直接放射照度
  5. solve_radiosity  相互反射 (波長ビンごとの線形系)
  6. gather_target    測定面の各セルで直接 + 間接を合成
  7. energy_balance   放射束の収支 (壁の吸収 / 群落の吸収)

エネルギー収支は検証の要。群落が無い閉キャビティでは、どの反射率でも
「壁の吸収 = 光源の放射束」が厳密に成り立つ (定常状態の熱力学的な要請で、
形状にはよらない)。ここからのずれがそのまま形態係数の求積誤差になる。
*/

#include "ppfd.h"

/* 測定面の各セルで直接光 + 壁からの間接光を合成する */
static void gather_target(ppfd_t *p)
{
	const int n = p->npatch;
	const int nl = p->nlam;
	double *kext = (double *)xcalloc((size_t)nl, sizeof(double));
	vec3_t  ctop[4], cbot[4];
	strans_t *str = (strans_t *)xmalloc((size_t)MAXSPECT * sizeof(strans_t));
	const int nstr = spec_transforms(p, str, MAXSPECT);
	int     it, i;

	if (p->canopy.on && p->canopy.scatter) {
		const double zt = p->canopy.ztop, zb = p->canopy.zbot;
		ctop[0] = v_make(0.0, 0.0, zt);   ctop[1] = v_make(p->Lx, 0.0, zt);
		ctop[2] = v_make(p->Lx, p->Ly, zt); ctop[3] = v_make(0.0, p->Ly, zt);
		cbot[0] = v_make(0.0, 0.0, zb);   cbot[1] = v_make(p->Lx, 0.0, zb);
		cbot[2] = v_make(p->Lx, p->Ly, zb); cbot[3] = v_make(0.0, p->Ly, zb);
	}
	else {
		memset(ctop, 0, sizeof(ctop));
		memset(cbot, 0, sizeof(cbot));
	}

	if (p->canopy.on) {
		for (i = 0; i < nl; i++) kext[i] = canopy_kext(p, i);
	}

	for (it = 0; it < p->ntarget; it++) {
		target_t *t = &p->target[it];
		const double dx = (t->x1 - t->x0) / t->nx;
		const double dy = (t->y1 - t->y0) / t->ny;
		const int ncell = t->nx * t->ny;
		int ic;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
		for (ic = 0; ic < ncell; ic++) {
			const int ix = ic % t->nx;
			const int iy = ic / t->nx;
			const vec3_t nz = {0.0, 0.0, 1.0};
			const vec3_t x = {t->x0 + ((ix + 0.5) * dx), t->y0 + ((iy + 0.5) * dy), t->z};
			double *E = &t->E[(size_t)ic * nl];
			int ip, il, i2;

			direct_point(p, x, nz, E);

			if (p->canopy.on && p->canopy.scatter) {
				/* 群落が返す拡散光 : 外なら上下面から、内部なら層の下向き流束 */
				if (x.z >= p->canopy.ztop) {
					const double f = ff_point_poly(x, nz, ctop, 4);
					for (il = 0; il < nl; il++) E[il] += f * p->cbtop[il];
				}
				else if (x.z <= p->canopy.zbot) {
					const double f = ff_point_poly(x, nz, cbot, 4);
					for (il = 0; il < nl; il++) E[il] += f * p->cbbot[il];
				}
				else {
					for (il = 0; il < nl; il++) E[il] += canopy_interior(p, il, x.z);
				}
			}

			for (ip = 0; ip < n; ip++) {
				const patch_t *q = &p->patch[ip];
				const double *B = &p->B[(size_t)ip * nl];
				double f, s;
				/* パッチの表側から見ているか (遮蔽物の裏面を拾わないため) */
				if (v_dot(q->n, v_sub(x, q->p[0])) <= 0.0) continue;
				f = ff_point_poly(x, nz, q->p, 4);
				if (f <= 0.0) continue;
				f *= vis_point_patch(p, x, q, NULL);
				if (f <= 0.0) continue;
				s = canopy_path(p, x, q->c);
				if (s > 0.0) {
					for (il = 0; il < nl; il++) {
						E[il] += B[il] * f * exp(-kext[il] * s);
					}
				}
				else {
					for (il = 0; il < nl; il++) {
						E[il] += B[il] * f;
					}
				}
			}

			/* 鏡面壁ごしの拡散光 (鏡像パッチからの寄与。setup_ff の拡張と同じ帳簿) */
			for (i2 = 0; i2 < nstr; i2++) {
				for (ip = 0; ip < n; ip++) {
					const patch_t *q = &p->patch[ip];
					const double *B = &p->B[(size_t)ip * nl];
					vec3_t q4[4], nj;
					double f, s;
					int    v;
					if (spec_on_mirror(p, &str[i2], q)) continue;
					for (v = 0; v < 4; v++) q4[v] = spec_apply(p, &str[i2], q->p[v]);
					nj = spec_apply_dir(&str[i2], q->n);
					if (v_dot(nj, v_sub(x, q4[0])) <= 0.0) continue;
					f = str[i2].w * ff_point_poly(x, nz, q4, 4);
					if (f <= 0.0) continue;
					/* 鏡像重心への直線から実経路の群落内長さを測る
					   (床・天井の鏡映ぶんは canopy_path が折り返す) */
					s = canopy_path(p, x, spec_apply(p, &str[i2], q->c));
					if (s > 0.0) {
						for (il = 0; il < nl; il++) {
							E[il] += B[il] * f * exp(-kext[il] * s);
						}
					}
					else {
						for (il = 0; il < nl; il++) {
							E[il] += B[il] * f;
						}
					}
				}
			}
		}
	}

	free(str);
	free(kext);
}

/* 放射束・光量子束の総量と収支 */
static void energy_balance(ppfd_t *p)
{
	const int nl = p->nlam;
	double *E = (double *)xcalloc((size_t)nl, sizeof(double));
	int     ie, ip, il;

	p->flux_total = 0.0;
	p->watt_total = 0.0;
	p->ppf_total = 0.0;
	/* 鏡像は実在の光源ではないので総量には数えない */
	for (ie = 0; ie < p->nemit - p->nimage; ie++) {
		const emitter_t *e = &p->emit[ie];
		p->flux_total += e->flux;
		p->watt_total += e->watt;
		for (il = 0; il < nl; il++) E[il] = e->flux * p->spec[e->ispec].w[il];
		p->ppf_total += calc_ppfd(p, E);
	}

	p->absorb_wall = 0.0;
	p->ppf_wall = 0.0;
	for (ip = 0; ip < p->npatch; ip++) {
		const pmat_t *m = &p->mat[p->patch[ip].imat];
		const double *Ei = &p->Einc[(size_t)ip * nl];
		double w = 0.0;
		for (il = 0; il < nl; il++) {
			/* 鏡面成分は鏡像が運ぶので、壁が吸収するのは 1 - rho_d - rho_s */
			double a = 1.0 - m->rho[il] - m->rhos;
			if (a < 0.0) a = 0.0;
			E[il] = Ei[il] * a;
			w += E[il];
		}
		p->absorb_wall += p->patch[ip].area * w;
		p->ppf_wall += p->patch[ip].area * calc_ppfd(p, E);
	}

	if (p->canopy.on && p->canopy.scatter) {
		/* 残差ではなく実計算 : 群落が吸収した分を層ごとに足し上げてある */
		double a = 0.0;
		for (il = 0; il < nl; il++) {
			E[il] = p->cabs[il];
			a += p->cabs[il];
		}
		p->absorb_canopy = a;
		p->canopy_abs = calc_ppfd(p, E);
	}
	else {
		p->absorb_canopy = p->flux_total - p->absorb_wall;
		p->canopy_abs = p->ppf_total - p->ppf_wall;
	}

	free(E);
}

void solve(ppfd_t *p)
{
	double t0;

	t0 = cputime();
	setup_patch(p);
	setup_target(p);
	plog(p, "patches = %d (chamber %d x %d x %d, occluders %d), spectral bins = %d\n",
		p->npatch, 6, p->ndivu, p->ndivv,
		p->npatch - (6 * p->ndivu * p->ndivv), p->nlam);

	setup_ff(p);
	plog(p, "form factor : row sum error = %.3e, reciprocity error = %.3e  (%.2f s)\n",
		p->ff_rowsum_err, p->ff_recip_err, cputime() - t0);
	if (p->canopy.on && p->canopy.scatter) {
		canopy_setup(p);
		plog(p, "canopy scatter : %d layers, plane closure up/down = %.4f / %.4f, "
			"buried patches = %d of %d\n",
			p->canopy.nlayer, p->cnorm_up, p->cnorm_dn, p->cnburied, p->npatch);
	}
	if (p->nocc > 0) {
		/* 遮蔽が部分的な対だけは標本化 (= 厳密でない) なので数を出す */
		plog(p, "occluders   : %d, partially shadowed patch pairs = %d, enclosed patches = %d\n",
			p->nocc, p->ff_npartial, p->ff_nblind);
	}

	setup_images(p);
	if (p->nimage > 0) {
		plog(p, "specular : %d mirror images (%d bounces), truncated flux = %.3e W\n",
			p->nimage, p->specbounce, p->spec_lost);
	}

	t0 = cputime();
	direct_patch(p);
	plog(p, "direct irradiance done  (%.2f s)\n", cputime() - t0);

	t0 = cputime();
	solve_radiosity(p);
	plog(p, "radiosity : %d iterations, residual = %.3e  (%.2f s)\n",
		p->niter, p->resid, cputime() - t0);

	t0 = cputime();
	gather_target(p);
	plog(p, "target planes done  (%.2f s)\n", cputime() - t0);

	energy_balance(p);
}
