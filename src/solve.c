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
	int     it, i;

	if (p->canopy.on) {
		for (i = 0; i < nl; i++) {
			double omega = p->mat[p->canopy.imat].rho[i] + p->mat[p->canopy.imat].tau[i];
			if (omega < 0.0) omega = 0.0;
			if (omega > 0.999) omega = 0.999;
			kext[i] = p->canopy.k0 * sqrt(1.0 - omega) * p->canopy.a;
		}
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
			int ip, il;

			direct_point(p, x, nz, E);

			for (ip = 0; ip < n; ip++) {
				const double f = ff_point_poly(x, nz, p->patch[ip].p, 4);
				const double *B = &p->B[(size_t)ip * nl];
				double s;
				if (f <= 0.0) continue;
				s = canopy_path(p, x, p->patch[ip].c);
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
	for (ie = 0; ie < p->nemit; ie++) {
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
			E[il] = Ei[il] * (1.0 - m->rho[il]);
			w += E[il];
		}
		p->absorb_wall += p->patch[ip].area * w;
		p->ppf_wall += p->patch[ip].area * calc_ppfd(p, E);
	}

	p->absorb_canopy = p->flux_total - p->absorb_wall;
	p->canopy_abs = p->ppf_total - p->ppf_wall;

	free(E);
}

void solve(ppfd_t *p)
{
	double t0;

	t0 = cputime();
	setup_patch(p);
	setup_target(p);
	plog(p, "patches = %d (%d faces x %d x %d), spectral bins = %d\n",
		p->npatch, 6, p->ndivu, p->ndivv, p->nlam);

	setup_ff(p);
	plog(p, "form factor : row sum error = %.3e, reciprocity error = %.3e  (%.2f s)\n",
		p->ff_rowsum_err, p->ff_recip_err, cputime() - t0);

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
