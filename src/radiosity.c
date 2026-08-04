/*
radiosity.c

壁面間の相互反射 (ラジオシティ)。波長ビンごとに独立な線形系

    B_i = ρ_i ( Ed_i + Σ_j F_ij τ_ij B_j )

を Jacobi 反復で解く。B は放射発散度 [W/m²]、Ed は直接放射照度、
τ_ij は群落スラブによるパッチ間の透過率 (パッチ重心どうしを結ぶ線分で
評価する近似)。ρ_i < 1 なので反復行列のスペクトル半径は 1 未満、
必ず収束する。

波長を外側ループにすると τ_ij の exp() が波長ごと 1 回で済み、反復の
内側では行列積だけになる。光源の放射束が 0 のビンは丸ごと飛ばす。
*/

#include "ppfd.h"

void solve_radiosity(ppfd_t *p)
{
	const int n = p->npatch;
	const int nl = p->nlam;
	double  *M = NULL;
	double  *Bold, *Bnew;
	int      il;
	int      itermax = 0;
	double   residmax = 0.0;

	p->B = (double *)xcalloc((size_t)n * nl, sizeof(double));
	p->Einc = (double *)xmalloc((size_t)n * nl * sizeof(double));
	memcpy(p->Einc, p->Ed, (size_t)n * nl * sizeof(double));

	if (p->canopy.on) {
		M = (double *)xmalloc((size_t)n * n * sizeof(double));
	}
	Bold = (double *)xcalloc((size_t)n, sizeof(double));
	Bnew = (double *)xcalloc((size_t)n, sizeof(double));

	for (il = 0; il < nl; il++) {
		const double *F;
		double kext = 0.0;
		double emax = 0.0;
		int    i, iter, active = 0;

		/* このビンに光が入っているか */
		for (i = 0; i < n; i++) {
			if (p->Ed[((size_t)i * nl) + il] > 0.0) {
				active = 1;
				break;
			}
		}
		if (!active) continue;

		/* 反射率が全パッチで 0 なら相互反射は無い */
		for (i = 0; i < n; i++) {
			const double r = p->mat[p->patch[i].imat].rho[il];
			if (r > emax) emax = r;
		}
		if (emax <= 0.0) continue;

		if (p->canopy.on) {
			double omega = p->mat[p->canopy.imat].rho[il] + p->mat[p->canopy.imat].tau[il];
			int j;
			if (omega < 0.0) omega = 0.0;
			if (omega > 0.999) omega = 0.999;
			kext = p->canopy.k0 * sqrt(1.0 - omega) * p->canopy.a;
#ifdef _OPENMP
#pragma omp parallel for private(j) schedule(static)
#endif
			for (i = 0; i < n; i++) {
				for (j = 0; j < n; j++) {
					const size_t k = ((size_t)i * n) + j;
					const double s = p->plen[k];
					M[k] = (s > 0.0) ? (p->ff[k] * exp(-kext * s)) : p->ff[k];
				}
			}
			F = M;
		}
		else {
			F = p->ff;
		}

		/* 初期値 : 直接光による 1 回反射 */
		for (i = 0; i < n; i++) {
			Bold[i] = p->mat[p->patch[i].imat].rho[il] * p->Ed[((size_t)i * nl) + il];
		}

		for (iter = 0; iter < p->maxiter; iter++) {
			double dmax = 0.0, bmax = 0.0;

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
			for (i = 0; i < n; i++) {
				const double *Fi = &F[(size_t)i * n];
				double g = p->Ed[((size_t)i * nl) + il];
				int    j;
				for (j = 0; j < n; j++) {
					g += Fi[j] * Bold[j];
				}
				Bnew[i] = p->mat[p->patch[i].imat].rho[il] * g;
			}

			for (i = 0; i < n; i++) {
				const double d = fabs(Bnew[i] - Bold[i]);
				if (d > dmax) dmax = d;
				if (Bnew[i] > bmax) bmax = Bnew[i];
				Bold[i] = Bnew[i];
			}

			if (iter + 1 > itermax) itermax = iter + 1;
			if (bmax <= 0.0) break;
			if ((dmax / bmax) < p->converg) {
				if ((dmax / bmax) > residmax) residmax = dmax / bmax;
				break;
			}
			if (iter == p->maxiter - 1) {
				if ((dmax / bmax) > residmax) residmax = dmax / bmax;
			}
		}

		/* 入射放射照度 (エネルギー収支の診断に使う) */
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
		for (i = 0; i < n; i++) {
			const double *Fi = &F[(size_t)i * n];
			double g = p->Ed[((size_t)i * nl) + il];
			int    j;
			for (j = 0; j < n; j++) {
				g += Fi[j] * Bold[j];
			}
			p->Einc[((size_t)i * nl) + il] = g;
			p->B[((size_t)i * nl) + il] = p->mat[p->patch[i].imat].rho[il] * g;
		}
	}

	p->niter = itermax;
	p->resid = residmax;

	free(Bold);
	free(Bnew);
	if (M != NULL) free(M);
}
