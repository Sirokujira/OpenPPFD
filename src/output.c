/*
output.c

ppfd.log (人間向け) と CSV (機械読み取り / GUI 向け) の出力。

一般照明 (lm / lx / lm/W) と植物用 (µmol/m²/s / µmol/J) を必ず並べて
出すのが本ソルバーの主眼。同じ器具でも視感度曲線 V(λ) と光合成の
光量子換算では効率の順位が入れ替わるため、両方を見ないと設計を誤る。
*/

#include "ppfd.h"

/* 測定面の 1 セル分のスペクトルへのポインタ */
#define TCELL(t, ic, nl) (&(t)->E[(size_t)(ic) * (nl)])

/* 測定面の統計 */
typedef struct {
	double avg, min, max, cv, u0;
	double ypfd, epar, lux, irr, rfr, dli;
	double *espec;                  /* 面平均スペクトル [W/m²] (呼び出し側で free) */
} tstat_t;

static void target_stat(const ppfd_t *p, const target_t *t, tstat_t *st)
{
	const int nl = p->nlam;
	const int ncell = t->nx * t->ny;
	double sum = 0.0, sum2 = 0.0, fr;
	int    ic, il;

	st->espec = (double *)xcalloc((size_t)nl, sizeof(double));
	st->min = 1e300;
	st->max = -1e300;

	for (ic = 0; ic < ncell; ic++) {
		const double *E = TCELL(t, ic, nl);
		const double q = calc_ppfd(p, E);
		sum += q;
		sum2 += q * q;
		if (q < st->min) st->min = q;
		if (q > st->max) st->max = q;
		for (il = 0; il < nl; il++) st->espec[il] += E[il];
	}
	for (il = 0; il < nl; il++) st->espec[il] /= ncell;

	st->avg = sum / ncell;
	{
		const double var = (sum2 / ncell) - (st->avg * st->avg);
		st->cv = (st->avg > 0.0) ? (sqrt((var > 0.0) ? var : 0.0) / st->avg) : 0.0;
	}
	st->u0 = (st->avg > 0.0) ? (st->min / st->avg) : 0.0;

	st->ypfd = calc_ypfd(p, st->espec);
	st->epar = calc_epar(p, st->espec);
	st->lux  = calc_lux(p, st->espec);
	st->irr  = calc_irradiance(p, st->espec);
	st->dli  = st->avg * p->photoperiod * 3600.0 * 1e-6;

	/* R:FR は狭帯域 (660±5 nm / 730±5 nm) の光量子束比 */
	fr = calc_bandppf(p, st->espec, 725.0, 735.0);
	st->rfr = (fr > 0.0) ? (calc_bandppf(p, st->espec, 655.0, 665.0) / fr) : 0.0;
}

void output_log(ppfd_t *p)
{
	const int nl = p->nlam;
	double  lm_total = 0.0;
	double  awall = 0.0, ewall = 0.0;
	int     ie, il, it, ib;

	plog(p, "\n");
	plog(p, "=== %s %d.%d.%d ===\n", PROGRAM, VERSION_MAJOR, VERSION_MINOR, VERSION_BUILD);
	if (p->title[0] != '\0') plog(p, "title :%s\n", p->title);
	plog(p, "chamber      : %g x %g x %g m\n", p->Lx, p->Ly, p->Lz);
	plog(p, "wavelength   : %g - %g nm, %g nm step (%d bins)\n",
		p->lam0, p->lam0 + ((nl - 1) * p->dlam), p->dlam, nl);
	plog(p, "sources      : %d\n", p->nemit);
	plog(p, "photoperiod  : %g h/day\n", p->photoperiod);
	for (ie = 0; ie < p->nocc; ie++) {
		const occluder_t *o = &p->occ[ie];
		plog(p, "occluder     : %-10s x %g..%g, y %g..%g, z %g..%g m, %s\n",
			o->name, o->x0, o->x1, o->y0, o->y1, o->z0, o->z1, p->mat[o->imat].name);
	}
	if (p->canopy.on) {
		plog(p, "canopy       : z = %g..%g m, LAI = %g, G = %g, leaf = %s\n",
			p->canopy.zbot, p->canopy.ztop, p->canopy.lai, p->canopy.k0,
			p->mat[p->canopy.imat].name);
	}
	for (ie = 0; ie < p->ndist; ie++) {
		const photdist_t *d = &p->dist[ie];
		plog(p, "photometry   : %s (%s, %d x %d angles",
			d->path, d->isldt ? "EULUMDAT" : "IES LM-63", d->nc, d->ng);
		if (d->lm > 0.0) plog(p, ", rated %.4g lm", d->lm);
		if (d->watt > 0.0) plog(p, ", %.4g W", d->watt);
		plog(p, ") -- shape only, flux taken from the input file\n");
		if (!d->isldt && (d->ptype != 1)) {
			plog(p, "  * warning : photometric type %d is not type C; treated as type C\n", d->ptype);
		}
	}

	/* 光源側の総量 */
	for (ie = 0; ie < p->nemit; ie++) {
		const emitter_t *e = &p->emit[ie];
		double v = 0.0;
		for (il = 0; il < nl; il++) v += p->spec[e->ispec].w[il] * p->vlambda[il];
		lm_total += e->flux * K_MAX * v;
	}

	plog(p, "\n--- source totals ---\n");
	plog(p, "radiant flux      = %12.6g W\n", p->flux_total);
	plog(p, "PPF (400-700nm)   = %12.6g umol/s\n", p->ppf_total);
	plog(p, "luminous flux     = %12.6g lm\n", lm_total);
	if (p->watt_total > 0.0) {
		plog(p, "input power       = %12.6g W\n", p->watt_total);
		plog(p, "photon efficacy   = %12.6g umol/J\n", p->ppf_total / p->watt_total);
		plog(p, "luminous efficacy = %12.6g lm/W\n", lm_total / p->watt_total);
	}
	else {
		plog(p, "input power       = (not given; add \"input <W>\" to a led/array key for efficacy)\n");
	}

	/* エネルギー収支 : 群落が無ければ壁の吸収 = 光源の放射束 (厳密) */
	for (it = 0; it < p->npatch; it++) {
		double s = 0.0;
		for (il = 0; il < nl; il++) s += p->Einc[((size_t)it * nl) + il];
		ewall += p->patch[it].area * s;
		awall += p->patch[it].area;
	}
	plog(p, "\n--- energy balance ---\n");
	plog(p, "emitted           = %12.6g W\n", p->flux_total);
	plog(p, "absorbed by walls = %12.6g W%s\n", p->absorb_wall,
		(p->nocc > 0) ? " (incl. occluders)" : "");
	plog(p, "mean wall irrad.  = %12.6g W/m2 (wall area %g m2)\n",
		(awall > 0.0) ? (ewall / awall) : 0.0, awall);
	if (p->canopy.on) {
		plog(p, "absorbed by canopy= %12.6g W\n", p->absorb_canopy);
		plog(p, "canopy PPF capture= %12.6g umol/s (%.2f %% of PPF)\n",
			p->canopy_abs, (p->ppf_total > 0.0) ? (100.0 * p->canopy_abs / p->ppf_total) : 0.0);
	}
	else {
		plog(p, "closure error     = %12.3e (should be 0; = form factor quadrature error)\n",
			(p->flux_total > 0.0) ? ((p->absorb_wall - p->flux_total) / p->flux_total) : 0.0);
	}

	/* 測定面 */
	for (it = 0; it < p->ntarget; it++) {
		target_t *t = &p->target[it];
		tstat_t   st;
		target_stat(p, t, &st);

		plog(p, "\n--- target \"%s\" (z = %g m, %d x %d cells) ---\n", t->name, t->z, t->nx, t->ny);
		plog(p, "%-22s %12s %12s %12s\n", "quantity", "average", "min", "max");
		plog(p, "%-22s %12.5g %12.5g %12.5g\n", "PPFD [umol/m2/s]", st.avg, st.min, st.max);
		plog(p, "%-22s %12.5g\n", "YPFD (McCree)", st.ypfd);
		plog(p, "%-22s %12.5g\n", "ePAR 400-750 [umol]", st.epar);
		plog(p, "%-22s %12.5g\n", "DLI [mol/m2/day]", st.dli);
		plog(p, "%-22s %12.5g\n", "irradiance [W/m2]", st.irr);
		plog(p, "%-22s %12.5g\n", "illuminance [lx]", st.lux);
		plog(p, "%-22s %12.5g\n", "uniformity min/avg", st.u0);
		plog(p, "%-22s %12.5g\n", "CV (sigma/avg)", st.cv);
		plog(p, "%-22s %12.5g\n", "R:FR (660/730)", st.rfr);
		if (st.avg > 0.0) {
			plog(p, "%-22s %12.5g\n", "YPFD/PPFD", st.ypfd / st.avg);
			plog(p, "%-22s %12.5g\n", "lx per umol/m2/s", st.lux / st.avg);
		}

		plog(p, "  spectral bands (photon flux fraction of 400-700nm PPFD):\n");
		for (ib = 0; ib < p->nband; ib++) {
			const double q = calc_bandppf(p, st.espec, p->band[ib].lam1, p->band[ib].lam2);
			plog(p, "    %-10s %7.1f-%7.1f nm : %10.5g umol/m2/s (%6.2f %%)\n",
				p->band[ib].name, p->band[ib].lam1, p->band[ib].lam2, q,
				(st.avg > 0.0) ? (100.0 * q / st.avg) : 0.0);
		}
		free(st.espec);
	}

	plog(p, "\n=== normal end ===\n");
}

void output_csv(ppfd_t *p)
{
	const int nl = p->nlam;
	FILE  *fp;
	int    it, il;

	/* --- 面ごとの分布 --- */
	for (it = 0; it < p->ntarget; it++) {
		target_t *t = &p->target[it];
		const double dx = (t->x1 - t->x0) / t->nx;
		const double dy = (t->y1 - t->y0) / t->ny;
		char   fn[BUFSIZ];
		int    ix, iy;

		snprintf(fn, sizeof(fn), "ppfd_map_%s.csv", t->name);
		fp = fopen(fn, "w");
		if (fp == NULL) continue;
		fprintf(fp, "ix,iy,x_m,y_m,ppfd_umol_m2_s,ypfd_umol_m2_s,epar_umol_m2_s,lux,irradiance_W_m2,r_fr\n");
		for (iy = 0; iy < t->ny; iy++) {
			for (ix = 0; ix < t->nx; ix++) {
				const int ic = (iy * t->nx) + ix;
				const double *E = TCELL(t, ic, nl);
				const double fr = calc_bandppf(p, E, 725.0, 735.0);
				fprintf(fp, "%d,%d,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g\n",
					ix, iy,
					t->x0 + ((ix + 0.5) * dx), t->y0 + ((iy + 0.5) * dy),
					calc_ppfd(p, E), calc_ypfd(p, E), calc_epar(p, E),
					calc_lux(p, E), calc_irradiance(p, E),
					(fr > 0.0) ? (calc_bandppf(p, E, 655.0, 665.0) / fr) : 0.0);
			}
		}
		fclose(fp);
	}

	/* --- 面ごとの要約 --- */
	fp = fopen(FN_SUMMARY, "w");
	if (fp != NULL) {
		fprintf(fp, "target,z_m,nx,ny,ppfd_avg,ppfd_min,ppfd_max,uniformity_min_avg,cv,"
		            "dli_mol_m2_day,ypfd,epar,lux,irradiance_W_m2,r_fr\n");
		for (it = 0; it < p->ntarget; it++) {
			target_t *t = &p->target[it];
			tstat_t   st;
			target_stat(p, t, &st);
			fprintf(fp, "%s,%.6g,%d,%d,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g\n",
				t->name, t->z, t->nx, t->ny,
				st.avg, st.min, st.max, st.u0, st.cv, st.dli,
				st.ypfd, st.epar, st.lux, st.irr, st.rfr);
			free(st.espec);
		}
		fclose(fp);
	}

	/* --- 面平均スペクトル --- */
	fp = fopen(FN_SPECTRUM, "w");
	if (fp != NULL) {
		double **sp = (double **)xcalloc((size_t)p->ntarget, sizeof(double *));
		fprintf(fp, "lambda_nm,action_mccree,v_lambda");
		for (it = 0; it < p->ntarget; it++) {
			tstat_t st;
			target_stat(p, p->target + it, &st);
			sp[it] = st.espec;
			fprintf(fp, ",%s_Ee_W_m2_nm,%s_Q_umol_m2_s_nm", p->target[it].name, p->target[it].name);
		}
		fprintf(fp, "\n");
		for (il = 0; il < nl; il++) {
			fprintf(fp, "%.6g,%.6g,%.6g", p->lam[il], p->action[il], p->vlambda[il]);
			for (it = 0; it < p->ntarget; it++) {
				const double ee = sp[it][il] / p->dlam;
				fprintf(fp, ",%.6g,%.6g", ee, ee * PHOTON_K(p->lam[il]));
			}
			fprintf(fp, "\n");
		}
		for (it = 0; it < p->ntarget; it++) free(sp[it]);
		free(sp);
		fclose(fp);
	}
}
