/*
main.c

OpenPPFD : 植物工場照明 (PPFD / スペクトル / 輻射伝達) ソルバー

    oppfd [-n <threads>] <input.ppfd>

出力は実行時のカレントディレクトリに ppfd.log / ppfd_summary.csv /
ppfd_spectrum.csv / ppfd_map_<target>.csv。
*/

#include "ppfd.h"

#ifdef _OPENMP
#include <omp.h>
#endif

static void usage(void)
{
	fprintf(stderr, "Usage: oppfd [-n <threads>] <input.ppfd>\n");
	fprintf(stderr, "  -n <threads>  number of OpenMP threads (default: system)\n");
	fprintf(stderr, "  --selftest    run internal algebraic identity checks and exit\n");
}

int main(int argc, char **argv)
{
	ppfd_t  p;
	FILE   *fp;
	const char *fname = NULL;
	int     i, nthread = 0;
	double  t0;

	for (i = 1; i < argc; i++) {
		if (!strcmp(argv[i], "-n") && (i + 1 < argc)) {
			nthread = atoi(argv[++i]);
		}
		else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
			usage();
			return 0;
		}
		else if (!strcmp(argv[i], "--selftest")) {
			/* 幾何にも入力にも依存しない代数的性質の検査 (3 OS の CI で回す) */
			return canopy_selftest() ? 1 : 0;
		}
		else if (argv[i][0] == '-') {
			fprintf(stderr, "*** unknown option : %s\n", argv[i]);
			usage();
			return 1;
		}
		else {
			fname = argv[i];
		}
	}
	if (fname == NULL) {
		usage();
		return 1;
	}

#ifdef _OPENMP
	if (nthread > 0) omp_set_num_threads(nthread);
#else
	(void)nthread;
#endif

	fp = fopen(fname, "r");
	if (fp == NULL) {
		fprintf(stderr, "*** cannot open : %s\n", fname);
		return 1;
	}

	if (input_data(fp, &p)) {
		fclose(fp);
		fprintf(stderr, "*** input error\n");
		return 1;
	}
	fclose(fp);

	if (p.nlam > MAXLAM) {
		fprintf(stderr, "*** too many spectral bins : %d > %d (increase the wavelength step)\n",
			p.nlam, MAXLAM);
		return 1;
	}

	p.fplog = fopen(FN_LOG, "w");
	if (p.fplog == NULL) {
		fprintf(stderr, "*** cannot open : %s\n", FN_LOG);
		return 1;
	}

	t0 = cputime();
	solve(&p);
	output_log(&p);
	output_csv(&p);
	plog(&p, "cpu time = %.2f sec\n", cputime() - t0);

	fclose(p.fplog);
	return 0;
}
