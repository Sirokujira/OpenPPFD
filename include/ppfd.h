/*
ppfd.h

OpenPPFD : 植物工場照明 (PPFD / スペクトル / 輻射伝達) ソルバー
定数・構造体・全プロトタイプ

OpenFDTD / OpenPEEC の姉妹プロジェクト。波動光学 (FDTD/RCWA) ではなく
光線・輻射伝達の問題を扱う。測光量が一般照明 (lm/lx, V(λ)) ではなく
光合成光量子束密度 (PPFD, µmol/m²/s, McCree 作用曲線) である点が本質。
*/

#ifndef _PPFD_H_
#define _PPFD_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define PROGRAM "OpenPPFD"
#define VERSION_MAJOR 1
#define VERSION_MINOR 0
#define VERSION_BUILD 0

#define FN_LOG     "ppfd.log"
#define FN_SUMMARY "ppfd_summary.csv"
#define FN_SPECTRUM "ppfd_spectrum.csv"

/* 数学・物理定数 (自前マクロ : <math.h> の M_PI には依存しない) */
#define PI      (4.0 * atan(1.0))
#define C0      (2.99792458e8)      /* 光速 [m/s] */
#define H_PLANCK (6.62607015e-34)   /* プランク定数 [J s] */
#define N_AVOG  (6.02214076e23)     /* アボガドロ定数 [1/mol] */
#define K_MAX   (683.0)             /* 最大視感効果度 [lm/W] @555nm */
#define EPS     (1e-12)

/*
光子換算係数 [µmol/(s·m²) per (W/m²)] @ λ [nm]
    1 mol の光子のエネルギー = h c N_A / λ [J/mol]  (λ [m])
    → µmol/J = 1e6 λ_nm 1e-9 / (h c N_A) = λ_nm * 1e-3 / 0.119626566
             = λ_nm * 8.359346e-3
h, c, N_A はいずれも SI 定義値なので、この係数は厳密量。
*/
#define PHOTON_K(lam) ((lam) * 1e-3 / (H_PLANCK * C0 * N_AVOG))

#define NAMELEN 32
#define MAXPATCH 2400               /* 形態係数行列が N^2 doubles */
#define MAXLAM 512                  /* 分光グリッド点数の上限 */

/* ---- ベクトル ---------------------------------------------------- */
typedef struct {double x, y, z;} vec3_t;

/* ---- 材料 (Lambert 拡散面 / 葉) ----------------------------------- */
typedef struct {
	char    name[NAMELEN];
	double *rho;                    /* [nlam] 分光反射率 0..1 */
	double *tau;                    /* [nlam] 分光透過率 0..1 (葉のみ使用) */
} pmat_t;

/* ---- スペクトル (波長ビンごとの放射束の重み、総和 1) --------------- */
typedef struct {
	char    name[NAMELEN];
	double *w;                      /* [nlam] Σw = 1 */
	int     i0, i1;                 /* w != 0 のビン範囲 (両端含む)。内側ループの短縮用 */
} spec_t;

/* ---- 実測配光 (IESNA LM-63 / EULUMDAT) ----------------------------- */
typedef struct {
	char    path[256];              /* ファイルパス (重複読み込みの判定にも使う) */
	int     isldt;                  /* 0 = IES, 1 = EULUMDAT */
	int     ptype;                  /* IES の光度分布型 (1 = type C) */
	int     ng, nc;                 /* 鉛直角 γ の数 / 水平角 C の数 (1 = 回転対称) */
	double *gam;                    /* [ng] γ [deg] 昇順、0 = 配光軸方向 */
	double *cpl;                    /* [nc] C [deg] 昇順 0..360 (対称性を展開済み) */
	double *I;                      /* [nc*ng] 相対光度 [1/sr]、∫I dΩ = 1 に正規化 */
	double  lm;                     /* ファイル記載の光度から求めた光束 [lm] (参考値) */
	double  watt;                   /* ファイル記載の消費電力 [W] (参考値) */
} photdist_t;

/* ---- 光源 --------------------------------------------------------- */
typedef struct {
	vec3_t  pos;                    /* 中心位置 [m] */
	vec3_t  dir;                    /* 配光軸 (単位ベクトル) */
	vec3_t  ax, ay;                 /* 面光源のローカル軸 (単位ベクトル) */
	vec3_t  cx, cy;                 /* 配光ファイルの C 面基準 (C=0 が cx、C=90 が cy) */
	double  flux;                   /* 放射束 [W] */
	double  watt;                   /* 消費電力 [W] (0 = 未指定) */
	int     ispec;                  /* スペクトル番号 */
	double  mexp;                   /* 配光 I(θ) ∝ cos^m θ (m<0 = 等方) */
	int     idist;                  /* 実測配光番号 (-1 = cos^m の解析形) */
	double  crot;                   /* 配光軸まわりの回転 [deg] */
	double  wsize, hsize;           /* 面光源の寸法 [m] (0 = 点光源) */
	int     nu, nv;                 /* 面光源の分割数 */
} emitter_t;

/* ---- 壁パッチ (平面四角形) ---------------------------------------- */
typedef struct {
	vec3_t  p[4];                   /* 頂点 (法線まわり CCW) */
	vec3_t  n;                      /* 単位法線 (キャビティ内向き) */
	vec3_t  c;                      /* 重心 */
	double  area;                   /* 面積 [m²] */
	int     imat;                   /* 材料番号 */
	int     iface;                  /* 所属面 0..5 */
} patch_t;

/* ---- 測定面 ------------------------------------------------------- */
typedef struct {
	char    name[NAMELEN];
	double  z;                      /* 高さ [m] (水平面のみ) */
	int     nx, ny;
	double  x0, x1, y0, y1;
	double *E;                      /* [nx*ny*nlam] 分光放射照度 [W/m²] */
} target_t;

/* ---- 帯域 (光量子束の割合を出す) ----------------------------------- */
typedef struct {
	char    name[NAMELEN];
	double  lam1, lam2;
} band_t;

/* ---- 群落 (Beer-Lambert 減衰スラブ) -------------------------------- */
typedef struct {
	int     on;
	double  ztop, zbot;             /* ztop > zbot */
	double  lai;                    /* 葉面積指数 [-] */
	double  k0;                     /* 葉群投影係数 G (球状分布 0.5 / 水平葉 1.0) */
	int     imat;                   /* 葉材料 (ω = rho + tau) */
	double  a;                      /* 葉面積密度 = lai / (ztop - zbot) [1/m] */
} canopy_t;

/* ---- 解析コンテキスト --------------------------------------------- */
typedef struct {
	char      title[BUFSIZ];

	/* 分光グリッド : lam[i] = lam0 + i*dlam、ビン幅 dlam */
	int       nlam;
	double    lam0, dlam;
	double   *lam;

	/* チャンバ */
	double    Lx, Ly, Lz;
	int       wallmat[6];           /* 面 0:xmin 1:xmax 2:ymin 3:ymax 4:zmin 5:zmax */
	int       ndivu, ndivv;         /* 1 面あたりの分割数 */

	/* 入力テーブル */
	int       nmat, nspec, nemit, ntarget, nband, ndist;
	pmat_t   *mat;
	spec_t   *spec;
	emitter_t *emit;
	target_t *target;
	band_t   *band;
	photdist_t *dist;
	canopy_t  canopy;

	/* 解析条件 */
	double    photoperiod;          /* [h/day] */
	int       maxiter;
	double    converg;
	int       msub;                 /* パッチ直接光の複合求積分割数 (0 = 自動) */
	double   *action;               /* [nlam] 光合成作用曲線 (McCree) */
	double   *vlambda;              /* [nlam] 標準比視感度 V(λ) */

	/* 導出量 */
	int       npatch;
	patch_t  *patch;
	double   *ff;                   /* [npatch*npatch] 形態係数 F_ij */
	float    *plen;                 /* [npatch*npatch] 群落内の経路長 [m] */
	double   *Ed;                   /* [npatch*nlam] 直接放射照度 [W/m²] */
	double   *B;                    /* [npatch*nlam] 放射発散度 [W/m²] */
	double   *Einc;                 /* [npatch*nlam] 入射放射照度 (直接+間接) */

	/* 診断 */
	double    ff_rowsum_err;        /* 行和 Σ_j F_ij の 1 からの最大偏差 */
	double    ff_recip_err;         /* 相反則 A_i F_ij = A_j F_ji の最大相対誤差 */
	int       niter;
	double    resid;
	double    flux_total;           /* 全光源の放射束 [W] */
	double    ppf_total;            /* 全光源の光量子束 PPF [µmol/s] */
	double    watt_total;           /* 全光源の消費電力 [W] */
	double    absorb_wall;          /* 壁が吸収した放射束 [W] */
	double    absorb_canopy;        /* 群落が吸収した放射束 [W] = 収支の残り */
	double    ppf_wall;             /* 壁が吸収した PPF [µmol/s] */
	double    canopy_abs;           /* 群落が吸収した PPF [µmol/s] */

	FILE     *fplog;
} ppfd_t;

/* ---- プロトタイプ -------------------------------------------------- */

/* utils.c */
int     tokenize(char *, const char *, char **, int);
int     streq_ci(const char *, const char *);
void    *xmalloc(size_t);
void    *xcalloc(size_t, size_t);
double  cputime(void);
void    plog(ppfd_t *, const char *, ...);

/* vector (utils.c) */
vec3_t  v_make(double, double, double);
vec3_t  v_add(vec3_t, vec3_t);
vec3_t  v_sub(vec3_t, vec3_t);
vec3_t  v_scale(vec3_t, double);
double  v_dot(vec3_t, vec3_t);
vec3_t  v_cross(vec3_t, vec3_t);
double  v_norm(vec3_t);
vec3_t  v_unit(vec3_t);

/* spectrum.c */
void    spectrum_grid(ppfd_t *);
void    spectrum_builtin(ppfd_t *);
double  interp_table(const double *, const double *, int, double);
void    spec_normalize(double *, int);
double  band_weight(const ppfd_t *, int, double, double);
double  calc_ppfd(const ppfd_t *, const double *);
double  calc_ypfd(const ppfd_t *, const double *);
double  calc_epar(const ppfd_t *, const double *);
double  calc_lux(const ppfd_t *, const double *);
double  calc_bandppf(const ppfd_t *, const double *, double, double);
double  calc_irradiance(const ppfd_t *, const double *);

/* photometry.c */
int     photdist_load(ppfd_t *, const char *, int);
double  photdist_value(const photdist_t *, double, double);

/* input_data.c */
int     input_data(FILE *, ppfd_t *);
int     find_mat(const ppfd_t *, const char *);
int     find_spec(const ppfd_t *, const char *);

/* geometry.c */
void    setup_patch(ppfd_t *);
void    setup_target(ppfd_t *);
double  canopy_path(const ppfd_t *, vec3_t, vec3_t);
void    canopy_trans(const ppfd_t *, double, double *);

/* formfactor.c */
double  ff_point_poly(vec3_t, vec3_t, const vec3_t *, int);
void    setup_ff(ppfd_t *);

/* direct.c */
void    direct_patch(ppfd_t *);
void    direct_point(const ppfd_t *, vec3_t, vec3_t, double *);

/* radiosity.c */
void    solve_radiosity(ppfd_t *);

/* solve.c */
void    solve(ppfd_t *);

/* output.c */
void    output_log(ppfd_t *);
void    output_csv(ppfd_t *);

#endif  /* _PPFD_H_ */
