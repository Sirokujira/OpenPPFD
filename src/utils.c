/*
utils.c

トークン化・ログ出力・ベクトル演算。
*/

#include "ppfd.h"
#include <stdarg.h>
#include <ctype.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* 空白区切りのトークン化 (strline は破壊される) */
int tokenize(char *str, const char *sep, char **token, int maxtoken)
{
	int n = 0;
	char *p = strtok(str, sep);
	while ((p != NULL) && (n < maxtoken)) {
		token[n++] = p;
		p = strtok(NULL, sep);
	}
	return n;
}

/* 大文字小文字を無視した比較 (strcasecmp は MSVC に無い) */
int streq_ci(const char *a, const char *b)
{
	while (*a && *b) {
		if (tolower((unsigned char)*a) != tolower((unsigned char)*b)) return 0;
		a++;
		b++;
	}
	return (*a == '\0') && (*b == '\0');
}

void *xmalloc(size_t n)
{
	void *p = malloc(n ? n : 1);
	if (p == NULL) {
		fprintf(stderr, "*** out of memory (%lu bytes)\n", (unsigned long)n);
		exit(1);
	}
	return p;
}

void *xcalloc(size_t n, size_t sz)
{
	void *p = calloc(n ? n : 1, sz);
	if (p == NULL) {
		fprintf(stderr, "*** out of memory (%lu x %lu bytes)\n", (unsigned long)n, (unsigned long)sz);
		exit(1);
	}
	return p;
}

double cputime(void)
{
#ifdef _OPENMP
	return omp_get_wtime();
#else
	return (double)clock() / CLOCKS_PER_SEC;
#endif
}

/* 標準出力とログファイルの両方へ */
void plog(ppfd_t *p, const char *fmt, ...)
{
	char buf[BUFSIZ];
	va_list ap;
	va_start(ap, fmt);
	vsnprintf(buf, sizeof(buf), fmt, ap);
	va_end(ap);

	fputs(buf, stdout);
	if (p->fplog != NULL) {
		fputs(buf, p->fplog);
		fflush(p->fplog);
	}
}

/* ---- ベクトル ---------------------------------------------------- */

vec3_t v_make(double x, double y, double z)
{
	vec3_t v;
	v.x = x;
	v.y = y;
	v.z = z;
	return v;
}

vec3_t v_add(vec3_t a, vec3_t b)
{
	return v_make(a.x + b.x, a.y + b.y, a.z + b.z);
}

vec3_t v_sub(vec3_t a, vec3_t b)
{
	return v_make(a.x - b.x, a.y - b.y, a.z - b.z);
}

vec3_t v_scale(vec3_t a, double s)
{
	return v_make(a.x * s, a.y * s, a.z * s);
}

double v_dot(vec3_t a, vec3_t b)
{
	return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
}

vec3_t v_cross(vec3_t a, vec3_t b)
{
	return v_make((a.y * b.z) - (a.z * b.y),
	              (a.z * b.x) - (a.x * b.z),
	              (a.x * b.y) - (a.y * b.x));
}

double v_norm(vec3_t a)
{
	return sqrt(v_dot(a, a));
}

vec3_t v_unit(vec3_t a)
{
	const double d = v_norm(a);
	return (d > EPS) ? v_scale(a, 1.0 / d) : v_make(0, 0, 1);
}
