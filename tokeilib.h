/**
* @file
* @brief 簡単な統計処理を実行するためのライブラリのヘッダーファイル
* @author Kenjiro Sugio
* @date 2025/11/19
* @copyright Kenjiro Sugio
*/

#ifndef _TOKEILIB_H
#define _TOKEILIB_H

#define TABLE_MAX 1000 /* 行の最大数　*/
#define VALS_MAX 50    /* 1行あたりのデータ最大数 */
#define BINS_MAX 100   /* ヒストグラムのビンの最大数 */

/* TABLE 構造体の定義 */
typedef struct _TABLE {
	char name[50];
	float vals[VALS_MAX];
} TABLE;

/* プロトタイプ宣言 */
int read_csv(char [], char **, TABLE [], int *);
void disp_table(int, int, char **, TABLE []);
int search_vid(char [], char **, int);
double calc_mean(int, TABLE [], int);
double calc_min(int, TABLE [], int);
double calc_max(int, TABLE [], int);
double calc_var(int, TABLE [], int, int);
double calc_std(int, TABLE [], int, int);
void bubble_sort(int, float []);
int qsort_comp(const void *, const void *);
double calc_median(int, TABLE [], int);
int summary_table(int, int, TABLE [], TABLE [], int);
void histogram(int, TABLE [], int, int);
void calc_score(int, int, TABLE [], TABLE [], int, float, float);
double pearson_corr(int, TABLE [], int, int, int);
int corr_table(int, int, char **, TABLE [], TABLE [], int);
void mat_trans(double *, int, int, double *);
void mat_multi(double *, int, int, double *, int, double *);
int mat_inverse(double *, int, int, double *);
int mat_regression(double *, int, int, double *, double *);
double regression(int, TABLE [], int, int [], double *);
void disp_coef(int, int [], char **, double *, double);

#endif
