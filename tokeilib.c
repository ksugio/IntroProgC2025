/**
* @file
* @brief 簡単な統計処理を実行するためのライブラリ
* @author Kenjiro Sugio
* @date 2025/11/19
* @copyright Kenjiro Sugio
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tokeilib.h"

/**
* @brief CSVファイルを読み込む関数
* @param[in] filename ファイル名
* @param[out] header 列名
* @param[out] table TABLE構造体の配列
* @param[out] ncol 列数
* @return bool 読み込み成功：1，読み込み失敗：0
* @details ファイル形式は1行目が列名（header），1列目が名前(name)とすること
*/
int read_csv(char filename[], char **header, TABLE table[], int *ncol) {
    FILE *fp;
    char s[256], *p;
    static char h[256];
    int i = 0, j;
    fp = fopen(filename, "r");
    if(fp == NULL) return 0;
    /* ヘッダーの読み込み */
    fgets(h, 256, fp);
    p = strtok(h, ",\n");
    j = 0;
    while (p != NULL) {
        header[j] = p;
        p = strtok(NULL, ",\n");
        j++;
    }
    *ncol = j;
    /* データの読み込み */
    while(1) {
        fgets(s, 256, fp);
        if(feof(fp)) break;
        p = strtok(s, ",\n");
        j = 0;
        while (p != NULL) {
            if (j == 0) strcpy(table[i].name, p);
            else table[i].vals[j-1] = atof(p);
            p = strtok(NULL, ",\n");
            j++;
        }
        i++;
    }
    fclose(fp);
    return i;
}

/**
* @brief テーブルを表示する関数
* @param[in] num 行数
* @param[in] ncol 列数
* @param[in] header 列名
* @param[in] table TABLE構造体の配列
*/
void disp_table(int num, int ncol, char **header, TABLE table[]) {
    int i, j;
    for (j = 0; j < ncol; j++)
        printf("%10s ", header[j]);
    printf("\n");
    for (i = 0; i < num; i++) {
        printf("%10s ", table[i].name);
        for (j = 0; j < ncol - 1; j++) {
            printf("%10.4f ", table[i].vals[j]);
        }
        printf("\n");
    }
}

/**
* @brief 列名を検索してTABLE構造体.valsのインデックスを返す関数
* @param[in] arg 検索名
* @param[in] header 列名
* @param[in] ncol 列数
* @return TABLE構造体.valsのインデックス，見つからない場合は-1を返す
*/
int search_vid(char arg[], char **header, int ncol) {
    int i;
    for (i = 1; i < ncol; i++)
        if (strcmp(header[i], arg) == 0) return i - 1;
    return -1;
}

/**
* @brief 指定した行の平均値を計算する関数
* @param[in] num 行数
* @param[in] table TABLE構造体の配列
* @param[in] vid TABLE構造体.valsのインデックス
* @return vid で指定した行の平均値
*/
double calc_mean(int num, TABLE table[], int vid) {
    int i;
    double sum = table[0].vals[vid];
    for (i = 1; i < num; i++) {
        sum += table[i].vals[vid];
    }
    return sum / num;
}

/**
* @brief 指定した行の最小値を計算する関数
* @param[in] num 行数
* @param[in] table TABLE構造体の配列
* @param[in] vid TABLE構造体.valsのインデックス
* @return vid で指定した行の最小値
*/
double calc_min(int num, TABLE table[], int vid) {
    int i;
    double v, min = table[0].vals[vid];
    for (i = 1; i < num; i++) {
        v = table[i].vals[vid];
        if (v < min) min = v;
    }
    return min;
}

/**
* @brief 指定した行の最大値を計算する関数
* @param[in] num 行数
* @param[in] table TABLE構造体の配列
* @param[in] vid TABLE構造体.valsのインデックス
* @return vid で指定した行の最大値
*/
double calc_max(int num, TABLE table[], int vid) {
    int i;
    double v, max = table[0].vals[vid];
    for (i = 1; i < num; i++) {
        v = table[i].vals[vid];
        if (v > max) max = v;
    }
    return max;
}

/**
* @brief 指定した行の分散を計算する関数
* @param[in] num 行数
* @param[in] table TABLE構造体の配列
* @param[in] vid TABLE構造体.valsのインデックス
* @param[in] unbias 0:標本分散，1:不偏分散
* @return vid で指定した行の分散
*/
double calc_var(int num, TABLE table[], int vid, int unbias) {
    int i;
    double mean = calc_mean(num, table, vid);
    double v, var = 0.0;
    for (i = 0; i < num; i++) {
        v = table[i].vals[vid];
        var += (v - mean) * (v - mean);
    }
    if (unbias) return var / (num - 1);
    else return var / num;
}

/**
* @brief 指定した行の標準偏差を計算する関数
* @param[in] num 行数
* @param[in] table TABLE構造体の配列
* @param[in] vid TABLE構造体.valsのインデックス
* @param[in] unbias 0:標準偏差，1:不偏標準偏差
* @return vid で指定した行の標準偏差
*/
double calc_std(int num, TABLE table[], int vid, int unbias) {
    double var = calc_var(num, table, vid, unbias);
    return sqrt(var);
}

/**
* @brief バブルソートで昇順に並べ替える関数
* @param[in] num 行数
* @param[in] data 配列
*/
void bubble_sort(int num, float data[]) {
    int i, j;
    float temp;
    for (i = 0; i < num - 1; i++) {
        for (j = num - 1; j > i; j--) {
            if (data[j-1] > data[j]) {
                temp = data[j-1];
                data[j-1] = data[j];
                data[j] = temp;
            }
        }
    }
}

/**
* @brief クイックソート用の比較関数（昇順に並べ替える）
* @param[in] a 比較する値のポインタ1
* @param[in] b 比較する値のポインタ2
* @return a > b ならば 1，a < b ならば -1，a = b ならば 0 を返す
* @details float型にキャストして比較する
*/
int qsort_comp(const void* a, const void* b)
{
    if (*(float*)a > *(float*)b) return 1;
    else if (*(float*)a < *(float*)b) return -1;
    else return 0;
}

/**
* @brief 指定した行の中央値を計算する関数
* @param[in] num 行数
* @param[in] table TABLE構造体の配列
* @param[in] vid TABLE構造体.valsのインデックス
* @return vid で指定した行の中央値
*/
double calc_median(int num, TABLE table[], int vid) {
    float data[TABLE_MAX];
    int i, id;
    for (i = 0; i < num; i++) {
        data[i] = table[i].vals[vid];
    }
    //bubble_sort(num, v);
    qsort(data, num, sizeof(float), qsort_comp);
    id = num / 2;
    if (num % 2 == 0)
        return (data[id - 1] + data[id]) / 2.0;
    else
        return data[id];
}

/**
* @brief 全ての行の統計値を計算してまとめテーブルを作成する関数
* @param[in] num 行数
* @param[in] ncol 列数
* @param[in] table TABLE構造体の配列
* @param[out] summary TABLE構造体の配列，まとめテーブル
* @param[in] unbias 0:標本分散および標準偏差，1:不偏分散および不偏標準偏差
* @return まとめテーブルの行数
*/
int summary_table(int num, int ncol, TABLE table [], TABLE summary [], int unbias)
{
    int j;
    strcpy(summary[0].name, "Mean");
    for (j = 0; j < ncol - 1; j++)
        summary[0].vals[j] = calc_mean(num, table, j);
    strcpy(summary[1].name, "Var");
    for (j = 0; j < ncol - 1; j++)
        summary[1].vals[j] = calc_var(num, table, j, unbias);
    strcpy(summary[2].name, "Std");
    for (j = 0; j < ncol - 1; j++)
        summary[2].vals[j] = calc_std(num, table, j, unbias);
    strcpy(summary[3].name, "Min");    
    for (j = 0; j < ncol - 1; j++)
        summary[3].vals[j] = calc_min(num, table, j);
    strcpy(summary[4].name, "Median");    
    for (j = 0; j < ncol - 1; j++)
        summary[4].vals[j] = calc_median(num, table, j);
    strcpy(summary[5].name, "Max");    
    for (j = 0; j < ncol - 1; j++)
        summary[5].vals[j] = calc_max(num, table, j);
    return 6;
}

/**
* @brief 指定した行のヒストグラムを表示する関数
* @param[in] num 行数
* @param[in] table TABLE構造体の配列
* @param[in] vid TABLE構造体.valsのインデックス
* @param[in] bins ビンの数
*/
void histogram(int num, TABLE table[], int vid, int bins) {
    double min = calc_min(num, table, vid);
    double max = calc_max(num, table, vid);
    double dv = (max - min) / bins;
    int count[BINS_MAX];
    int i, j, id;    
    for (i = 0; i < bins; i++) count[i] = 0;
    for (i = 0; i < num; i++) {
        id = (table[i].vals[vid] - min) / dv;
        if (id >= bins) id = bins - 1; 
        count[id]++;
    }
    for (i = 0; i < bins; i++) {
        printf("%10.4f %4d ", i * dv + min, count[i]);
        for (j = 0; j < count[i]; j++) printf("*");
        printf("\n");
    }
    printf("%10.4f\n", max);
}

/**
* @brief 全ての行を標準化する関数
* @param[in] num 行数
* @param[in] ncol 列数
* @param[in] table TABLE構造体の配列
* @param[out] score TABLE構造体の配列，標準化されたテーブル
* @param[in] unbias 0:標本分散および標準偏差，1:不偏分散および不偏標準偏差
* @param[in] mu 標準化されたデータの平均
* @param[in] sigma 標準化されたデータの標準偏差
* @details muを50，sigmaを10とした場合，一般的に使われる偏差値が計算される
*/
void calc_score(int num, int ncol, TABLE table[], TABLE score[], int unbias, float mu, float sigma) {
    double mean, std;
    int i, j;
    for (i = 0; i < num; i++) {
        strcpy(score[i].name, table[i].name);
    }
    for (j = 0; j < ncol - 1; j++) {
        mean = calc_mean(num, table, j);
        std = calc_std(num, table, j, unbias);
        for (i = 0; i < num; i++) {
            score[i].vals[j] = (table[i].vals[j] - mean) / std * sigma + mu;
        }
    }
}

/**
* @brief 指定した２つの行の共分散を計算する関数
* @param[in] num 行数
* @param[in] table TABLE構造体の配列
* @param[in] vid1 TABLE構造体.valsのインデックス１
* @param[in] vid2 TABLE構造体.valsのインデックス２
* @param[in] unbias 0:標本分散，1:不偏分散
* @return vid1, vid2 で指定した行の共分散
*/
double calc_covar(int num, TABLE table[], int vid1, int vid2, int unbias) {
    int i;
    double mean1 = calc_mean(num, table, vid1);
    double mean2 = calc_mean(num, table, vid2);
    double var = 0.0;
    for (i = 0; i < num; i++) {
        var += (table[i].vals[vid1] - mean1)
            * (table[i].vals[vid2] - mean2);
    }
    if (unbias) return var / (num - 1);
    else return var / num;
}

/**
* @brief 指定した２つの行のPearson相関係数を計算する関数
* @param[in] num 行数
* @param[in] table TABLE構造体の配列
* @param[in] vid1 TABLE構造体.valsのインデックス１
* @param[in] vid2 TABLE構造体.valsのインデックス２
* @param[in] unbias 0:標本分散，1:不偏分散
* @return vid1, vid2 で指定した行のPearson相関係数
*/
double pearson_corr(int num, TABLE table[], int vid1, int vid2, int unbias) {
    double std1 = calc_std(num, table, vid1, unbias);
    double std2 = calc_std(num, table, vid2, unbias);
    double covar12 = calc_covar(num, table, vid1, vid2, unbias);
    return covar12 / (std1 * std2);
}

/**
* @brief 全ての行の組み合わせの相関係数を計算して相関テーブルを作成する関数
* @param[in] num 行数
* @param[in] ncol 列数
* @param[in] header 列名
* @param[in] table TABLE構造体の配列
* @param[out] corre TABLE構造体の配列，相関テーブル
* @param[in] unbias 0:標本分散，1:不偏分散
* @return 相関テーブルの行数
*/
int corr_table(int num, int ncol, char **header, TABLE table [], TABLE corre [], int unbias) {
    int i, j;
    for (i = 0; i < ncol - 1; i++) {
        strcpy(corre[i].name, header[i + 1]);
        for (j = 0; j < ncol - 1; j++) {
            corre[i].vals[j] = pearson_corr(num, table, i, j, unbias);
        }
    }
    return i;
}

/**
* @brief 1次元配列を2次元配列として取り扱うためのマクロ
* @param[in] mat 配列の先頭ポインタ
* @param[in] ncol 列数
* @param[in] i 行インデックス
* @param[in] j 列インデックス
* @return 配列のポインタが示す実体
*/
#define MATRIX(mat, ncol, i, j) *((mat)+(i)*(ncol)+(j))

/**
* @brief 転置行列を作成する関数
* @param[in] mat 行列の先頭ポインタ
* @param[in] m 行数
* @param[in] n 列数
* @param[out] matt 転置行列の先頭ポインタ
*/
void mat_trans(double *mat, int m, int n, double *matt) {
    int i, j;
    for (i = 0; i < n; i++)
        for (j = 0;j < m; j++)
            MATRIX(matt, m, i, j) = MATRIX(mat, n, j, i);
}

/**
* @brief 行列A，行列Bのかけ算（A×B）を計算する関数
* @param[in] mata 行列Aの先頭ポインタ
* @param[in] m 行数
* @param[in] n 行列Aの列数
* @param[in] matb 行列Bの先頭ポインタ
* @param[in] p 行列Bの列数
* @param[out] matm かけ算された行列の先頭ポインタ
*/
void mat_multi(double *mata, int m, int n, double *matb, int p, double *matm) {
    int i, j, k;
    for (i = 0; i < m; i++) {
        for (j = 0; j < p; j++) {
            MATRIX(matm, p, i, j) = 0.0;
            for (k = 0; k < n; k++)
                MATRIX(matm, p, i, j) +=
                    MATRIX(mata, n, i, k) *
                    MATRIX(matb, p, k, j);
        }
    }
}

/**
* @brief ガウス・ジョルダン法により逆行列を計算する関数
* @param[in] mat 行列の先頭ポインタ
* @param[in] m 行数
* @param[in] n 列数
* @param[out] mati 逆行列の先頭ポインタ
* @return 0:対角要素に0が存在，1:成功
*/
int mat_inverse(double *mat, int m, int n, double *mati) {
    int i, j, k;
    int ne = 2 * n;
    double pivot, factor;
    // 拡大行列 [A | I] を作成
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            MATRIX(mati, ne, i, j) = MATRIX(mat, n, i, j);
            MATRIX(mati, ne, i, j + n) = (i == j) ? 1.0 : 0.0;
        }
    }
    // ガウス・ジョルダン法
    for (i = 0; i < m; i++) {
        pivot = MATRIX(mati, ne, i, i);
        if (pivot == 0.0) return 0;
        for (j = 0; j < ne; j++) 
            MATRIX(mati, ne, i, j) /= pivot;
        for (k = 0; k < n; k++) {
            if (k != i) {
                factor = MATRIX(mati, ne, k, i);
                for (j = 0; j < ne; j++) {
                    MATRIX(mati, ne, k, j) -= factor * MATRIX(mati, ne, i, j);
                }
            }
        }
    }
    // 逆行列を抽出
    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
            MATRIX(mati, n, i, j) = MATRIX(mati, ne, i, j + n);
    return 1;
}

/**
* @brief 重回帰分析を実行する関数
* @param[in] feat 特徴行列の先頭ポインタ
* @param[in] m 行数
* @param[in] n 列数
* @param[in] obj 目的配列の先頭ポインタ
* @param[in] coef 回帰係数配列の先頭ポインタ
* @return 0:失敗，1:成功
*/
int mat_regression(double *feat, int m, int n, double *obj, double *coef) {
    double *matt, *matm, *mati;
    int res;
    matt = (double *)malloc(sizeof(double) * m * n);
    matm = (double *)malloc(sizeof(double) * m * n);
    mati = (double *)malloc(sizeof(double) * 2 * m * n);
    mat_trans(feat, m, n, matt);
    mat_multi(matt, n, m, feat, n, matm);
    res = mat_inverse(matm, n, n, mati);
    if (res) {
        mat_multi(mati, n, n, matt, m, matm);
        mat_multi(matm, n, m, obj, 1, coef);
    }
    free(matt);
    free(matm);
    free(mati);
    return res;
}

/**
* @brief 決定係数を計算する関数
* @param[in] num データ数
* @param[in] y 真値配列の先頭ポインタ
* @param[in] y_pred 予測値配列の先頭ポインタ
* @return 決定係数
*/
double r2_score(int num, double *y, double *y_pred) {
    int i;
    double mean = 0.0;
    double tss = 0.0;
    double rss = 0.0;
    for (i = 0; i < num; i++) mean += y[i];
    mean /= num;
    for (i = 0; i < num; i++) 
        tss += (y[i] - mean) * (y[i] - mean);
    for (i = 0; i < num; i++) 
        rss += (y[i] - y_pred[i]) * (y[i] - y_pred[i]);
    return 1.0 - rss / tss;
}

/**
* @brief 目的変数および説明変数を指定して重回帰分析を実行する
* @param[in] num 行数
* @param[in] table TABLE構造体の配列
* @param[in] nvids TABLE構造体.valsのインデックス配列の個数
* @param[in] vids TABLE構造体.valsのインデックス配列
* @param[in] coef 回帰係数配列の先頭ポインタ
* @return 決定係数
* @details vids[0]が目的変数，それ以外が説明変数として使用される
*/
double regression(int num, TABLE table[], int nvids, int vids[], double *coef) {
    double *feat, *obj, *pred;
    int m = num, n = nvids - 1;
    int i, j;
    double r2 = 0.0;
    feat = (double *)malloc(sizeof(double) * m * n);
    obj = (double *)malloc(sizeof(double) * m);
    pred = (double *)malloc(sizeof(double) * m);
    for (i = 0; i < num; i++) {
        obj[i] = table[i].vals[vids[0]];
        for (j = 1; j < nvids; j++)
            MATRIX(feat, n, i, j - 1) = table[i].vals[vids[j]];
    }
    if (mat_regression(feat, m, n, obj, coef)) {
        for (i = 0; i < m; i++) {
            pred[i] = 0.0;
            for (j = 0; j < n; j++)
                pred[i] += MATRIX(feat, n, i, j) * coef[j];
        }
        r2 = r2_score(m, obj, pred);
    }
    free(feat);
    free(obj);
    free(pred);
    return r2;
}

/**
* @brief 回帰係数および決定係数を表示する関数
* @param[in] nvids TABLE構造体.valsのインデックス配列の個数
* @param[in] vids TABLE構造体.valsのインデックス配列
* @param[in] header 列名
* @param[in] coef 回帰係数配列の先頭ポインタ
* @param[in] r2 決定係数
*/
void disp_coef(int nvids, int vids[], char **header, double *coef, double r2) {
    int i;
    for (i = 1; i < nvids; i++)
        printf("%10s", header[vids[i] + 1]);
    printf(" |         R2\n");
    for (i = 0; i < nvids - 1; i++)
        printf("%10.4f", coef[i]);
    printf(" | %10.4f", r2);
}