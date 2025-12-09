/**
* @file
* @brief 簡単な統計処理を実行するためのコマンドライン・ツール
* @author Kenjiro Sugio
* @date 2025/11/19
* @copyright Kenjiro Sugio
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tokeilib.h"

char Help[] = "\
tokei [filename] [command] [op1] [op2] [op3] ...\n\
[command]\n\
    disp : display table\n\
    mean : calculate mean [op1:header]\n\
    min : calculate minimum [op1:header]\n\
    max : calculate maximum [op1:header]\n\
    med : calculate median [op1:header]\n\
    var : calculate variance [op1:header] [op2:unbias(0|1)]\n\
    std : calculate standard deviation [op1:header] [op2:unbias(0|1)]\n\
    summ : display summary [op1:unbias(0|1)]\n\
    hist : display histgram [op1:header] [op2:bins]\n\
    score : display score [op1:unbias(0|1)] [op2:mu] [op3:sigma]\n\
    corr : calculate pearson correlation [op1:header1] [op2:header2]\n\
    cmap : display pearson correlation map\n\
    regr : multiple regression analysis [op1:obj_header] [op2:feat_header1] [op3:feat_header2] ...\n\
";

int main(int argc, char* argv[])
{
    char *header[VALS_MAX+1];
    TABLE table[TABLE_MAX];
    TABLE summary[10];
    TABLE score[TABLE_MAX];
    TABLE corre[VALS_MAX];
    double coef[VALS_MAX], r2;
    int num, ncol, nnum;
    int vid, vid1, vid2;
    int nvids, vids[VALS_MAX];
    int i;
    if (argc >= 3) {
        num = read_csv(argv[1], header, table, &ncol);
        if (num == 0) printf("Can't read %s.", argv[1]);
    }
    if (argc == 3 && strcmp(argv[2], "disp") == 0) {
        disp_table(num, ncol, header, table);
    }
    else if (argc == 4 && strcmp(argv[2], "mean") == 0) {
        vid = search_vid(argv[3], header, ncol);
        if (vid >= 0) printf("Mean %f\n",
            calc_mean(num, table, vid));
    }
    else if (argc == 4 && strcmp(argv[2], "min") == 0) {
        vid = search_vid(argv[3], header, ncol);
        if (vid >= 0) printf("Min %f\n",
            calc_min(num, table, vid));
    }
    else if (argc == 4 && strcmp(argv[2], "max") == 0) {
        vid = search_vid(argv[3], header, ncol);
        if (vid >= 0) printf("Max %f\n",
            calc_max(num, table, vid));
    }
    else if (argc == 4 && strcmp(argv[2], "med") == 0) {
        vid = search_vid(argv[3], header, ncol);
        if (vid >= 0) printf("Median %f\n",
            calc_median(num, table, vid));
    }
    else if (argc == 5 && strcmp(argv[2], "var") == 0) {
        vid = search_vid(argv[3], header, ncol);
        if (vid >= 0) printf("Var %f\n",
            calc_var(num, table, vid, atoi(argv[4])));
    }
    else if (argc == 5 && strcmp(argv[2], "std") == 0) {
        vid = search_vid(argv[3], header, ncol);
        if (vid >= 0) printf("Std %f\n",
            calc_std(num, table, vid, atoi(argv[4])));
    }
    else if (argc == 4 && strcmp(argv[2], "summ") == 0) {
        nnum = summary_table(num, ncol, table, summary, atoi(argv[3]));
        header[0] = "";
        disp_table(nnum, ncol, header, summary);
    }
    else if (argc == 5 && strcmp(argv[2], "hist") == 0) {
        vid = search_vid(argv[3], header, ncol);
        if (vid >= 0)
            histogram(num, table, vid, atoi(argv[4]));
    }
    else if (argc == 6 && strcmp(argv[2], "score") == 0) {
        calc_score(num, ncol, table, score, atoi(argv[3]), atof(argv[4]), atof(argv[5]));
        disp_table(num, ncol, header, score);
    }
    else if (argc == 5 && strcmp(argv[2], "corr") == 0) {
        vid1 = search_vid(argv[3], header, ncol);
        vid2 = search_vid(argv[4], header, ncol);
        if (vid1 >= 0 && vid2 >= 0) printf("Corr %f\n",
            pearson_corr(num, table, vid1, vid2, 0));
    }
    else if (argc == 3 && strcmp(argv[2], "cmap") == 0) {
        nnum = corr_table(num, ncol, header, table, corre, 0);
        header[0] = "";
        disp_table(nnum, ncol, header, corre);
    }
    else if (argc >= 5 && strcmp(argv[2], "regr") == 0) {
        nvids = argc - 3;
        for (i = 0; i < nvids; i++) {
            vids[i] = search_vid(argv[i + 3], header, ncol);
            if (vids[i] == -1) {
                printf("No such header : %s\n", argv[i + 3]);
                return 0;
            }
        }
        calc_score(num, ncol, table, score, 0, 1.0, 1.0);
        r2 = regression(num, score, nvids, vids, coef);
        disp_coef(nvids, vids, header, coef, r2);
    }
    else {
        printf("%s", Help);
    }
}
