/* -------------------------------------------------------------------------
* SUPREME_GraphPlot.cpp : Plot Output Results
*
*           Copyright (C) by C. Zhao 2016-2018 All rights reserved
*           Contact zhaoC.@whigg.ac.cn
*
* Update date: 2017.9.16
* Update date: 2017.12.4
* ------------------------------------------------------------------------- */
#ifndef SUPREME_GRAPH_PLOT_H_HH
#define SUPREME_GRAPH_PLOT_H_HH

#include "SUPREME_Options.h"

using namespace std;

#define GNULINELENGTH            100         /* gnu file path length                  */
#define GNUFILELENGTH            200         /* gnu input and output file paht length */
#define GNULINETYPE_POINT        0           /* gnu line type: points                 */
#define GNULINETYPE_LINE         1           /* gnu line type: lines                  */
#define GNULINETYPE_LINEPOINT    2           /* gnu line type: linepoints             */
#define GNUGRAPHTYPE_LINE        0           /* gnu graph type: broken line graph     */
#define GNUGRAPHTYPE_HISTOGRAM   1           /* gnu graph type: histogram graph       */
#define GNUOFF                   0           /* gnu off mode                          */
#define GNUON                    1           /* gnu on mode                           */

#define GNU_ORBIT                0
#define GNU_CLOCK                1


static int PLOT_SAT_FLAG[100];

const char GnuplotColor[50][20]
{
	"black", "red", "dark-red", "green", "blue", "slateblue1", "cyan", "magenta", "gold", "goldenrod",
		"navy", "purple", "mediumpurple3", "pink", "brown", "dark-khaki", "olive", "orchid", "orchid4", "plum",
		"dark-plum", "dark-green", "web-blue", "yello4", "orange", "salmon", "plum", "coral", "steelblue", "dark-blue",
		"dark-olivegreen", "spring-green", "dark-orange", "dark-pink", "bisque", "tan1", "lignt-goldenrod", "dark-violet", "royablue", "wine",
		"", "", "", "", "", "", "", "", "", "",
};

typedef struct    /* gnuplot config type */
{
	int data_num;                 /* data number needed to plot */
	int graph_type;               /* graph type: lines,histogram */

	int width;                    /* graph width */
	int length;                   /* graph length */

	int font_size;                /* graph font size */
	int tfont_size;               /* title font size */

	int xfont_size;               /* x label font size */
	int yfont_size;               /* y label font size */

	int xgrid;                    /* x label grid mode 0:off 1:on */
	int ygrid;                    /* y label grid mode 0:off 1:on */
	int grid;                     /* grid mode 0:off 1:on */

	int xdata_time;               /* x data time mode 0:off 1:on */
	int ydata_time;               /* y data time mode 0:off 1:on */

	int xrotate;                  /* x rotate degree */

	double xtics;                 /* x label tics */
	double ytics;                 /* y label tics */

	double ymin;                  /* y min value */
	double ymax;                  /* y max value */
	double xmin;                  /* x min value */
	double xmax;                  /* y max value */

	int data1_column;             /* first data column index */
	int data2_column;             /* second data column index */
	int data3_column;             /* third data column index */

	int data1_lt;                  /* data 1 line type */
	int data2_lt;                  /* data 2 line type */
	int data3_lt;                  /* data 3 line type */
	int data1_lw;                  /* data 1 line width */
	int data2_lw;                  /* data 2 line width */
	int data3_lw;                  /* data 3 line width */

	char data1_inputf[PATH_LENGTH];    /* x input file path */
	char data2_inputf[PATH_LENGTH];    /* y input file path */
	char data3_inputf[PATH_LENGTH];    /* z input file path */
	char outputf[PATH_LENGTH];         /* output graph path */

	char title[GNULINELENGTH];      /* graph title */
	char xtitle[GNULINELENGTH];     /* x title */
	char ytitle[GNULINELENGTH];     /* y title */

	char time_format[GNULINELENGTH]; /* time format */
	char xformat[GNULINELENGTH];     /* x label format */
	char yformat[GNULINELENGTH];     /* y label format */

	char data1_title[GNULINELENGTH]; /* data 1 title */
	char data2_title[GNULINELENGTH]; /* data 2 title */
	char data3_title[GNULINELENGTH]; /* data 3 title */

}gnu_cfg_t;

extern const gnu_cfg_t gnucfg_default;
extern const gnu_cfg_t gnucfg_ppp;
extern const gnu_cfg_t gnucfg_isb;
extern const gnu_cfg_t gnucfg_orb_acc;
extern const gnu_cfg_t gnucfg_clk_acc;

int create_gnuscrp(gnu_cfg_t *gnucfg, char* scrpname);

int Plot_SUPREME_Post_Results(ppp_option_t *popt);


#endif
