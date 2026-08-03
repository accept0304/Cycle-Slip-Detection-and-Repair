#include <stdio.h>
#include <string.h>

#include "SUPREME_GraphPlot.h"
#include "SUPREME_CommonFunction.h"

const gnu_cfg_t gnucfg_default = { /* Defalult settings of gnuplot configuration */
	3, GNUGRAPHTYPE_LINE,          /* data number, gnu graph type */
	600, 1200,                     /* graph width,graph length */
	12, 16,                        /* graph font size, title font size */
	14, 14,                        /* x font size, y font size */
	GNUON, GNUON, GNUON,           /* xgrid mode, ygrid mode, grid mode */
	GNUON, GNUOFF,                 /* x data time mode, y data time mode */
	0, 0, 0,                       /* x rotate degree, xtics, ytics */
	0, 0, 0, 0,                    /* ymin, ymax, xmin, xmax */
	0, 0, 0,                       /* data1 column, data2 column, data3 column */
	GNULINETYPE_LINE,              /* data1 line type: lines */
	GNULINETYPE_LINE,              /* data2 line type: lines */
	GNULINETYPE_LINE,              /* data3 line type: lines */
	2, 2, 2,                       /* data1 line width, data2 line width, data3 line width */
	"", "", "",                    /* data1,data2,data3 input file */
	"",                            /* graph output file path */
	"", "", "",                    /* graph title, x title, y title */
	"", "", "",                    /* time format, x format, y format */
	"", "", ""                     /* data1 label, data2 label, data3 label */
};

/// Orbit accuracy gnuplot config
const gnu_cfg_t gnucfg_orb_acc = { /* PPP default settings of gnuplot configuration */
	3, GNUGRAPHTYPE_HISTOGRAM,     /* data number, gnu graph type */
	600, 1200,                     /* graph width,graph length */
	12, 16,                        /* graph font size, title font size */
	14, 14,                        /* x font size, y font size */
	GNUOFF, GNUON, GNUON,          /* xgrid mode, ygrid mode, grid mode */
	GNUOFF, GNUOFF,                /* x data time mode, y data time mode */
	45, 3600, 0,                   /* x rotate degree, xtics, ytics */
	0, 10, 0, 0,                   /* ymin, ymax, xmin, xmax */
	5, 6, 7,                       /* data1 column, data2 column, data3 column */
	GNULINETYPE_LINE,              /* data1 line type: points */
	GNULINETYPE_LINE,              /* data2 line type: points */
	GNULINETYPE_LINE,              /* data3 line type: points */
	1, 1, 1,                       /* data1 line width, data2 line width, data3 line width */
	"", "", "",                    /* data1,data2,data3 input file */
	"",                            /* graph output file path */
	"",                            /* graph title */
	"SAT PRN",                     /* x title */
	"Accuracy / cm",               /* y title */
	"%Y:%m:%d:%H:%M:%S",           /* time format */
	"%H", "",                      /* x format, y format */
	"", "", ""                     /* data1 label, data2 label, data3 label */
};

/// Clock accuracy gnuplot config
const gnu_cfg_t gnucfg_clk_acc = { /* PPP default settings of gnuplot configuration */
	2, GNUGRAPHTYPE_HISTOGRAM,     /* data number, gnu graph type */
	600, 1200,                     /* graph width,graph length */
	12, 16,                        /* graph font size, title font size */
	14, 14,                        /* x font size, y font size */
	GNUOFF, GNUON, GNUON,          /* xgrid mode, ygrid mode, grid mode */
	GNUOFF, GNUOFF,                /* x data time mode, y data time mode */
	45, 3600, 0,                   /* x rotate degree, xtics, ytics */
	0, 0, 0, 0,                    /* ymin, ymax, xmin, xmax */
	2, 3, 4,                       /* data1 column, data2 column, data3 column */
	GNULINETYPE_LINE,              /* data1 line type: points */
	GNULINETYPE_LINE,              /* data2 line type: points */
	GNULINETYPE_LINE,              /* data3 line type: points */
	1, 1, 1,                       /* data1 line width, data2 line width, data3 line width */
	"", "", "",                    /* data1,data2,data3 input file */
	"",                            /* graph output file path */
	"",                            /* graph title */
	"SAT PRN",                     /* x title */
	"Accuracy / ns",               /* y title */
	"%Y:%m:%d:%H:%M:%S",           /* time format */
	"%H", "",                      /* x format, y format */
	"", "", ""                     /* data1 label, data2 label, data3 label */
};

/// PPP gnuplot config
const gnu_cfg_t gnucfg_ppp = {     /* PPP default settings of gnuplot configuration */
	3, GNUGRAPHTYPE_LINE,          /* data number, gnu graph type */
	600, 1200,                     /* graph width,graph length */
	12, 16,                        /* graph font size, title font size */
	14, 14,                        /* x font size, y font size */
	GNUOFF, GNUON, GNUON,          /* xgrid mode, ygrid mode, grid mode */
	GNUON, GNUOFF,                 /* x data time mode, y data time mode */
	0, 3600, 0.1,                  /* x rotate degree, xtics, ytics */
	-0.5, 0.5, 0, 0,               /* ymin, ymax, xmin, xmax */
	5, 3, 4,                       /* data1 column, data2 column, data3 column */
	//GNULINETYPE_POINT,             /* data1 line type: points */
	//GNULINETYPE_POINT,             /* data2 line type: points */
	//GNULINETYPE_POINT,             /* data3 line type: points */
	GNULINETYPE_LINE,              /* data1 line type: points */
	GNULINETYPE_LINE,              /* data2 line type: points */
	GNULINETYPE_LINE,              /* data3 line type: points */
	1.5, 1.5, 1.5,                 /* data1 line width, data2 line width, data3 line width */
	"", "", "",                    /* data1,data2,data3 input file */
	"",                            /* graph output file path */
	"",                            /* graph title */
	"Epoch-Time(GPST) / Hour",     /* x title */
	"Positioning Errors / m",      /* y title */
	"%Y:%m:%d:%H:%M:%S",           /* time format */
	"%H", "",                      /* x format, y format */
	"U", "N", "E"                  /* data1 label, data2 label, data3 label */
};

/// Orbit compare gnuplot config
const gnu_cfg_t gnucfg_orbit = {   /* PPP default settings of gnuplot configuration */
	3, GNUGRAPHTYPE_LINE,          /* data number, gnu graph type */
	600, 1200,                     /* graph width,graph length */
	12, 16,                        /* graph font size, title font size */
	14, 14,                        /* x font size, y font size */
	GNUOFF, GNUON, GNUON,          /* xgrid mode, ygrid mode, grid mode */
	GNUON, GNUOFF,                 /* x data time mode, y data time mode */
	0, 3600, 0,                    /* x rotate degree, xtics, ytics */
	0, 0, 0, 0,                    /* ymin, ymax, xmin, xmax */
	2, 3, 4,                       /* data1 column, data2 column, data3 column */
	GNULINETYPE_LINEPOINT,         /* data1 line type: points */
	GNULINETYPE_LINEPOINT,         /* data2 line type: points */
	GNULINETYPE_LINEPOINT,         /* data3 line type: points */
	1, 1, 1,                       /* data1 line width, data2 line width, data3 line width */
	"", "", "",                    /* data1,data2,data3 input file */
	"",                            /* graph output file path */
	"",                            /* graph title */
	"Epoch-Time(GPST) / Hour",     /* x title */
	"Orbit Difference / cm",       /* y title */
	"%Y:%m:%d:%H:%M:%S",           /* time format */
	"%H", "",                      /* x format, y format */
	"X", "Y", "Z"                  /* data1 label, data2 label, data3 label */
};

/// Clock compare gnuplot config
const gnu_cfg_t gnucfg_clock = {   /* PPP default settings of gnuplot configuration */
	3, GNUGRAPHTYPE_LINE,          /* data number, gnu graph type */
	600, 1200,                     /* graph width,graph length */
	12, 16,                        /* graph font size, title font size */
	14, 14,                        /* x font size, y font size */
	GNUOFF, GNUON, GNUON,          /* xgrid mode, ygrid mode, grid mode */
	GNUON, GNUOFF,                 /* x data time mode, y data time mode */
	0, 3600, 0,                    /* x rotate degree, xtics, ytics */
	0, 0, 0, 0,                    /* ymin, ymax, xmin, xmax */
	2, 3, 4,                       /* data1 column, data2 column, data3 column */
	GNULINETYPE_LINEPOINT,         /* data1 line type: points */
	GNULINETYPE_LINEPOINT,         /* data2 line type: points */
	GNULINETYPE_LINEPOINT,         /* data3 line type: points */
	1, 1, 1,                       /* data1 line width, data2 line width, data3 line width */
	"", "", "",                    /* data1,data2,data3 input file */
	"",                            /* graph output file path */
	"",                            /* graph title */
	"Epoch-Time(GPST) / Hour",     /* x title */
	"Clock Diff / ns",             /* y title */
	"%Y:%m:%d:%H:%M:%S",           /* time format */
	"%H", "",                      /* x format, y format */
	"X", "Y", "Z"                  /* data1 label, data2 label, data3 label */
};

/// Ionosphere slant delay gnuplot config
const gnu_cfg_t gnucfg_ion = {   /* PPP default settings of gnuplot configuration */
	3, GNUGRAPHTYPE_LINE,          /* data number, gnu graph type */
	600, 1200,                     /* graph width,graph length */
	12, 16,                        /* graph font size, title font size */
	14, 14,                        /* x font size, y font size */
	GNUOFF, GNUON, GNUON,          /* xgrid mode, ygrid mode, grid mode */
	GNUON, GNUOFF,                 /* x data time mode, y data time mode */
	0, 3600, 0,                    /* x rotate degree, xtics, ytics */
	0, 0, 0, 0,                    /* ymin, ymax, xmin, xmax */
	2, 3, 4,                       /* data1 column, data2 column, data3 column */
	GNULINETYPE_LINEPOINT,         /* data1 line type: points */
	GNULINETYPE_LINEPOINT,         /* data2 line type: points */
	GNULINETYPE_LINEPOINT,         /* data3 line type: points */
	1, 1, 1,                       /* data1 line width, data2 line width, data3 line width */
	"", "", "",                    /* data1,data2,data3 input file */
	"",                            /* graph output file path */
	"",                            /* graph title */
	"Epoch-Time(GPST) / Hour",     /* x title */
	"Ionosphere Slant Delay / TECu", /* y title */
	"%Y:%m:%d:%H:%M:%S",           /* time format */
	"%H", "",                      /* x format, y format */
	"", "", ""                     /* data1 label, data2 label, data3 label */
};

const gnu_cfg_t gnucfg_spp = {     /* SPP default settings of gnuplot configuration */
	3, GNUGRAPHTYPE_LINE,          /* data number, gnu graph type */
	400, 800,                      /* graph width,graph length */
	12, 16,                        /* graph font size, title font size */
	14, 14,                        /* x font size, y font size */
	GNUOFF, GNUON, GNUON,          /* xgrid mode, ygrid mode, grid mode */
	GNUON, GNUOFF,                 /* x data time mode, y data time mode */
	0, 7200, 1,                    /* x rotate degree, xtics, ytics */
	-10, 10, 0, 0,                 /* ymin, ymax, xmin, xmax */
	2, 3, 4,                       /* data1 column, data2 column, data3 column */
	GNULINETYPE_LINE,              /* data1 line type: lines */
	GNULINETYPE_LINE,              /* data2 line type: lines */
	GNULINETYPE_LINE,              /* data3 line type: lines */
	2, 2, 2,                       /* data1 line width, data2 line width, data3 line width */
	"", "", "",                    /* data1,data2,data3 input file */
	"",                            /* graph output file path */
	"SPP_Precision",               /* graph title */
	"Epoch-Time",                  /* x title */
	"Positon Errors /m",           /* y title */
	"%Y:%m:%d:%H:%M:%S",           /* time format */
	"%H:%M", "",                   /* x format, y format */
	"North", "East", "Up"          /* data1 label, data2 label, data3 label */
};

const gnu_cfg_t gnucfg_isb = {     /* ISB default settings of gnuplot configuration */
	3, GNUGRAPHTYPE_LINE,          /* data number, gnu graph type */
	600, 1200,                     /* graph width,graph length */
	12, 16,                        /* graph font size, title font size */
	14, 14,                        /* x font size, y font size */
	GNUOFF, GNUON, GNUON,          /* xgrid mode, ygrid mode, grid mode */
	GNUON, GNUOFF,                 /* x data time mode, y data time mode */
	0, 7200, 0.2,                  /* x rotate degree, xtics, ytics */
	-1, 1, 0, 0,                   /* ymin, ymax, xmin, xmax */
	2, 3, 4,                       /* data1 column, data2 column, data3 column */
	GNULINETYPE_LINE,              /* data1 line type: lines */
	GNULINETYPE_LINE,              /* data2 line type: lines */
	GNULINETYPE_LINE,              /* data3 line type: lines */
	1, 1, 1,                       /* data1 line width, data2 line width, data3 line width */
	"", "", "",                    /* data1,data2,data3 input file */
	"",                            /* graph output file path */
	"",                            /* graph title */
	"Epoch-Time",                  /* x title */
	"Positon Errors /m",           /* y title */
	"%Y:%m:%d:%H:%M:%S",           /* time format */
	"%H:%M", "",                   /* x format, y format */
	"p1-c1", "rms", "Up"           /* data1 label, data2 label, data3 label */
};

/// Create graph title
int Create_GNU_Title(ppp_option_t& popt, char* title)
{
	char sta_str[6] = { 0 };
	char sys_str[10] = { 0 };
	char mode_str[15] = { 0 };
	char mathmode_str[15] = { 0 };
	char freq_flag = 0;

	if (popt.process_time == MODE_POST)
		memcpy(sta_str, popt.StaName, sizeof(char) * 4);
	else if (popt.process_time == MODE_REALTIME)
		memcpy(sta_str, popt.StaName, sizeof(char) * 6);

	if (popt.System & PRO_SYS_GPS)
		strcat(sys_str, "G");
	if (popt.System & PRO_SYS_GLO)
		strcat(sys_str, "R");
	if (popt.System & PRO_SYS_GAL)
		strcat(sys_str, "E");
	if (popt.System & PRO_SYS_BDS)
		strcat(sys_str, "C");

	switch (popt.mode_ppp)
	{
	case MODE_STATIC: { sprintf(mode_str, "%s", "Static"); break; }
	case MODE_KINEMATIC: { sprintf(mode_str, "%s", "Kinematic"); break; }
	default: { printf("Error:Mode setting error...\n"); return 0; }
	}

	switch (popt.math_model)
	{
	case MODE_UNCOMBINE: { sprintf(mathmode_str, "%s", "Uncombined"); break; }
	case MODE_IF_COMBINE: { sprintf(mathmode_str, "%s", "Combined"); break; }
	default: { printf("Error:Math model setting error...\n"); return 0; }
	}

	switch (popt.freqn)
	{
	case 1: { freq_flag = 'S'; break; }
	case 2: { freq_flag = 'D'; break; }
	case 3: { freq_flag = 'T'; break; }
	default: { printf("Error:Frequency number setting error...\n"); return 0; }
	}

	if (title)
		sprintf(title, "S5.0_%s_%s-%cF-%s-%s-PPP_%s", sta_str, sys_str, freq_flag, mode_str, mathmode_str, popt.IOfile.out_f.time_str);

	return 1;
}

/* Create gnuplot script -------------------------------------------------
* gnu_cfg_t    *gnucfg       I        gnuplot confifuration and settings
* Function:
*           1. Line graph
*           2. Histogram graph
* ----------------------------------------------------------------------- */
int create_gnuscrp(gnu_cfg_t* gnucfg, char* scrpname)
{
	if (strlen(scrpname) == 0 || scrpname == NULL) return 0;
	FILE* fp = fopen(scrpname, "w");

	if (fp)
	{
		fprintf(fp, "set output '%s' \n", gnucfg->outputf);

		/* font size and title */
		fprintf(fp, "set term gif font \",%d\" size %d,%d \n", gnucfg->font_size, gnucfg->length, gnucfg->width);
		fprintf(fp, "set title \"%s\" font \",%d\" \n", gnucfg->title, gnucfg->tfont_size);
		fprintf(fp, "set ylabel font \",%d\" \"%s\" \n", gnucfg->yfont_size, gnucfg->ytitle);
		fprintf(fp, "set xlabel font \",%d\" \"%s\" \n", gnucfg->xfont_size, gnucfg->xtitle);

		/* x rotate degree */
		fprintf(fp, "set xtics rotate by -%d \n", gnucfg->xrotate);

		/* y range */
		if (gnucfg->ymin != gnucfg->ymax)
			fprintf(fp, "set yrange [%.1f:%.1f] \n", gnucfg->ymin, gnucfg->ymax);

		/* x,y grid mode */
		if (gnucfg->ygrid)
			fprintf(fp, "set grid ytics \n");
		if (gnucfg->xgrid)
			fprintf(fp, "set grid xtics \n");

		if (gnucfg->graph_type == GNUGRAPHTYPE_HISTOGRAM)
		{
			fprintf(fp, "set style data histogram\n");
			fprintf(fp, "set style histogram clustered\n");
			fprintf(fp, "set style fill solid\n");
		}

		if (gnucfg->xdata_time)
		{
			fprintf(fp, "set xdata time \n");
			fprintf(fp, "set timefmt '%s' \n", gnucfg->time_format);
			fprintf(fp, "set format x '%s' \n", gnucfg->xformat);
		}
		if (gnucfg->ydata_time)
		{
			fprintf(fp, "set ydata time \n");
			fprintf(fp, "set timefmt '%s' \n", gnucfg->time_format);
			fprintf(fp, "set format y '%s' \n", gnucfg->yformat);
		}

		/* x and y tics */
		if (gnucfg->ytics != 0)
			fprintf(fp, "set ytics %.1f \n", gnucfg->ytics);
		if (gnucfg->xtics != 0)
			fprintf(fp, "set xtics %.1f \n", gnucfg->xtics);

		/* LINE graph */
		if (gnucfg->graph_type == GNUGRAPHTYPE_LINE)
		{
			/* 1st data */
			fprintf(fp, "plot '%s' using 1:%d t '%s' with ", gnucfg->data1_inputf, gnucfg->data1_column, gnucfg->data1_title);
			if (gnucfg->data1_lt == GNULINETYPE_POINT)
			{
				fprintf(fp, "points pt 1 ps 0.1 ");
			}
			if (gnucfg->data1_lt == GNULINETYPE_LINE)
			{
				fprintf(fp, "lines lw %d ", gnucfg->data1_lw);
			}
			if (gnucfg->data1_lt == GNULINETYPE_LINEPOINT)
			{
				fprintf(fp, "linespoints lw %d pt 4", gnucfg->data1_lw);
			}

			if (gnucfg->data_num == 1)
			{
				fprintf(fp, "\n");
				fclose(fp);
				return 1;
			}

			/* 2nd data */
			fprintf(fp, ", '%s' using 1:%d t '%s' with ", gnucfg->data2_inputf, gnucfg->data2_column, gnucfg->data2_title);
			if (gnucfg->data2_lt == GNULINETYPE_POINT)
			{
				fprintf(fp, "points pt 2 ps 0.1 ");
			}
			if (gnucfg->data2_lt == GNULINETYPE_LINE)
			{
				fprintf(fp, "lines lw %d ", gnucfg->data2_lw);
			}
			if (gnucfg->data2_lt == GNULINETYPE_LINEPOINT)
			{
				fprintf(fp, "linespoints lw %d pt 8", gnucfg->data2_lw);
			}

			if (gnucfg->data_num == 2)
			{
				fprintf(fp, "\n");
				fclose(fp);
				return 1;
			}

			/* 3rd data */
			fprintf(fp, ", '%s' using 1:%d t '%s' with ", gnucfg->data3_inputf, gnucfg->data3_column, gnucfg->data3_title);
			if (gnucfg->data3_lt == GNULINETYPE_POINT)
			{
				fprintf(fp, "points pt 3 ps 0.1 ");
			}
			if (gnucfg->data3_lt == GNULINETYPE_LINE)
			{
				fprintf(fp, "lines lw %d ", gnucfg->data3_lw);
			}
			if (gnucfg->data3_lt == GNULINETYPE_LINEPOINT)
			{
				fprintf(fp, "linespoints lw %d pt 14", gnucfg->data3_lw);
			}

			if (gnucfg->data_num == 3)
			{
				fprintf(fp, "\n");
				fclose(fp);
				return 1;
			}
		}
		/* Histogram graph */
		else if (gnucfg->graph_type == GNUGRAPHTYPE_HISTOGRAM)
		{
			fprintf(fp, "set key reverse\n");

			/* 1st data */
			fprintf(fp, "plot '%s' using %d:xticlabels(1) lc rgb\"navy\" t '%s' ", gnucfg->data1_inputf, gnucfg->data1_column, gnucfg->data1_title);

			if (gnucfg->data_num == 1)
			{
				fprintf(fp, "\n");
				fclose(fp);
				return 1;
			}

			/* 2nd data */
			fprintf(fp, ", '%s' using %d:xticlabels(1) lc rgb\"dark-turquoise\" t '%s' ", gnucfg->data2_inputf, gnucfg->data2_column, gnucfg->data2_title);

			if (gnucfg->data_num == 2)
			{
				fprintf(fp, "\n");
				fclose(fp);
				return 1;
			}

			/* 3rd data */
			fprintf(fp, ", '%s' using %d:xticlabels(1) lc rgb\"gold\" t '%s' ", gnucfg->data3_inputf, gnucfg->data3_column, gnucfg->data3_title);

			if (gnucfg->data_num == 3)
			{
				fprintf(fp, "\n");
				fclose(fp);
				return 1;
			}
		}
	}

	if (fp) fclose(fp);

	return 1;
}

int create_Ionosphere_gnuscrp(gnu_cfg_t* gnucfg, char* scrpname, unsigned int column, char sys, char* timestr, char* data_src)
{
	if (strlen(scrpname) == 0 || scrpname == NULL) return 0;
	FILE* fp = fopen(scrpname, "w");

	if (fp)
	{
		fprintf(fp, "set output '%s' \n", gnucfg->outputf);

		//% font size and title
		fprintf(fp, "set term gif font \",%d\" size %d,%d \n", gnucfg->font_size, gnucfg->length, gnucfg->width);
		fprintf(fp, "set title \"%s\" font \",%d\" \n", gnucfg->title, gnucfg->tfont_size);
		fprintf(fp, "set ylabel font \",%d\" \"%s\" \n", gnucfg->yfont_size, gnucfg->ytitle);
		fprintf(fp, "set xlabel font \",%d\" \"%s\" \n", gnucfg->xfont_size, gnucfg->xtitle);

		//% Put key outside
		fprintf(fp, "set key outside \n");

		//% x rotate degree
		fprintf(fp, "set xtics rotate by -%d \n", gnucfg->xrotate);

		//% y range
		if (gnucfg->xmin != gnucfg->xmax)
			fprintf(fp, "set xrange [%.1f:%.1f] \n", gnucfg->xmin, gnucfg->xmax);

		//% y range
		if (gnucfg->ymin != gnucfg->ymax)
			fprintf(fp, "set yrange [%.1f:%.1f] \n", gnucfg->ymin, gnucfg->ymax);

		//% x,y grid mode
		if (gnucfg->ygrid)
			fprintf(fp, "set grid ytics \n");
		if (gnucfg->xgrid)
			fprintf(fp, "set grid xtics \n");


		if (gnucfg->xdata_time)
		{
			fprintf(fp, "set xdata time \n");
			fprintf(fp, "set timefmt '%s' \n", gnucfg->time_format);
			fprintf(fp, "set format x '%s' \n", gnucfg->xformat);
		}
		if (gnucfg->ydata_time)
		{
			fprintf(fp, "set ydata time \n");
			fprintf(fp, "set timefmt '%s' \n", gnucfg->time_format);
			fprintf(fp, "set format y '%s' \n", gnucfg->yformat);
		}

		/* x and y tics */
		if (gnucfg->ytics != 0)
			fprintf(fp, "set ytics %.1f \n", gnucfg->ytics);
		if (gnucfg->xtics != 0)
			fprintf(fp, "set xtics %.1f \n", gnucfg->xtics);

		unsigned int HeadPlot = 0;
		for (int i = 0; i < gnucfg->data_num; ++i)
		{
			unsigned int prn = i + 1;

			if (HeadPlot == 0)
			{
				fprintf(fp, "plot ");
				HeadPlot = 1;
			}
			else
				fprintf(fp, ",");

			fprintf(fp, "'%s%c%02d_%s.ion' using 1:%d t '%c%02d' with ",
				data_src, sys, prn, timestr, 2, sys, prn);

			fprintf(fp, "points ps 1 pt %d  lw %d lc rgb\"%s\" ", i, 1, GnuplotColor[i]);
		}
	}

	if (fp) fclose(fp);

	return 1;
}

/// Plot Trop results
int Plot_Trop_PPP_Results(ppp_option_t& popt)
{
	char trop_fig_out[200] = { 0 };
	sprintf(trop_fig_out, "%s%s/%s_%s_Trop.jpg", popt.IOfile.out_f.out_fp, popt.ProjName, popt.StaName, popt.IOfile.out_f.time_str);

	gnu_cfg_t gnuplot = gnucfg_ppp;

	strcpy(gnuplot.outputf, trop_fig_out);

	gnuplot.data_num = 1;
	gnuplot.data1_column = 14;
	gnuplot.ymin = 1.5;
	gnuplot.ymax = 3.5;
	gnuplot.ytics = 0;
	gnuplot.data1_lt = GNULINETYPE_LINEPOINT;
	sprintf(gnuplot.title, "%s Troposphere Zenith Delay", popt.StaName);
	memset(gnuplot.ytitle, 0, sizeof(char) * 100);
	sprintf(gnuplot.ytitle, "Delay Value / (m)");
	memset(gnuplot.data1_title, 0, sizeof(char) * 100);
	strcpy(gnuplot.data1_title, "TropDelay");
	strcpy(gnuplot.data1_inputf, popt.IOfile.out_f.out_fp_fpath);

	char scrp_name[200] = { 0 };
	sprintf(scrp_name, "%s", "Temp/trop_gnuplot.scrp");
	create_gnuscrp(&gnuplot, scrp_name);

	char gnu_commond[100] = { 0 };
	sprintf(gnu_commond, "gnuplot %s", scrp_name);
	system(gnu_commond);

	return 1;
}

int Plot_SUPREME_ION_Results(char sys, char* title, unsigned int column, char* timestr, char* data_src, char* outputfile)
{
	char gnu_commond[100] = { 0 };    // system commond line for plot
	char out_graph[100] = { 0 };      // output graph file
	char scrp_name[100] = { 0 };      // gnuplot scrip file

	unsigned int SatN = 0;

	//% Get satellite number
	switch (sys)
	{
	case 'G': { SatN = GPS_SATNUM; break; }
	case 'R': { SatN = GLO_SATNUM; break; }
	case 'E': { SatN = GAL_SATNUM; break; }
	case 'C': { SatN = BDS_SATNUM; break; }
	default:return 0;
	}

	sprintf(out_graph, outputfile);

	if (!FileExist("Temp/")) CreatePath("Temp/");
	sprintf(scrp_name, "Temp/ion_gnuplot.scrp");

	gnu_cfg_t gnuplot = gnucfg_ion;
	strcpy(gnuplot.outputf, out_graph);

	strcpy(gnuplot.title, title);

	gnuplot.data_num = SatN;

	create_Ionosphere_gnuscrp(&gnuplot, scrp_name, column, sys, timestr, data_src);

	sprintf(gnu_commond, "gnuplot %s", scrp_name);
	system(gnu_commond);

	return 1;
}

/// Plot ppp results of SUPREME ppp
int Plot_SUPREME_PPP_Results(ppp_option_t& popt)
{
	char gnu_commond[100] = { 0 };            // system commond line for plot
	char scrp_name[100] = { 0 };              // gnuplot scrip file
	char out_graph[PATH_LENGTH] = { 0 };      // output graph file

	sprintf(out_graph, "%s%s/%s.jpg", popt.IOfile.out_f.out_fp, popt.ProjName, popt.IOfile.out_f.OutName);

	if (!FileExist("Temp/")) CreatePath("Temp/");
	sprintf(scrp_name, "Temp/%s_gnuplot.scrp", popt.StaName);

	gnu_cfg_t gnuplot = gnucfg_ppp;
	strcpy(gnuplot.data1_inputf, popt.IOfile.out_f.out_fp_fpath);
	strcpy(gnuplot.data2_inputf, popt.IOfile.out_f.out_fp_fpath);
	strcpy(gnuplot.data3_inputf, popt.IOfile.out_f.out_fp_fpath);
	strcpy(gnuplot.outputf, out_graph);

	if (popt.freqn == 1)
	{
		if (popt.mode_ppp == MODE_KINEMATIC)
		{
			gnuplot.ymin = -2;
			gnuplot.ymax = 2;
			gnuplot.ytics = 0.4;

			//gnuplot.ymin = -10;
			//gnuplot.ymax = 10;
			//gnuplot.ytics = 2;
		}

		else if (popt.process_mode == MODE_STATIC)
		{
			gnuplot.ymin = -1;
			gnuplot.ymax = 1;
			gnuplot.ytics = 0.1;
			if (popt.math_model == MODE_IF_COMBINE)
			{
				/*gnuplot.ymin = -0.5;
				gnuplot.ymax = 0.5;
				gnuplot.ytics = 0.1;*/

				gnuplot.ymin = -1;
				gnuplot.ymax = 1;
				gnuplot.ytics = 0.5;
			}

		}

	}

	/*if ((popt.System == PRO_SYS_GAL) || (popt.System == PRO_SYS_BDS) || (popt.System == PRO_SYS_GPS && popt.freqn == 2))
	{
		gnuplot.ymin = -1;
		gnuplot.ymax = 1;
		gnuplot.ytics = 0.2;

		gnuplot.ymin = -0.5;
		gnuplot.ymax = 0.5;
		gnuplot.ytics = 0.1;
	}*/

	if (popt.process_mode == MODE_SPP)
	{
		gnuplot.ymin = -15;
		gnuplot.ymax = 15;
		gnuplot.ytics = 3;
	}

	//Create_GNU_Title(popt, gnuplot.title);

	create_gnuscrp(&gnuplot, scrp_name);

	sprintf(gnu_commond, "gnuplot %s", scrp_name);
	system(gnu_commond);

	return 1;
}

/// Plot Post PPP Results
int Plot_SUPREME_Post_Results(ppp_option_t* popt)
{
	if (popt->Gnuplot)
	{
		if (popt->estpara.xyz)
		{
			if (Norm(popt->StaPos_CRD, 3) != 0)
				Plot_SUPREME_PPP_Results(*popt);
		}

		/*if (popt->estpara.trop)
		{
			Plot_Trop_PPP_Results(*popt);
		}*/

		if (!popt->trace) return 1;

#if 0
		if (popt->math_model == MODE_UNCOMBINE)
		{
			if (popt->System & PRO_SYS_GPS)
			{
				char ion_out[200] = { 0 };
				sprintf(ion_out, "%s%s_%s_%c_Ion.jpg",
					popt->IOfile.out_f.out_ion, popt->StaName, popt->IOfile.out_f.time_str, 'G');
				Plot_SUPREME_ION_Results('G', "", 2, popt->IOfile.out_f.time_str, popt->IOfile.out_f.out_ion, ion_out);
			}
			if (popt->System & PRO_SYS_GLO)
			{
				char ion_out[200] = { 0 };
				sprintf(ion_out, "%s%s_%s_%c_Ion.jpg",
					popt->IOfile.out_f.out_ion, popt->StaName, popt->IOfile.out_f.time_str, 'R');
				Plot_SUPREME_ION_Results('R', "", 2, popt->IOfile.out_f.time_str, popt->IOfile.out_f.out_ion, ion_out);
			}
			if (popt->System & PRO_SYS_GAL)
			{
				char ion_out[200] = { 0 };
				sprintf(ion_out, "%s%s_%s_%c_Ion.jpg",
					popt->IOfile.out_f.out_ion, popt->StaName, popt->IOfile.out_f.time_str, 'E');
				Plot_SUPREME_ION_Results('E', "", 2, popt->IOfile.out_f.time_str, popt->IOfile.out_f.out_ion, ion_out);
			}
			if (popt->System & PRO_SYS_BDS)
			{
				char ion_out[200] = { 0 };
				sprintf(ion_out, "%s%s_%s_%c_Ion.jpg",
					popt->IOfile.out_f.out_ion, popt->StaName, popt->IOfile.out_f.time_str, 'C');
				Plot_SUPREME_ION_Results('C', "", 2, popt->IOfile.out_f.time_str, popt->IOfile.out_f.out_ion, ion_out);
			}
		}
#endif
	}

	return 1;
}
