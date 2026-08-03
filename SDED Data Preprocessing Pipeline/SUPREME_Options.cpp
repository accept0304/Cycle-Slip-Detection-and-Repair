#include <string.h>

#include "SUPREME_Options.h"
#include "SUPREME_CommonFunction.h"

#define CONFIG_LINE      400        // Config line length
#define MAX_PATH_LEN     400        // Max line length
#define OPTION_FLAG_LEN  13         // Option flag length

ppp_config::ppp_config()
{
	StaNum = 0;
	popt = NULL;
}
ppp_config::~ppp_config()
{
	ppp_option_t *opt = NULL, *tail = NULL;
	opt = popt; tail = opt->next;
	while (opt != NULL)
	{
		if (opt) free(opt);
		opt = tail;
		if (opt != NULL)
			tail = opt->next;
	}
}

/// Add option to config
int ppp_config::add_opt(ppp_option_t *p, ppp_option_t **tail)
{
	if (p == NULL) return 0;

	if (popt == NULL)
	{
		popt = p;
		*tail = p;
	}
	else
	{
		(*tail)->next = p;
		(*tail) = p;
	}

	StaNum++;

	return 1;
}

/// Read multi-gnss ppp configuration
int ppp_config::ReadConfig(const char *cfg_path)
{
	int i = 0, j = 0, k = 0, pl = 0;
	int lenl = CONFIG_LINE + 1;                 // config string length
	char strl[CONFIG_LINE + 1] = { 0 };         // config string line
	char str_flg[OPTION_FLAG_LEN + 1] = { 0 };  // config option flag

	char project_name[MAX_PATH_LEN] = { 0 };    // project name

	/* Nav, Obs, Eph, Clk file path */
	char navpath[MAX_PATH_LEN] = { 0 }, obspath[MAX_PATH_LEN] = { 0 };
	char pre_eph_path[MAX_PATH_LEN] = { 0 }, pre_clk_path[MAX_PATH_LEN] = { 0 };

	/* General data path */
	char atx_path[MAX_PATH_LEN] = { 0 }, pco_path[MAX_PATH_LEN] = { 0 }, pcv_path[MAX_PATH_LEN] = { 0 };
	char dcb_p1c1_path[MAX_PATH_LEN] = { 0 }, dcb_p1p2_path[MAX_PATH_LEN] = { 0 }, dcb_p2c2_path[MAX_PATH_LEN] = { 0 };
	char erp_path[MAX_PATH_LEN] = { 0 }, blq_path[MAX_PATH_LEN] = { 0 };
	char trop_path[MAX_PATH_LEN] = { 0 }, ion_path[MAX_PATH_LEN] = { 0 };
	char crd_path[MAX_PATH_LEN] = { 0 };
	char ntrip_path[MAX_PATH_LEN] = { 0 };

	/* output file path */
	char out_path[MAX_PATH_LEN] = { 0 };            // result output directory
	char res_path[MAX_PATH_LEN] = { 0 };            // residual file directory
	char vbl_path[MAX_PATH_LEN] = { 0 };            // process variable directory
	char log_path[MAX_PATH_LEN] = { 0 };            // log file directory
	char trace_path[MAX_PATH_LEN] = { 0 };          // trace file directory
	char sat_path[MAX_PATH_LEN] = { 0 };            // satellite information file directory
	char cycleslip_path[MAX_PATH_LEN] = { 0 };      // cycle slip directory

	/* Option 1 */
	char psystem[8] = { 0 };                        // process system  : GREC
	char pmode[6] = { 0 };                          // process mode    : spp/ppp
	char tmode[6] = { 0 };                          // post/realtime
	char gnuplot[4] = { 0 };                        // gnuplot option  : 1/0: on/off
	char trace[4] = { 0 };                          // trace option  : 1/0: on/off

	/* Option 2 */
	char trop_model[4] = { 0 };                     // troposphere model
	char iono_model[4] = { 0 };                     // ionosphere model
	char rclk_model[4] = { 0 };                     // receiver clock model
	char process_mode[4] = { 0 };                   // 0:static  1:kinematic
	char mathmode[4] = { 0 };                       // combine or not
	char gim_p1p2[4] = { 0 };
	char is_Dcb[4] = { 0 };
	char is_Isb[4] = { 0 };
	char isbmode[4] = { 0 };
	char estxyz[4] = { 0 };                         // estimate xyz
	char estrclk[4] = { 0 };                        // estimate receiver clock
	char esttrop[4] = { 0 };                        // estimate troposphere delay
	char estiono[4] = { 0 };                        // estimate ionosphere delay
	char estamb[4] = { 0 };                         // estimate ambiguity
	char estcbias[4] = { 0 };                       // estimate code bias
	char estpbias[4] = { 0 };                       // estimate phase bias
	char est_mode[4] = { 0 };                       // parameter estimation mode
	char filter_mode[4] = { 0 };                    // filting mode

	/* Option 3 */
	char interval[4] = { 0 };                       // sampling interval
	char freqn[4] = { 0 };                          // frequency number
	char minelev[4] = { 0 };                        // minimum satellite elevation
	char timedelay[4] = { 0 };                      // ppp time delay
	char sigemap[10] = { 0 };
	char sigemal[10] = { 0 };
	char trop_noise[10] = { 0 }, trop_gn_noise[10] = { 0 }, trop_ge_noise[10] = { 0 };
	char ion_noise[10] = { 0 };
	char rclk_noise[10] = { 0 };
	char ant_delta_n[10] = { 0 };                   // receiver antenna offset: N
	char ant_delta_e[10] = { 0 };                   // receiver antenna offset: E
	char ant_delta_u[10] = { 0 };                   // receiver antenna offset: U
	char antenna_type[21] = { 0 };                  // antenna type
	char receiver_type[21] = { 0 };                 // receiver type

	/* Option 4 */
	char ssr_mount_point[50] = { 0 };               // ssr mount point
	char nav_mount_point[50] = { 0 };               // nav mount point
	char ssr_id[5] = { 0 }, nav_id[5] = { 0 };      // ssr and nav id (mainly for ntrip in linux)
	char station_name[50] = { 0 };                  // station name
	char use_orbclk[4] = { 0 };                     // use orbit and clock product
	char ssr_ref_point[4] = { 0 };                  // orbit reference point

	/* Option 5 */
	char threadnum[4] = { 0 };
	char batch[4] = { 0 };
	char center[10] = { 0 };
	char batch_continue[4] = { 0 };
	char batch_obspath[MAX_PATH_LEN] = { 0 };
	char batch_prepath[MAX_PATH_LEN] = { 0 };
	char batch_genpath[MAX_PATH_LEN] = { 0 };
	char stationlist[MAX_PATH_LEN] = { 0 };
	char startt[20] = { 0 }, endt[20] = { 0 };

	char sta_id[4] = { 0 }, comapc[4] = { 0 };
	char nmea[2] = "", lat[15] = "", lon[15] = "";

	ppp_option_t *opt = NULL, *tail = NULL;

	FILE *fp = fopen(cfg_path, "r");
	if (fp == NULL)
	{
		printf("Error:Open config file %s failed.\n", cfg_path);
#ifdef _WIN32
		system("pause");
#endif
		return 0;
	}

	while (!feof(fp))
	{
		if (GetFileLine(fp, strl, lenl) == 0) continue;
		if (strl[0] == '#' || strl[0] == ' ' || strl[0] == '!' || strl[0] == '%') continue;

		if (strl[0] == '*')
		{
			GetSubStr(strl, lenl, str_flg, OPTION_FLAG_LEN + 1, 0, OPTION_FLAG_LEN);

			j = 0;

			/* Project Name */
			if (strcmp(str_flg, "*PROJ_NAME  :") == 0)
			{
				memset((void*)(project_name), 0, sizeof(project_name));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j < MAX_PATH_LEN - 1 && strl[i] != ' ')  project_name[j++] = strl[i];
				}
				continue;
			}

			/* navigation file path */
			if (strcmp(str_flg, "*NAV_PATH   :") == 0)
			{
				memset((void*)(navpath), 0, sizeof(navpath));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j < MAX_PATH_LEN - 1 && strl[i] != ' ')  navpath[j++] = strl[i];
				}
				continue;
			}
			/* precise ephemeris path */
			if (strcmp(str_flg, "*PRE_PATH   :") == 0)
			{
				memset((void*)(pre_eph_path), 0, sizeof(pre_eph_path));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') pre_eph_path[j++] = strl[i];
				}
				continue;
			}
			/* precise clock file path */
			if (strcmp(str_flg, "*CLK_PATH   :") == 0)
			{
				memset((void*)(pre_clk_path), 0, sizeof(pre_clk_path));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') pre_clk_path[j++] = strl[i];
				}
				continue;
			}
			/* observation file path */
			if (strcmp(str_flg, "*OBS_PATH   :") == 0)
			{
				memset((void*)(obspath), 0, sizeof(obspath));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') obspath[j++] = strl[i];
				}
				continue;
			}

			/* ionosphere tec file path */
			if (strcmp(str_flg, "*ION_PATH   :") == 0)
			{
				memset((void*)(ion_path), 0, sizeof(ion_path));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') ion_path[j++] = strl[i];
				}
				continue;
			}
			/* earth rotation parameter file path */
			if (strcmp(str_flg, "*ERP_PATH   :") == 0)
			{
				memset((void*)(erp_path), 0, sizeof(erp_path));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') erp_path[j++] = strl[i];
				}
				continue;
			}
			if (strcmp(str_flg, "*BLQ_PATH   :") == 0)
			{
				memset((void*)(blq_path), 0, sizeof(blq_path));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') blq_path[j++] = strl[i];
				}
				continue;
			}

#if 0
			/* JPL planet ephemeris path */
			if (strcmp(str_flg, "*PLA_PATH   :") == 0)
			{
				memset((void*)(pla_path), 0, sizeof(pla_path));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') pla_path[j++] = strl[i];
				}
				continue;
			}
			/* iau file path */
			if (strcmp(str_flg, "*IAU_PATH   :") == 0)
			{
				memset((void*)(iau_path), 0, sizeof(iau_path));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') iau_path[j++] = strl[i];
				}
				continue;
			}
			if (strcmp(str_flg, "*SEC_PATH   :") == 0)
			{
				memset((void*)(sec_path), 0, sizeof(sec_path));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') sec_path[j++] = strl[i];
				}
				continue;
			}
#endif

			/* ATX file path */
			if (strcmp(str_flg, "*ATX_PATH   :") == 0)
			{
				memset((void*)(atx_path), 0, sizeof(atx_path));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') atx_path[j++] = strl[i];
				}
				continue;
			}
			/* pco file path */
			if (strcmp(str_flg, "*PCO_PATH   :") == 0)
			{
				memset((void*)(pco_path), 0, sizeof(pco_path));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') pco_path[j++] = strl[i];
				}
				continue;
			}
			/* pcv file path */
			if (strcmp(str_flg, "*PCV_PATH   :") == 0)
			{
				memset((void*)(pcv_path), 0, sizeof(pcv_path));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') pcv_path[j++] = strl[i];
				}
				continue;
			}
			/* differential code bias file path */
			if (strcmp(str_flg, "*P1C1_PATH  :") == 0)
			{
				memset((void*)(dcb_p1c1_path), 0, sizeof(dcb_p1c1_path));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') dcb_p1c1_path[j++] = strl[i];
				}
				continue;
			}
			if (strcmp(str_flg, "*P1P2_PATH  :") == 0)
			{
				memset((void*)(dcb_p1p2_path), 0, sizeof(dcb_p1p2_path));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') dcb_p1p2_path[j++] = strl[i];
				}
				continue;
			}
			if (strcmp(str_flg, "*P2C2_PATH  :") == 0)
			{
				memset((void*)(dcb_p2c2_path), 0, sizeof(dcb_p2c2_path));
				for (i = OPTION_FLAG_LEN; i < lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j < MAX_PATH_LEN - 1 && strl[i] != ' ') dcb_p2c2_path[j++] = strl[i];
				}
				continue;
			}
			/* station coordinate file path */
			if (strcmp(str_flg, "*CRD_PATH   :") == 0)
			{
				memset((void*)(crd_path), 0, sizeof(crd_path));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') crd_path[j++] = strl[i];
				}
				continue;
			}

			/* output file path */
			if (strcmp(str_flg, "*OUT_PATH   :") == 0)
			{
				memset((void*)(out_path), 0, sizeof(out_path));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') out_path[j++] = strl[i];
				}
				continue;
			}
			/* output file path */
			if (strcmp(str_flg, "*VBL_PATH   :") == 0)
			{
				memset((void*)(vbl_path), 0, sizeof(vbl_path));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') vbl_path[j++] = strl[i];
				}
				continue;
			}
			/* log file path */
			if (strcmp(str_flg, "*LOG_PATH   :") == 0)
			{
				memset((void*)(log_path), 0, sizeof(log_path));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') log_path[j++] = strl[i];
				}
				continue;
			}
			/* detail information file path */
			if (strcmp(str_flg, "*TRA_PATH   :") == 0)
			{
				memset((void*)(trace_path), 0, sizeof(trace_path));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') trace_path[j++] = strl[i];
				}
				continue;
			}
			/* satellite cycle slip file path */
			if (strcmp(str_flg, "*SAT_PATH   :") == 0)
			{
				memset((void*)(sat_path), 0, sizeof(sat_path));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') sat_path[j++] = strl[i];
				}
				continue;
			}
			/* residual file path */
			if (strcmp(str_flg, "*RES_PATH   :") == 0)
			{
				memset((void*)(res_path), 0, sizeof(res_path));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') res_path[j++] = strl[i];
				}
				continue;
			}
			/* gf and mw statistic file path */
			if (strcmp(str_flg, "*CYS_PATH   :") == 0)
			{
				memset((void*)(cycleslip_path), 0, sizeof(cycleslip_path));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') cycleslip_path[j++] = strl[i];
				}
				continue;
			}

			/* Process Option ---------------------------------------------------------- */
			/* proess system */
			if (strstr(str_flg, "*SYSTEM     :"))
			{
				memset((void*)(psystem), 0, sizeof(psystem));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') psystem[j++] = strl[i];
				}
				continue;
			}
			/* read process mode (spp ppp) */
			if (strstr(str_flg, "*PMODE      :"))
			{
				memset((void*)(pmode), 0, sizeof(pmode));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') pmode[j++] = strl[i];
				}
				continue;
			}
			if (strstr(str_flg, "*TMODE      :"))
			{
				memset((void*)(tmode), 0, sizeof(tmode));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') tmode[j++] = strl[i];
				}
				continue;
			}
			/* gnuplot option */
			if (strstr(str_flg, "*GNUPLOT    :"))
			{
				memset((void*)(gnuplot), 0, sizeof(gnuplot));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') gnuplot[j++] = strl[i];
				}
				continue;
			}
			/* trace option */
			if (strstr(str_flg, "*TRACE      :"))
			{
				memset((void*)(trace), 0, sizeof(trace));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') trace[j++] = strl[i];
				}
				continue;
			}

			/* ppp math mode combined/uncombined mode */
			if (strstr(str_flg, "*MATHMODE   :"))
			{
				memset((void*)(mathmode), 0, sizeof(mathmode));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') mathmode[j++] = strl[i];
				}
				continue;
			}
			if (strstr(str_flg, "*GIM_P1P2   :"))
			{
				memset((void*)(gim_p1p2), 0, sizeof(gim_p1p2));
				for (i = OPTION_FLAG_LEN; i < lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j < MAX_PATH_LEN - 1 && strl[i] != ' ') gim_p1p2[j++] = strl[i];
				}
				continue;
			}
			if (strstr(str_flg, "*DCB        :"))
			{
				memset((void*)(is_Dcb), 0, sizeof(is_Dcb));
				for (i = OPTION_FLAG_LEN; i < lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j < MAX_PATH_LEN - 1 && strl[i] != ' ') is_Dcb[j++] = strl[i];
				}
				continue;
			}
			if (strstr(str_flg, "*ISB        :"))
			{
				memset((void*)(is_Isb), 0, sizeof(is_Isb));
				for (i = OPTION_FLAG_LEN; i < lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j < MAX_PATH_LEN - 1 && strl[i] != ' ') is_Isb[j++] = strl[i];
				}
				continue;
			}
			if (strstr(str_flg, "*ISBMODE    :"))
			{
				memset((void*)(isbmode), 0, sizeof(isbmode));
				for (i = OPTION_FLAG_LEN; i < lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j < MAX_PATH_LEN - 1 && strl[i] != ' ') isbmode[j++] = strl[i];
				}
				continue;
			}
			/* filting mode */
			if (strstr(str_flg, "*FILTERMODE :"))
			{
				memset((void*)(filter_mode), 0, sizeof(filter_mode));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') filter_mode[j++] = strl[i];
				}
				continue;
			}
			/* static/kinematic mode */
			if (strstr(str_flg, "*MODE       :"))
			{
				memset((void*)(process_mode), 0, sizeof(process_mode));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') process_mode[j++] = strl[i];
				}
				continue;
			}
			/* troposphere model */
			if (strstr(str_flg, "*TROP_MODEL :"))
			{
				memset((void*)(trop_model), 0, sizeof(trop_model));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') trop_model[j++] = strl[i];
				}
				continue;
			}
			if (strstr(str_flg, "*IONO_MODEL :"))
			{
				memset((void*)(iono_model), 0, sizeof(iono_model));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') iono_model[j++] = strl[i];
				}
				continue;
			}
			if (strstr(str_flg, "*RCLK_MODEL :"))
			{
				memset((void*)(rclk_model), 0, sizeof(rclk_model));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') rclk_model[j++] = strl[i];
				}
				continue;
			}
			/* Parameter Estimation Mode */
			if (strstr(str_flg, "*PARAESTMODE:"))
			{
				memset((void*)(est_mode), 0, sizeof(est_mode));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') est_mode[j++] = strl[i];
				}
				continue;
			}
			if (strstr(str_flg, "*ESTXYZ     :"))
			{
				memset((void*)(estxyz), 0, sizeof(estxyz));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') estxyz[j++] = strl[i];
				}
				continue;
			}
			if (strstr(str_flg, "*ESTRCLK    :"))
			{
				memset((void*)(estrclk), 0, sizeof(estrclk));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') estrclk[j++] = strl[i];
				}
				continue;
			}
			if (strstr(str_flg, "*ESTTROP    :"))
			{
				memset((void*)(esttrop), 0, sizeof(esttrop));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') esttrop[j++] = strl[i];
				}
				continue;
			}
			if (strstr(str_flg, "*ESTIONO    :"))
			{
				memset((void*)(estiono), 0, sizeof(estiono));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') estiono[j++] = strl[i];
				}
				continue;
			}
			if (strstr(str_flg, "*ESTAMB     :"))
			{
				memset((void*)(estamb), 0, sizeof(estamb));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') estamb[j++] = strl[i];
				}
				continue;
			}
			if (strstr(str_flg, "*CBIAS      :"))
			{
				memset((void*)(estcbias), 0, sizeof(estcbias));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') estcbias[j++] = strl[i];
				}
				continue;
			}
			if (strstr(str_flg, "*PBIAS      :"))
			{
				memset((void*)(estpbias), 0, sizeof(estpbias));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') estpbias[j++] = strl[i];
				}
				continue;
			}

			/* ppp process interval */
			if (strstr(str_flg, "*INTERVAL   :"))
			{
				memset((void*)(interval), 0, sizeof(interval));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') interval[j++] = strl[i];
				}
				continue;
			}
			/* read frequence number */
			if (strstr(str_flg, "*FRNUM      :"))
			{
				memset((void*)(freqn), 0, sizeof(freqn));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') freqn[j++] = strl[i];
				}
				continue;
			}
			/* read elevation mask */
			if (strstr(str_flg, "*EMASK      :"))
			{
				memset((void*)(minelev), 0, sizeof(minelev));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') minelev[j++] = strl[i];
				}
				continue;
			}
			/* Sigema P */
			if (strstr(str_flg, "*SIGEMAP    :"))
			{
				memset((void*)(sigemap), 0, sizeof(sigemap));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') sigemap[j++] = strl[i];
				}
				continue;
			}
			/* Sigema L */
			if (strstr(str_flg, "*SIGEMAL    :"))
			{
				memset((void*)(sigemal), 0, sizeof(sigemal));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') sigemal[j++] = strl[i];
				}
				continue;
			}

			/* Trop Noise */
			if (strstr(str_flg, "*TROPNOISE  :"))
			{
				memset((void*)(trop_noise), 0, sizeof(trop_noise));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') trop_noise[j++] = strl[i];
				}
				continue;
			}
			/* Trop Gradient N Noise */
			if (strstr(str_flg, "*TROPGNNOISE:"))
			{
				memset((void*)(trop_gn_noise), 0, sizeof(trop_gn_noise));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') trop_gn_noise[j++] = strl[i];
				}
				continue;
			}
			/* Trop Gradient E Noise */
			if (strstr(str_flg, "*TROPGENOISE:"))
			{
				memset((void*)(trop_ge_noise), 0, sizeof(trop_ge_noise));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') trop_ge_noise[j++] = strl[i];
				}
				continue;
			}
			/* Iono Noise */
			if (strstr(str_flg, "*IONONOISE  :"))
			{
				memset((void*)(ion_noise), 0, sizeof(ion_noise));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') ion_noise[j++] = strl[i];
				}
				continue;
			}
			/* Receiver clock Noise */
			if (strstr(str_flg, "*RCLKNOISE  :"))
			{
				memset((void*)(rclk_noise), 0, sizeof(rclk_noise));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') rclk_noise[j++] = strl[i];
				}
				continue;
			}

			/* receiver type */
			if (strstr(str_flg, "*RECTYPE    :"))
			{
				memset((void*)(receiver_type), 0, sizeof(receiver_type));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<21 - 1) receiver_type[j++] = strl[i];
				}
				continue;
			}
			/* antenna type */
			if (strstr(str_flg, "*ANTTYPE    :"))
			{
				memset((void*)(antenna_type), 0, sizeof(antenna_type));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<21 - 1) antenna_type[j++] = strl[i];
				}
				continue;
			}
			/* receiver antenna delta N */
			if (strstr(str_flg, "*ANTDELTA_N :"))
			{
				memset((void*)(ant_delta_n), 0, sizeof(ant_delta_n));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') ant_delta_n[j++] = strl[i];
				}
				continue;
			}
			/* receiver antenna delta N */
			if (strstr(str_flg, "*ANTDELTA_E :"))
			{
				memset((void*)(ant_delta_e), 0, sizeof(ant_delta_e));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') ant_delta_e[j++] = strl[i];
				}
				continue;
			}
			/* receiver antenna delta N */
			if (strstr(str_flg, "*ANTDELTA_U :"))
			{
				memset((void*)(ant_delta_u), 0, sizeof(ant_delta_u));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') ant_delta_u[j++] = strl[i];
				}
				continue;
			}

			/* thread number */
			if (strstr(str_flg, "*THREAD_NUM :"))
			{
				memset((void*)(threadnum), 0, sizeof(threadnum));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') threadnum[j++] = strl[i];
				}
				continue;
			}
			/* analysis center */
			if (strstr(str_flg, "*CENTRE     :"))
			{
				memset((void*)(center), 0, sizeof(center));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') center[j++] = strl[i];
				}
				continue;
			}
			/* batch process mode */
			if (strstr(str_flg, "*BATCH_MODE :"))
			{
				memset((void*)(batch), 0, sizeof(batch));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') batch[j++] = strl[i];
				}
				continue;
			}
			/* batch continue process mode */
			if (strstr(str_flg, "*CONTINUE   :"))
			{
				memset((void*)(batch_continue), 0, sizeof(batch_continue));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') batch_continue[j++] = strl[i];
				}
				continue;
			}
			/* station list path */
			if (strcmp(str_flg, "*STATIONLIST:") == 0)
			{
				memset((void*)(stationlist), 0, sizeof(stationlist));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') stationlist[j++] = strl[i];
				}
				continue;
			}
			/* batch process obs path */
			if (strstr(str_flg, "*OBSPATH    :"))
			{
				memset((void*)(batch_obspath), 0, sizeof(batch_obspath));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') batch_obspath[j++] = strl[i];
				}
				continue;
			}
			/* batch process pre path */
			if (strstr(str_flg, "*PREPATH    :"))
			{
				memset((void*)(batch_prepath), 0, sizeof(batch_prepath));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') batch_prepath[j++] = strl[i];
				}
				continue;
			}
			/* batch process gen path */
			if (strstr(str_flg, "*GENBSPATH    :"))
			{
				memset((void*)(batch_genpath), 0, sizeof(batch_genpath));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') batch_genpath[j++] = strl[i];
				}
				continue;
			}
			/* batch process start time */
			if (strcmp(str_flg, "*START_TIME :") == 0)
			{
				memset((void*)(startt), 0, sizeof(startt));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') startt[j++] = strl[i];
				}
				continue;
			}
			/* batch process end time */
			if (strcmp(str_flg, "*END_TIME   :") == 0)
			{
				memset((void*)(endt), 0, sizeof(endt));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') endt[j++] = strl[i];
				}
				continue;
			}

			//
			/* IGG Ntrip path */
			if (strcmp(str_flg, "*NTRIP_PATH :") == 0)
			{
				memset((void*)(ntrip_path), 0, sizeof(ntrip_path));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') ntrip_path[j++] = strl[i];
				}
				continue;
			}

			/* ppp time delay */
			if (strstr(str_flg, "*TIME_DELAY :"))
			{
				memset((void*)(timedelay), 0, sizeof(timedelay));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') timedelay[j++] = strl[i];
				}
				continue;
			}

			/* ssr mount point */
			if (strstr(str_flg, "*SSRMOUNTP  :"))
			{
				memset((void*)(ssr_mount_point), 0, sizeof(ssr_mount_point));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') ssr_mount_point[j++] = strl[i];
				}
				continue;
			}

			/* ssr orbit reference point */
			if (strstr(str_flg, "*SSRREFPOINT:"))
			{
				memset((void*)(ssr_ref_point), 0, sizeof(ssr_ref_point));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') ssr_ref_point[j++] = strl[i];
				}
				continue;
			}

			/* navigation mount point */
			if (strstr(str_flg, "*NAVMOUNTP  :"))
			{
				memset((void*)(nav_mount_point), 0, sizeof(nav_mount_point));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') nav_mount_point[j++] = strl[i];
				}
				continue;
			}
			/* ssr shared memory id */
			if (strstr(str_flg, "*SSRSHMID   :"))
			{
				memset((void*)(ssr_id), 0, sizeof(ssr_id));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') ssr_id[j++] = strl[i];
				}
				continue;
			}
			/* navigation shared memory id */
			if (strstr(str_flg, "*NAVSHMID   :"))
			{
				memset((void*)(nav_id), 0, sizeof(nav_id));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') nav_id[j++] = strl[i];
				}
				continue;
			}
			/* use orbit/clock product */
			if (strstr(str_flg, "*USEORBCLK  :"))
			{
				memset((void*)(use_orbclk), 0, sizeof(use_orbclk));
				for (i = OPTION_FLAG_LEN; i<lenl; ++i)
				{
					if (strl[i] == '#') break;
					if (j<MAX_PATH_LEN - 1 && strl[i] != ' ') use_orbclk[j++] = strl[i];
				}
				continue;
			}
			
		}
		else
		{
			opt = (ppp_option_t*)(calloc(1, sizeof(ppp_option_t)));

			if (opt == NULL)
			{
				printf("Error:PPP Option Memory calloc failed.\n");
#ifdef _WIN32
				system("pause");
#endif
				continue;
			}
			else
			{
				strcpy(opt->ProjName, project_name);                   // Project name

				strcpy(opt->IOfile.rinex_f.navf, navpath);             // Navigation file path
				strcpy(opt->IOfile.rinex_f.obsf, obspath);             // Observation file path
				strcpy(opt->IOfile.rinex_f.prephf, pre_eph_path);      // Precise ephemeris file path
				strcpy(opt->IOfile.rinex_f.preclkf, pre_clk_path);     // Precise clock file path

				strcpy(opt->IOfile.gen_f.trp_fl, trop_path);           // ION file path
				strcpy(opt->IOfile.gen_f.ion_fl, ion_path);            // ION file path

				strcpy(opt->IOfile.gen_f.erp_fl, erp_path);            // Erp file path
				strcpy(opt->IOfile.gen_f.blq_fl, blq_path);            // Blq file path

#if 0
				strcpy(opt->IOfile.gen_f.pla_fl, pla_path);            // Planet ephemeris file path
				strcpy(opt->IOfile.gen_f.iau_fl, iau_path);
				strcpy(opt->IOfile.gen_f.sec_fl, sec_path);
#endif

				strcpy(opt->IOfile.gen_f.atx_fl, atx_path);            // ATX file path
				strcpy(opt->IOfile.gen_f.pco_fl, pco_path);            // PCO file path
				strcpy(opt->IOfile.gen_f.pcv_fl, pcv_path);            // PCO file path

				strcpy(opt->IOfile.gen_f.dcb_p1c1_fl, dcb_p1c1_path);  // DCB file path
				strcpy(opt->IOfile.gen_f.dcb_p1p2_fl, dcb_p1p2_path);  // DCB file path
				strcpy(opt->IOfile.gen_f.dcb_p2c2_fl, dcb_p2c2_path);
				strcpy(opt->IOfile.gen_f.crd_fl, crd_path);            // CRD file path

				strcpy(opt->IGGNtripPath, ntrip_path);                 // IggNtrip run path (for linux)

				strcpy(opt->IOfile.out_f.out_fp, out_path);
				strcpy(opt->IOfile.out_f.out_vbl, vbl_path);
				strcpy(opt->IOfile.out_f.out_log, log_path);
				strcpy(opt->IOfile.out_f.out_trace, trace_path);
				strcpy(opt->IOfile.out_f.out_sat, sat_path);
				strcpy(opt->IOfile.out_f.out_res, res_path);
				strcpy(opt->IOfile.out_f.out_cycle_slip, cycleslip_path);

				pl = strlen(opt->IOfile.rinex_f.navf);
				if (pl && (opt->IOfile.rinex_f.navf[pl - 1] == '\\' || opt->IOfile.rinex_f.navf[pl - 1] == '/'))
					opt->IOfile.rinex_f.navf[pl - 1] = '\0';
				pl = strlen(opt->IOfile.rinex_f.obsf);
				if (pl && (opt->IOfile.rinex_f.obsf[pl - 1] == '\\' || opt->IOfile.rinex_f.obsf[pl - 1] == '/'))
					opt->IOfile.rinex_f.obsf[pl - 1] = '\0';
				pl = strlen(opt->IOfile.rinex_f.prephf);
				if (pl && (opt->IOfile.rinex_f.prephf[pl - 1] == '\\' || opt->IOfile.rinex_f.prephf[pl - 1] == '/'))
					opt->IOfile.rinex_f.prephf[pl - 1] = '\0';
				pl = strlen(opt->IOfile.rinex_f.preclkf);
				if (pl && (opt->IOfile.rinex_f.preclkf[pl - 1] == '\\' || opt->IOfile.rinex_f.preclkf[pl - 1] == '/'))
					opt->IOfile.rinex_f.preclkf[pl - 1] = '\0';

				/* Process system ---------------- */
				if (strstr(psystem, "G"))
				{
					opt->System |= PRO_SYS_GPS;
				}
				if (strstr(psystem, "R"))
				{
					opt->System |= PRO_SYS_GLO;
				}
				if (strstr(psystem, "C"))
				{
					opt->System |= PRO_SYS_BDS;
				}
				if (strstr(psystem, "E"))
				{
					opt->System |= PRO_SYS_GAL;
				}
				/* ------------------------------- */

				opt->process_mode = str_to_i(pmode);
				opt->process_time = str_to_i(tmode);
				opt->Gnuplot = str_to_i(gnuplot);
				opt->trace = str_to_i(trace);

				opt->trop_model = str_to_i(trop_model);
				opt->ion_model = str_to_i(iono_model);
				opt->clk_model = str_to_i(rclk_model);
				opt->math_model = str_to_i(mathmode);
				opt->GIM_P1P2 = str_to_i(gim_p1p2);
				opt->IS_DCB = str_to_i(is_Dcb);
				opt->IS_ISB = str_to_i(is_Isb);
				opt->ISB_MODE = str_to_i(isbmode);


				opt->mode_ppp = str_to_i(process_mode);
				opt->mode_filter = str_to_i(filter_mode);              /* filting mode                      */

				opt->estpara.xyz = str_to_i(estxyz);                   /* estimate xyz parameter            */
				opt->estpara.rclk = str_to_i(estrclk);                 /* estimate receiver clock parameter */
				opt->estpara.trop = str_to_i(esttrop);                 /* estimate troposphere delay(wet)   */
				opt->estpara.ion = str_to_i(estiono);                  /* estimate ionosphere delay         */
				opt->estpara.amb = str_to_i(estamb);                   /* estimate ambiguity                */
				
				opt->interval = str_to_i(interval);                    /* process interval                  */
				opt->Freqn = str_to_i(freqn);                          /* frequence number                  */
				opt->freqn = 1;
				opt->ele_mask = str_to_f(minelev);                     /* minimum elevation angle           */
				opt->time_delay = str_to_i(timedelay);                 /* ppp time delay                    */
				opt->sigmaP = str_to_f(sigemap);                       /* Sigema P                          */
				opt->sigmaL = str_to_f(sigemal);                       /* Sigema L                          */
				opt->trop_noise = str2num(trop_noise,0,10);
				opt->trop_gradientNnoise = str2num(trop_gn_noise, 0, 10);
				opt->trop_gradientEnoise = str2num(trop_ge_noise, 0, 10);
				opt->ion_noise = str2num(ion_noise, 0, 10);
				opt->rclk_noise = str2num(rclk_noise, 0, 10);
				opt->antenna_offset[0] = str_to_f(ant_delta_n);        /* receiver antenna delta N          */
				opt->antenna_offset[1] = str_to_f(ant_delta_e);        /* receiver antenna delta E          */
				opt->antenna_offset[2] = str_to_f(ant_delta_u);        /* receiver antenna delta U          */
				strcpy(opt->antenna_type, antenna_type);               /* antenna type                      */
				strcpy(opt->receiver_type, receiver_type);             /* receiver type                     */

				strcpy(opt->ssr_mountpoint, ssr_mount_point);          /* ssr mountain point                */
				opt->ssr_com = str_to_i(ssr_ref_point);
				strcpy(opt->nav_mountpoint, nav_mount_point);          /* navigation mount point            */
				opt->ssr_id = str_to_i(ssr_id);                        /* ssr shared memory id              */
				opt->nav_id = str_to_i(nav_id);                        /* navigation shared memory id       */
				opt->use_orbclk = str_to_i(use_orbclk);                /* use orbit and clock product flag  */

				opt->para_est_mode = str_to_i(est_mode);               // parameter estimation mode
				opt->threadn = str_to_i(threadnum);
				batchmode = str_to_i(batch);
				opt->ContinueMode = str_to_i(batch_continue);
				strcpy(opt->CentreName, center);
				strcpy(opt->stationlist, stationlist);
				strcpy(opt->ObsPath, batch_obspath);
				strcpy(opt->PrePath, batch_prepath);
				strcpy(opt->GenPath, batch_genpath);

				int ep_s[3] = { 0 }, ep_e[3] = { 0 }, count = 0;
				char *r = strtok(startt, ":");
				while (r != NULL)
				{
					ep_s[count++] = str_to_i(r);
					r = strtok(NULL, ":");
				}
				count = 0;
				r = strtok(endt, ":");
				while (r != NULL)
				{
					ep_e[count++] = str_to_i(r);
					r = strtok(NULL, ":");
				}
				// batch process start time
				opt->starttime.m_Year = ep_s[0];
				opt->starttime.m_Month = ep_s[1];
				opt->starttime.m_Day = ep_s[2];
				// batch process end time
				opt->endtime.m_Year = ep_e[0];
				opt->endtime.m_Month = ep_e[1];
				opt->endtime.m_Day = ep_e[2];

				opt->starttime._date2gpst(); opt->endtime._date2gpst();

				if (strl[0] != '>')
					continue;
				else
				{
					memset((void*)(sta_id), 0, sizeof(sta_id));        // station shared memory id (for linux)
					memset((void*)(comapc), 0, sizeof(comapc));
					memset((void*)(nmea), 0, sizeof(nmea));
					memset((void*)(lat), 0, sizeof(lat));              // station latitude
					memset((void*)(lon), 0, sizeof(lon));              // station lontitude

					j = 0; k = 0;
					for (i = 0; i < lenl; ++i)
					{
						if (strl[i] == '>' || strl[i] == ' ') continue;

						if (strl[i] == ',') { j++; k = 0; continue; }

						if (j == 0 && k < 4){ sta_id[k++] = strl[i]; }
						else if (j == 1 && k <50){ opt->StaName[k++] = strl[i]; }
						else if (j == 2 && k < 4){ comapc[k++] = strl[i]; }
						else if (j == 3 && k < 1){ nmea[k++] = strl[i]; }
						else if (j == 4 && k <14){ lon[k++] = strl[i]; }
						else if (j == 5 && k <14){ lat[k++] = strl[i]; }
						else                    { free(opt);  continue; }
					}

					opt->nmea = str_to_i(nmea);
					opt->Sta_id = str_to_i(sta_id);     // station id
					opt->lat = str_to_f(lat);           // station latitude
					opt->lon = str_to_f(lon);           // station lontitude

					add_opt(opt, &tail);
				}
			}
			//
		}
	}

	if (fp != NULL) fclose(fp);

	return 1;
}

/// Write multi-gnss ppp configuration
int ppp_config::WriteConfig(ppp_option_t *popt, char *cfg_path)
{
	FILE *fp = fopen(cfg_path, "w");

	if (fp)
	{
		fprintf(fp, "# SUPREME_Ver1.0 Configuration File ================================\n");
		fprintf(fp, "!>\n");
		fprintf(fp, "#===============================================================\n");
		fprintf(fp, "*PROJ_NAME  :%s\n", popt->ProjName);
		fprintf(fp, "*PRE_PATH   :%s\n", popt->IOfile.rinex_f.prephf);
		fprintf(fp, "*CLK_PATH   :%s\n", popt->IOfile.rinex_f.preclkf);
		fprintf(fp, "*OBS_PATH   :%s\n", popt->IOfile.rinex_f.obsf);
		fprintf(fp, "*NAV_PATH   :%s\n", popt->IOfile.rinex_f.navf);
		fprintf(fp, "#-----------------------------------------\n");
		fprintf(fp, "#General Input File Path\n");
		fprintf(fp, "*ION_PATH   :%s\n", popt->IOfile.gen_f.ion_fl);
		fprintf(fp, "*TRO_PATH   :%s\n", popt->IOfile.gen_f.trp_fl);
		fprintf(fp, "*ERP_PATH   :%s\n", popt->IOfile.gen_f.erp_fl);
		fprintf(fp, "*BLQ_PATH   :%s\n", popt->IOfile.gen_f.blq_fl);
		//fprintf(fp, "*PCO_PATH   :%s\n", popt->IOfile.gen_f.pco_fl);
		//fprintf(fp, "*PCV_PATH   :%s\n", popt->IOfile.gen_f.pcv_fl);
		fprintf(fp, "*ATX_PATH   :%s\n", popt->IOfile.gen_f.atx_fl);
		fprintf(fp, "*P1C1_PATH  :%s\n", popt->IOfile.gen_f.dcb_p1c1_fl);
		fprintf(fp, "*P1P2_PATH  :%s\n", popt->IOfile.gen_f.dcb_p1p2_fl);
		fprintf(fp, "*P2C2_PATH  :%s\n", popt->IOfile.gen_f.dcb_p2c2_fl);
		fprintf(fp, "*CRD_PATH   :%s\n", popt->IOfile.gen_f.crd_fl);
		fprintf(fp, "*NTRIP_PATH :%s\n", popt->IGGNtripPath);
		fprintf(fp, "#-----------------------------------------\n");
		fprintf(fp, "#Output Information Directory\n");
		fprintf(fp, "*OUT_PATH   :%s\n", popt->IOfile.out_f.out_fp);
		fprintf(fp, "*VBL_PATH   :%s\n", popt->IOfile.out_f.out_vbl);
		fprintf(fp, "*LOG_PATH   :%s\n", popt->IOfile.out_f.out_log);
		fprintf(fp, "*TRA_PATH   :%s\n", popt->IOfile.out_f.out_trace);
		fprintf(fp, "*SAT_PATH   :%s\n", popt->IOfile.out_f.out_sat);
		fprintf(fp, "*RES_PATH   :%s\n", popt->IOfile.out_f.out_res);
		fprintf(fp, "*CYS_PATH   :%s\n", popt->IOfile.out_f.out_cycle_slip);
		fprintf(fp, "#=========================================\n");
		fprintf(fp, "!<\n\n");

		fprintf(fp, "!>\n");
		fprintf(fp, "#=========================================\n");
		fprintf(fp, "# Option1 --------------------------------\n");

		char sys = 0;
		unsigned int mutli_sys = 1;

		if (popt->System == PRO_SYS_GPS)
		{
			fprintf(fp, "*SYSTEM     :%c", 'G');
			mutli_sys = 0;
		}
		if (popt->System == PRO_SYS_GLO)
		{
			fprintf(fp, "*SYSTEM     :%c", 'R');
			mutli_sys = 0;
		}
		if (popt->System == PRO_SYS_GAL)
		{
			fprintf(fp, "*SYSTEM     :%c", 'E');
			mutli_sys = 0;
		}
		if (popt->System == PRO_SYS_BDS)
		{
			fprintf(fp, "*SYSTEM     :%c", 'C');
			mutli_sys = 0;
		}

		if (mutli_sys)
		{
			if (popt->System&PRO_SYS_GPS)
				fprintf(fp, "*SYSTEM     :%c", 'G');
			if (popt->System&PRO_SYS_GLO)
				fprintf(fp, "%c", 'R');
			if (popt->System&PRO_SYS_GAL)
				fprintf(fp, "%c", 'E');
			if (popt->System&PRO_SYS_BDS)
				fprintf(fp, "%c", 'C');
		}

		fprintf(fp, "\n");
		//fprintf(fp, "*SYSTEM     :%c\n", sys);

		fprintf(fp, "*PMODE      :%d\n", popt->process_mode);
		fprintf(fp, "*TMODE      :%d\n", popt->process_time);
		fprintf(fp, "*GNUPLOT    :%d\n", popt->Gnuplot);
		fprintf(fp, "*TRACE      :%d\n", popt->trace);
		fprintf(fp, "#-----------------------------------------\n");
		fprintf(fp, "# Option2 --------------------------------\n");
		fprintf(fp, "*TROP_MODEL :%d\n", popt->trop_model);
		fprintf(fp, "*IONO_MODEL :%d\n", popt->ion_model);
		fprintf(fp, "*RCLK_MODEL :%d\n", popt->clk_model);
		fprintf(fp, "*MODE       :%d\n", popt->mode_ppp);
		fprintf(fp, "*MATHMODE   :%d\n", popt->math_model);
		fprintf(fp, "*FILTERMODE :%d\n", popt->mode_filter);
		fprintf(fp, "#Parameter estimation\n");
		fprintf(fp, "*PARAESTMODE:%d\n", 0);
		fprintf(fp, "*ESTXYZ     :%d\n", popt->estpara.xyz);
		fprintf(fp, "*ESTRCLK    :%d\n", popt->estpara.rclk);
		fprintf(fp, "*ESTTROP    :%d\n", popt->estpara.trop);
		fprintf(fp, "*ESTIONO    :%d\n", popt->estpara.ion);
		fprintf(fp, "*ESTAMB     :%d\n", popt->estpara.amb);
		fprintf(fp, "#-----------------------------------------\n");
		fprintf(fp, "# Option3 --------------------------------\n");
		fprintf(fp, "*INTERVAL   :%d\n", popt->interval);
		fprintf(fp, "*FRNUM      :%d\n", popt->freqn);
		fprintf(fp, "*EMASK      :%d\n", (int)popt->ele_mask);
		fprintf(fp, "*SIGEMAP    :%.4f\n", popt->sigmaP);
		fprintf(fp, "*SIGEMAL    :%.4f\n", popt->sigmaL);
		fprintf(fp, "*TROPNOISE  :%.2e\n", popt->trop_noise);
		fprintf(fp, "*TROPGNNOISE:%.2e\n", popt->trop_gradientNnoise);
		fprintf(fp, "*TROPGENOISE:%.2e\n", popt->trop_gradientEnoise);
		fprintf(fp, "*IONONOISE  :%.2e\n", popt->ion_noise);
		fprintf(fp, "*RCLKNOISE  :%.2e\n", popt->rclk_noise);
		fprintf(fp, "#Station Infomation ----------------------\n");
		fprintf(fp, "*RECVANT    :%s\n", popt->receiver_type);
		fprintf(fp, "*ANTTYPE    :%s\n", popt->antenna_type);
		fprintf(fp, "*ANTDELTA_N :%d\n", popt->antenna_offset[0]);
		fprintf(fp, "*ANTDELTA_E :%d\n", popt->antenna_offset[1]);
		fprintf(fp, "*ANTDELTA_U :%d\n", popt->antenna_offset[2]);
		fprintf(fp, "#-----------------------------------------\n");
		fprintf(fp, "# Option4 --------------------------------\n");
		fprintf(fp, "*TIME_DELAY :%d\n", popt->time_delay);
		fprintf(fp, "*SSRMOUNTP  :%s\n", popt->ssr_mountpoint);
		fprintf(fp, "*NAVMOUNTP  :%s\n", popt->nav_mountpoint);
		fprintf(fp, "*SSRSHMID   :%d\n", popt->ssr_id);
		fprintf(fp, "*NAVSHMID   :%d\n", popt->nav_id);
		fprintf(fp, "*USEORBCLK  :%d\n", popt->use_orbclk);
		fprintf(fp, "#-----------------------------------------\n");
		fprintf(fp, "# Option5 --------------------------------\n");
		fprintf(fp, "*THREAD_NUM :%d\n", 1);
		fprintf(fp, "*BATCH_MODE :%d\n", 0);
		fprintf(fp, "*CENTRE     :\n");
		fprintf(fp, "*OBSCENTRE  :\n");
		fprintf(fp, "*STATIONLIST:\n");
		fprintf(fp, "*START_TIME :\n");
		fprintf(fp, "*END_TIME   :\n");
		fprintf(fp, "#-----------------------------------------\n");
		fprintf(fp, "> 001, %s\n", popt->StaName);
		fprintf(fp, "#-----------------------------------------\n");
		fprintf(fp, "!<\n");
	}

	if (fp) fclose(fp);

	return 1;
}
