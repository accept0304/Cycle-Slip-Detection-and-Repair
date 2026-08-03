#ifndef SUPREME_PPP_H_HH
#define SUPREME_PPP_H_HH
#include <deque>
#include <vector>
#include <fstream>
#include "SUPREME_Options.h"


int SUPREME_PPP_LS_Post(ppp_config &cfg);

int SUPREME_PPP_LS_RealTime(ppp_config &cfg);

int SUPREME_PPP_Processing(char *cfg_path);

#endif