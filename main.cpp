/*
 * main.cpp
 *
 *  Created on: 2021/07/22
 *      Author: matsu
 */

#include <stdio.h>
#include <string>
#include "./lib/filter_param.hpp"

using namespace std;

void test_BandParam_new();
void test_Band_generator();
void test_analyze_edges();
void test_FilterParam_read_csv();

int main(void)
{
	printf("example run\n");


	return 0;
}

void test_BandParam_new()
{
	auto bp = BandParam(BandType::Pass, 0.0, 0.2);
	printf("%s\n", bp.sprint().c_str());
}

void test_Band_generator()
{
	auto bands = FilterParam::gen_bands(FilterType::LPF, 0.2, 0.3);
	for(auto bp :bands)
	{
		printf("%s\n", bp.sprint().c_str());
	}
}

void test_analyze_edges()
{
	string type("LPF(0.2, 0.3)");
	auto ftype = FilterParam::analyze_type(type);
	auto edges = FilterParam::analyze_edges(type);

	if(ftype == FilterType::LPF)
	{
		printf("LPF\n");
	}
	for(auto v :edges)
	{
		printf("%f ", v);
	}
	printf("\n");

}

void test_FilterParam_read_csv()
{
	string filename("./desire_filter.csv");
	auto params = FilterParam::read_csv(filename);
	for(auto param :params)
	{
		printf("order(zero/pole) : %d/%d\n", param.zero_order(), param.pole_order());
		printf("optimization order : %d\n", param.opt_order());
		printf("nsplit(approx-transition) : %d-%d\n", param.partition_approx(), param.partition_transition());
		printf("group delay : %f\n\n", param.gd());

		for(auto band :param.fbands())
		{
			printf("%s\n", band.sprint().c_str());
		}
		printf("---------------------------\n");
	}
}
