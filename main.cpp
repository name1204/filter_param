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
void test_FilterParam_csw();
void test_FilterParam_freq_res_se();

int main(void)
{
	printf("example run\n");

	test_FilterParam_freq_res_se();

	return 0;
}

/* 周波数帯域の構造体
 *   生成と表示のテスト
 *
 */
void test_BandParam_new()
{
	auto bp = BandParam(BandType::Pass, 0.0, 0.2);
	printf("%s\n", bp.sprint().c_str());
}

/* フィルタ構造体
 *   フィルタタイプから周波数帯域を生成するテスト
 *   大抵、複数の周波数帯域からフィルタが成るため、
 *   vectorで周波数帯域を返却する
 */
void test_Band_generator()
{
	auto bands = FilterParam::gen_bands(FilterType::LPF, 0.2, 0.3);
	for(auto bp :bands)
	{
		printf("%s\n", bp.sprint().c_str());
	}
}

/* フィルタ構造体
 *   文字列(string)から、フィルタタイプと
 *   周波数帯域端を分離するテスト
 */
void test_analyze_edges()
{
	string type("LPF(0.2 : 0.3)");
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

/* フィルタ構造体
 *   CSVファイルから所望特性を読み取るテスト
 *   書式は以下の通り
 *
 *     No,Numerator,Denominator,State,GroupDelay,NsplitApprox,NspritTransition
 *
 *   vectorのインデックス = Noなので
 *   "No"は読みこまない
 */
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

void test_FilterParam_csw()
{
	auto band = BandParam(BandType::Pass, 0.0, 0.2);
	auto csw = FilterParam::gen_csw(band, 100);
	auto csw2 = FilterParam::gen_csw2(band, 100);

	for(auto z :csw)
	{
		printf("%6f %6f\n", z.real(), z.imag());
	}
	printf("\n\n\n\n");
	for(auto z :csw2)
	{
		printf("%6f %6f\n", z.real(), z.imag());
	}

}

void test_FilterParam_freq_res_se()
{
	vector<double> coef
	{
		0.018656458,
		1.969338828,
		1.120102082,
		0.388717952,
		0.996398946,
		1.048137529,
		1.037079725,
		-4.535575709,
		6.381429398,
		-0.139429968,
		0.763426685
	};
	auto bands = FilterParam::gen_bands(FilterType::LPF, 0.2, 0.3);
	FilterParam fparam(8, 2, bands, 200, 50, 5.0);

	auto freq_res = fparam.freq_res(coef);

	for(auto band_res :freq_res)
	{
		for(auto res :band_res)
		{
			printf("%f\n", abs(res));
		}
	}
}
