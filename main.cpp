/*
 * main.cpp
 *
 *  Created on: 2021/07/22
 *      Author: matsu
 */

#include "./lib/filter_param.hpp"

#include <stdio.h>
#include <string>
#include <chrono>

using namespace std;

void test_BandParam_new();
void test_Band_generator();
void test_analyze_edges();
void test_FilterParam_read_csv();
void test_FilterParam_csw();
void test_FilterParam_desire_res();
void test_FilterParam_freq_res_speed();
void test_FilterParam_freq_res_se();
void test_FilterParam_freq_res_so();
void test_FilterParam_freq_res_no();
void test_FilterParam_freq_res_mo();
void test_FilterParam_judge_stability_even();
void test_FilterParam_judge_stability_odd();
void test_FilterParam_evaluate_objective_function();
void test_FilterParam_init_coef();
void test_FilterParam_init_stable_coef();

int main(void)
{
	printf("example run\n");

	test_FilterParam_freq_res_so();

	return 0;
}

/* 周波数帯域の構造体
 *   生成と表示のテスト
 *
 */
void test_BandParam_new()
{
	auto bp = BandParam(BandType::Pass, 0.0, 0.2175);
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
	for (auto bp : bands)
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

	if (ftype == FilterType::LPF)
	{
		printf("LPF\n");
	}
	for (auto v : edges)
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
	for (auto param : params)
	{
		printf("order(zero/pole) : %d/%d\n", param.zero_order(), param.pole_order());
		printf("optimization order : %d\n", param.opt_order());
		printf("nsplit(approx-transition) : %d-%d\n", param.partition_approx(), param.partition_transition());
		printf("group delay : %f\n\n", param.gd());

		for (auto band : param.fbands())
		{
			printf("%s\n", band.sprint().c_str());
		}
		printf("---------------------------\n");
	}
}

/* フィルタ構造体
 *   複素正弦波の１次・２次の確認用
 */
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

/* # フィルタ構造体
 *   所望特性の周波数特性についてのテスト関数
 *   通過域でe^jωτ(τ= group delay)，
 *   阻止域で０，遷移域で要素なしの出力
 */
void test_FilterParam_desire_res()
{
	double gd = 5.0;

	auto pass_band = BandParam(BandType::Pass, 0.0, 0.2);
	auto desire_pass = FilterParam::gen_desire_res(pass_band, 100, gd);
	printf("-----------pass band-----------------\n");
	for(auto z :desire_pass)
	{
		printf("%6f %6f\n", z.real(), z.imag());
	}
	printf("----------------------------\n\n");

	auto stop_band = BandParam(BandType::Stop, 0.0, 0.2);
	auto desire_stop = FilterParam::gen_desire_res(stop_band, 100, gd);
	printf("-----------stop band-----------------\n");
	for(auto z :desire_stop)
	{
		printf("%6f %6f\n", z.real(), z.imag());
	}
	printf("----------------------------\n\n");

	auto trans_band = BandParam(BandType::Transition, 0.0, 0.2);
	auto desire_trans = FilterParam::gen_desire_res(trans_band, 100, gd);
	printf("-----------transition band-----------------\n");
	for(auto z :desire_trans)
	{
		printf("%6f %6f\n", z.real(), z.imag());
	}
	printf("----------------------------\n\n");

}

/* # フィルタ構造体
 *   周波数特性計算関数の実行速度を計算する
 *
 *   trialについての平均で判断する
 *   また各１回の実行時間は微小のため，
 *   repeat分繰り返して割ることで
 *   精度よく測定することを試みる
 */
void test_FilterParam_freq_res_speed()
{
	printf("thread will ce locked about 2 minutes.\n");

// 時間計測関連
	int trial = 100;
	int exp_ = 5;
	int repeat = pow(10, exp_);
	double ave = 0.0;
	double ave_all = 0.0;

// 計測雑利用
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
    vector<vector<complex<double>>> freq_res;

	for(int i = 0 ; i < trial ; ++i)
	{
		auto start1 = chrono::system_clock::now();      // 計測スタート時刻を保存

		for(int j = 0 ; j < repeat ; ++j)
		{
		    freq_res = fparam.freq_res(coef);
		}

		auto end1 = chrono::system_clock::now();       // 計測終了時刻を保存
		double msec1 = chrono::duration_cast<chrono::milliseconds>(end1 - start1).count();
		double once_time1 = msec1 / (double)repeat;
		ave += once_time1;
		ave_all += msec1;
	}

	printf("\n------------------------------------\n\n\n\n");
	printf("using functional : Average %15.15f[ns]\n", 1000*1000*ave / (double)trial);
	printf("using functional : All(10^%d) %15.15f[ms]\n", exp_, ave_all / (double)trial);
	printf("Size : %lld\n", sizeof(fparam));
}

/* フィルタ構造体
 *   偶数次/偶数次の場合の周波数特性確認用
 *
 */
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

void test_FilterParam_freq_res_so()
{
	vector<double> coef
	{
		-0.040659737,
		-2.372311969,

		-2.144646171,
		4.343497453,
		1.359348897,
		0.984834163,

		-0.710147059,
		-0.696696684,
		0.514853197,
		0.503697311,
		0.70680348
	};

	auto bands = FilterParam::gen_bands(FilterType::LPF, 0.3, 0.345);
	FilterParam fparam(5, 5, bands, 200, 50, 5.0);

	auto freq_res = fparam.freq_res(coef);

	for (auto band_res : freq_res)
	{
		for (auto res : band_res)
		{
			printf("%f\n", abs(res));
		}
	}
}

void test_FilterParam_freq_res_no()
{
	vector<double> coef
	{
		0.025247504683641238,

		0.8885952985540255,
		-4.097963802039866,
		5.496940685423355,
		0.3983519261092186,
		0.9723236917140877,
		1.1168784833810899,
		0.8492039597182939,

		-0.686114259307724,
		0.22008381076439384,
		-0.22066728558327908,
		0.7668032045079851
	};
	auto bands = FilterParam::gen_bands(FilterType::LPF, 0.2, 0.275);
	FilterParam Fparam(7,4,bands,200,50,5.0);

	vector<vector<complex<double>>> freq=Fparam.freq_res(coef);

	for(auto band:freq)
	{
		for(auto amp:band)
		{
			printf("%f\n",abs(amp));
		}
	}
}

void test_FilterParam_freq_res_mo()
{
	vector<double> coef	
	{
		-0.040404875,

		0.957674103,
		0.765466003,
		-1.585891794,
		-1.903482473,
		-0.441904071,
		0.79143639,
		-1.149627531,
		0.965348065,
		
		-0.434908839,
		-1.332562129,
		0.838349784
	};
	auto bands = FilterParam::gen_bands(FilterType::LPF, 0.1, 0.145);
	FilterParam fparam(8, 3, bands, 200, 50, 5.0);

	auto freq_res = fparam.freq_res(coef);

	for(auto band_res :freq_res)
	{
		for(auto res :band_res)
		{
			printf("%f\n", abs(res));
		}
	}
}

void test_FilterParam_judge_stability_even()
{
	auto bands = FilterParam::gen_bands(FilterType::LPF, 0.2, 0.275);
	FilterParam fparam(7, 4, bands, 200, 50, 5.0);

	vector<double> coef_1	//安定テスト
	{
		0.025247504683641238,

		0.8885952985540255,
		-4.097963802039866,
		5.496940685423355,
		0.3983519261092186,
		0.9723236917140877,
		1.1168784833810899,
		0.8492039597182939,

		0.686114259307724,
		0.22008381076439384,
		-0.22066728558327908,
		0.7668032045079851
	};
	double penalty = fparam.judge_stability(coef_1);
	printf("stable %f\n", penalty);

	vector<double> coef_2	//b_2>1の時
	{
		0.025247504683641238,

		0.8885952985540255,
		-4.097963802039866,
		5.496940685423355,
		0.3983519261092186,
		0.9723236917140877,
		1.1168784833810899,
		0.8492039597182939,

		0.686114259307724,
		1.22008381076439384,
		-0.22066728558327908,
		0.7668032045079851
	};
	penalty = fparam.judge_stability(coef_2);
	printf("unstable(b_2>1) %f\n", penalty);

	vector<double> coef_3	//b_1-1>b_2の時
	{
		0.025247504683641238,

		0.8885952985540255,
		-4.097963802039866,
		5.496940685423355,
		0.3983519261092186,
		0.9723236917140877,
		1.1168784833810899,
		0.8492039597182939,

		1.686114259307724,
		0.22008381076439384,
		-0.22066728558327908,
		0.7668032045079851
	};
	penalty = fparam.judge_stability(coef_3);
	printf("unstable(b_1-1>b_2) %f\n", penalty);

	vector<double> coef_4	//両方不安定の場合
	{
		0.025247504683641238,

		0.8885952985540255,
		-4.097963802039866,
		5.496940685423355,
		0.3983519261092186,
		0.9723236917140877,
		1.1168784833810899,
		0.8492039597182939,

		2.686114259307724,
		1.22008381076439384,
		-0.22066728558327908,
		0.7668032045079851
	};
	
	penalty = fparam.judge_stability(coef_4);
	printf("unstable(b_2>1,b_1-1>b_2) %f\n", penalty);
}

void test_FilterParam_judge_stability_odd()
{
	auto bands = FilterParam::gen_bands(FilterType::LPF, 0.1, 0.145);
	FilterParam fparam(8, 3, bands, 200, 50, 5.0);
	double penalty = 0.0;

	vector<double> coef_test1
	{
		-0.040404875,

		0.957674103,
		0.765466003,
		-1.585891794,
		-1.903482473,
		-0.441904071,
		0.79143639,
		-1.149627531,
		0.965348065,
		
		-0.434908839,
		-1.332562129,
		0.838349784
	};

	penalty = fparam.judge_stability(coef_test1);
	printf("stability %f\n", penalty);

	vector<double> coef_test2
	{
		-0.040404875,

		0.957674103,
		0.765466003,
		-1.585891794,
		-1.903482473,
		-0.441904071,
		0.79143639,
		-1.149627531,
		0.965348065,
		
		-1.434908839,
		-1.332562129,
		0.838349784
	};

	penalty = fparam.judge_stability(coef_test2);
	printf("instability %f\n", penalty);

	vector<double> coef_test3
	{
		-0.040404875,

		0.957674103,
		0.765466003,
		-1.585891794,
		-1.903482473,
		-0.441904071,
		0.79143639,
		-1.149627531,
		0.965348065,
		
		-0.434908839,
		-1.332562129,
		1.838349784
	};

	penalty = fparam.judge_stability(coef_test3);
	printf("instability %f\n", penalty);

	vector<double> coef_test4
	{
		-0.040404875,

		0.957674103,
		0.765466003,
		-1.585891794,
		-1.903482473,
		-0.441904071,
		0.79143639,
		-1.149627531,
		0.965348065,
		
		-0.434908839,
		-1.332562129,
		0.238349784
	};

	penalty = fparam.judge_stability(coef_test4);
	printf("instability %f\n", penalty);

	vector<double> coef_test5
	{
		-0.040404875,

		0.957674103,
		0.765466003,
		-1.585891794,
		-1.903482473,
		-0.441904071,
		0.79143639,
		-1.149627531,
		0.965348065,
		
		-0.434908839,
		-2.332562129,
		1.238349784
	};

	penalty = fparam.judge_stability(coef_test5);
	printf("instability %f\n", penalty);
}

void test_FilterParam_evaluate_objective_function()
{
	vector<double> coef
	{
		0.025247504683641238,

		0.8885952985540255,
		-4.097963802039866,
		5.496940685423355,
		0.3983519261092186,
		0.9723236917140877,
		1.1168784833810899,
		0.8492039597182939,

		-0.686114259307724,
		0.22008381076439384,
		-0.22066728558327908,
		0.7668032045079851
	};

	auto bands = FilterParam::gen_bands(FilterType::LPF, 0.2, 0.275);
	FilterParam Fparam(7,4,bands,200,50,5.0);

	auto objective_function_value = Fparam.evaluate(coef);

	printf("objective_function_value %f\n",objective_function_value);
}

void test_FilterParam_init_coef()
{
	auto bands = FilterParam::gen_bands(FilterType::LPF, 0.2, 0.3);
	FilterParam fparam(8, 9, bands, 200, 50, 5.0);
	double a0 = 0.5;
	double a = 3.0;
	double b = 3.0;

	for (unsigned int i = 0; i < 20; ++i)
	{
		auto coef = fparam.init_coef(a0, a, b);

		for (auto c :coef)
		{
			printf("% 3.3f ", c);
		}
		printf("\n");
	}
}

void test_FilterParam_init_stable_coef()
{
	auto bands = FilterParam::gen_bands(FilterType::LPF, 0.2, 0.3);
	FilterParam fparam(2, 9, bands, 200, 50, 5.0);
	double a0 = 0.5;
	double a = 3.0;

	for (unsigned int i = 0; i < 20; ++i)
	{
		auto coef = fparam.init_stable_coef(a0, a);

		printf("[ %3.3f]", fparam.judge_stability(coef));
		for (auto c :coef)
		{
			printf("% 3.3f ", c);
		}
		printf("\n");
	}
}
