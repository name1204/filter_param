/*
 * filter_param.cpp
 *
 *  Created on: 2021/07/22
 *      Author: matsu
 */

#include "filter_param.hpp"

using namespace std;

FILE *fileopen(const string &filename, const char mode, const string &call_file, const int call_line)
{
	FILE *fp = fopen(filename.c_str(), &mode);
	if (fp == NULL)
	{
		fprintf(stderr,
				"Error: [%s l.%d]Can't open file.(file name : %s, mode : %c)\n",
				call_file.c_str(), call_line, filename.c_str(), mode);
		exit(EXIT_FAILURE);
	}
	return fp;
}

vector<string> split_char(string &str, char spliter)
{
	stringstream ss{str};
	vector<std::string> strs;
	string buf;
	while (std::getline(ss, buf, spliter))
	{
		strs.emplace_back(buf);
	}
	return strs;
}

string BandParam::sprint()
{
	string str;
	switch (band_type)
	{
		case BandType::Pass:
		{
			str = format("PassBand(%6.3f, %6.3f)", left_side, right_side);
			break;
		}
		case BandType::Stop:
		{
			str = format("StopBand(%6.3f, %6.3f)", left_side, right_side);
			break;
		}
		case BandType::Transition:
		{
			str = format("TransitionBand(%6.3f, %6.3f)", left_side, right_side);
			break;
		}
	}
	return str;
}

FilterParam::FilterParam
(unsigned int zero, unsigned int pole, vector<BandParam> input_bands,
	unsigned int input_nsplit_approx, unsigned int input_nsplit_transition,
	double gd)
:n_order(zero), m_order(pole), bands(input_bands),
 nsplit_approx(input_nsplit_approx), nsplit_transition(input_nsplit_transition),
 group_delay(gd),
 threshold_riple(1.0)
{
	// 周波数帯域の整合性チェック
	double band_left = 0.0;
	for (auto bp : bands)
	{
		if (bp.left() != band_left)
		{
			fprintf(stderr, "Error: [%s l.%d]Adjacent band edge is necessary same.\n"
							"Or first band of left side must be 0.0.\n",
					__FILE__, __LINE__);
			for (auto bp : bands)
			{
				fprintf(stderr, "%s\n", bp.sprint().c_str());
			}
			exit(EXIT_FAILURE);
		}
		band_left = bp.right();
	}
	if (bands.back().right() != 0.5)
	{
		fprintf(stderr,
				"Error: [%s l.%d]Last band of right side must be 0.5.\n",
				__FILE__, __LINE__);
		for (auto bp : bands)
		{
			fprintf(stderr, "%s\n", bp.sprint().c_str());
		}
		exit(EXIT_FAILURE);
	}

	// 帯域ごとの分割数算出
	double approx_range = 0.0;
	double transition_range = 0.0;

	for (auto bp : bands)
	{
		switch (bp.type())
		{
		case BandType::Pass:
		{
			approx_range += bp.width();
			break;
		}
		case BandType::Stop:
		{
			approx_range += bp.width();
			break;
		}
		case BandType::Transition:
		{
			transition_range += bp.width();
			break;
		}
		}
	}

	vector<unsigned int> split;
	split.reserve(bands.size());
	for (auto bp : bands)
	{
		switch (bp.type())
		{
			case BandType::Pass:
			{
				split.emplace_back(
					(unsigned int)((double)nsplit_approx * bp.width() / approx_range));
				break;
			}
			case BandType::Stop:
			{
				split.emplace_back(
					(unsigned int)((double)nsplit_approx * bp.width() / approx_range));
				break;
			}
			case BandType::Transition:
			{
				split.emplace_back(
					(unsigned int)((double)nsplit_transition * bp.width() / transition_range));
				break;
			}
		}
	}
	split.at(0) += 1;

	// generate complex sin wave(e^-jω)
	// desire frequency response
	csw.reserve(bands.size());
	csw2.reserve(bands.size());
	desire_res.reserve(bands.size());
	for (unsigned int i = 0; i < bands.size(); ++i)
	{
		csw.emplace_back(gen_csw(bands.at(i), split.at(i)));
		csw2.emplace_back(gen_csw2(bands.at(i), split.at(i)));
		desire_res.emplace_back(gen_desire_res(bands.at(i), split.at(i), group_delay));
	}

	// decide using function
	if ((n_order % 2) == 0)
	{
		if ((m_order % 2) == 0)
		{
			freq_res_func = bind(&FilterParam::freq_res_se, this, placeholders::_1);
			stability_func = bind(&FilterParam::judge_stability_even, this, placeholders::_1);
		}
		else
		{
			freq_res_func = bind(&FilterParam::freq_res_mo, this, placeholders::_1);
		}
	}
	else
	{
		if ((m_order % 2) == 0)
		{
			freq_res_func = bind(&FilterParam::freq_res_no, this, placeholders::_1);
			stability_func = bind(&FilterParam::judge_stability_even, this, placeholders::_1);
		}
		else
		{
			printf("so is unimplemented.\n");
		}
	}
}

/* # フィルタ構造体
 *   CSVファイルから所望特性を読み取る関数
 *   複数の所望特性を読み取り，フィルタ構造体の配列で返却
 *   CSVファイルの書式は以下の通り
 *
 *     No, Numerator, Denominator, State, GroupDelay, NsplitApprox, NsplitTransition
 *         No : 所望特性のナンバリング
 *         Numerator : 分子次数
 *         Denominator : 分母次数
 *         State : フィルタの特性情報(L.P.F.など)
 *         GroupDelay : 所望群遅延
 *         NsplitApprox : 近似帯域の分割数
 *         NsplitTransition : 遷移域の分割数
 *
 *   vectorのインデックス = Noなので
 *   "No"は読みこまない
 *   また，１行目はヘッダなので読み飛ばす
 *
 * # 引数
 * string& filename : CSVファイルのパス(exeからの相対パスでも可能)
 * # 返り値
 * vector<FilterParam> params : CSVファイルにある分，全部の所望特性を格納した配列
 */
vector<FilterParam> FilterParam::read_csv(string& filename)
{
	vector<FilterParam> filter_params;
	ifstream ifs(filename);
	if (!ifs)
	{
		fprintf(stderr,
				"Error: [%s l.%d]Can't open file.(file name : %s, mode : read)\n",
				__FILE__, __LINE__, filename.c_str());
		exit(EXIT_FAILURE);
	}

	string buf;
	getline(ifs, buf); // ヘッダ読み飛ばし
	while (getline(ifs, buf))
	{
		auto vals = split_char(buf, ',');
		auto type = FilterParam::analyze_type(vals.at(3));
		vector<BandParam> bands;

		switch (type)
		{
			case FilterType::LPF:
			{
				auto edges = FilterParam::analyze_edges(vals.at(3));
				if (edges.size() != 2)
				{
					fprintf(stderr,
							"Error: [%s l.%d]Format of filter state is illegal.(input : \"%s\")\n"
									"If you assign filter type L.P.F. , length of edges is only 2.\n",
							__FILE__, __LINE__, vals.at(3).c_str());
					exit(EXIT_FAILURE);
				}
				bands = FilterParam::gen_bands(FilterType::LPF, edges.at(0), edges.at(1));
				break;
			}
			default:
			{
				fprintf(stderr, 
					"Error: [%s l.%d]It has not been implement yet.(input : \"%s\")\n",
					__FILE__, __LINE__, vals.at(3).c_str());
				exit(EXIT_FAILURE);
			}
		}

		filter_params.emplace_back(
			FilterParam(atoi(vals.at(1).c_str()), atoi(vals.at(2).c_str()),
							bands, atoi(vals.at(5).c_str()),
							atoi(vals.at(6).c_str()), atof(vals.at(4).c_str())));
	}

	return filter_params;
}

/* フィルタのタイプを簡易入力(文字列)する
 *   LPF(0.2, 0.3) : L.P.F.で通過域端 0.2, 阻止域端 0.3のフィルタ
 *
 *   Other(path) : その他の特性のフィルタ。pathは周波数帯域の情報を記したファイルへの相対パス
 */
FilterType FilterParam::analyze_type(const string& input)
{
	FilterType type;
	string str = input;
	str.erase(remove(str.begin(), str.end(), ' '), str.end());
	string csv_space = regex_replace(str, regex("[(),:]+"), " ");
	auto input_type = split_char(csv_space, ' ');

	if (input_type.size() <= 1)
	{
		fprintf(stderr,
				"Error: [%s l.%d]Format of filter type is illegal.(input : \"%s\")\n",
				__FILE__, __LINE__, input.c_str());
		exit(EXIT_FAILURE);
	}
	else if (input_type.at(0) == "LPF")
	{
		type = FilterType::LPF;
	}
	else if (true)	// hpf, bpf, bef バリエーション
	{
		fprintf(stderr,
				"Error: [%s l.%d]It has not been implement yet.(input : \"%s\")\n",
				__FILE__, __LINE__, input.c_str());
		exit(EXIT_FAILURE);
	}
	else
	{
		fprintf(stderr,
				"Error: [%s l.%d]Filter Type is undefined.(input : \"%s\")\n",
				__FILE__, __LINE__, input.c_str());
		exit(EXIT_FAILURE);
	}

	return type;
}

/* # フィルタ構造体
 *   フィルタの周波数帯域端を文字列から分離する
 *   LPF(0.2, 0.3) : 帯域域端 0.2, 帯域端 0.3
 *
 * # 引数
 * string& input : 解析する文字列
 * # 返り値
 * vector<double> edge : 解析・分離した周波数帯域端の配列
 */
vector<double> FilterParam::analyze_edges(const string& input)
{
	vector<double> edge;
	string str = input;
	str.erase(remove(str.begin(), str.end(), ' '), str.end());
	string csv_space = regex_replace(str, regex("[(),:]+"), " ");
	auto input_type = split_char(csv_space, ' ');
	input_type.erase(input_type.begin());

	if (input_type.size() <= 0)
	{
		fprintf(stderr,
				"Error: [%s l.%d]Format of filter state is illegal.(input : \"%s\")\n",
				__FILE__, __LINE__, input.c_str());
		exit(EXIT_FAILURE);
	}
	for (auto v : input_type)
	{
		edge.emplace_back(atof(v.c_str()));
	}
	return edge;
}

/* # フィルタ構造体
 *   複素正弦波の基本波(e^-jω)を生成する関数
 *   刻みは引数の周波数帯域幅と分割数に応じる
 *
 * # 引数
 * BandParam& bp : 周波数帯域情報をもった構造体
 * unsigned int nsplit : 周波数帯域の分割数
 * # 返り値
 * vector<complex<double>> csw : 引数の周波数帯域幅と分割数に応じた
 *                                        複素正弦波の配列
 *                                        csw = Complex Sin Wave
 *
 */
vector<complex<double>> FilterParam::gen_csw
(const BandParam& bp, const unsigned int nsplit)
{
	vector<complex<double>> csw;
		csw.reserve(nsplit);
	const double step_size = bp.width() / (double)nsplit;
	const double left = bp.left();
	constexpr double dpi = -2.0 * M_PI;			// double pi

	for (unsigned int i = 0; i < nsplit; ++i)
	{
		csw.emplace_back(polar(1.0, dpi * (left + step_size * (double) i)));
	}

	return csw;
}

/* # フィルタ構造体
 *   複素正弦波の第２次高調波(e^-j2ω)を生成する関数
 *   刻みは引数の周波数帯域幅と分割数に応じる
 *
 * # 引数
 * BandParam& bp : 周波数帯域情報をもった構造体
 * unsigned int nsplit : 周波数帯域の分割数
 * # 返り値
 * vector<complex<double>> csw2 : 引数の周波数帯域幅と分割数に応じた
 *                                        複素正弦波の配列
 *                                        csw2 = Complex Sin Wave 2(second)
 *
 */
vector<complex<double>> FilterParam::gen_csw2
(const BandParam& bp, const unsigned int nsplit)
{
	vector<complex<double>> csw2;
		csw2.reserve(nsplit);
	const double step_size = bp.width() / (double)nsplit;
	const double left = bp.left();
	constexpr double dpi = -4.0 * M_PI;			// quadrical pi

	for (unsigned int i = 0; i < nsplit; ++i)
	{
		csw2.emplace_back(polar(1.0, dpi * (left + step_size * (double) i)));
	}

	return csw2;
}

vector<complex<double>> FilterParam::gen_desire_res
(const BandParam& bp, const unsigned int nsplit, const double group_delay)
{
	vector<complex<double>> desire;

	switch (bp.type())
	{
		case BandType::Pass:
		{
			desire.reserve(nsplit);
			const double step_size = bp.width() / (double)nsplit;
			const double left = bp.left();
			constexpr double dpi = 2.0 * M_PI;			// double pi
			const double ang = dpi*group_delay;

			for (unsigned int i = 0; i < nsplit; ++i)
			{
				desire.emplace_back(polar(1.0, ang * (left + step_size * (double) i)));
			}
			break;
		}
		case BandType::Stop:
		{
			desire.resize(nsplit, 0.0);
			break;
		}
		case BandType::Transition:
			break;
	}

	return desire;
}

vector<vector<complex<double>>> FilterParam::freq_res_se(const vector<double>& coef) const
{
	vector<vector<complex<double>>> res;
		res.reserve(bands.size());

	for (unsigned int i = 0; i < bands.size(); ++i)    // 周波数帯域のループ
	{
		vector<complex<double>> band_res;
			band_res.reserve(csw.at(i).size());

		for (unsigned int j = 0; j < csw.at(i).size(); ++j)  // 周波数帯域内の分割数によるループ
		{
			complex<double> frac_over(1.0, 1.0);
			complex<double> frac_under(1.0, 1.0);

			for (unsigned int n = 1; n < n_order; n += 2)
			{
				frac_over *= 1.0 + coef.at(n)*csw.at(i).at(j) + coef.at(n + 1)*csw2.at(i).at(j);
			}
			for (unsigned int m = n_order + 1; m < opt_order(); m += 2)
			{
				frac_under *= 1.0 + coef.at(m)*csw.at(i).at(j) + coef.at(m + 1)*csw2.at(i).at(j);
			}
			band_res.emplace_back( coef.at(0)*(frac_over / frac_under) );
		}
		res.emplace_back(band_res);
	}

	return res;
}

vector<vector<complex<double>>> FilterParam::freq_res_no(const vector<double>& coef) const
{
	vector<vector<complex<double>>> freq;
		freq.reserve(bands.size());

	for (unsigned int i = 0; i < bands.size(); ++i)    // 周波数帯域のループ
	{
		vector<complex<double>> freq_band;
			freq_band.reserve(csw.at(i).size());

		for (unsigned int j = 0; j < csw.at(i).size(); ++j)  // 周波数帯域内の分割数によるループ
		{
			complex<double> freq_denominator(1.0, 1.0);
			complex<double> freq_numerator(1.0, 1.0);

			freq_numerator *= 1.0 + coef.at(1)*csw.at(i).at(j);
			for (unsigned int n = 2; n < n_order; n += 2)		//分子の総乗ループ
			{
				freq_numerator *= 1.0 + coef.at(n)*csw.at(i).at(j) + coef.at(n + 1)*csw2.at(i).at(j);
			}
			for (unsigned int m = n_order + 1; m < opt_order(); m += 2)		//分母の総乗ループ
			{
				freq_denominator *= 1.0 + coef.at(m)*csw.at(i).at(j) + coef.at(m + 1)*csw2.at(i).at(j);
			}

			freq_band.emplace_back( coef.at(0)*(freq_numerator / freq_denominator));
		}
		freq.emplace_back(freq_band);
	}
	
	return freq;
}

vector<vector<complex<double>>> FilterParam::freq_res_mo(const vector<double>& coef) const
{
	vector<vector<complex<double>>> freq;
		freq.reserve(bands.size());

	for (unsigned int i = 0; i < bands.size(); ++i) // 周波数帯域のループ
	{
		vector<complex<double>> freq_band;
			freq_band.reserve(csw.at(i).size());

		for (unsigned int j = 0; j < csw.at(i).size(); ++j) // 周波数帯域内の分割数によるループ
		{
			complex<double> nume(1.0, 1.0);
			complex<double> deno(1.0, 1.0);

			for (unsigned int n = 1; n < n_order; n += 2)
			{
				nume *= 1.0 + coef.at(n)*csw.at(i).at(j) + coef.at(n + 1)*csw2.at(i).at(j);
			}
			deno *= 1.0 + coef.at(n_order + 1)*csw.at(i).at(j);
			for (unsigned int m = n_order + 2; m < opt_order(); m += 2)
			{
				deno *= 1.0 + coef.at(m)*csw.at(i).at(j) + coef.at(m + 1)*csw2.at(i).at(j);
				freq_band.emplace_back( coef.at(0)*(nume / deno) );
		}
		freq.emplace_back(freq_band);

	}
	return freq;
}

double FilterParam::judge_stability_even(const vector<double>& coef) const
{
	double penalty = 0.0;

	for (unsigned int i = n_order + 1; i < this->opt_order(); i += 2)
	{
		if(abs(coef.at(i+1)) >= 1 || coef.at(i+1) <= abs(coef.at(i)) - 1)
		{
			penalty += coef.at(i)*coef.at(i) + coef.at(i+1)*coef.at(i+1);
		}
	}
	return penalty;
	
double FilterParam::evaluate(const vector<double> &coef) // 目的関数地計算関数
{
	constexpr double cs = 100;	//安定性のペナルティの重み
	constexpr double ct = 100;	//振幅隆起のペナルティの重み

	double max_error = 0.0;	//最大誤差
	double max_riple = 0.0;	//振幅隆起のペナルティの値

	auto bands = FilterParam::gen_bands(FilterType::LPF, 0.2, 0.275);
	FilterParam Fparam(7,4,bands,200,50,5.0);
	//const BandParam& bp;	//初期化子が必要と言われますが、どのようにすればよいのかわかりません

	//double penalty_stability = penalty_stability(coef);
	vector<vector<complex<double>>> freq=Fparam.freq_res(coef);

	for (unsigned int i = 0; i < bands.size(); ++i)    // 周波数帯域のループ(L.P.F.なら３つ)
	{
		for (unsigned int j = 0; j < csw.at(i).size(); ++j)  // 周波数帯域内の分割数によるループ
		{
			/*switch (bp.type())	
			{
			case BandType::Pass:
			case BandType::Stop:
				{
					if(max_error < abs(1.0 - freq.at(i).at(j)))
						{
							max_error = abs(1.0 - freq.at(i).at(j));
						}  
				}
			case BandType::Transition:
				{
					if(max_riple > 1.0 )
						{
							max_riple = abs(freq.at(i).at(j));
						}
				}
			}*/
		}
	}
	return(max_error + ct*max_riple*max_riple /*+ cs*penalty_stability*/);
}
